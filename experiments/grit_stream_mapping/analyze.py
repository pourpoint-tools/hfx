"""Recover the historical directed-cycle population and measure native-cell overlap."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
from typing import TypedDict

import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
import rasterio
from matplotlib.patches import Rectangle
from pyproj import Transformer
from rasterio.windows import Window
from shapely import from_wkb, transform
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry
from shapely.strtree import STRtree

D8 = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}


class Member(TypedDict):
    cell: list[int]
    raw: int


class MeasuredMember(Member):
    inside_vector_mask: bool
    distance_m: float
    segment_distance_m: float
    reach_distance_m: float
    nearest_river_id: int
    nearest_source_id: str
    nearest_level: int


class Cycle(TypedDict):
    cycle_id: int
    members: list[MeasuredMember]
    length: int
    touching_members: int
    classification: str
    min_distance_m: float
    max_member_distance_m: float


def recover_cycles(records: list[list[Member]], expected: int = 294) -> list[list[Member]]:
    unique: dict[tuple[tuple[int, int], ...], list[Member]] = {}
    for members in records:
        cells = [tuple(m["cell"]) for m in members]
        if len(set(cells)) != len(cells):
            raise ValueError("repeated member")
        for i, member in enumerate(members):
            r, c = member["cell"]
            dr, dc = D8[member["raw"]]
            if (r + dr, c + dc) != cells[(i + 1) % len(cells)]:
                raise ValueError("directed closure mismatch")
        start = cells.index(min(cells))
        canonical = members[start:] + members[:start]
        unique[tuple((m["cell"][0], m["cell"][1]) for m in canonical)] = canonical
    if len(unique) != expected:
        raise ValueError(f"population mismatch: {len(unique)} != {expected}")
    return [unique[key] for key in sorted(unique)]


def overlap(cells: list[BaseGeometry], rivers: list[BaseGeometry]) -> list[float]:
    if not rivers:
        raise ValueError("no river coverage")
    tree = STRtree(rivers)
    distances: list[float] = []
    for cell in cells:
        nearest = tree.nearest(cell)
        assert nearest is not None
        distances.append(float(cell.distance(rivers[int(nearest)])))
    return distances


def run(evidence: Path, prior: Path, output: Path) -> None:
    initial_path = evidence / "hfx269-coverage/exhaustive_summary.json"
    continuation_path = evidence / "hfx269-coverage/edge_continuations.json"
    initial = json.loads(initial_path.read_text())
    continued = json.loads(continuation_path.read_text())
    cycles = recover_cycles(
        [e["cycle"] for e in initial["cycles"]] + [e["cycle"] for e in continued["results"] if e["stop"] == "cycle"]
    )
    tile = evidence / "hfx269-source/GRITv1.0_drainage_direction_W008100000N003900000_EPSG8857.tif"
    cache = prior / "independent_trace_arrays.npz"
    z = np.load(cache)
    mask = z["vector_mask"]
    r0, c0 = int(z["request_row"]), int(z["request_col"])
    assert str(z["crs"]) == "EPSG:8857"
    assert list(z["transform"][:6]) == [30, 0, -15000000 + 30 * c0, 0, -30, 8400000 - 30 * r0]
    projected = Transformer.from_crs(4326, 8857, always_xy=True)
    rows = pq.read_table(output / "rivers.parquet").to_pylist()
    rivers = [transform(from_wkb(row["geometry"]), projected.transform, interleaved=False) for row in rows]
    river_tree = STRtree(rivers)
    segments = [g for g, row in zip(rivers, rows, strict=True) if row["level"] == 0]
    reaches = [g for g, row in zip(rivers, rows, strict=True) if row["level"] == 1]

    def inside(r: int, c: int) -> bool:
        return bool(0 <= r - r0 < mask.shape[0] and 0 <= c - c0 < mask.shape[1] and mask[r - r0, c - c0])

    records: list[Cycle] = []
    with rasterio.open(tile) as ds:
        source_r0 = round((8400000 - ds.transform.f) / 30)
        source_c0 = round((ds.transform.c + 15000000) / 30)
        assert ds.crs == rasterio.crs.CRS.from_epsg(8857)
        assert ds.transform.a == 30 and ds.transform.e == -30
        assert ds.transform.b == ds.transform.d == 0
        assert ds.transform.c == -15000000 + 30 * source_c0
        assert ds.transform.f == 8400000 - 30 * source_r0
        for number, cycle in enumerate(cycles, 1):
            cells = [footprint(*m["cell"]) for m in cycle]
            distances = overlap(cells, rivers)
            seg_distances = overlap(cells, segments)
            reach_distances = overlap(cells, reaches)
            for member in cycle:
                r, c = member["cell"]
                assert 0 <= r - source_r0 < ds.height and 0 <= c - source_c0 < ds.width
                window = Window(c - source_c0, r - source_r0, 1, 1)  # ty: ignore[too-many-positional-arguments] -- attrs-generated signature
                raw = int(ds.read(1, window=window)[0, 0])
                assert raw == member["raw"], (r, c, raw, member)
            members: list[MeasuredMember] = []
            for m, cell, d, sd, rd in zip(cycle, cells, distances, seg_distances, reach_distances, strict=True):
                nearest = river_tree.nearest(cell)
                assert nearest is not None
                river = rows[int(nearest)]
                members.append(
                    MeasuredMember(
                        cell=m["cell"],
                        raw=m["raw"],
                        inside_vector_mask=inside(*m["cell"]),
                        distance_m=d,
                        segment_distance_m=sd,
                        reach_distance_m=rd,
                        nearest_river_id=river["id"],
                        nearest_source_id=river["source_id"],
                        nearest_level=river["level"],
                    )
                )
            touched = sum(d == 0 for d in distances)
            records.append(
                {
                    "cycle_id": number,
                    "members": members,
                    "length": len(cycle),
                    "touching_members": touched,
                    "classification": "all" if touched == len(cycle) else "some" if touched else "none",
                    "min_distance_m": min(distances),
                    "max_member_distance_m": max(distances),
                }
            )
        representatives = [
            next(r for r in records if r["classification"] == "all"),
            next(r for r in records if r["classification"] == "some"),
            max(records, key=lambda r: r["min_distance_m"]),
        ]
        # A junction candidate is selected by maximum distinct published segment overlap.
        junction = max(
            records,
            key=lambda rec: sum(any(footprint(*m["cell"]).intersects(g) for m in rec["members"]) for g in segments),
        )
        representatives.append(junction)
        farthest = representatives[2]
        center = footprint(*farthest["members"][0]["cell"]).centroid
        fig, ax = plt.subplots(figsize=(7, 7))
        for river, row in zip(rivers, rows, strict=True):
            if river.distance(center) < 14000:
                parts = river.geoms if river.geom_type == "MultiLineString" else [river]
                for part in parts:
                    x, y = part.xy
                    ax.plot(x, y, color="#0072b2" if row["level"] == 0 else "#e69f00", lw=0.6)
        ax.scatter([center.x], [center.y], color="red", marker="s", s=25)
        ax.set_xlim(center.x - 14000, center.x + 14000)
        ax.set_ylim(center.y - 14000, center.y + 14000)
        ax.set_aspect("equal")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_title("Cycle 4 (red, enlarged): published river context")
        ax.set_xlabel("EPSG:8857 easting (m); blue segments, orange reaches")
        ax.set_ylabel("EPSG:8857 northing (m)")
        fig.text(0.1, 0.01, "GRIT / Wortmann et al., CC BY-NC 4.0", fontsize=8)
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        fig.savefig(output / "farthest-context.png", dpi=160)
        plt.close(fig)
        for record, label in zip(representatives, ["all-members", "some-members", "farthest", "junction"], strict=True):
            draw(record, ds, source_r0, source_c0, rivers, rows, output / f"{label}.png")

    def digest(path: Path) -> str:
        with path.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()

    nonoverlap = [r["min_distance_m"] for r in records if r["classification"] == "none"]
    receipt = json.loads((output / "river-acquisition.json").read_text())
    bbox = receipt["query_bounds_epsg8857"]
    edge_clearance = min(
        min(
            footprint(*m["cell"]).bounds[0] - bbox[0],
            footprint(*m["cell"]).bounds[1] - bbox[1],
            bbox[2] - footprint(*m["cell"]).bounds[2],
            bbox[3] - footprint(*m["cell"]).bounds[3],
        )
        for cycle in cycles
        for m in cycle
    )
    assert max(r["max_member_distance_m"] for r in records) < edge_clearance
    result = {
        "cycles": records,
        "summary": {
            "cycles": len(records),
            "initial_cycles": len(initial["cycles"]),
            "continuation_records": sum(e["stop"] == "cycle" for e in continued["results"]),
            "members": sum(r["length"] for r in records),
            "classifications": {k: sum(r["classification"] == k for r in records) for k in ["all", "some", "none"]},
            "members_outside_mask": sum(not m["inside_vector_mask"] for r in records for m in r["members"]),
            "cycles_with_outside_members": sum(any(not m["inside_vector_mask"] for m in r["members"]) for r in records),
            "nonoverlap_min_m": min(nonoverlap) if nonoverlap else None,
            "nonoverlap_median_m": float(np.median(nonoverlap)) if nonoverlap else None,
            "nonoverlap_max_m": max(nonoverlap) if nonoverlap else None,
            "river_records": len(rows),
            "segment_records": len(segments),
            "reach_records": len(reaches),
            "query_edge_clearance_m": edge_clearance,
            "representative_cycle_ids": [r["cycle_id"] for r in representatives],
        },
        "inputs": {
            str(p): digest(p)
            for p in [initial_path, continuation_path, tile, cache, output / "manifest.json", output / "rivers.parquet"]
        },
    }
    # One directed cycle per line keeps the complete inventory reviewable without expanding every scalar.
    text = '{"cycles":[\n' + ",\n".join(json.dumps(cycle) for cycle in records)
    text += '\n],"summary":' + json.dumps(result["summary"]) + ',"inputs":' + json.dumps(result["inputs"]) + "}\n"
    (output / "inventory.json").write_text(text)
    print(json.dumps(result["summary"], indent=2))


def footprint(r: int, c: int) -> BaseGeometry:
    x, y = -15000000 + 30 * c, 8400000 - 30 * r
    return box(x, y - 30, x + 30, y)


def draw(
    record: Cycle,
    ds: Any,
    r0: int,
    c0: int,
    rivers: list[BaseGeometry],
    rows: list[dict[str, Any]],
    destination: Path,
) -> None:
    members = {tuple(m["cell"]) for m in record["members"]}
    minr, maxr = min(r for r, c in members) - 3, max(r for r, c in members) + 3
    minc, maxc = min(c for r, c in members) - 3, max(c for r, c in members) + 3
    fig, ax = plt.subplots(figsize=(8, 7))
    for r in range(minr, maxr + 1):
        for c in range(minc, maxc + 1):
            cell = footprint(r, c)
            x, y = cell.centroid.coords[0]
            ax.add_patch(
                Rectangle(
                    (x - 15, y - 15),
                    30,
                    30,
                    facecolor="#ffc4c4" if (r, c) in members else "white",
                    edgecolor="#cccccc",
                    lw=0.5,
                )
            )
            raw = int(ds.read(1, window=Window(col_off=c - c0, row_off=r - r0, width=1, height=1))[0, 0])  # ty: ignore[unknown-argument] -- attrs-generated Window signature
            if raw in D8:
                dr, dc = D8[raw]
                ax.arrow(
                    x,
                    y,
                    dc * 11,
                    -dr * 11,
                    width=0.35,
                    head_width=3,
                    length_includes_head=True,
                    color="#222222",
                    zorder=4,
                )
            ax.text(x - 12, y + 8, str(raw), fontsize=7, zorder=5)
    bounds = box(*footprint(maxr, minc).bounds[:2], *footprint(minr, maxc).bounds[2:])
    for river, row in zip(rivers, rows, strict=True):
        if river.intersects(bounds):
            parts = river.geoms if river.geom_type == "MultiLineString" else [river]
            for part in parts:
                x, y = part.xy
                ax.plot(
                    x,
                    y,
                    color="#0072b2" if row["level"] == 0 else "#e69f00",
                    lw=2 if row["level"] == 0 else 1,
                    alpha=0.8,
                )
    ax.set_xlim(bounds.bounds[0], bounds.bounds[2])
    ax.set_ylim(bounds.bounds[1], bounds.bounds[3])
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_xlabel("EPSG:8857 easting (m)")
    ax.set_ylabel("EPSG:8857 northing (m)")
    ax.set_title(
        f"Cycle {record['cycle_id']}: {record['classification']} members touch; min {record['min_distance_m']:.2f} m\nRows {minr}–{maxr}; columns {minc}–{maxc}"
    )
    fig.text(
        0.1,
        0.02,
        "Red: cycle members • black: raw GRASS codes/arrows • blue: segments • orange: reaches\nGRIT / Wortmann et al., CC BY-NC 4.0",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(destination, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path)
    parser.add_argument("prior", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.evidence, args.prior, args.output)
