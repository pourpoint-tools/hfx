"""Independent D8 graph checks and spatial comparison (never warp pointers)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, TypedDict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from affine import Affine
from matplotlib.colors import SymLogNorm
from numba import njit
from numpy.typing import NDArray
from process import basin_geometry, write_raster
from pyproj import Transformer
from rasterio.features import geometry_mask
from shapely.ops import transform

OFFSETS = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}


class GraphResult(TypedDict):
    accumulation: NDArray[np.int64]
    connected: NDArray[np.bool_]
    cycle_mask: NDArray[np.bool_]
    cycle_cells: int
    invalid_cells: int
    downstream: NDArray[np.int64]


@njit(cache=True)
def _walk(
    ds: NDArray[np.int64], valid: NDArray[np.bool_], outlet: int
) -> tuple[NDArray[np.int64], NDArray[np.bool_], NDArray[np.bool_]]:
    size = ds.size
    indegree = np.zeros(size, dtype=np.int32)
    accumulation = valid.astype(np.int64)
    for i in range(size):
        if ds[i] >= 0:
            indegree[ds[i]] += 1
    queue = np.empty(size, dtype=np.int64)
    end = 0
    for i in range(size):
        if valid[i] and indegree[i] == 0:
            queue[end] = i
            end += 1
    start = 0
    while start < end:
        i = queue[start]
        start += 1
        j = ds[i]
        if j >= 0:
            accumulation[j] += accumulation[i]
            indegree[j] -= 1
            if indegree[j] == 0:
                queue[end] = j
                end += 1
    connected = np.zeros(size, dtype=np.bool_)
    connected[outlet] = valid[outlet]
    for k in range(end - 1, -1, -1):
        i = queue[k]
        if ds[i] >= 0 and connected[ds[i]]:
            connected[i] = True
    return accumulation, connected, indegree > 0


def trace_d8(directions: NDArray[np.uint8], outlet: int) -> GraphResult:
    if not 0 <= outlet < directions.size:
        raise ValueError("Outlet index outside raster")
    height, width = directions.shape
    valid = directions != 247
    invalid = valid & ~np.isin(directions, [0, *OFFSETS])
    ds = np.full(directions.size, -1, dtype=np.int64)
    for code, (dr, dc) in OFFSETS.items():
        rows, cols = np.where(directions == code)
        rr, cc = rows + dr, cols + dc
        inside = (rr >= 0) & (rr < height) & (cc >= 0) & (cc < width)
        rows, cols, rr, cc = rows[inside], cols[inside], rr[inside], cc[inside]
        target_valid = valid[rr, cc] & ~invalid[rr, cc]
        ds[rows[target_valid] * width + cols[target_valid]] = rr[target_valid] * width + cc[target_valid]
    accumulation, connected, cycles = _walk(ds, (valid & ~invalid).ravel(), outlet)
    return {
        "accumulation": accumulation.reshape(directions.shape),
        "connected": connected.reshape(directions.shape),
        "cycle_mask": cycles.reshape(directions.shape),
        "cycle_cells": int(cycles.sum()),
        "invalid_cells": int(invalid.sum()),
        "downstream": ds,
    }


REFERENCE = (-89.32299563991558, 31.639357487809487)


def evaluate(base: Path, buffer_km: float, grit_npz: Path) -> None:
    output = base / "outputs" / f"buffer-{buffer_km:g}km"
    with rasterio.open(output / "direction.tif") as src:
        directions, affine, crs = src.read(1), src.transform, src.crs
    with rasterio.open(output / "accumulation_cells.tif") as src:
        counts = src.read(1)
    geometry = transform(
        Transformer.from_crs(4326, crs, always_xy=True).transform, basin_geometry(base / "inputs/coarse_union.geojson")
    )
    mask = geometry_mask([geometry], directions.shape, affine, invert=True)
    x, y = Transformer.from_crs(4326, crs, always_xy=True).transform(*REFERENCE)
    row, col = rasterio.transform.rowcol(affine, x, y)
    index = int(row * directions.shape[1] + col)
    result = trace_d8(directions, index)
    edge = np.zeros(directions.shape, dtype=bool)
    edge[[0, -1], :] = True
    edge[:, [0, -1]] = True
    # Additional terrain-only comparison: largest accumulation within 300 m.
    rr, cc = np.indices(directions.shape)
    xx, yy = affine.c + (cc + 0.5) * 30, affine.f - (rr + 0.5) * 30
    nearby = (xx - x) ** 2 + (yy - y) ** 2 <= 300**2
    adjusted_index = int(np.argmax(np.where(nearby, counts, -1)))
    adjusted = trace_d8(directions, adjusted_index)

    def outlet_metrics(idx: int, graph: GraphResult) -> dict[str, Any]:
        r, c = divmod(idx, directions.shape[1])
        center = rasterio.transform.xy(affine, r, c)
        lonlat = Transformer.from_crs(crs, 4326, always_xy=True).transform(*center)
        connected = graph["connected"]
        return {
            "row": r,
            "column": c,
            "center_lonlat": list(lonlat),
            "displacement_m": float(np.hypot(center[0] - x, center[1] - y)),
            "accumulation_cells": int(counts[r, c]),
            "area_km2": float(counts[r, c] * 0.0009),
            "vector_connected_cells": int((connected & mask).sum()),
            "vector_coverage_percent": float(100 * (connected & mask).sum() / mask.sum()),
            "outside_vector_cells": int((connected & ~mask).sum()),
            "boundary_connected_cells": int((connected & edge).sum()),
        }

    metrics = {
        "buffer_km": buffer_km,
        "vector_cells": int(mask.sum()),
        "cycle_cells": result["cycle_cells"],
        "invalid_cells": result["invalid_cells"],
        "accumulation_mismatch_cells": int((result["accumulation"] != counts).sum()),
        "reference": outlet_metrics(index, result),
        "adjusted_max_accumulation_within_300m": outlet_metrics(adjusted_index, adjusted),
    }
    write_raster(output / "reference_connected.tif", result["connected"].astype(np.uint8), affine, 255)
    write_raster(output / "adjusted_connected.tif", adjusted["connected"].astype(np.uint8), affine, 255)
    # Sample terrain values at GRIT native cell centers, without warping pointers.
    with np.load(grit_npz) as grit:
        gt = Affine(*grit["transform"][:6])
        gm = grit["vector_mask"]
        gr, gc = np.indices(gm.shape)
        gx, gy = gt.c + (gc + 0.5) * gt.a, gt.f + (gr + 0.5) * gt.e
        tx, ty = Transformer.from_crs(str(grit["crs"]), crs, always_xy=True).transform(gx, gy)
        tc = np.floor((tx - affine.c) / 30).astype(int)
        tr = np.floor((affine.f - ty) / 30).astype(int)
        inside = (tr >= 0) & (tr < directions.shape[0]) & (tc >= 0) & (tc < directions.shape[1])
        if not inside[gm].all():
            raise ValueError("Terrain does not contain full vector comparison mask")
        sample = np.full(gm.shape, np.nan)
        sample[inside] = counts[tr[inside], tc[inside]] * 0.0009
        comparable = gm & grit["accumulation_valid"] & inside
        difference = sample[comparable] - grit["accumulation"][comparable]
        gcodes = np.full(gm.shape, 247, dtype=np.uint8)
        for code, converted in enumerate([128, 64, 32, 16, 8, 4, 2, 1], start=1):
            gcodes[grit["direction_valid"] & (grit["raw_direction"] == code)] = converted
        ox, oy = Transformer.from_crs(4326, str(grit["crs"]), always_xy=True).transform(*REFERENCE)
        orow, ocol = rasterio.transform.rowcol(gt, ox, oy)
        gi = int(orow * gm.shape[1] + ocol)
        graph = trace_d8(gcodes, gi)
        metrics["grit_public_window"] = {
            "cycle_cells": graph["cycle_cells"],
            "invalid_cells": graph["invalid_cells"],
            "vector_cells": int(gm.sum()),
            "vector_connected_cells": int((graph["connected"] & gm).sum()),
            "outlet_accumulation_km2": float(grit["accumulation"].ravel()[gi]),
            "sampled_accumulation_difference_km2_percentiles_0_50_95_100": np.percentile(
                difference, [0, 50, 95, 100]
            ).tolist(),
        }
        fig, axes = plt.subplots(1, 3, figsize=(15, 8), layout="constrained")
        extent = (gt.c, gt.c + gm.shape[1] * 30, gt.f - gm.shape[0] * 30, gt.f)
        for ax, values, title in zip(
            axes,
            [grit["accumulation"], sample, sample - grit["accumulation"]],
            ["GRIT width-partitioned km²", "Terrain D8 inclusive km²", "Terrain minus GRIT km²"],
            strict=True,
        ):
            image = np.where(gm, values, np.nan)
            plot = ax.imshow(
                image,
                extent=extent,
                norm=SymLogNorm(linthresh=0.01, vmin=-2000 if "minus" in title else 0, vmax=2000),
                cmap="coolwarm" if "minus" in title else "viridis",
            )
            ax.set_title(title)
            ax.set_xlabel("EPSG:8857 easting (m)")
            fig.colorbar(plot, ax=ax, shrink=0.6)
        fig.savefig(base / "report" / f"accumulation-{buffer_km:g}km.png", dpi=160)
        plt.close(fig)
        terrain_connected = np.zeros(gm.shape, dtype=bool)
        terrain_connected[inside] = adjusted["connected"][tr[inside], tc[inside]]
        comparison = graph["connected"].astype(np.uint8) + 2 * terrain_connected
        fig, ax = plt.subplots(figsize=(8, 10), layout="constrained")
        plot = ax.imshow(
            np.where(gm, comparison, np.nan), extent=extent, vmin=0, vmax=3, cmap="viridis", interpolation="nearest"
        )
        ax.set_title("Outlet-connected vector cells: GRIT vs terrain adjusted outlet")
        bar = fig.colorbar(plot, ax=ax, ticks=[0, 1, 2, 3])
        bar.ax.set_yticklabels(["Neither", "GRIT only", "Terrain only", "Both"])
        fig.savefig(base / "report" / f"connectivity-{buffer_km:g}km.png", dpi=160)
        plt.close(fig)
        # Native-grid arrows near the diagnostic outlet; coordinates transformed, never pointer codes.
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
        for ax, codes, areas, native, native_crs, title in [
            (axes[0], gcodes, grit["accumulation"], gt, str(grit["crs"]), "GRIT native D8"),
            (axes[1], directions, counts * 0.0009, affine, str(crs), "Terrain native D8"),
        ]:
            px, py = Transformer.from_crs(4326, native_crs, always_xy=True).transform(*REFERENCE)
            r0, c0 = rasterio.transform.rowcol(native, px, py)
            rows, cols = np.mgrid[r0 - 12 : r0 + 13, c0 - 12 : c0 + 13]
            sx = native.c + (cols + 0.5) * native.a
            sy = native.f + (rows + 0.5) * native.e
            dx = np.zeros(rows.shape)
            dy = np.zeros(rows.shape)
            for code, (dr, dc) in OFFSETS.items():
                select = codes[rows, cols] == code
                dx[select] = dc * native.a
                dy[select] = dr * native.e
            conv = Transformer.from_crs(native_crs, crs, always_xy=True)
            axx, ayy = conv.transform(sx, sy)
            bxx, byy = conv.transform(sx + dx, sy + dy)
            ax.scatter(axx, ayy, c=np.log10(areas[rows, cols] + 0.0009), s=20, cmap="viridis")
            ax.quiver(axx, ayy, bxx - axx, byy - ayy, angles="xy", scale_units="xy", scale=1, width=0.002)
            ax.plot(x, y, "r+", markersize=12, label="reference coordinate")
            ar, ac = divmod(adjusted_index, directions.shape[1])
            ax.plot(*rasterio.transform.xy(affine, ar, ac), "m*", markersize=12, label="adjusted (300m search)")
            for boundary in (
                [geometry.exterior] if geometry.geom_type == "Polygon" else [g.exterior for g in geometry.geoms]
            ):
                bx, by = boundary.xy
                ax.plot(bx, by, color="orange", linewidth=1)
            ax.set_xlim(x - 400, x + 400)
            ax.set_ylim(y - 400, y + 400)
            ax.set_title(title)
            ax.set_aspect("equal")
            ax.legend()
            ax.set_xlabel("EPSG:32616 easting (m)")
        fig.savefig(base / "report" / f"outlet-directions-{buffer_km:g}km.png", dpi=160)
        plt.close(fig)
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("grit_npz", type=Path)
    parser.add_argument("--buffer-km", type=float, required=True)
    args = parser.parse_args()
    evaluate(args.base, args.buffer_km, args.grit_npz)
