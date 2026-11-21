"""Normalize original TDX vectors and derive evidence-backed catchment connectivity."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sqlite3
from collections import Counter, defaultdict, deque
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import UTC, datetime
from itertools import chain, groupby
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq
from pyproj import Geod
from shapely import MultiPolygon, Polygon, from_wkb, get_coordinates, to_wkb
from shapely.validation import explain_validity

from hfx.candidates import fingerprint
from hfx.tdx_source import read_layer
from hfx.vector_storage import write_vectors

BOUNDS = pa.struct([pa.field(name, pa.float32(), nullable=False) for name in ("xmin", "ymin", "xmax", "ymax")])
SOURCE_FIDS_TYPE = pa.list_(pa.field("item", pa.int64(), nullable=False))
CATCHMENT_DERIVATION = "river_local_associations_with_zero_area_virtual_stub_contraction"
ADAPTER_VERSION = "tdx-connectivity-1"
GRAPH_SCHEMA = pa.schema(
    [
        pa.field("id", pa.int64(), nullable=False),
        pa.field("level", pa.int16(), nullable=False),
        pa.field("upstream_ids", pa.list_(pa.int64()), nullable=False),
    ]
    + [pa.field("bbox_" + name, pa.float32(), nullable=False) for name in ("minx", "miny", "maxx", "maxy")]
)
ASSOCIATION_SCHEMA = pa.schema(
    [pa.field("river_id", pa.int64(), nullable=False), pa.field("catchment_id", pa.int64(), nullable=False)]
)
_GEOD = Geod(ellps="WGS84")


@dataclass(frozen=True)
class BasinConnectivityResult:
    """Derived catchment connectivity and exhaustive audit stats for one basin."""

    catchment_upstream: dict[int, list[int]]
    catchment_cuts: dict[int, int]
    contracted_stub_count: int
    max_stub_chain_depth: int
    contracted_confluence_degrees: dict[int, int]
    non_contracted_linkno0_count: int
    non_contracted_linkno0_headwater_count: int
    non_contracted_linkno0_isolated_count: int
    non_contracted_linkno0_internal_count: int
    non_contracted_other_count: int
    non_contracted_other_samples: list[dict[str, Any]]
    catchments_with_cuts_count: int
    affected_downstream_catchments_count: int


def _outward_bbox(bounds: tuple[float, float, float, float]) -> dict[str, float]:
    raw = np.asarray(bounds, dtype=np.float64)
    narrowed = raw.astype(np.float32)
    for index in range(4):
        if (index < 2 and narrowed[index] > raw[index]) or (index >= 2 and narrowed[index] < raw[index]):
            narrowed[index] = np.nextafter(narrowed[index], np.float32(-np.inf if index < 2 else np.inf))
    return dict(zip(("xmin", "ymin", "xmax", "ymax"), narrowed.tolist(), strict=True))


def _geodesic_polygon_area_km2(shape: Polygon | MultiPolygon) -> float:
    parts = list(shape.geoms) if isinstance(shape, MultiPolygon) else [shape]
    area = 0.0
    for part in parts:
        exterior = abs(_GEOD.polygon_area_perimeter(*part.exterior.xy)[0])
        holes = sum(abs(_GEOD.polygon_area_perimeter(*ring.xy)[0]) for ring in part.interiors)
        area += (exterior - holes) / 1e6
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Original polygon has no positive finite geodesic area")
    return area


def _validate_polygon_shape(shape: Any) -> Polygon | MultiPolygon:
    if shape.geom_type not in {"Polygon", "MultiPolygon"} or shape.is_empty or shape.has_z or shape.has_m:
        raise ValueError("Unsupported original geometry type, dimensions or empty geometry")
    coordinates = get_coordinates(shape)
    if (
        not np.all(np.isfinite(coordinates))
        or np.any(np.abs(coordinates[:, 0]) > 180)
        or np.any(np.abs(coordinates[:, 1]) > 90)
    ):
        raise ValueError("Original coordinate outside finite longitude/latitude domain")
    if not shape.is_valid:
        raise ValueError("Invalid original polygon requires explicit source decision")
    return cast(Polygon | MultiPolygon, shape)


def collect_polygon_fragments(
    basin: str,
    stream_id: int,
    fid_wkb_pairs: list[tuple[int, bytes]],
) -> tuple[bytes, float, dict[str, float], list[int]]:
    """Collect single-part POLYGON fragments for one (basin, streamID) into one MultiPolygon.

    Single-row groups preserve exact original WKB bytes untouched. Multi-row groups
    pack original polygon parts in ascending fid order without dissolving, simplifying,
    or repairing, and exhaustively reject any OGC-invalid MultiPolygon.
    """
    if not fid_wkb_pairs:
        raise ValueError(f"Empty polygon fragment group for basin={basin} streamID={stream_id}")
    ordered = sorted(fid_wkb_pairs, key=lambda pair: pair[0])
    fids = [fid for fid, _ in ordered]
    if len(set(fids)) != len(fids):
        raise ValueError(f"Duplicate fid in polygon group for basin={basin} streamID={stream_id}: {fids}")
    if len(ordered) == 1:
        wkb = ordered[0][1]
        shape = _validate_polygon_shape(from_wkb(wkb))
        return wkb, _geodesic_polygon_area_km2(shape), _outward_bbox(shape.bounds), fids

    parts: list[Polygon] = []
    for _, wkb in ordered:
        shape = _validate_polygon_shape(from_wkb(wkb))
        if isinstance(shape, MultiPolygon):
            parts.extend(shape.geoms)
        else:
            parts.append(shape)
    grouped = MultiPolygon(parts)
    if not grouped.is_valid:
        reason = explain_validity(grouped)
        raise ValueError(f"Invalid grouped MultiPolygon for basin={basin} streamID={stream_id} fids={fids}: {reason}")
    grouped_wkb = to_wkb(grouped, output_dimension=2, byte_order=1)
    return grouped_wkb, _geodesic_polygon_area_km2(grouped), _outward_bbox(grouped.bounds), fids


def is_contractible_virtual_stub(
    *,
    linkno: int,
    length: float | None,
    geom_length: float | None,
    ds_cont_area: float | None,
    us_cont_area: float | None,
    upstream_degree: int,
    downstream_degree: int,
) -> bool:
    """Return True iff all five zero-length/zero-area internal virtual stub checks pass.

    Fails closed on LINKNO == 0 or any null / non-finite attribute.
    """
    if linkno == 0:
        return False
    if upstream_degree < 1 or downstream_degree != 1:
        return False
    for val in (length, geom_length, ds_cont_area, us_cont_area):
        if val is None or isinstance(val, bool) or not isinstance(val, (int, float)) or not math.isfinite(val):
            return False
    return length == 0.0 and geom_length == 0.0 and ds_cont_area == us_cont_area


def derive_basin_catchment_connectivity(
    river_flow: dict[int, list[int]],
    river_to_catchment: dict[int, int],
    polygonless_attrs: dict[int, tuple[int, float | None, float | None, float | None, float | None]],
) -> BasinConnectivityResult:
    """Derive catchment flow graph and upstream boundary cuts for one verified river basin.

    `polygonless_attrs` maps each polygonless river_id to
    `(linkno, length, geom_length, ds_cont_area, us_cont_area)`.
    """
    downstream: dict[int, list[int]] = {rid: [] for rid in river_flow}
    in_degree: dict[int, int] = {}
    for target, origins in river_flow.items():
        in_degree[target] = len(origins)
        for origin in origins:
            downstream[origin].append(target)

    contracted_stubs: set[int] = set()
    linkno0_count = 0
    linkno0_hw = 0
    linkno0_iso = 0
    linkno0_int = 0
    other_count = 0
    other_samples: list[dict[str, Any]] = []

    for rid in river_flow:
        if rid in river_to_catchment:
            continue
        if rid not in polygonless_attrs:
            raise ValueError(f"Missing polygonless attributes for river_id={rid}")
        linkno, length, geom_length, ds_area, us_area = polygonless_attrs[rid]
        up_deg = len(river_flow[rid])
        down_deg = len(downstream[rid])
        if is_contractible_virtual_stub(
            linkno=linkno,
            length=length,
            geom_length=geom_length,
            ds_cont_area=ds_area,
            us_cont_area=us_area,
            upstream_degree=up_deg,
            downstream_degree=down_deg,
        ):
            contracted_stubs.add(rid)
        elif linkno == 0:
            linkno0_count += 1
            if up_deg == 0 and down_deg == 0:
                linkno0_iso += 1
            elif up_deg == 0:
                linkno0_hw += 1
            else:
                linkno0_int += 1
        else:
            other_count += 1
            if len(other_samples) < 25:
                other_samples.append(
                    {
                        "river_id": rid,
                        "LINKNO": linkno,
                        "Length": length,
                        "geom_length": geom_length,
                        "DSContArea": ds_area,
                        "USContArea": us_area,
                        "upstream_degree": up_deg,
                        "downstream_degree": down_deg,
                    }
                )

    queue: deque[int] = deque(rid for rid, deg in in_degree.items() if deg == 0)
    visited = 0
    catchment_upstream: dict[int, list[int]] = {}
    catchment_cuts: dict[int, int] = {}
    forwarded_catchments: dict[int, list[int]] = {}
    forwarded_cuts: dict[int, int] = {}
    stub_chain_depth: dict[int, int] = {}
    has_cut_in_closure: dict[int, bool] = {}
    confluence_degrees: Counter[int] = Counter()
    max_chain_depth = 0

    while queue:
        rid = queue.popleft()
        visited += 1
        inc_cats: list[int] = []
        inc_cuts = 0
        max_up_stub_depth = 0
        had_contracted_upstream = False
        for up in river_flow[rid]:
            if up in river_to_catchment:
                inc_cats.append(river_to_catchment[up])
            elif up in contracted_stubs:
                had_contracted_upstream = True
                inc_cats.extend(forwarded_catchments[up])
                inc_cuts += forwarded_cuts[up]
                max_up_stub_depth = max(max_up_stub_depth, stub_chain_depth[up])
            else:
                # Non-contracted polygonless river: never bridge; record 1 incoming cut.
                inc_cuts += 1

        if rid in river_to_catchment:
            cid = river_to_catchment[rid]
            if cid in inc_cats or len(set(inc_cats)) != len(inc_cats):
                raise ValueError(f"Duplicate or self-linked derived catchment edge at catchment_id={cid}")
            ordered_ups = sorted(inc_cats)
            catchment_upstream[cid] = ordered_ups
            catchment_cuts[cid] = inc_cuts
            if had_contracted_upstream:
                confluence_degrees[len(ordered_ups)] += 1
            has_cut_in_closure[cid] = inc_cuts > 0 or any(has_cut_in_closure[u] for u in ordered_ups)
        elif rid in contracted_stubs:
            depth = max_up_stub_depth + 1
            stub_chain_depth[rid] = depth
            max_chain_depth = max(max_chain_depth, depth)
            forwarded_catchments[rid] = inc_cats
            forwarded_cuts[rid] = inc_cuts
        # Non-contracted polygonless rivers do not forward upstream catchments or cuts.

        for down in downstream[rid]:
            in_degree[down] -= 1
            if in_degree[down] == 0:
                queue.append(down)

    if visited != len(river_flow):
        raise ValueError("Supplied river graph contains a cycle")

    return BasinConnectivityResult(
        catchment_upstream=catchment_upstream,
        catchment_cuts=catchment_cuts,
        contracted_stub_count=len(contracted_stubs),
        max_stub_chain_depth=max_chain_depth,
        contracted_confluence_degrees=dict(sorted(confluence_degrees.items())),
        non_contracted_linkno0_count=linkno0_count,
        non_contracted_linkno0_headwater_count=linkno0_hw,
        non_contracted_linkno0_isolated_count=linkno0_iso,
        non_contracted_linkno0_internal_count=linkno0_int,
        non_contracted_other_count=other_count,
        non_contracted_other_samples=other_samples,
        catchments_with_cuts_count=sum(1 for cuts in catchment_cuts.values() if cuts > 0),
        affected_downstream_catchments_count=sum(1 for affected in has_cut_in_closure.values() if affected),
    )


def normalized_batch(
    source: pa.RecordBatch,
    domain: Literal["catchments", "rivers"],
    basin: str,
    *,
    first_id: int,
    local_catchment_coverage: Literal["supplied", "unavailable"] = "supplied",
) -> pa.RecordBatch:
    """Map one fid-ordered native batch; numbering is dataset-local, not a source ID."""
    if not basin.isdecimal() or first_id <= 0:
        raise ValueError("Invalid basin identity or first dataset ID")
    common = [
        pa.field("id", pa.int64(), nullable=False),
        pa.field("level", pa.int16(), nullable=False),
        pa.field("parent_id", pa.int64()),
        pa.field("source_namespace", pa.string(), nullable=False),
        pa.field("source_id", pa.string(), nullable=False),
        pa.field("processing_basin", pa.string(), nullable=False),
        pa.field("bbox", BOUNDS, nullable=False),
        pa.field("geometry", pa.binary(), nullable=False),
    ]
    extra = (
        [
            pa.field("area_km2", pa.float32(), nullable=False),
            pa.field("up_area_km2", pa.float32()),
            pa.field("outlet_lon", pa.float64()),
            pa.field("outlet_lat", pa.float64()),
        ]
        if domain == "catchments"
        else [
            pa.field("local_catchment_coverage", pa.string(), nullable=False),
            pa.field("component_catchment_id", pa.int64()),
            pa.field("upstream_node_id", pa.int64()),
            pa.field("downstream_node_id", pa.int64()),
            pa.field("geometry_orientation", pa.int8(), nullable=False),
        ]
    )
    native = [field for field in source.schema if field.name != "geom"]
    if {field.name for field in native} & {field.name for field in common + extra}:
        raise ValueError("Native column collides with normalized field; explicit rename required")
    if domain == "catchments":
        if source.schema.field("streamID").type not in (pa.int32(), pa.int64()):
            raise ValueError("streamID must use a supported integer type")
        native = [
            pa.field(
                "streamID",
                pa.int64(),
                nullable=field.nullable,
                metadata={"hfx:source_mapping": "streamID INTEGER/MEDIUMINT to int64, lossless"},
            )
            if field.name == "streamID"
            else field
            for field in native
        ]
    schema = pa.schema(common + extra + native)
    rows: list[dict[str, object]] = []
    allowed = {"Polygon", "MultiPolygon"} if domain == "catchments" else {"LineString", "MultiLineString"}
    for offset, row in enumerate(source.to_pylist()):
        wkb = row.pop("geom")
        shape = from_wkb(wkb)
        if shape.geom_type not in allowed or shape.is_empty or shape.has_z or shape.has_m:
            raise ValueError("Unsupported original geometry type, dimensions or empty geometry")
        coordinates = get_coordinates(shape)
        if (
            not np.all(np.isfinite(coordinates))
            or np.any(np.abs(coordinates[:, 0]) > 180)
            or np.any(np.abs(coordinates[:, 1]) > 90)
        ):
            raise ValueError("Original coordinate outside finite longitude/latitude domain")
        row.update(
            id=first_id + offset,
            level=0,
            parent_id=None,
            source_namespace=f"NGA-TDX-Hydro-20230126/{basin}/{domain}",
            source_id=str(row["fid"]),
            processing_basin=basin,
            bbox=_outward_bbox(shape.bounds),
            geometry=wkb,
        )
        if domain == "catchments":
            if not shape.is_valid:
                raise ValueError("Invalid original polygon requires explicit source decision")
            row.update(area_km2=_geodesic_polygon_area_km2(shape), up_area_km2=None, outlet_lon=None, outlet_lat=None)
        else:
            row.update(
                component_catchment_id=None,
                upstream_node_id=None,
                downstream_node_id=None,
                geometry_orientation=0,
                local_catchment_coverage=local_catchment_coverage,
            )
        rows.append(row)
    return pa.RecordBatch.from_pylist(rows, schema=schema)


def source_relationships(
    rivers: list[tuple[int, int, int, int, int]], polygons: list[tuple[int, int]]
) -> tuple[dict[int, list[int]], list[tuple[int, int]]]:
    """Validate both supplied directions before mapping one processing basin.

    River tuples contain HFX ID, LINKNO, DSLINKNO, USLINKNO1, USLINKNO2.
    Polygon tuples contain HFX ID and native association streamID.
    No graph contraction through polygonless rivers and no cross-basin inference.
    """
    native = {link: identity for identity, link, _, _, _ in rivers}
    if len(native) != len(rivers) or len({row[0] for row in rivers}) != len(rivers):
        raise ValueError("Duplicate native river or dataset identity")
    downstream: set[tuple[int, int]] = set()
    upstream: set[tuple[int, int]] = set()
    for _, link, target, first, second in rivers:
        if link < 0:
            raise ValueError("Invalid native river identity")
        for other in (target, first, second):
            if other != -1 and other not in native:
                raise ValueError("Unresolved supplied river relationship")
            if other == link:
                raise ValueError("Self-linked supplied river relationship")
        if target != -1:
            downstream.add((link, target))
        for origin in (first, second):
            if origin != -1:
                edge = (origin, link)
                if edge in upstream:
                    raise ValueError("Duplicate supplied upstream relationship")
                upstream.add(edge)
    if upstream != downstream:
        raise ValueError("Supplied upstream and downstream relationships disagree")
    incoming: dict[int, list[int]] = {identity: [] for identity in native.values()}
    successors: dict[int, list[int]] = {link: [] for link in native}
    degree = dict.fromkeys(native, 0)
    for origin, target in downstream:
        incoming[native[target]].append(native[origin])
        successors[origin].append(target)
        degree[target] += 1
    pending = [link for link, count in degree.items() if count == 0]
    visited = 0
    while pending:
        origin = pending.pop()
        visited += 1
        for target in successors[origin]:
            degree[target] -= 1
            if degree[target] == 0:
                pending.append(target)
    if visited != len(native):
        raise ValueError("Supplied river graph contains a cycle")
    associations: list[tuple[int, int]] = []
    if len({identity for identity, _ in polygons}) != len(polygons):
        raise ValueError("Duplicate polygon feature identity")
    for identity, stream in polygons:
        if stream not in native:
            raise ValueError("Unmatched polygon streamID")
        associations.append((native[stream], identity))
    return {identity: sorted(values) for identity, values in incoming.items()}, sorted(associations)


def _catchment_schema_with_connectivity(base_schema: pa.Schema) -> pa.Schema:
    """Insert upstream_boundary_cuts and source_fids into a normalized TDX catchment schema."""
    fields = list(base_schema)
    insert_at = base_schema.get_field_index("outlet_lat") + 1
    extra = [
        pa.field("upstream_boundary_cuts", pa.int32(), nullable=True),
        pa.field("source_fids", SOURCE_FIDS_TYPE, nullable=False),
    ]
    return pa.schema(fields[:insert_at] + extra + fields[insert_at:])


def _tdx_notice(audit: dict[str, Any]) -> str:
    return (
        "Original NGA TDX-Hydro, CC BY-SA 4.0. Source release NGA-TDX-Hydro-20230126. "
        "GeoPackage headers removed; single-row catchment WKB and all river WKB retained unchanged. "
        f"In basins exported with single-part POLYGON geometry ({', '.join(sorted(audit['grouped_basins'])) or 'none'}), "
        f"{audit['grouped_multipolygon_count']} multi-part catchments ({audit['extra_fragment_rows_collapsed']} extra "
        "fragment rows) are collected into one OGC-valid MultiPolygon per (processing_basin, streamID) without "
        "dissolving, simplifying, or repairing coordinates; fid and source_id retain min(source_fids), and "
        "source_fids preserves every original fid in ascending order. "
        "Native fields, including raw area values, retained. Derived local area is WGS84 geodesic km2; "
        "upstream area and outlets unknown. Native zero-length rivers retained; river flow graph untouched. "
        "Catchment connectivity is derived from river connectivity and 1-to-1 local associations, contracting "
        f"{audit['contracted_stub_count']} exhaustively verified zero-length, zero-area internal virtual confluence "
        "stubs (Length == 0, geom.length == 0, DSContArea == USContArea, upstream_degree >= 1, downstream_degree == 1) "
        "and never bridging non-contracted polygonless rivers "
        f"({audit['non_contracted_linkno0_count']} LINKNO == 0 omissions and "
        f"{audit['non_contracted_other_count']} other nonzero polygonless rivers; "
        f"{audit['catchments_with_cuts_count']} catchments with upstream_boundary_cuts > 0, "
        f"{audit['affected_downstream_catchments_count']} affected downstream closures). "
        "Missing links do not classify terminals; processing-boundary truncation, coastal/dateline cuts, and "
        "unverified vector-raster consistency remain source limitations. "
        "Catchment streamID INTEGER/MEDIUMINT is stored losslessly as int64; original declarations remain in "
        "source-metadata.json. River local_catchment_coverage distinguishes supplied polygon sources from "
        "unavailable sources.\n"
    )


def _catchment_connectivity_declaration(completeness: Literal["complete", "partial"]) -> list[dict[str, Any]]:
    return [
        {
            "level": 0,
            "status": "supplied",
            "completeness": completeness,
            "upstream_cuts_column": "upstream_boundary_cuts",
            "derivation": CATCHMENT_DERIVATION,
        }
    ]


def _projected_rows(path: Path, names: list[str]) -> Iterator[dict[str, Any]]:
    """Project IPC fields before decompression/decoding; never load geometry here."""
    with ipc.open_file(str(path)) as reader:
        fields = [reader.schema.get_field_index(name) for name in names]
    with ipc.open_file(str(path), options=ipc.IpcReadOptions(included_fields=fields)) as reader:
        for index in range(reader.num_record_batches):
            yield from reader.get_batch(index).to_pylist()


def _aggregate_audit(
    basin_audits: dict[str, dict[str, Any]],
    *,
    grouped_basins: dict[str, dict[str, int]],
) -> dict[str, Any]:
    confluence_hist: Counter[int] = Counter()
    other_samples: list[dict[str, Any]] = []
    for entry in basin_audits.values():
        for deg, count in entry["contracted_confluence_degrees"].items():
            confluence_hist[int(deg)] += int(count)
        if len(other_samples) < 25:
            other_samples.extend(entry["non_contracted_other_samples"][: 25 - len(other_samples)])
    return {
        "derivation": CATCHMENT_DERIVATION,
        "grouped_basins": dict(sorted(grouped_basins.items())),
        "grouped_multipolygon_count": sum(v["grouped_multipolygon_count"] for v in grouped_basins.values()),
        "extra_fragment_rows_collapsed": sum(v["extra_fragment_rows_collapsed"] for v in grouped_basins.values()),
        "contracted_stub_count": sum(e["contracted_stub_count"] for e in basin_audits.values()),
        "max_stub_chain_depth": max((e["max_stub_chain_depth"] for e in basin_audits.values()), default=0),
        "contracted_confluence_degrees": {str(k): v for k, v in sorted(confluence_hist.items())},
        "non_contracted_linkno0_count": sum(e["non_contracted_linkno0_count"] for e in basin_audits.values()),
        "non_contracted_linkno0_headwater_count": sum(
            e["non_contracted_linkno0_headwater_count"] for e in basin_audits.values()
        ),
        "non_contracted_linkno0_isolated_count": sum(
            e["non_contracted_linkno0_isolated_count"] for e in basin_audits.values()
        ),
        "non_contracted_linkno0_internal_count": sum(
            e["non_contracted_linkno0_internal_count"] for e in basin_audits.values()
        ),
        "non_contracted_other_count": sum(e["non_contracted_other_count"] for e in basin_audits.values()),
        "non_contracted_other_samples": other_samples,
        "catchments_with_cuts_count": sum(e["catchments_with_cuts_count"] for e in basin_audits.values()),
        "affected_downstream_catchments_count": sum(
            e["affected_downstream_catchments_count"] for e in basin_audits.values()
        ),
        "per_basin": basin_audits,
    }


def _basin_result_to_dict(res: BasinConnectivityResult) -> dict[str, Any]:
    return {
        "contracted_stub_count": res.contracted_stub_count,
        "max_stub_chain_depth": res.max_stub_chain_depth,
        "contracted_confluence_degrees": {str(k): v for k, v in res.contracted_confluence_degrees.items()},
        "non_contracted_linkno0_count": res.non_contracted_linkno0_count,
        "non_contracted_linkno0_headwater_count": res.non_contracted_linkno0_headwater_count,
        "non_contracted_linkno0_isolated_count": res.non_contracted_linkno0_isolated_count,
        "non_contracted_linkno0_internal_count": res.non_contracted_linkno0_internal_count,
        "non_contracted_other_count": res.non_contracted_other_count,
        "non_contracted_other_samples": res.non_contracted_other_samples,
        "catchments_with_cuts_count": res.catchments_with_cuts_count,
        "affected_downstream_catchments_count": res.affected_downstream_catchments_count,
    }


def _write_graph_rows(
    writer: pq.ParquetWriter,
    rows: list[dict[str, Any]],
) -> None:
    if rows:
        writer.write_table(pa.Table.from_pylist(rows, schema=GRAPH_SCHEMA))


def produce(inventory_path: Path, output: Path) -> None:
    """Create private normalized staging from raw GeoPackages with derived catchment connectivity."""
    inventory = json.loads(inventory_path.read_bytes())
    if inventory["release"] != "NGA-TDX-Hydro-20230126":
        raise ValueError("Unsupported source release")
    sources = sorted(inventory["sources"], key=lambda item: item["basin"])
    basins = [item["basin"] for item in sources]
    if not basins or len(set(basins)) != len(basins):
        raise ValueError("Source processing basins must be nonempty and unique")
    for item in sources:
        for domain in ("catchments", "rivers"):
            entry = item.get(domain)
            if not isinstance(entry, dict):
                raise ValueError("Every basin must declare both source domains")
            if "status" in entry:
                if (
                    domain != "catchments"
                    or entry.get("status") != "unavailable"
                    or set(entry) != {"status", "reason"}
                    or not isinstance(entry["reason"], str)
                    or not entry["reason"].strip()
                ):
                    raise ValueError("Unavailable catchments require an explicit reason and no source path")
            elif not {"path", "sha256"}.issubset(entry):
                raise ValueError("Supplied source requires path and sha256")
    if all("status" in item["catchments"] for item in sources):
        raise ValueError("At least one supplied catchment source is required")

    before: dict[Path, str] = {}
    source_metadata: list[dict[str, object]] = []
    for item in sources:
        for domain in ("catchments", "rivers"):
            if item[domain].get("status") == "unavailable":
                source_metadata.append({"basin": item["basin"], "domain": domain, **item[domain]})
                continue
            path = Path(item[domain]["path"])
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            if digest != item[domain]["sha256"]:
                raise ValueError("Original source digest differs from acquisition inventory")
            before[path] = digest
            connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
            try:
                connection.row_factory = sqlite3.Row
                contents = [dict(row) for row in connection.execute("SELECT * FROM gpkg_contents")]
                layers = [row["table_name"] for row in contents if row["data_type"] == "features"]
                schemas = {}
                for layer in layers:
                    quoted = '"' + layer.replace('"', '""') + '"'
                    schemas[layer] = [dict(row) for row in connection.execute(f"PRAGMA table_info({quoted})")]
                source_metadata.append(
                    {
                        "basin": item["basin"],
                        "domain": domain,
                        "status": "supplied",
                        "sha256": digest,
                        "contents": contents,
                        "geometry_columns": [
                            dict(row) for row in connection.execute("SELECT * FROM gpkg_geometry_columns")
                        ],
                        "schemas": schemas,
                    }
                )
            finally:
                connection.close()

    output.mkdir(parents=True, exist_ok=False)
    counts = {"catchments": 0, "rivers": 0}
    overall = [float("inf"), float("inf"), -float("inf"), -float("inf")]

    # First pass: write rivers.arrow (preserving exact river IDs and attributes).
    river_schema_ref: pa.Schema | None = None
    river_geom_lengths: dict[int, float] = {}

    def river_batches() -> Iterator[pa.RecordBatch]:
        nonlocal river_schema_ref
        for item in sources:
            basin = item["basin"]
            for incoming in read_layer(Path(item["rivers"]["path"]), "rivers"):
                mapped = normalized_batch(
                    incoming,
                    "rivers",
                    basin,
                    first_id=counts["rivers"] + 1,
                    local_catchment_coverage="unavailable" if "status" in item["catchments"] else "supplied",
                )
                river_schema_ref = mapped.schema
                counts["rivers"] += len(mapped)
                for row in mapped.select(
                    ["id", "bbox", "geometry", "LINKNO", "DSLINKNO", "USLINKNO1", "USLINKNO2"]
                ).to_pylist():
                    bbox = row["bbox"]
                    for index, name in enumerate(("xmin", "ymin", "xmax", "ymax")):
                        overall[index] = (
                            min(overall[index], bbox[name]) if index < 2 else max(overall[index], bbox[name])
                        )
                    values = [row[name] for name in ("LINKNO", "DSLINKNO", "USLINKNO1", "USLINKNO2")]
                    if any(type(value) is not int for value in values):
                        raise ValueError("Null or noninteger native connectivity")
                    river_geom_lengths[row["id"]] = float(from_wkb(row["geometry"]).length)
                yield mapped

    riv_incoming = river_batches()
    riv_first = next(riv_incoming)
    riv_directory = write_vectors(
        riv_first.schema,
        chain([riv_first], riv_incoming),
        output / "rivers.arrow",
        output / "river-lookup.arrow",
        lookup_binding="hfx:data_sha256",
    )
    (output / "river-lookup-directory.json").write_text(json.dumps(riv_directory) + "\n")

    # Second pass: per-basin polygon fragment collection + river/catchment connectivity derivation.
    grouped_basins: dict[str, dict[str, int]] = {}
    basin_audits: dict[str, dict[str, Any]] = {}
    has_unavailable_catchments = any("status" in item["catchments"] for item in sources)

    # Read projected river attributes grouped by basin.
    riv_proj_cols = ["processing_basin", "id", "bbox", "LINKNO", "DSLINKNO", "USLINKNO1", "USLINKNO2"]
    assert river_schema_ref is not None
    has_stub_attrs = all(col in river_schema_ref.names for col in ("Length", "DSContArea", "USContArea"))
    if has_stub_attrs:
        riv_proj_cols.extend(["Length", "DSContArea", "USContArea"])

    river_groups = groupby(
        _projected_rows(output / "rivers.arrow", riv_proj_cols),
        key=lambda row: row["processing_basin"],
    )

    cat_schema_final: pa.Schema | None = None

    with (
        pq.ParquetWriter(output / "graph.parquet", GRAPH_SCHEMA) as cat_graph,
        pq.ParquetWriter(output / "river-graph.parquet", GRAPH_SCHEMA) as riv_graph,
        pq.ParquetWriter(output / "local-associations.parquet", ASSOCIATION_SCHEMA) as local_assoc,
    ):

        def catchment_batches() -> Iterator[pa.RecordBatch]:
            nonlocal cat_schema_final
            for item in sources:
                basin = item["basin"]
                river_basin, basin_river_iter = next(river_groups, (None, iter(())))
                if river_basin != basin:
                    raise ValueError("Staging processing basin coverage differs")
                basin_rivers = list(basin_river_iter)

                riv_tuples = [
                    (r["id"], r["LINKNO"], r["DSLINKNO"], r["USLINKNO1"], r["USLINKNO2"]) for r in basin_rivers
                ]
                riv_bounds = {r["id"]: r["bbox"] for r in basin_rivers}

                if "status" in item["catchments"]:
                    rflow, _ = source_relationships(riv_tuples, [])
                    for start in range(0, len(basin_rivers), 4096):
                        chunk = basin_rivers[start : start + 4096]
                        _write_graph_rows(
                            riv_graph,
                            [
                                {
                                    "id": r["id"],
                                    "level": 0,
                                    "upstream_ids": rflow[r["id"]],
                                    **dict(
                                        zip(
                                            ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                            riv_bounds[r["id"]].values(),
                                            strict=True,
                                        )
                                    ),
                                }
                                for r in chunk
                            ],
                        )
                    continue

                # Collect raw polygons for this basin, grouped by streamID.
                by_stream: dict[int, list[tuple[int, bytes]]] = defaultdict(list)
                base_schema: pa.Schema | None = None
                for incoming in read_layer(Path(item["catchments"]["path"]), "catchments"):
                    if base_schema is None:
                        # Build canonical catchment schema from first batch's native schema.
                        sample = normalized_batch(incoming.slice(0, 1), "catchments", basin, first_id=1)
                        base_schema = _catchment_schema_with_connectivity(sample.schema)
                    for row in incoming.to_pylist():
                        sid = row["streamID"]
                        if type(sid) is not int:
                            raise ValueError("Null or noninteger streamID")
                        by_stream[sid].append((int(row["fid"]), bytes(row["geom"])))
                if base_schema is None or not by_stream:
                    raise ValueError(f"Empty catchment source for basin {basin}")
                if cat_schema_final is None:
                    cat_schema_final = base_schema

                # Order groups deterministically by min(fid) ascending.
                ordered_groups = sorted(by_stream.items(), key=lambda kv: min(f for f, _ in kv[1]))
                grouped_count = 0
                extra_rows = 0
                cat_rows_assembled: list[dict[str, Any]] = []
                polygon_pairs: list[tuple[int, int]] = []
                cat_bounds: dict[int, dict[str, float]] = {}

                for stream_id, fid_wkbs in ordered_groups:
                    counts["catchments"] += 1
                    cid = counts["catchments"]
                    if len(fid_wkbs) > 1:
                        grouped_count += 1
                        extra_rows += len(fid_wkbs) - 1
                    wkb, area, bbox, fids = collect_polygon_fragments(basin, stream_id, fid_wkbs)
                    for index, name in enumerate(("xmin", "ymin", "xmax", "ymax")):
                        overall[index] = (
                            min(overall[index], bbox[name]) if index < 2 else max(overall[index], bbox[name])
                        )
                    cat_bounds[cid] = bbox
                    polygon_pairs.append((cid, stream_id))
                    cat_rows_assembled.append(
                        {
                            "id": cid,
                            "level": 0,
                            "parent_id": None,
                            "source_namespace": f"NGA-TDX-Hydro-20230126/{basin}/catchments",
                            "source_id": str(fids[0]),
                            "processing_basin": basin,
                            "bbox": bbox,
                            "geometry": wkb,
                            "area_km2": area,
                            "up_area_km2": None,
                            "outlet_lon": None,
                            "outlet_lat": None,
                            "upstream_boundary_cuts": 0,
                            "source_fids": fids,
                            "fid": fids[0],
                            "streamID": stream_id,
                        }
                    )

                if grouped_count > 0:
                    grouped_basins[basin] = {
                        "grouped_multipolygon_count": grouped_count,
                        "extra_fragment_rows_collapsed": extra_rows,
                    }

                rflow, assocs = source_relationships(riv_tuples, polygon_pairs)
                river_to_cat: dict[int, int] = {}
                for rid, cid in assocs:
                    if rid in river_to_cat:
                        raise ValueError(f"River {rid} still has multiple catchments after grouping")
                    river_to_cat[rid] = cid

                polygonless_attrs: dict[int, tuple[int, float | None, float | None, float | None, float | None]] = {}
                for r in basin_rivers:
                    rid = r["id"]
                    if rid not in river_to_cat:
                        polygonless_attrs[rid] = (
                            r["LINKNO"],
                            r.get("Length") if has_stub_attrs else None,
                            river_geom_lengths[rid],
                            r.get("DSContArea") if has_stub_attrs else None,
                            r.get("USContArea") if has_stub_attrs else None,
                        )

                conn_res = derive_basin_catchment_connectivity(rflow, river_to_cat, polygonless_attrs)
                basin_audits[basin] = _basin_result_to_dict(conn_res)

                for row in cat_rows_assembled:
                    row["upstream_boundary_cuts"] = conn_res.catchment_cuts[row["id"]]

                # Write river-graph.parquet, local-associations.parquet, and graph.parquet.
                for start in range(0, len(basin_rivers), 4096):
                    chunk = basin_rivers[start : start + 4096]
                    _write_graph_rows(
                        riv_graph,
                        [
                            {
                                "id": r["id"],
                                "level": 0,
                                "upstream_ids": rflow[r["id"]],
                                **dict(
                                    zip(
                                        ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                        riv_bounds[r["id"]].values(),
                                        strict=True,
                                    )
                                ),
                            }
                            for r in chunk
                        ],
                    )
                for start in range(0, len(assocs), 4096):
                    local_assoc.write_table(
                        pa.Table.from_pylist(
                            [{"river_id": rid, "catchment_id": cid} for rid, cid in assocs[start : start + 4096]],
                            schema=ASSOCIATION_SCHEMA,
                        )
                    )
                for start in range(0, len(cat_rows_assembled), 4096):
                    chunk_rows = cat_rows_assembled[start : start + 4096]
                    _write_graph_rows(
                        cat_graph,
                        [
                            {
                                "id": row["id"],
                                "level": 0,
                                "upstream_ids": conn_res.catchment_upstream[row["id"]],
                                **dict(
                                    zip(
                                        ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                        cat_bounds[row["id"]].values(),
                                        strict=True,
                                    )
                                ),
                            }
                            for row in chunk_rows
                        ],
                    )
                    yield pa.RecordBatch.from_pylist(chunk_rows, schema=cat_schema_final)

            if next(river_groups, None) is not None:
                raise ValueError("Extra staging processing basin")

        cat_incoming = catchment_batches()
        cat_first = next(cat_incoming)
        assert cat_schema_final is not None
        cat_directory = write_vectors(
            cat_schema_final,
            chain([cat_first], cat_incoming),
            output / "catchments.arrow",
            output / "lookup.arrow",
            lookup_binding="hfx:catchments_sha256",
        )
        (output / "lookup-directory.json").write_text(json.dumps(cat_directory) + "\n")

    for path, digest in before.items():
        with path.open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != digest:
                raise ValueError("Original source changed during normalization")

    audit = _aggregate_audit(basin_audits, grouped_basins=grouped_basins)
    completeness: Literal["complete", "partial"] = (
        "partial" if (has_unavailable_catchments or audit["catchments_with_cuts_count"] > 0) else "complete"
    )
    (output / "original-inventory.json").write_bytes(inventory_path.read_bytes())
    (output / "source-metadata.json").write_text(
        json.dumps({"sources": source_metadata, "catchment_connectivity": audit}, indent=2) + "\n"
    )
    (output / "NOTICE").write_text(_tdx_notice(audit))

    manifest = {
        "format_version": "0.4.0",
        "fabric_name": "tdx",
        "fabric_version": inventory["release"],
        "crs": "EPSG:4326",
        "has_up_area": False,
        "topology": "tree",
        "bbox": overall,
        "unit_count": counts["catchments"],
        "created_at": datetime.now(UTC).isoformat(),
        "adapter_version": ADAPTER_VERSION,
        "region": ",".join(basins),
        "catchments": "catchments.arrow",
        "graph": "graph.parquet",
        "lookup": "lookup.arrow",
        "directory": "lookup-directory.json",
        "index_batch_rows": 65536,
        "catchment_batch_bytes": 8388608,
        "auxiliary": [],
        "categorical_columns": [
            {
                "artifact": "rivers.arrow",
                "column": "local_catchment_coverage",
                "physical_type": "string",
                "codes": [
                    {"code": "supplied", "meaning": "Polygon source supplied; local associations exhaustively checked"},
                    {"code": "unavailable", "meaning": "Polygon source unavailable; local associations unknown"},
                ],
                "null_meaning": "not permitted",
            }
        ],
        "catchment_connectivity": _catchment_connectivity_declaration(completeness),
        "rivers": {
            "data": "rivers.arrow",
            "graph": "river-graph.parquet",
            "lookup": "river-lookup.arrow",
            "directory": "river-lookup-directory.json",
            "count": counts["rivers"],
            "topology": "tree",
            "batch_bytes": 8388608,
            "index_batch_rows": 65536,
        },
        "local_associations": {
            "data": "local-associations.parquet",
            "river_column": "river_id",
            "catchment_column": "catchment_id",
        },
        "artifacts": {
            p.relative_to(output).as_posix(): fingerprint(p) for p in sorted(output.rglob("*")) if p.is_file()
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def _iter_parquet_batches(
    pf: pq.ParquetFile,
    *,
    columns: list[str] | None = None,
) -> Iterator[pa.RecordBatch]:
    """Stream row groups via read_row_group (avoids PyArrow iter_batches buffer retention)."""
    pool = pa.default_memory_pool()
    for rg in range(pf.num_row_groups):
        table = pf.read_row_group(rg, columns=columns)
        yield from table.to_batches()
        del table
        if (rg + 1) % 500 == 0:
            pool.release_unused()


def _read_polygonless_river_geom_lengths(
    rivers_parquet: Path,
    rivers_lookup: Path,
    polygonless_ids: set[int],
) -> dict[int, float]:
    """Read geometry ONLY for polygonless rivers using the river lookup index."""
    if not polygonless_ids:
        return {}
    by_block: dict[int, list[tuple[int, int]]] = defaultdict(list)
    with ipc.open_file(str(rivers_lookup)) as reader:
        for batch_idx in range(reader.num_record_batches):
            batch = reader.get_batch(batch_idx)
            ids = batch.column("id").to_numpy()
            if ids[-1] < min(polygonless_ids) or ids[0] > max(polygonless_ids):
                continue
            blocks = batch.column("block").to_numpy()
            rows = batch.column("row").to_numpy()
            for idx, rid in enumerate(ids.tolist()):
                if rid in polygonless_ids:
                    by_block[int(blocks[idx])].append((int(rows[idx]), rid))
    lengths: dict[int, float] = {}
    with pq.ParquetFile(rivers_parquet) as pf:
        for block, entries in sorted(by_block.items()):
            table = pf.read_row_group(block, columns=["id", "geometry"]).combine_chunks()
            id_col = table.column("id")
            geom_col = table.column("geometry")
            for row_pos, expected_rid in entries:
                actual_rid = id_col[row_pos].as_py()
                if actual_rid != expected_rid:
                    raise ValueError(f"River lookup mismatch for river_id={expected_rid}: got {actual_rid}")
                lengths[expected_rid] = float(from_wkb(geom_col[row_pos].as_py()).length)
    if len(lengths) != len(polygonless_ids):
        raise ValueError("Missing polygonless river geometries from rivers-lookup.arrow")
    return lengths


def produce_from_bundle(
    base_spatial_dir: Path,
    output: Path,
    *,
    verify_base_artifacts: bool = True,
) -> dict[str, Any]:
    """Derive a fresh normalized staging bundle with catchment connectivity from a verified TDX 0.4.0 spatial bundle."""
    base_manifest_path = base_spatial_dir / "manifest.json"
    base_manifest_bytes = base_manifest_path.read_bytes()
    base_manifest_sha = hashlib.sha256(base_manifest_bytes).hexdigest()
    base_manifest = json.loads(base_manifest_bytes)
    if (
        base_manifest.get("format_version") != "0.4.0"
        or base_manifest.get("fabric_name") != "tdx"
        or base_manifest.get("representation") != "spatial-geoparquet-adjacency"
        or base_manifest.get("representation_version") != 1
    ):
        raise ValueError("Expected verified HFX 0.4.0 TDX spatial-geoparquet-adjacency bundle")

    if verify_base_artifacts:
        for rel, identity in base_manifest["artifacts"].items():
            path = base_spatial_dir / rel
            if (
                Path(rel).is_absolute()
                or not path.resolve().is_relative_to(base_spatial_dir.resolve())
                or not path.is_file()
                or fingerprint(path) != identity
            ):
                raise ValueError(f"Base bundle artifact differs or is missing: {rel}")

    output.mkdir(parents=True, exist_ok=False)

    base_cat_path = base_spatial_dir / str(base_manifest["families"]["catchments"]["data"])
    base_riv_path = base_spatial_dir / str(base_manifest["families"]["rivers"]["data"])
    base_riv_lookup = base_spatial_dir / str(base_manifest["families"]["rivers"]["lookup"])

    # Pass 1a: scan catchment metadata (no geometry) to plan grouping and assign new catchment IDs.
    min_fid_by_key: dict[str, dict[int, int]] = defaultdict(dict)
    count_by_key: dict[str, dict[int, int]] = defaultdict(dict)
    with pq.ParquetFile(base_cat_path) as cat_pf:
        base_cat_schema_no_geo = cat_pf.schema_arrow.with_metadata(
            {k: v for k, v in (cat_pf.schema_arrow.metadata or {}).items() if k != b"geo"}
        )
        for batch in _iter_parquet_batches(cat_pf, columns=["processing_basin", "streamID", "fid"]):
            for basin, sid, fid in zip(
                cast(list[str], batch.column("processing_basin").to_pylist()),
                cast(list[int], batch.column("streamID").to_pylist()),
                cast(list[int], batch.column("fid").to_pylist()),
                strict=True,
            ):
                basin_min = min_fid_by_key[basin]
                basin_cnt = count_by_key[basin]
                if sid not in basin_min:
                    basin_min[sid] = fid
                    basin_cnt[sid] = 1
                else:
                    if fid < basin_min[sid]:
                        basin_min[sid] = fid
                    basin_cnt[sid] += 1

    stream_to_cid: dict[str, dict[int, int]] = {}
    multi_group_size: dict[int, int] = {}
    grouped_basins: dict[str, dict[str, int]] = {}
    next_cid = 0
    for basin in sorted(min_fid_by_key):
        basin_min = min_fid_by_key[basin]
        basin_cnt = count_by_key[basin]
        mapping: dict[int, int] = {}
        g_count = 0
        g_extra = 0
        for sid, _ in sorted(basin_min.items(), key=lambda kv: (kv[1], kv[0])):
            next_cid += 1
            mapping[sid] = next_cid
            cnt = basin_cnt[sid]
            if cnt > 1:
                multi_group_size[next_cid] = cnt
                g_count += 1
                g_extra += cnt - 1
        stream_to_cid[basin] = mapping
        if g_count > 0:
            grouped_basins[basin] = {
                "grouped_multipolygon_count": g_count,
                "extra_fragment_rows_collapsed": g_extra,
            }
    del min_fid_by_key, count_by_key
    new_cat_count = next_cid

    # Pass 1b: scan river metadata (no geometry) per basin.
    rivers_by_basin: dict[str, list[tuple[int, int, int, int, int, float | None, float | None, float | None]]] = (
        defaultdict(list)
    )
    basin_coverage: dict[str, str] = {}
    polygonless_ids: set[int] = set()
    overall = [float("inf"), float("inf"), -float("inf"), -float("inf")]

    with pq.ParquetFile(base_riv_path) as riv_pf:
        base_riv_schema_no_geo = riv_pf.schema_arrow.with_metadata(
            {k: v for k, v in (riv_pf.schema_arrow.metadata or {}).items() if k != b"geo"}
        )
        riv_count = riv_pf.metadata.num_rows
        riv_cols = [
            "id",
            "processing_basin",
            "LINKNO",
            "DSLINKNO",
            "USLINKNO1",
            "USLINKNO2",
            "Length",
            "DSContArea",
            "USContArea",
            "local_catchment_coverage",
            "bbox",
        ]
        for batch in _iter_parquet_batches(riv_pf, columns=riv_cols):
            for rid, basin, linkno, ds, u1, u2, length, ds_area, us_area, cov, bbox in zip(
                cast(list[int], batch.column("id").to_pylist()),
                cast(list[str], batch.column("processing_basin").to_pylist()),
                cast(list[int], batch.column("LINKNO").to_pylist()),
                cast(list[int], batch.column("DSLINKNO").to_pylist()),
                cast(list[int], batch.column("USLINKNO1").to_pylist()),
                cast(list[int], batch.column("USLINKNO2").to_pylist()),
                cast(list[float | None], batch.column("Length").to_pylist()),
                cast(list[float | None], batch.column("DSContArea").to_pylist()),
                cast(list[float | None], batch.column("USContArea").to_pylist()),
                cast(list[str], batch.column("local_catchment_coverage").to_pylist()),
                cast(list[dict[str, float]], batch.column("bbox").to_pylist()),
                strict=True,
            ):
                for index, name in enumerate(("xmin", "ymin", "xmax", "ymax")):
                    overall[index] = min(overall[index], bbox[name]) if index < 2 else max(overall[index], bbox[name])
                if basin not in basin_coverage:
                    basin_coverage[basin] = cov
                elif basin_coverage[basin] != cov:
                    raise ValueError(f"Inconsistent local_catchment_coverage in basin {basin}")
                rivers_by_basin[basin].append((rid, linkno, ds, u1, u2, length, ds_area, us_area))
                if cov == "supplied" and linkno not in stream_to_cid.get(basin, {}):
                    polygonless_ids.add(rid)

    # Pass 1c: read geometry ONLY for the ~4,077 polygonless rivers.
    poly_geom_lengths = _read_polygonless_river_geom_lengths(base_riv_path, base_riv_lookup, polygonless_ids)

    # Pass 1d: derive catchment connectivity per basin.
    catchment_upstream: list[list[int]] = [[] for _ in range(new_cat_count + 1)]
    catchment_cuts = np.zeros(new_cat_count + 1, dtype=np.int32)
    basin_link_to_rid: dict[str, dict[int, int]] = {}
    basin_audits: dict[str, dict[str, Any]] = {}
    has_unavailable_catchments = any(cov == "unavailable" for cov in basin_coverage.values())

    for basin in sorted(rivers_by_basin):
        b_rivers = rivers_by_basin[basin]
        riv_tuples = [(r[0], r[1], r[2], r[3], r[4]) for r in b_rivers]
        basin_link_to_rid[basin] = {r[1]: r[0] for r in b_rivers}
        b_stream_to_cid = stream_to_cid.get(basin, {})
        poly_pairs = [(cid, sid) for sid, cid in b_stream_to_cid.items()]
        rflow, assocs = source_relationships(riv_tuples, poly_pairs)
        if basin_coverage[basin] == "unavailable":
            if poly_pairs:
                raise ValueError(f"Unexpected catchments in unavailable basin {basin}")
            continue
        river_to_cat = dict(assocs)
        polygonless_attrs: dict[int, tuple[int, float | None, float | None, float | None, float | None]] = {}
        for rid, linkno, _, _, _, length, ds_area, us_area in b_rivers:
            if rid not in river_to_cat:
                polygonless_attrs[rid] = (linkno, length, poly_geom_lengths[rid], ds_area, us_area)
        res = derive_basin_catchment_connectivity(rflow, river_to_cat, polygonless_attrs)
        basin_audits[basin] = _basin_result_to_dict(res)
        for cid, ups in res.catchment_upstream.items():
            catchment_upstream[cid] = ups
        for cid, cuts in res.catchment_cuts.items():
            catchment_cuts[cid] = cuts
    del rivers_by_basin, poly_geom_lengths

    # Pass 2: stream catchments.parquet -> catchments.arrow, lookup.arrow, graph.parquet.
    new_cat_schema = _catchment_schema_with_connectivity(base_cat_schema_no_geo)
    in_flight: dict[int, list[tuple[int, bytes]]] = {}
    emitted_cats = 0
    base_cat_names = list(base_cat_schema_no_geo.names)
    insert_at = base_cat_names.index("outlet_lat") + 1
    prefix_indices = [base_cat_names.index(n) for n in base_cat_names[1:insert_at]]
    native_indices = [base_cat_names.index(n) for n in base_cat_names[insert_at:]]

    def _write_graph_arrays(
        writer: pq.ParquetWriter,
        cids: list[int],
        ups_list: list[list[int]],
        bbox_struct: pa.StructArray,
    ) -> None:
        n = len(cids)
        if not n:
            return
        table = pa.table(
            {
                "id": pa.array(cids, type=pa.int64()),
                "level": pa.array([0] * n, type=pa.int16()),
                "upstream_ids": pa.array(ups_list, type=pa.list_(pa.int64())),
                "bbox_minx": bbox_struct.field("xmin"),
                "bbox_miny": bbox_struct.field("ymin"),
                "bbox_maxx": bbox_struct.field("xmax"),
                "bbox_maxy": bbox_struct.field("ymax"),
            },
            schema=GRAPH_SCHEMA,
        )
        writer.write_table(table)

    with (
        pq.ParquetFile(base_cat_path) as cat_pf,
        pq.ParquetWriter(output / "graph.parquet", GRAPH_SCHEMA) as cat_graph,
    ):

        def derived_catchment_batches() -> Iterator[pa.RecordBatch]:
            nonlocal emitted_cats
            for batch in _iter_parquet_batches(cat_pf):
                basins_col = cast(list[str], batch.column("processing_basin").to_pylist())
                sids_col = cast(list[int], batch.column("streamID").to_pylist())
                fids_col = cast(list[int], batch.column("fid").to_pylist())
                ungrouped_idx: list[int] = []
                ungrouped_cids: list[int] = []
                ungrouped_cuts: list[int] = []
                ungrouped_fids: list[list[int]] = []
                grouped_idx: list[tuple[int, str, int, int, int]] = []
                for pos, (basin, sid, fid) in enumerate(zip(basins_col, sids_col, fids_col, strict=True)):
                    cid = stream_to_cid[basin][sid]
                    expected_size = multi_group_size.get(cid, 1)
                    if expected_size == 1:
                        ungrouped_idx.append(pos)
                        ungrouped_cids.append(cid)
                        ungrouped_cuts.append(int(catchment_cuts[cid]))
                        ungrouped_fids.append([fid])
                    else:
                        grouped_idx.append((pos, basin, sid, fid, cid))

                if ungrouped_idx:
                    sub = batch if len(ungrouped_idx) == batch.num_rows else batch.take(pa.array(ungrouped_idx))
                    bbox_col = cast(pa.StructArray, sub.column("bbox"))
                    for idx_b, name in enumerate(("xmin", "ymin", "xmax", "ymax")):
                        vals = bbox_col.field(name).to_numpy()
                        if len(vals):
                            overall[idx_b] = (
                                min(overall[idx_b], float(vals.min()))
                                if idx_b < 2
                                else max(overall[idx_b], float(vals.max()))
                            )
                    _write_graph_arrays(
                        cat_graph,
                        ungrouped_cids,
                        [catchment_upstream[c] for c in ungrouped_cids],
                        bbox_col,
                    )
                    cols = [
                        pa.array(ungrouped_cids, type=pa.int64()),
                        *[sub.column(i) for i in prefix_indices],
                        pa.array(ungrouped_cuts, type=pa.int32()),
                        pa.array(ungrouped_fids, type=SOURCE_FIDS_TYPE),
                        *[sub.column(i) for i in native_indices],
                    ]
                    emitted_cats += len(ungrouped_cids)
                    yield pa.RecordBatch.from_arrays(cols, schema=new_cat_schema)

                if grouped_idx:
                    geom_col = batch.column("geometry")
                    out_rows: list[dict[str, Any]] = []
                    graph_rows: list[dict[str, Any]] = []
                    for pos, basin, sid, fid, cid in grouped_idx:
                        expected_size = multi_group_size[cid]
                        buf = in_flight.setdefault(cid, [])
                        buf.append((fid, bytes(geom_col[pos].as_py())))
                        if len(buf) == expected_size:
                            del in_flight[cid]
                            wkb, area, bbox, fids = collect_polygon_fragments(basin, sid, buf)
                            for index, name in enumerate(("xmin", "ymin", "xmax", "ymax")):
                                overall[index] = (
                                    min(overall[index], bbox[name]) if index < 2 else max(overall[index], bbox[name])
                                )
                            out_rows.append(
                                {
                                    "id": cid,
                                    "level": 0,
                                    "parent_id": None,
                                    "source_namespace": f"NGA-TDX-Hydro-20230126/{basin}/catchments",
                                    "source_id": str(fids[0]),
                                    "processing_basin": basin,
                                    "bbox": bbox,
                                    "geometry": wkb,
                                    "area_km2": area,
                                    "up_area_km2": None,
                                    "outlet_lon": None,
                                    "outlet_lat": None,
                                    "upstream_boundary_cuts": int(catchment_cuts[cid]),
                                    "source_fids": fids,
                                    "fid": fids[0],
                                    "streamID": sid,
                                }
                            )
                            graph_rows.append(
                                {
                                    "id": cid,
                                    "level": 0,
                                    "upstream_ids": catchment_upstream[cid],
                                    **dict(
                                        zip(
                                            ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                            bbox.values(),
                                            strict=True,
                                        )
                                    ),
                                }
                            )
                            emitted_cats += 1
                    if out_rows:
                        _write_graph_rows(cat_graph, graph_rows)
                        yield pa.RecordBatch.from_pylist(out_rows, schema=new_cat_schema)
            if in_flight or emitted_cats != new_cat_count:
                raise ValueError(
                    f"Incomplete catchment grouping: emitted={emitted_cats}, expected={new_cat_count}, pending={len(in_flight)}"
                )

        cat_directory = write_vectors(
            new_cat_schema,
            derived_catchment_batches(),
            output / "catchments.arrow",
            output / "lookup.arrow",
            lookup_binding="hfx:catchments_sha256",
        )
        (output / "lookup-directory.json").write_text(json.dumps(cat_directory) + "\n")

    catchment_upstream.clear()
    multi_group_size.clear()

    # Pass 3: stream rivers.parquet -> rivers.arrow, river-lookup.arrow, river-graph.parquet, local-associations.parquet.
    with (
        pq.ParquetFile(base_riv_path) as riv_pf,
        pq.ParquetWriter(output / "river-graph.parquet", GRAPH_SCHEMA) as riv_graph,
        pq.ParquetWriter(output / "local-associations.parquet", ASSOCIATION_SCHEMA) as local_assoc,
    ):

        def derived_river_batches() -> Iterator[pa.RecordBatch]:
            for batch in _iter_parquet_batches(riv_pf):
                ids = cast(list[int], batch.column("id").to_pylist())
                basins_col = cast(list[str], batch.column("processing_basin").to_pylist())
                linknos = cast(list[int], batch.column("LINKNO").to_pylist())
                u1s = cast(list[int], batch.column("USLINKNO1").to_pylist())
                u2s = cast(list[int], batch.column("USLINKNO2").to_pylist())
                bboxes = cast(list[dict[str, float]], batch.column("bbox").to_pylist())
                graph_rows: list[dict[str, Any]] = []
                assoc_rows: list[dict[str, int]] = []
                for rid, basin, linkno, u1, u2, bbox in zip(ids, basins_col, linknos, u1s, u2s, bboxes, strict=True):
                    link_map = basin_link_to_rid[basin]
                    ups = sorted(link_map[u] for u in (u1, u2) if u != -1)
                    graph_rows.append(
                        {
                            "id": rid,
                            "level": 0,
                            "upstream_ids": ups,
                            **dict(
                                zip(
                                    ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                    bbox.values(),
                                    strict=True,
                                )
                            ),
                        }
                    )
                    cid = stream_to_cid.get(basin, {}).get(linkno)
                    if cid is not None:
                        assoc_rows.append({"river_id": rid, "catchment_id": cid})
                _write_graph_rows(riv_graph, graph_rows)
                if assoc_rows:
                    local_assoc.write_table(pa.Table.from_pylist(assoc_rows, schema=ASSOCIATION_SCHEMA))
                yield batch.replace_schema_metadata(base_riv_schema_no_geo.metadata)

        riv_directory = write_vectors(
            base_riv_schema_no_geo,
            derived_river_batches(),
            output / "rivers.arrow",
            output / "river-lookup.arrow",
            lookup_binding="hfx:data_sha256",
        )
        (output / "river-lookup-directory.json").write_text(json.dumps(riv_directory) + "\n")

    # Pass 4: preserve base provenance files and write updated source-metadata.json, NOTICE, manifest.json.
    for rel in base_manifest["artifacts"]:
        if rel in {"manifest.json", "source-manifest.json", "source-metadata.json", "NOTICE"}:
            continue
        suffix = Path(rel).suffix.lower()
        if suffix in {".arrow", ".parquet", ".tif", ".tiff"} or "directory" in rel:
            continue
        src = base_spatial_dir / rel
        dest = output / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)

    base_source_meta = json.loads((base_spatial_dir / "source-metadata.json").read_bytes())
    sources_list = base_source_meta["sources"] if isinstance(base_source_meta, dict) else base_source_meta
    audit = _aggregate_audit(basin_audits, grouped_basins=grouped_basins)
    completeness: Literal["complete", "partial"] = (
        "partial" if (has_unavailable_catchments or audit["catchments_with_cuts_count"] > 0) else "complete"
    )
    (output / "source-metadata.json").write_text(
        json.dumps(
            {
                "derived_from_verified_bundle_manifest_sha256": base_manifest_sha,
                "sources": sources_list,
                "catchment_connectivity": audit,
            },
            indent=2,
        )
        + "\n"
    )
    (output / "NOTICE").write_text(_tdx_notice(audit))

    manifest = {
        "format_version": "0.4.0",
        "fabric_name": "tdx",
        "fabric_version": base_manifest["fabric_version"],
        "crs": "EPSG:4326",
        "has_up_area": False,
        "topology": "tree",
        "bbox": overall,
        "unit_count": new_cat_count,
        "created_at": datetime.now(UTC).isoformat(),
        "adapter_version": ADAPTER_VERSION,
        "region": base_manifest["region"],
        "catchments": "catchments.arrow",
        "graph": "graph.parquet",
        "lookup": "lookup.arrow",
        "directory": "lookup-directory.json",
        "index_batch_rows": 65536,
        "catchment_batch_bytes": 8388608,
        "auxiliary": list(base_manifest.get("auxiliary", [])),
        "categorical_columns": [
            {
                "artifact": "rivers.arrow",
                "column": "local_catchment_coverage",
                "physical_type": "string",
                "codes": [
                    {"code": "supplied", "meaning": "Polygon source supplied; local associations exhaustively checked"},
                    {"code": "unavailable", "meaning": "Polygon source unavailable; local associations unknown"},
                ],
                "null_meaning": "not permitted",
            }
        ],
        "catchment_connectivity": _catchment_connectivity_declaration(completeness),
        "rivers": {
            "data": "rivers.arrow",
            "graph": "river-graph.parquet",
            "lookup": "river-lookup.arrow",
            "directory": "river-lookup-directory.json",
            "count": riv_count,
            "topology": "tree",
            "batch_bytes": 8388608,
            "index_batch_rows": 65536,
        },
        "local_associations": {
            "data": "local-associations.parquet",
            "river_column": "river_id",
            "catchment_column": "catchment_id",
        },
        "artifacts": {
            p.relative_to(output).as_posix(): fingerprint(p) for p in sorted(output.rglob("*")) if p.is_file()
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return {"manifest": manifest, "audit": audit}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Inventory JSON path, or base spatial bundle dir with --from-bundle")
    parser.add_argument("output", type=Path)
    parser.add_argument("--from-bundle", action="store_true", help="Derive staging from a verified spatial bundle")
    parser.add_argument("--skip-base-verify", action="store_true", help="Skip base artifact SHA-256 verification")
    args = parser.parse_args()
    if args.from_bundle:
        produce_from_bundle(args.source, args.output, verify_base_artifacts=not args.skip_base_verify)
    else:
        produce(args.source, args.output)


if __name__ == "__main__":
    main()
