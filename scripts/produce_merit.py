"""Normalize MERIT Basins vectors (+ optional MERIT Hydro rasters) with narrow ring normalization."""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import os
import platform
import shutil
import sqlite3
import sys
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import fiona
import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq
import shapely
from pyproj import Geod
from shapely import STRtree, from_wkb, get_coordinates, make_valid, prepare
from shapely import box as shapely_box
from shapely.geometry import shape as shapely_shape
from shapely.validation import explain_validity

from hfx.candidates import fingerprint
from hfx.vector_storage import write_vectors

BOUNDS = pa.struct([pa.field(name, pa.float32(), nullable=False) for name in ("xmin", "ymin", "xmax", "ymax")])
MIN_WGS84_LON = -180.0
MAX_UNWRAPPED_ANTIMERIDIAN_LON = 195.0
ANTIMERIDIAN_QUERY_MIN_X = -165.0
ANTIMERIDIAN_MAX_X = 180.0

GEOMETRY_NORMALIZATION_RULE = "merit_ring_self_intersection_lossless_split_v1"
GEOMETRY_NORMALIZATION_ALGORITHM = "shapely.make_valid(method='structure', keep_collapsed=False)"
GEODESIC_AREA_ABS_TOL_KM2 = 1e-9
GEODESIC_AREA_REL_TOL = 1e-12
REPAIRED_FEATURES_ARTIFACT = "repaired-features.jsonl.gz"

CATCHMENT_GEOMETRY_NORMALIZATION_CODES: list[dict[str, str]] = [
    {"code": "original", "meaning": "Source polygon was OGC-valid and preserved unchanged"},
    {
        "code": "ring_self_intersection_split",
        "meaning": (
            "Invalid MERIT Ring Self-intersection polygon normalized into valid rings preserving exact float64 "
            "vertices, undirected edges, geodesic area, and zero sibling interior overlap"
        ),
    },
]

STRUCTURAL_FIELDS = {
    "id",
    "level",
    "parent_id",
    "source_namespace",
    "source_id",
    "processing_partition",
    "catchment_kind",
    "geometry_normalization",
    "bbox",
    "geometry",
    "area_km2",
    "up_area_km2",
    "outlet_lon",
    "outlet_lat",
    "local_catchment_id",
    "local_catchment_status",
    "downstream_boundary_status",
    "upstream_boundary_cuts",
    "component_catchment_id",
    "upstream_node_id",
    "downstream_node_id",
    "geometry_orientation",
    "FID",
    "fid",
}

CATCHMENT_RAW_FIELDS: list[pa.Field] = [
    pa.field("id", pa.int64(), nullable=False),
    pa.field("level", pa.int16(), nullable=False),
    pa.field("parent_id", pa.int64()),
    pa.field("source_namespace", pa.string(), nullable=False),
    pa.field("source_id", pa.string(), nullable=False),
    pa.field("processing_partition", pa.string(), nullable=False),
    pa.field("catchment_kind", pa.string(), nullable=False),
    pa.field("geometry_normalization", pa.string(), nullable=False),
    pa.field("bbox", BOUNDS, nullable=False),
    pa.field("geometry", pa.binary(), nullable=False),
    pa.field("area_km2", pa.float32(), nullable=False),
    pa.field("up_area_km2", pa.float32()),
    pa.field("outlet_lon", pa.float64()),
    pa.field("outlet_lat", pa.float64()),
]

CATCHMENT_FIELDS: list[pa.Field] = [
    *CATCHMENT_RAW_FIELDS,
    pa.field("downstream_boundary_status", pa.string(), nullable=False),
    pa.field("upstream_boundary_cuts", pa.int32()),
]

CATCHMENT_DOWNSTREAM_STATUS_CODES: list[dict[str, str]] = [
    {"code": "internal", "meaning": "Paired river NextDownID resolves to a delivered river with a delivered catchment"},
    {"code": "source_zero_sentinel", "meaning": "Paired river carries native NextDownID=0 absent-link sentinel"},
    {
        "code": "cut_missing_downstream_catchment",
        "meaning": "Paired river NextDownID is nonzero, but the downstream river or its catchment was not delivered",
    },
    {
        "code": "omitted_degenerate_river",
        "meaning": "Paired source river was omitted for degenerate self-referencing geometry",
    },
    {"code": "river_not_in_source", "meaning": "Supplied river shapefile had no matching COMID"},
    {"code": "river_source_unavailable", "meaning": "River partition source was not available"},
    {
        "code": "unassociated_coastal_hillslope",
        "meaning": "Coastal hillslope polygon lacks a COMID and river-network association",
    },
]

RIVER_FIELDS: list[pa.Field] = [
    pa.field("id", pa.int64(), nullable=False),
    pa.field("level", pa.int16(), nullable=False),
    pa.field("parent_id", pa.int64()),
    pa.field("source_namespace", pa.string(), nullable=False),
    pa.field("source_id", pa.string(), nullable=False),
    pa.field("processing_partition", pa.string(), nullable=False),
    pa.field("bbox", BOUNDS, nullable=False),
    pa.field("geometry", pa.binary(), nullable=False),
    pa.field("local_catchment_status", pa.string(), nullable=False),
    pa.field("downstream_boundary_status", pa.string(), nullable=False),
    pa.field("upstream_boundary_cuts", pa.int8(), nullable=False),
    pa.field("component_catchment_id", pa.int64()),
    pa.field("upstream_node_id", pa.int64()),
    pa.field("downstream_node_id", pa.int64()),
    pa.field("geometry_orientation", pa.int8(), nullable=False),
]


def _sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _arrow_type(declared: str) -> pa.DataType:
    kind = declared.split(":", 1)[0].lower()
    if kind in ("int", "int32", "int64"):
        return pa.int64()
    if kind in ("float", "float32", "float64"):
        return pa.float64()
    if kind in ("str", "date", "time", "datetime"):
        return pa.string()
    raise ValueError(f"Unsupported shapefile property declaration: {declared}")


def _rename_column(name: str) -> str:
    return f"native_{name}" if name in STRUCTURAL_FIELDS else name


def _verify_sidecars(entry: dict[str, Any]) -> dict[str, dict[str, Any]]:
    sidecars = cast(dict[str, dict[str, Any]], entry["sidecars"])
    shp_path = Path(entry["shp"])
    hashes: dict[str, dict[str, Any]] = {}
    for side_name, expected in sidecars.items():
        path = Path(expected["path"])
        if path.parent.resolve() != shp_path.parent.resolve() or path.name != side_name:
            raise ValueError(f"Sidecar path mismatch: {path} vs {side_name}")
        digest = _sha256(path)
        size = path.stat().st_size
        if digest != expected["sha256"] or size != expected["size_bytes"]:
            raise ValueError(f"Source digest or size mismatch for {path}")
        hashes[str(path)] = {"sha256": digest, "size_bytes": size}
    return hashes


def _check_crs(src: fiona.Collection, entry: dict[str, Any]) -> str:
    has_prj = Path(entry["shp"]).with_suffix(".prj").exists()
    if has_prj:
        crs = src.crs
        if not crs:
            raise ValueError(f"Empty CRS despite .prj in {entry['shp']}")
        from pyproj import CRS

        parsed = CRS.from_user_input(crs)
        if parsed != CRS.from_epsg(4326):
            raise ValueError(f"Expected EPSG:4326 CRS in {entry['shp']}, got {parsed}")
        return "supplied_prj_epsg4326"
    evidence = entry.get("crs_evidence")
    if not isinstance(evidence, str) or not evidence.strip():
        raise ValueError(f"Missing .prj without documented crs_evidence in {entry['shp']}")
    return evidence


def _narrowed_bbox(raw_bounds: tuple[float, float, float, float]) -> dict[str, float]:
    bounds = np.asarray(raw_bounds, dtype=np.float64)
    narrowed = bounds.astype(np.float32)
    for idx in range(4):
        if (idx < 2 and narrowed[idx] > bounds[idx]) or (idx >= 2 and narrowed[idx] < bounds[idx]):
            narrowed[idx] = np.nextafter(narrowed[idx], np.float32(-np.inf if idx < 2 else np.inf))
    return dict(zip(("xmin", "ymin", "xmax", "ymax"), narrowed.tolist(), strict=True))


def _geodesic_area_km2(geom: Any, geod: Geod) -> float:
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    area = 0.0
    for poly in parts:
        ext = abs(geod.polygon_area_perimeter(*poly.exterior.xy)[0])
        holes = sum(abs(geod.polygon_area_perimeter(*ring.xy)[0]) for ring in poly.interiors)
        area += (ext - holes) / 1e6
    return area


def _validate_coordinates(geom: Any) -> None:
    coords = get_coordinates(geom)
    if (
        not len(coords)
        or not np.all(np.isfinite(coords))
        or np.any(coords[:, 0] < MIN_WGS84_LON)
        or np.any(coords[:, 0] > MAX_UNWRAPPED_ANTIMERIDIAN_LON)
        or np.any(np.abs(coords[:, 1]) > 90.0)
    ):
        raise ValueError("Coordinates outside finite WGS84 longitude/latitude bounds")


def _ring_edges_and_vertices(
    geom: Any,
) -> tuple[Counter[tuple[tuple[float, float], tuple[float, float]]], set[tuple[float, float]], int]:
    """Return exact IEEE 754 float64 undirected edge multiset, vertex set, and zero-length edge count."""
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    edges: Counter[tuple[tuple[float, float], tuple[float, float]]] = Counter()
    verts: set[tuple[float, float]] = set()
    zero_len = 0
    for poly in parts:
        for ring in [poly.exterior, *poly.interiors]:
            pts = [(float(x), float(y)) for x, y in ring.coords]
            verts.update(pts)
            for a, b in zip(pts[:-1], pts[1:], strict=True):
                if a == b:
                    zero_len += 1
                edges[(a, b) if a <= b else (b, a)] += 1
    return edges, verts, zero_len


def _part_and_hole_counts(geom: Any) -> tuple[int, int]:
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    return len(parts), sum(len(poly.interiors) for poly in parts)


def _shift_x(geom: Any, dx: float) -> Any:
    return shapely.transform(geom, lambda coords: coords + np.array([dx, 0.0], dtype=np.float64))


def _attempt_ring_normalization(
    raw_geom: Any, raw_reason: str, geod: Geod
) -> tuple[Any | None, dict[str, Any] | None, str | None]:
    """Attempt narrow lossless MERIT Ring Self-intersection split; fail closed on any invariant breach."""
    if not raw_reason.startswith("Ring Self-intersection["):
        return None, None, "ineligible_validity_reason"

    repaired = make_valid(raw_geom, method="structure", keep_collapsed=False)
    if repaired is None or repaired.is_empty:
        return None, None, "empty_after_repair"
    if repaired.geom_type not in {"Polygon", "MultiPolygon"}:
        return None, None, f"unsupported_repaired_type:{repaired.geom_type}"
    if repaired.has_z or repaired.has_m:
        return None, None, "repaired_has_zm"
    if not repaired.is_valid:
        return None, None, f"invalid_after_repair:{explain_validity(repaired)}"
    try:
        _validate_coordinates(repaired)
    except ValueError:
        return None, None, "repaired_coords_out_of_bounds"

    raw_edges, raw_verts, raw_zero = _ring_edges_and_vertices(raw_geom)
    rep_edges, rep_verts, rep_zero = _ring_edges_and_vertices(repaired)
    if raw_verts != rep_verts:
        return None, None, "vertex_set_mismatch"
    if raw_edges != rep_edges or raw_zero > 0 or rep_zero > 0:
        return None, None, "edge_multiset_mismatch"

    raw_area = _geodesic_area_km2(raw_geom, geod)
    rep_area = _geodesic_area_km2(repaired, geod)
    if not (np.isfinite(raw_area) and np.isfinite(rep_area)) or raw_area <= 0 or rep_area <= 0:
        return None, None, "non_positive_geodesic_area"

    area_diff = abs(rep_area - raw_area)
    area_rel_diff = area_diff / raw_area
    if area_diff > max(GEODESIC_AREA_ABS_TOL_KM2, GEODESIC_AREA_REL_TOL * raw_area):
        return None, None, "geodesic_area_mismatch"

    raw_parts, raw_holes = _part_and_hole_counts(raw_geom)
    rep_parts, rep_holes = _part_and_hole_counts(repaired)
    metrics = {
        "raw_geom_type": raw_geom.geom_type,
        "raw_part_count": raw_parts,
        "raw_hole_count": raw_holes,
        "repaired_geom_type": repaired.geom_type,
        "repaired_part_count": rep_parts,
        "repaired_hole_count": rep_holes,
        "vertex_count": len(rep_verts),
        "edge_count": sum(rep_edges.values()),
        "raw_geodesic_area_km2": raw_area,
        "repaired_geodesic_area_km2": rep_area,
        "geodesic_area_diff_km2": rep_area - raw_area,
        "geodesic_area_abs_diff_km2": area_diff,
        "geodesic_area_rel_diff": area_rel_diff,
    }
    return repaired, metrics, None


@dataclass(slots=True)
class _RepairCandidate:
    partition: str
    ns_domain: str
    ns: str
    raw_sid: str
    comid: int | None
    is_hillslope: bool
    raw_reason: str
    repaired_geom: Any
    original_wkb: bytes
    metrics: dict[str, Any]
    native_row_values: dict[str, Any]


@dataclass(frozen=True)
class _CatchmentConnectivity:
    upstream_ids: list[list[int]]
    downstream_status: list[str]
    upstream_cuts: list[int | None]
    counts: dict[str, Any]
    not_in_source_rivers: list[dict[str, Any]]
    not_in_source_rivers_truncated: bool


def _derive_catchment_connectivity(
    cat_records: list[tuple[bool, str, int | None]],
    *,
    db: sqlite3.Connection,
    unavailable_catchments: dict[str, str],
    unavailable_rivers: dict[str, str],
    omitted_river_comids: dict[int, str],
    riv_comid_to_hfx: dict[int, int],
    riv_partition_by_comid: dict[int, str],
    riv_next_down_by_comid: dict[int, int],
    riv_incoming_by_comid: dict[int, list[int]],
    riv_upstream_cuts_by_id: dict[int, int],
) -> _CatchmentConnectivity:
    """Derive immediate catchment edges from river topology and 1-to-1 COMID associations."""
    chan_rows = db.execute("SELECT comid, hfx_id, partition FROM catchments").fetchall()
    chan_comid_to_cat: dict[int, int] = {int(r[0]): int(r[1]) for r in chan_rows}
    chan_part_by_comid: dict[int, str] = {int(r[0]): sys.intern(str(r[2])) for r in chan_rows}

    for comid, riv_part in riv_partition_by_comid.items():
        if comid in chan_part_by_comid and chan_part_by_comid[comid] != riv_part:
            raise ValueError(
                f"Cross-partition COMID collision between river {riv_part}/{comid} "
                f"and catchment {chan_part_by_comid[comid]}/{comid}"
            )

    omitted_catchments_set = {
        (str(r[0]), int(r[1])) for r in db.execute("SELECT partition, comid FROM omitted_catchments").fetchall()
    }
    rivers_omitted_invalid_polygon = 0
    rivers_catchment_unavailable = 0
    not_in_source_rivers: list[dict[str, Any]] = []
    for comid, rid in riv_comid_to_hfx.items():
        if comid in chan_comid_to_cat:
            continue
        pname = riv_partition_by_comid[comid]
        if (pname, comid) in omitted_catchments_set:
            rivers_omitted_invalid_polygon += 1
        elif pname in unavailable_catchments:
            rivers_catchment_unavailable += 1
        else:
            not_in_source_rivers.append(
                {
                    "COMID": comid,
                    "partition": pname,
                    "hfx_river_id": rid,
                    "NextDownID": riv_next_down_by_comid[comid],
                }
            )

    cat_count = len(cat_records)
    cat_up_ids: list[list[int]] = [[] for _ in range(cat_count)]
    cat_down_status: list[str] = ["" for _ in range(cat_count)]
    cat_up_cuts: list[int | None] = [None for _ in range(cat_count)]

    paired_count = 0
    unpaired_omitted_riv = 0
    unpaired_riv_not_in_source = 0
    unpaired_riv_unavailable = 0
    hillslope_count = 0
    derived_edges = 0
    cuts_cat_count = 0
    total_up_cut_edges = 0
    downstream_cuts = 0

    for idx, (is_hillslope, pname, comid) in enumerate(cat_records):
        cid = idx + 1
        if is_hillslope:
            hillslope_count += 1
            cat_down_status[idx] = sys.intern("unassociated_coastal_hillslope")
            continue
        assert comid is not None
        if comid not in riv_comid_to_hfx:
            if comid in omitted_river_comids:
                cat_down_status[idx] = sys.intern("omitted_degenerate_river")
                unpaired_omitted_riv += 1
            elif pname in unavailable_rivers:
                cat_down_status[idx] = sys.intern("river_source_unavailable")
                unpaired_riv_unavailable += 1
            else:
                cat_down_status[idx] = sys.intern("river_not_in_source")
                unpaired_riv_not_in_source += 1
            continue

        paired_count += 1
        rid = riv_comid_to_hfx[comid]
        cuts = riv_upstream_cuts_by_id[rid]
        ups: list[int] = []
        for u_comid in riv_incoming_by_comid[comid]:
            u_cid = chan_comid_to_cat.get(u_comid)
            if u_cid is not None:
                if u_cid == cid:
                    raise ValueError(f"Self-loop derived for catchment {cid} (COMID {comid})")
                ups.append(u_cid)
            else:
                cuts += 1
        ups.sort()
        if len(set(ups)) != len(ups):
            raise ValueError(f"Duplicate derived upstream edge for catchment {cid} (COMID {comid})")
        cat_up_ids[idx] = ups
        cat_up_cuts[idx] = cuts
        derived_edges += len(ups)
        if cuts > 0:
            cuts_cat_count += 1
            total_up_cut_edges += cuts

        nd = riv_next_down_by_comid[comid]
        if nd == 0:
            cat_down_status[idx] = sys.intern("source_zero_sentinel")
        elif nd in chan_comid_to_cat and nd in riv_comid_to_hfx:
            cat_down_status[idx] = sys.intern("internal")
        else:
            cat_down_status[idx] = sys.intern("cut_missing_downstream_catchment")
            downstream_cuts += 1

    # Structural validation: Kahn topological sort on derived catchment graph + closure completeness.
    cat_in_deg = [len(ups) for ups in cat_up_ids]
    cat_down_adj: list[list[int]] = [[] for _ in range(cat_count)]
    for idx, ups in enumerate(cat_up_ids):
        for u_cid in ups:
            if not 1 <= u_cid <= cat_count:
                raise ValueError(f"Unresolved derived catchment edge {u_cid} -> {idx + 1}")
            cat_down_adj[u_cid - 1].append(idx)

    closure_clean = [cuts == 0 for cuts in cat_up_cuts]
    cat_queue = [i for i, deg in enumerate(cat_in_deg) if deg == 0]
    cat_visited = 0
    while cat_queue:
        cur = cat_queue.pop()
        cat_visited += 1
        for nxt in cat_down_adj[cur]:
            if not closure_clean[cur]:
                closure_clean[nxt] = False
            cat_in_deg[nxt] -= 1
            if cat_in_deg[nxt] == 0:
                cat_queue.append(nxt)
    if cat_visited != cat_count:
        raise ValueError("Derived catchment graph contains a cycle")

    # Cross-check against full river-network closure oracle.
    riv_in_deg = {comid: len(riv_incoming_by_comid[comid]) for comid in riv_comid_to_hfx}
    anc_has_gap = {
        comid: (comid not in chan_comid_to_cat) or (riv_upstream_cuts_by_id[riv_comid_to_hfx[comid]] > 0)
        for comid in riv_comid_to_hfx
    }
    anc_has_delivered = {comid: (comid in chan_comid_to_cat) for comid in riv_comid_to_hfx}
    anc_has_stranded_delivered = dict.fromkeys(riv_comid_to_hfx, False)
    riv_queue = [comid for comid, deg in riv_in_deg.items() if deg == 0]
    while riv_queue:
        cur = riv_queue.pop()
        nxt = riv_next_down_by_comid[cur]
        if nxt not in riv_comid_to_hfx:
            continue
        if anc_has_gap[cur]:
            anc_has_gap[nxt] = True
        if anc_has_delivered[cur]:
            anc_has_delivered[nxt] = True
        if (nxt not in chan_comid_to_cat) and anc_has_delivered[cur]:
            anc_has_stranded_delivered[nxt] = True
        if anc_has_stranded_delivered[cur]:
            anc_has_stranded_delivered[nxt] = True
        riv_in_deg[nxt] -= 1
        if riv_in_deg[nxt] == 0:
            riv_queue.append(nxt)

    stranded_count = 0
    for idx, (is_hillslope, _pname, comid) in enumerate(cat_records):
        if is_hillslope or comid is None or comid not in riv_comid_to_hfx:
            continue
        if closure_clean[idx] != (not anc_has_gap[comid]):
            raise ValueError(f"Completeness theorem mismatch on catchment {idx + 1} (COMID {comid})")
        if anc_has_stranded_delivered[comid]:
            stranded_count += 1

    complete_count = sum(closure_clean)
    unverified_count = sum(1 for cuts in cat_up_cuts if cuts is None)
    incomplete_paired_count = cat_count - complete_count - unverified_count
    leaf_only_count = incomplete_paired_count - stranded_count
    status_counts = dict(sorted(Counter(cat_down_status).items()))

    return _CatchmentConnectivity(
        upstream_ids=cat_up_ids,
        downstream_status=cat_down_status,
        upstream_cuts=cat_up_cuts,
        counts={
            "derived_flow_edges": derived_edges,
            "paired_channelized_catchments": paired_count,
            "unpaired_channelized_catchments_omitted_river": unpaired_omitted_riv,
            "unpaired_channelized_catchments_river_not_in_source": unpaired_riv_not_in_source,
            "unpaired_channelized_catchments_river_unavailable": unpaired_riv_unavailable,
            "unassociated_coastal_hillslopes": hillslope_count,
            "rivers_omitted_invalid_polygon": rivers_omitted_invalid_polygon,
            "rivers_catchment_not_in_source": len(not_in_source_rivers),
            "rivers_catchment_unavailable": rivers_catchment_unavailable,
            "downstream_boundary_status_counts": status_counts,
            "catchments_with_upstream_cuts": cuts_cat_count,
            "total_upstream_catchment_cut_edges": total_up_cut_edges,
            "catchments_with_downstream_cuts": downstream_cuts,
            "upstream_complete_catchments": complete_count,
            "upstream_incomplete_paired_catchments": incomplete_paired_count,
            "upstream_incomplete_leaf_only_catchments": leaf_only_count,
            "upstream_incomplete_stranded_delivered_catchments": stranded_count,
            "upstream_unverified_catchments": unverified_count,
        },
        not_in_source_rivers=not_in_source_rivers[:100],
        not_in_source_rivers_truncated=len(not_in_source_rivers) > 100,
    )


GRAPH_SCHEMA = pa.schema(
    [
        pa.field("id", pa.int64(), nullable=False),
        pa.field("level", pa.int16(), nullable=False),
        pa.field("upstream_ids", pa.list_(pa.int64()), nullable=False),
    ]
    + [pa.field("bbox_" + name, pa.float32(), nullable=False) for name in ("minx", "miny", "maxx", "maxy")]
)
ASSOC_SCHEMA = pa.schema(
    [pa.field("river_id", pa.int64(), nullable=False), pa.field("catchment_id", pa.int64(), nullable=False)]
)


@dataclass(slots=True)
class _TopologyResult:
    cat_conn: _CatchmentConnectivity
    cut_edges: list[dict[str, Any]]
    overflow_confluences: list[dict[str, Any]]
    downstream_status_by_id: dict[int, str]
    upstream_cuts_by_id: dict[int, int]
    comid_to_hfx: dict[int, int]
    incoming_by_comid: dict[int, list[int]]


def _analyze_topology(
    db: sqlite3.Connection,
    cat_records: list[tuple[bool, str, int | None]],
    *,
    unavailable_catchments: dict[str, str],
    unavailable_rivers: dict[str, str],
    omitted_river_comids: dict[int, str],
) -> _TopologyResult:
    """Validate river reciprocity and acyclicity, then derive catchment connectivity."""
    cut_edges: list[dict[str, Any]] = []
    downstream_status_by_id: dict[int, str] = {}
    upstream_cuts_by_id: dict[int, int] = {}

    all_rows = db.execute("SELECT comid, hfx_id, partition, next_down, up1, up2, up3, up4 FROM rivers").fetchall()
    comid_to_hfx: dict[int, int] = {r[0]: r[1] for r in all_rows}
    riv_partition_by_comid: dict[int, str] = {r[0]: sys.intern(str(r[2])) for r in all_rows}
    comid_to_ups: dict[int, tuple[int, ...]] = {r[0]: tuple(u for u in r[4:8] if u != 0) for r in all_rows}
    next_down_by_comid: dict[int, int] = {r[0]: r[3] for r in all_rows}

    for comid, _hfx_id, _pname, _next_down, up1, up2, up3, up4 in all_rows:
        for u in (up1, up2, up3, up4):
            if u != 0 and u in comid_to_hfx and next_down_by_comid[u] != comid:
                raise ValueError(
                    f"Supplied river edges disagree: {comid} lists upstream {u} with mismatched NextDownID"
                )

    in_degree: dict[int, int] = dict.fromkeys(comid_to_hfx, 0)
    successors: dict[int, list[int]] = {c: [] for c in comid_to_hfx}
    incoming_by_comid: dict[int, list[int]] = {c: [] for c in comid_to_hfx}
    overflow_confluences: list[dict[str, Any]] = []

    for comid, hfx_id, pname, next_down, up1, up2, up3, up4 in all_rows:
        if next_down == 0:
            downstream_status_by_id[hfx_id] = "source_zero_sentinel"
        elif next_down in comid_to_hfx:
            if comid not in comid_to_ups[next_down]:
                if len(comid_to_ups[next_down]) == 4:
                    overflow_confluences.append(
                        {
                            "upstream_comid": comid,
                            "downstream_comid": next_down,
                            "partition": pname,
                            "target_up1_to_up4": list(comid_to_ups[next_down]),
                            "reason": "dbf_up1_to_up4_four_slot_overflow",
                        }
                    )
                else:
                    raise ValueError(
                        f"Supplied river edges disagree: {comid} NextDownID={next_down} absent from up1..4"
                    )
            downstream_status_by_id[hfx_id] = "internal"
            incoming_by_comid[next_down].append(comid)
            successors[comid].append(next_down)
            in_degree[next_down] += 1
        else:
            downstream_status_by_id[hfx_id] = "cut_omitted_target"
            reason = (
                f"omitted_target:{omitted_river_comids[next_down]}"
                if next_down in omitted_river_comids
                else "target_not_in_delivered_partitions"
            )
            cut_edges.append(
                {
                    "direction": "downstream",
                    "origin_comid": comid,
                    "target_comid": next_down,
                    "origin_partition": pname,
                    "reason": reason,
                }
            )

        up_cuts = 0
        for u in (up1, up2, up3, up4):
            if u == 0 or u in comid_to_hfx:
                continue
            up_cuts += 1
            reason = (
                f"omitted_origin:{omitted_river_comids[u]}"
                if u in omitted_river_comids
                else "origin_not_in_delivered_partitions"
            )
            cut_edges.append(
                {
                    "direction": "upstream",
                    "origin_comid": u,
                    "target_comid": comid,
                    "target_partition": pname,
                    "reason": reason,
                }
            )
        upstream_cuts_by_id[hfx_id] = up_cuts

    pending = [c for c, deg in in_degree.items() if deg == 0]
    visited = 0
    while pending:
        cur = pending.pop()
        visited += 1
        for nxt in successors[cur]:
            in_degree[nxt] -= 1
            if in_degree[nxt] == 0:
                pending.append(nxt)
    if visited != len(comid_to_hfx):
        raise ValueError("Supplied river graph contains a cycle")
    del all_rows, successors, in_degree, comid_to_ups

    cat_conn = _derive_catchment_connectivity(
        cat_records,
        db=db,
        unavailable_catchments=unavailable_catchments,
        unavailable_rivers=unavailable_rivers,
        omitted_river_comids=omitted_river_comids,
        riv_comid_to_hfx=comid_to_hfx,
        riv_partition_by_comid=riv_partition_by_comid,
        riv_next_down_by_comid=next_down_by_comid,
        riv_incoming_by_comid=incoming_by_comid,
        riv_upstream_cuts_by_id=upstream_cuts_by_id,
    )
    return _TopologyResult(
        cat_conn=cat_conn,
        cut_edges=cut_edges,
        overflow_confluences=overflow_confluences,
        downstream_status_by_id=downstream_status_by_id,
        upstream_cuts_by_id=upstream_cuts_by_id,
        comid_to_hfx=comid_to_hfx,
        incoming_by_comid=incoming_by_comid,
    )


def _write_staging_outputs(
    output: Path,
    *,
    release: str,
    organization: Any,
    region: Any,
    original_inventory_bytes: bytes,
    omitted_features: list[dict[str, Any]],
    cut_edges: list[dict[str, Any]],
    crs_records: dict[str, str],
    cat_orig_schemas: dict[str, Any],
    riv_orig_schemas: dict[str, Any],
    cat_renames: dict[str, Any],
    riv_renames: dict[str, Any],
    unavailable_catchments: dict[str, str],
    unavailable_rivers: dict[str, str],
    cat_count: int,
    channelized_count: int,
    hillslope_count: int,
    riv_count: int,
    omitted_cat_count: int,
    omitted_riv_count: int,
    cat_conn_counts: dict[str, Any],
    not_in_source_rivers: list[dict[str, Any]],
    not_in_source_rivers_truncated: bool,
    overall_bbox: list[float],
    overflow_confluences: list[dict[str, Any]],
    aux_declarations: list[dict[str, Any]],
    partial_reason: str | None,
    geometry_normalization_meta: dict[str, Any],
    derived_from_verified_bundle_manifest_sha256: str | None = None,
    expected_raster_sha256: dict[str, tuple[Path, str]] | None = None,
) -> dict[str, Any]:
    (output / "original-inventory.json").write_bytes(original_inventory_bytes)
    (output / "omitted-features.json").write_text(json.dumps(omitted_features, indent=2) + "\n")
    (output / "cut-edges.json").write_text(json.dumps(cut_edges, indent=2) + "\n")
    source_metadata: dict[str, Any] = {
        "release": release,
        "organization": organization,
        "crs_evidence": crs_records,
        "original_schemas": {"catchments": cat_orig_schemas, "rivers": riv_orig_schemas},
        "column_renames": {"catchments": cat_renames, "rivers": riv_renames},
        "unavailable_partitions": {"catchments": unavailable_catchments, "rivers": unavailable_rivers},
        "counts": {
            "delivered_catchments": cat_count,
            "delivered_channelized_catchments": channelized_count,
            "delivered_coastal_hillslopes": hillslope_count,
            "delivered_rivers": riv_count,
            "omitted_catchments": omitted_cat_count,
            "omitted_rivers": omitted_riv_count,
            "cut_edges": len(cut_edges),
        },
        "geometry_normalization": geometry_normalization_meta,
        "catchment_connectivity": {
            "derivation_rule": (
                "Joined channelized catchments and delivered rivers on globally unique COMID within matching "
                "processing_partition; mapped each river's NextDownID!=0 to a catchment edge iff both upstream and "
                "downstream rivers have delivered associated catchments. Omitted polygons, unmatched river COMIDs, "
                "cut river boundaries, omitted degenerate self-loop rivers, and FID-only coastal hillslopes are "
                "never bridged across gaps or coerced into natural outlets."
            ),
            "counts": cat_conn_counts,
            "not_in_source_rivers": not_in_source_rivers,
            "not_in_source_rivers_truncated": not_in_source_rivers_truncated,
            "limitations": (
                "Partial connectivity: missing catchment polygons remain missing area, and raster refinement "
                "within the input polygon footprint cannot recover unrepresented area. Incomplete upstream "
                "closures (upstream_boundary_cuts > 0) and unverified units (upstream_boundary_cuts is null: "
                "coastal hillslopes and catchments whose degenerate self-loop river was omitted) must not be "
                "presented as ordinary complete watersheds without explicit consumer opt-in."
            ),
        },
        "unwrapped_antimeridian_max_lon": float(overall_bbox[2]) if overall_bbox[2] > 180.0 else None,
        "up1_to_up4_overflow_confluences": overflow_confluences,
        "transformations": (
            "Valid source catchments and rivers preserve original coordinates, part order, and native attributes "
            "without clipping, simplification, reversal, or dissolving. Narrow MERIT ring normalization admits "
            "invalid source catchment polygons whose raw OGC validity failure begins with 'Ring Self-intersection[' "
            "strictly via shapely.make_valid(method='structure', keep_collapsed=False) only when output is valid "
            "Polygon/MultiPolygon, exact float64 vertex set and undirected edge multiset match, WGS84 geodesic area "
            "changes by <= max(1e-9, 1e-12 * raw_area_km2) km2, and global sibling interior non-overlap holds across "
            "all partitions and coastal hillslopes; original WKBs and per-feature audit records are preserved in "
            "repaired-features.jsonl.gz, with geometry_normalization flagging 'original' vs "
            "'ring_self_intersection_split'. Ineligible or failed source polygons and degenerate zero-length "
            "self-referencing river records (NextDownID=up1=COMID, maxup=0) are omitted with explicit provenance. "
            "River graph upstream_ids derived from verified NextDownID links so >4-tributary confluences keep all incoming "
            "edges while native maxup and up1..up4 attributes remain intact. Catchment connectivity in graph.parquet "
            "derived from verified NextDownID links and 1-to-1 COMID associations without bridging omitted polygons or "
            "unmatched rivers; completeness is declared 'partial' via upstream_boundary_cuts and "
            "downstream_boundary_status. Eastern Chukotka features (pfaf_35 and hillslope_35) retain native unwrapped "
            "antimeridian longitudes (>180 deg E) without clipping or splitting, while MERIT Hydro rasters use wrapped "
            "~[-180, +180) longitudes (consumers querying rasters across >180 deg E must wrap longitudes by -360 deg and "
            "split reads across the antimeridian). Geodesic WGS84 ellipsoid area_km2 computed over native parts. Native "
            "COMIDs preserved; dataset-local HFX IDs assigned independently per domain."
        ),
    }
    if derived_from_verified_bundle_manifest_sha256 is not None:
        source_metadata["derived_from_verified_bundle_manifest_sha256"] = derived_from_verified_bundle_manifest_sha256
    (output / "source-metadata.json").write_text(json.dumps(source_metadata, indent=2) + "\n")
    _default_notices(output, partial_reason=partial_reason)

    manifest: dict[str, Any] = {
        "format_version": "0.4.0",
        "fabric_name": "merit",
        "fabric_version": release,
        "crs": "EPSG:4326",
        "has_up_area": False,
        "topology": "tree",
        "bbox": overall_bbox,
        "unit_count": cat_count,
        "created_at": datetime.now(UTC).isoformat(),
        "adapter_version": "merit-original-3",
        "region": region,
        "catchments": "catchments.arrow",
        "graph": "graph.parquet",
        "lookup": "lookup.arrow",
        "directory": "lookup-directory.json",
        "index_batch_rows": 65536,
        "catchment_batch_bytes": 8388608,
        "auxiliary": aux_declarations,
        "categorical_columns": [
            {
                "artifact": "catchments.arrow",
                "column": "catchment_kind",
                "physical_type": "string",
                "codes": [
                    {"code": "channelized", "meaning": "MERIT-Basins channelized unit catchment (25 km2 threshold)"},
                    {"code": "coastal_hillslope", "meaning": "MERIT-Basins coastal hillslope below 25 km2 threshold"},
                ],
                "null_meaning": "not permitted",
            },
            {
                "artifact": "catchments.arrow",
                "column": "geometry_normalization",
                "physical_type": "string",
                "codes": CATCHMENT_GEOMETRY_NORMALIZATION_CODES,
                "null_meaning": "not permitted",
            },
            {
                "artifact": "catchments.arrow",
                "column": "downstream_boundary_status",
                "physical_type": "string",
                "codes": CATCHMENT_DOWNSTREAM_STATUS_CODES,
                "null_meaning": "not permitted",
            },
            {
                "artifact": "rivers.arrow",
                "column": "local_catchment_status",
                "physical_type": "string",
                "codes": [
                    {"code": "supplied", "meaning": "Paired unit catchment delivered in local_catchment_ids"},
                    {
                        "code": "omitted_invalid_polygon",
                        "meaning": "Source polygon existed but was omitted due to invalid geometry",
                    },
                    {"code": "not_in_source", "meaning": "Supplied catchment shapefile had no matching COMID"},
                    {"code": "source_unavailable", "meaning": "Catchment partition source was not available"},
                ],
                "null_meaning": "not permitted",
            },
            {
                "artifact": "rivers.arrow",
                "column": "downstream_boundary_status",
                "physical_type": "string",
                "codes": [
                    {"code": "internal", "meaning": "Downstream NextDownID resolves to a delivered river"},
                    {"code": "source_zero_sentinel", "meaning": "Native NextDownID=0 absent-link sentinel"},
                    {"code": "cut_omitted_target", "meaning": "Downstream target river was omitted or unavailable"},
                ],
                "null_meaning": "not permitted",
            },
        ],
        "catchment_connectivity": [
            {
                "level": 0,
                "status": "supplied",
                "derivation": "river_local_associations",
                "completeness": "partial",
                "upstream_cuts_column": "upstream_boundary_cuts",
            }
        ],
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
    }
    artifacts = {p.relative_to(output).as_posix(): fingerprint(p) for p in sorted(output.rglob("*")) if p.is_file()}
    for rel_ref, (src_path, expected_digest) in (expected_raster_sha256 or {}).items():
        if rel_ref not in artifacts or artifacts[rel_ref]["sha256"] != expected_digest:
            raise ValueError(f"Raster source changed during normalization: {src_path}")
        dest = output / rel_ref
        try:
            same_file = os.path.samefile(src_path, dest)
        except OSError:
            same_file = False
        if not same_file and (not src_path.is_file() or _sha256(src_path) != expected_digest):
            raise ValueError(f"Raster source changed during normalization: {src_path}")
    manifest["artifacts"] = artifacts
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    receipt: dict[str, Any] = {
        "status": "passed",
        "delivered_catchments": cat_count,
        "delivered_channelized_catchments": channelized_count,
        "delivered_coastal_hillslopes": hillslope_count,
        "repaired_catchments": geometry_normalization_meta["global_counts"]["accepted"],
        "rejected_repair_catchments": geometry_normalization_meta["global_counts"]["rejected"],
        "delivered_rivers": riv_count,
        "omitted_catchments": omitted_cat_count,
        "omitted_rivers": omitted_riv_count,
        "cut_edges": len(cut_edges),
        "catchment_flow_edges": cat_conn_counts["derived_flow_edges"],
        "upstream_complete_catchments": cat_conn_counts["upstream_complete_catchments"],
        "manifest_sha256": _sha256(output / "manifest.json"),
    }
    if derived_from_verified_bundle_manifest_sha256 is not None:
        receipt["derived_from_verified_bundle_manifest_sha256"] = derived_from_verified_bundle_manifest_sha256
    return receipt


@contextmanager
def _workspace_db(scratch_parent: Path | None = None) -> Iterator[tuple[sqlite3.Connection, Path]]:
    with TemporaryDirectory(prefix="hfx-merit-work-", dir=scratch_parent) as scratch:
        scratch_dir = Path(scratch)
        conn = sqlite3.connect(scratch_dir / "merit.sqlite")
        try:
            conn.execute("PRAGMA cache_size = -32768")
            conn.execute("PRAGMA temp_store = FILE")
            conn.execute("PRAGMA synchronous = OFF")
            conn.execute(
                "CREATE TABLE catchments ("
                "partition TEXT NOT NULL, comid INTEGER NOT NULL, hfx_id INTEGER NOT NULL, "
                "PRIMARY KEY(partition, comid))"
            )
            conn.execute("CREATE INDEX idx_cat_comid ON catchments(comid)")
            conn.execute(
                "CREATE TABLE omitted_catchments ("
                "partition TEXT NOT NULL, comid INTEGER NOT NULL, reason TEXT NOT NULL, "
                "PRIMARY KEY(partition, comid))"
            )
            conn.execute(
                "CREATE TABLE rivers ("
                "comid INTEGER PRIMARY KEY, hfx_id INTEGER NOT NULL UNIQUE, partition TEXT NOT NULL, "
                "next_down INTEGER NOT NULL, maxup INTEGER NOT NULL, "
                "up1 INTEGER NOT NULL, up2 INTEGER NOT NULL, up3 INTEGER NOT NULL, up4 INTEGER NOT NULL)"
            )
            conn.execute("CREATE INDEX idx_riv_next ON rivers(next_down)")
            yield conn, scratch_dir
        finally:
            conn.close()


def _unify_native_schema(
    sources: list[tuple[str, dict[str, Any]]],
) -> tuple[list[pa.Field], dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Collect unified nullable native Arrow fields and track column renames/types per partition."""
    field_types: dict[str, pa.DataType] = {}
    field_order: list[str] = []
    renames: dict[str, dict[str, str]] = {}
    original_schemas: dict[str, dict[str, str]] = {}

    for key, entry in sources:
        with fiona.open(entry["shp"]) as src:
            props = dict(src.schema["properties"])
        original_schemas[key] = {k: str(v) for k, v in props.items()}
        rename_map: dict[str, str] = {}
        for orig_name, decl in props.items():
            mapped_name = _rename_column(orig_name)
            if mapped_name != orig_name:
                rename_map[orig_name] = mapped_name
            dtype = _arrow_type(str(decl))
            if mapped_name in field_types:
                if field_types[mapped_name] != dtype:
                    raise ValueError(f"Conflicting type for {mapped_name}: {field_types[mapped_name]} vs {dtype}")
            else:
                field_types[mapped_name] = dtype
                field_order.append(mapped_name)
        if rename_map:
            renames[key] = rename_map

    fields: list[pa.Field] = []
    for name in field_order:
        fields.append(pa.field(name, field_types[name], nullable=True))
    return fields, renames, original_schemas


def _coerce_value(value: object, dtype: pa.DataType) -> object:
    if value is None:
        return None
    if dtype == pa.int64():
        if isinstance(value, (int, str)):
            return int(value)
        raise ValueError(f"Expected integer-compatible value, got {type(value).__name__}")
    if dtype == pa.float64():
        if isinstance(value, (int, float, str)):
            return float(value)
        raise ValueError(f"Expected float-compatible value, got {type(value).__name__}")
    return str(value)


def _normalize_bundle_rel_path(rel_path: str, root: Path) -> str:
    rel = Path(rel_path)
    dest = root / rel
    if rel.is_absolute() or ".." in rel.parts or not dest.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Raster path escapes bundle: {rel_path}")
    return dest.resolve().relative_to(root.resolve()).as_posix()


def _check_and_stage_rasters(
    auxiliary_spec: dict[str, Any] | None, output: Path
) -> tuple[list[dict[str, Any]], dict[str, tuple[Path, str]]]:
    if auxiliary_spec is None:
        return [], {}
    entries = auxiliary_spec.get("auxiliary", [])
    source_paths = cast(dict[str, str], auxiliary_spec.get("files", {}))
    raster_hashes: dict[str, tuple[Path, str]] = {}
    for rel_path, expected_src in source_paths.items():
        src_path = Path(expected_src)
        digest = _sha256(src_path)
        norm_rel = _normalize_bundle_rel_path(rel_path, output)
        dest = output / norm_rel
        raster_hashes[norm_rel] = (src_path, digest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(src_path, dest)
        except OSError:
            shutil.copyfile(src_path, dest)
    return list(entries), raster_hashes


def _default_notices(output: Path, *, partial_reason: str | None) -> None:
    caveat = f" Coverage note: {partial_reason}" if partial_reason else ""
    (output / "NOTICE").write_text(
        "Private HFX 0.4.0 MERIT Hydro + MERIT Basins dataset.\n\n"
        "Sources:\n"
        "- MERIT Hydro v1.0.1 (Yamazaki et al., 2019), dual-licensed CC BY-NC 4.0 / ODbL 1.0; "
        "used under CC BY-NC 4.0 (NonCommercial). Original-format complete redistribution requires explicit permission.\n"
        "- MERIT-Basins v0.7/v1.0_bugfix1 and Coastal Hillslopes (Lin et al., 2019), distributed via reachhydro.org "
        "under CC BY-NC-SA 4.0 and UCSD non-commercial/research license terms.\n\n"
        "Transformations: native raster grids and direction codes preserved with ocean code 247 explicitly tagged as "
        "GeoTIFF nodata, river mouth code 0 and inland depression code 255 retained distinct under flow_dir_encoding='merit'. "
        "Vectors keep original coordinates, part order, and native attributes on valid features. Narrow MERIT ring "
        "normalization admits invalid source catchment polygons with raw OGC 'Ring Self-intersection[' only when "
        "shapely.make_valid(method='structure', keep_collapsed=False) produces a valid Polygon/MultiPolygon preserving "
        "exact float64 vertices, undirected edge multiset, WGS84 geodesic area within max(1e-9, 1e-12*raw_area_km2), "
        "and zero sibling interior overlap globally (audit records and original WKBs preserved in "
        "repaired-features.jsonl.gz; features marked via geometry_normalization). Ineligible or failed source polygons "
        "and degenerate zero-length self-referencing river records (NextDownID=up1=COMID, maxup=0) are omitted with "
        "explicit provenance. Upstream IDs in river-graph.parquet derive from verified incoming NextDownID links so "
        "confluences exceeding the shapefile's four up1..up4 slots retain all incoming edges while native maxup and up1..up4 "
        "columns stay intact. Catchment connectivity in graph.parquet derives from verified river NextDownID links and "
        "1-to-1 COMID associations without bridging omitted polygons or unmatched rivers; completeness is declared "
        "'partial' via upstream_boundary_cuts and downstream_boundary_status (coastal hillslopes and catchments whose "
        "degenerate self-loop river was omitted carry null upstream_boundary_cuts). Eastern Chukotka vectors (pfaf_35 and "
        "hillslope_35) preserve native unwrapped antimeridian longitudes (>180 deg E, up to ~191.003 deg E) without clipping "
        "or splitting, while MERIT Hydro rasters span wrapped longitudes (~-180.0004 to +179.9996 deg), requiring -360 deg "
        "longitude wrapping and split antimeridian raster reads. Derived geodesic area_km2 uses the WGS84 ellipsoid and is "
        f"distinct from native unitarea and uparea; outlets, nodes, and upstream polygon area remain null/absent.{caveat}\n"
    )
    (output / "CITATION.txt").write_text(
        "Yamazaki, D., Ikeshima, D., Sosa, J., Bates, P. D., Allen, G. H., & Pavelsky, T. M. (2019). "
        "MERIT Hydro: A high-resolution global hydrography map based on latest topography datasets. "
        "Water Resources Research, 55, 5053-5073. https://doi.org/10.1029/2019WR024873\n\n"
        "Lin, P., Pan, M., Beck, H. E., Yang, Y., Yamazaki, D., Frasson, R., David, C. H., Durand, M., "
        "Pavelsky, T. M., Allen, G. H., Gleason, C. J., & Wood, E. F. (2019). "
        "Global reconstruction of naturalized river flows at 2.94 million reaches. "
        "Water Resources Research, 55, 6499-6516. https://doi.org/10.1029/2019WR025287\n"
    )


def _open_verified_base_bundle(base_spatial_dir: Path) -> tuple[dict[str, Any], str]:
    """Verify and return (manifest_dict, manifest_sha256) for a base HFX 0.4.0 MERIT bundle."""
    base_manifest_path = base_spatial_dir / "manifest.json"
    base_manifest = json.loads(base_manifest_path.read_bytes())
    base_manifest_sha = _sha256(base_manifest_path)
    if (
        base_manifest.get("format_version") != "0.4.0"
        or base_manifest.get("fabric_name") != "merit"
        or base_manifest.get("representation") != "spatial-geoparquet-adjacency"
        or base_manifest.get("representation_version") != 1
    ):
        raise ValueError("Expected verified HFX 0.4.0 MERIT spatial-geoparquet-adjacency bundle")

    for rel, identity in base_manifest["artifacts"].items():
        path = base_spatial_dir / rel
        if (
            Path(rel).is_absolute()
            or not path.resolve().is_relative_to(base_spatial_dir.resolve())
            or not path.is_file()
            or fingerprint(path) != identity
        ):
            raise ValueError(f"Base bundle artifact differs or is missing: {rel}")
    return base_manifest, base_manifest_sha


def _validate_and_open_base_bundle(
    base_spatial_dir: Path,
    inventory: dict[str, Any],
    inherited_hillslopes: set[str],
) -> tuple[dict[str, Any], str, dict[str, Any]]:
    """Verify base spatial bundle integrity and sidecar equivalence for hybrid inheritance."""
    base_manifest, base_manifest_sha = _open_verified_base_bundle(base_spatial_dir)

    base_orig = json.loads((base_spatial_dir / "original-inventory.json").read_bytes())
    if base_orig.get("release") != inventory.get("release"):
        raise ValueError("Release differs from base original-inventory.json")

    base_parts = {str(p["partition"]): p for p in base_orig["partitions"]}
    inv_parts = {str(p["partition"]): p for p in inventory["partitions"]}
    if set(inv_parts) != set(base_parts):
        raise ValueError("pfaf_* partition set differs from base original-inventory.json")
    for pname, inv_p in inv_parts.items():
        base_p = base_parts[pname]
        for domain in ("catchments", "rivers"):
            if inv_p[domain].get("status") != "supplied" or base_p[domain].get("status") != "supplied":
                raise ValueError(f"Partition {pname}/{domain} must be supplied when using base_spatial_dir")
            inv_sides = {n: (int(m["size_bytes"]), str(m["sha256"])) for n, m in inv_p[domain]["sidecars"].items()}
            base_sides = {n: (int(m["size_bytes"]), str(m["sha256"])) for n, m in base_p[domain]["sidecars"].items()}
            if inv_sides != base_sides:
                raise ValueError(f"Supplied sidecars differ from base original-inventory.json for {pname}/{domain}")

    base_hills = {str(h["partition"]): h for h in base_orig.get("coastal_hillslopes", [])}
    inv_hills = {str(h["partition"]): h for h in inventory.get("coastal_hillslopes", [])}
    if set(inv_hills) != set(base_hills):
        raise ValueError("Coastal hillslope partition set differs from base original-inventory.json")
    for pname, inv_h in inv_hills.items():
        base_h = base_hills[pname]
        if pname in inherited_hillslopes:
            if base_h.get("status") != "supplied":
                raise ValueError(f"Inherited hillslope {pname} was not supplied in base original-inventory.json")
        elif inv_h.get("status") == "supplied" and base_h.get("status") == "supplied":
            inv_sides = {n: (int(m["size_bytes"]), str(m["sha256"])) for n, m in inv_h["sidecars"].items()}
            base_sides = {n: (int(m["size_bytes"]), str(m["sha256"])) for n, m in base_h["sidecars"].items()}
            if inv_sides != base_sides:
                raise ValueError(f"Supplied hillslope sidecars differ from base original-inventory.json for {pname}")

    base_meta = json.loads((base_spatial_dir / "source-metadata.json").read_bytes())
    return base_manifest, base_manifest_sha, base_meta


def _stage_base_rasters(
    base_spatial_dir: Path, base_manifest: dict[str, Any], output: Path
) -> tuple[list[dict[str, Any]], dict[str, tuple[Path, str]]]:
    aux_declarations = list(base_manifest.get("auxiliary", []))
    base_artifacts = base_manifest.get("artifacts", {})
    raster_before: dict[str, tuple[Path, str]] = {}
    for entry in aux_declarations:
        if entry.get("schema") != "hfx.aux.d8_raster.v2":
            continue
        for rel_ref in entry.get("artifacts", {}).values():
            if not isinstance(rel_ref, str):
                continue
            norm_rel = _normalize_bundle_rel_path(rel_ref, output)
            src_path = base_spatial_dir / norm_rel
            dest = output / norm_rel
            identity = base_artifacts.get(rel_ref) or base_artifacts.get(norm_rel)
            if not isinstance(identity, dict) or not isinstance(identity.get("sha256"), str):
                raise ValueError(f"Base manifest missing verified sha256 for raster: {rel_ref}")
            raster_before[norm_rel] = (src_path, str(identity["sha256"]))
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(src_path, dest)
            except OSError:
                shutil.copyfile(src_path, dest)
    return aux_declarations, raster_before


def _validate_candidate_overlaps(
    candidates: list[_RepairCandidate],
    valid_cat_path: Path,
) -> tuple[dict[int, str], dict[int, int], dict[str, Any]]:
    """Verify sibling interior non-overlap across all candidate repairs and all valid catchments globally."""
    overlap_stats: dict[str, int] = {
        "same_partition_pairs": 0,
        "cross_partition_pairs": 0,
        "coastal_hillslope_pairs": 0,
        "candidate_candidate_pairs": 0,
        "candidate_interior_overlaps_detected": 0,
        "accepted_interior_overlaps": 0,
    }
    reject_reasons: dict[int, str] = {}
    neighbors_checked: dict[int, int] = dict.fromkeys(range(len(candidates)), 0)
    if not candidates:
        return reject_reasons, neighbors_checked, overlap_stats

    cand_geoms = [c.repaired_geom for c in candidates]
    for g in cand_geoms:
        prepare(g)
    cand_bounds = [g.bounds for g in cand_geoms]
    has_antimeridian = any(b[0] < ANTIMERIDIAN_QUERY_MIN_X or b[2] > ANTIMERIDIAN_MAX_X for b in cand_bounds)
    tree = STRtree(cand_geoms)

    def _record_pair(p_a: str, p_b: str, involves_hillslope: bool) -> None:
        if involves_hillslope:
            overlap_stats["coastal_hillslope_pairs"] += 1
        elif p_a == p_b:
            overlap_stats["same_partition_pairs"] += 1
        else:
            overlap_stats["cross_partition_pairs"] += 1

    # 1. Candidate vs candidate (all partitions globally, including antimeridian wrap).
    for i, cand in enumerate(candidates):
        g_i = cand_geoms[i]
        b_i = cand_bounds[i]
        queries: list[tuple[float, Any]] = [(0.0, g_i)]
        if has_antimeridian:
            if b_i[2] > ANTIMERIDIAN_MAX_X:
                queries.append((-360.0, _shift_x(g_i, -360.0)))
            if b_i[0] < ANTIMERIDIAN_QUERY_MIN_X:
                queries.append((360.0, _shift_x(g_i, 360.0)))
        for shift_dx, q_geom in queries:
            for j in tree.query(q_geom, predicate="intersects"):
                j_idx = int(j)
                if j_idx <= i:
                    continue
                cand_j = candidates[j_idx]
                g_j = _shift_x(cand_geoms[j_idx], -shift_dx) if shift_dx else cand_geoms[j_idx]
                if shift_dx != 0.0 and not g_i.intersects(g_j):
                    continue
                overlap_stats["candidate_candidate_pairs"] += 1
                neighbors_checked[i] += 1
                neighbors_checked[j_idx] += 1
                _record_pair(cand.partition, cand_j.partition, cand.is_hillslope or cand_j.is_hillslope)
                if g_i.relate_pattern(g_j, "2********"):
                    overlap_stats["candidate_interior_overlaps_detected"] += 1
                    reject_reasons.setdefault(i, "candidate_interior_overlap")
                    reject_reasons.setdefault(j_idx, "candidate_interior_overlap")

    # 2. Candidate vs all valid catchments + coastal hillslopes (streamed in Arrow batches).
    with pa.memory_map(str(valid_cat_path), "r") as valid_stream:
        valid_reader = ipc.open_file(valid_stream)
        for batch_idx in range(valid_reader.num_record_batches):
            batch = valid_reader.get_batch(batch_idx)
            bounds_arr = cast(pa.StructArray, batch.column("bbox"))
            xmins = bounds_arr.field("xmin").to_numpy().astype(np.float64)
            ymins = bounds_arr.field("ymin").to_numpy().astype(np.float64)
            xmaxs = bounds_arr.field("xmax").to_numpy().astype(np.float64)
            ymaxs = bounds_arr.field("ymax").to_numpy().astype(np.float64)
            pnames = cast(list[str], batch.column("processing_partition").to_pylist())
            kinds = cast(list[str], batch.column("catchment_kind").to_pylist())
            wkbs = batch.column("geometry")

            hit_batches: list[tuple[float, np.ndarray, np.ndarray]] = []
            rows_0, cands_0 = tree.query(shapely_box(xmins, ymins, xmaxs, ymaxs))
            if len(rows_0):
                hit_batches.append((0.0, rows_0, cands_0))
            if has_antimeridian:
                west_rows = np.flatnonzero(xmins < ANTIMERIDIAN_QUERY_MIN_X)
                if len(west_rows):
                    wr, wc = tree.query(
                        shapely_box(
                            xmins[west_rows] + 360.0, ymins[west_rows], xmaxs[west_rows] + 360.0, ymaxs[west_rows]
                        )
                    )
                    if len(wr):
                        hit_batches.append((360.0, west_rows[wr], wc))
                east_rows = np.flatnonzero(xmaxs > ANTIMERIDIAN_MAX_X)
                if len(east_rows):
                    er, ec = tree.query(
                        shapely_box(
                            xmins[east_rows] - 360.0, ymins[east_rows], xmaxs[east_rows] - 360.0, ymaxs[east_rows]
                        )
                    )
                    if len(er):
                        hit_batches.append((-360.0, east_rows[er], ec))

            decoded_cache: dict[int, Any] = {}
            for shift_dx, row_arr, cand_arr in hit_batches:
                for row_pos_raw, cand_idx_raw in zip(row_arr.tolist(), cand_arr.tolist(), strict=True):
                    row_pos = int(row_pos_raw)
                    cand_idx = int(cand_idx_raw)
                    v_base = decoded_cache.get(row_pos)
                    if v_base is None:
                        v_base = from_wkb(wkbs[row_pos].as_py())
                        decoded_cache[row_pos] = v_base
                    v_geom = _shift_x(v_base, shift_dx) if shift_dx else v_base
                    g_cand = cand_geoms[cand_idx]
                    if not g_cand.intersects(v_geom):
                        continue
                    neighbors_checked[cand_idx] += 1
                    _record_pair(
                        candidates[cand_idx].partition,
                        pnames[row_pos],
                        candidates[cand_idx].is_hillslope or kinds[row_pos] == "coastal_hillslope",
                    )
                    if g_cand.relate_pattern(v_geom, "2********"):
                        overlap_stats["candidate_interior_overlaps_detected"] += 1
                        reject_reasons.setdefault(cand_idx, "sibling_interior_overlap")

    del tree
    return reject_reasons, neighbors_checked, overlap_stats


def _build_geometry_normalization_meta(
    *,
    enabled: bool,
    part_stats: dict[str, dict[str, Any]],
    overlap_stats: dict[str, int],
) -> dict[str, Any]:
    global_counts = {
        "total_source": sum(s["total_source"] for s in part_stats.values()),
        "already_valid": sum(s["already_valid"] for s in part_stats.values()),
        "attempted": sum(s["attempted"] for s in part_stats.values()),
        "accepted": sum(s["accepted"] for s in part_stats.values()),
        "rejected": sum(s["rejected"] for s in part_stats.values()),
        "omitted_other": sum(s["omitted_other"] for s in part_stats.values()),
        "max_abs_area_diff_km2": max((s["max_abs_area_diff_km2"] for s in part_stats.values()), default=0.0),
        "max_rel_area_diff": max((s["max_rel_area_diff"] for s in part_stats.values()), default=0.0),
    }
    return {
        "artifact": REPAIRED_FEATURES_ARTIFACT,
        "enabled": enabled,
        "rule": GEOMETRY_NORMALIZATION_RULE if enabled else "disabled",
        "algorithm": GEOMETRY_NORMALIZATION_ALGORITHM,
        "versions": {
            "python": platform.python_version(),
            "shapely": shapely.__version__,
            "geos": ".".join(str(x) for x in shapely.geos_version),
            "pyproj": sys.modules["pyproj"].__version__,
            "fiona": fiona.__version__,
        },
        "tolerances": {
            "geodesic_area_abs_tol_km2": GEODESIC_AREA_ABS_TOL_KM2,
            "geodesic_area_rel_tol": GEODESIC_AREA_REL_TOL,
            "interior_overlap_pattern": "2********",
        },
        "global_counts": global_counts,
        "overlap_validation": {
            **overlap_stats,
            "total_pairs_checked": (
                overlap_stats["same_partition_pairs"]
                + overlap_stats["cross_partition_pairs"]
                + overlap_stats["coastal_hillslope_pairs"]
            ),
        },
        "by_partition": {k: part_stats[k] for k in sorted(part_stats)},
    }


def produce(
    inventory_path: Path,
    output: Path,
    *,
    auxiliary_spec_path: Path | None = None,
    base_spatial_dir: Path | None = None,
    normalize_rings: bool = True,
    batch_rows: int = 4096,
    partial_reason: str | None = None,
) -> dict[str, Any]:
    """Normalize MERIT Basins vectors (+ optional rasters) into HFX 0.4.0 staging."""
    inventory = json.loads(inventory_path.read_bytes())
    release = str(inventory["release"])
    partitions = sorted(inventory["partitions"], key=lambda item: item["partition"])
    hillslopes = sorted(inventory.get("coastal_hillslopes", []), key=lambda item: item["partition"])

    if not partitions:
        raise ValueError("Inventory must include at least one partition")

    before_hashes: dict[str, dict[str, Any]] = {}
    supplied_cat_sources: list[tuple[str, dict[str, Any]]] = []
    supplied_riv_sources: list[tuple[str, dict[str, Any]]] = []
    unavailable_catchments: dict[str, str] = {}
    unavailable_rivers: dict[str, str] = {}
    inherited_hillslopes: set[str] = set()

    for part in partitions:
        pname = str(part["partition"])
        c_entry = part["catchments"]
        r_entry = part["rivers"]
        if c_entry.get("status") == "supplied":
            before_hashes.update(_verify_sidecars(c_entry))
            supplied_cat_sources.append((pname, c_entry))
        else:
            unavailable_catchments[pname] = str(c_entry.get("reason", "unavailable"))
        if r_entry.get("status") == "supplied":
            before_hashes.update(_verify_sidecars(r_entry))
            supplied_riv_sources.append((pname, r_entry))
        else:
            unavailable_rivers[pname] = str(r_entry.get("reason", "unavailable"))

    for h_entry in hillslopes:
        pname = str(h_entry["partition"])
        status = h_entry.get("status")
        if status == "supplied":
            before_hashes.update(_verify_sidecars(h_entry))
            supplied_cat_sources.append((pname, h_entry))
        elif status == "inherited_from_verified_bundle":
            inherited_hillslopes.add(pname)
        else:
            raise ValueError(f"Unsupported coastal_hillslopes status {status!r} for {pname}")

    if inherited_hillslopes and base_spatial_dir is None:
        raise ValueError("Inherited coastal_hillslopes require base_spatial_dir")
    if (not supplied_cat_sources and not inherited_hillslopes) or not supplied_riv_sources:
        raise ValueError("Both supplied catchments and rivers are required")

    base_manifest: dict[str, Any] | None = None
    base_manifest_sha: str | None = None
    base_meta: dict[str, Any] | None = None
    if base_spatial_dir is not None:
        base_manifest, base_manifest_sha, base_meta = _validate_and_open_base_bundle(
            base_spatial_dir, inventory, inherited_hillslopes
        )

    auxiliary_spec = json.loads(auxiliary_spec_path.read_bytes()) if auxiliary_spec_path else None

    output.mkdir(parents=True, exist_ok=False)
    if auxiliary_spec is not None:
        aux_declarations, raster_before = _check_and_stage_rasters(auxiliary_spec, output)
    elif base_spatial_dir is not None and base_manifest is not None:
        aux_declarations, raster_before = _stage_base_rasters(base_spatial_dir, base_manifest, output)
    else:
        aux_declarations, raster_before = [], {}

    cat_native_fields, cat_renames, cat_orig_schemas = _unify_native_schema(supplied_cat_sources)
    if inherited_hillslopes and base_spatial_dir is not None and base_manifest is not None and base_meta is not None:
        base_cat_pf = pq.ParquetFile(base_spatial_dir / str(base_manifest["families"]["catchments"]["data"]))
        cat_field_types = {f.name: f.type for f in cat_native_fields}
        base_structural = {f.name for f in CATCHMENT_FIELDS}
        for field in base_cat_pf.schema_arrow.with_metadata({}):
            if field.name in base_structural:
                continue
            if field.name in cat_field_types:
                if cat_field_types[field.name] != field.type:
                    raise ValueError(f"Conflicting inherited type for {field.name}")
            else:
                cat_native_fields.append(pa.field(field.name, field.type, nullable=True))
                cat_field_types[field.name] = field.type
        for pname in sorted(inherited_hillslopes):
            if pname in base_meta.get("original_schemas", {}).get("catchments", {}):
                cat_orig_schemas[pname] = base_meta["original_schemas"]["catchments"][pname]
            if pname in base_meta.get("column_renames", {}).get("catchments", {}):
                cat_renames[pname] = base_meta["column_renames"]["catchments"][pname]

    riv_native_fields, riv_renames, riv_orig_schemas = _unify_native_schema(supplied_riv_sources)

    raw_cat_schema = pa.schema(CATCHMENT_RAW_FIELDS + cat_native_fields)
    cat_schema = pa.schema(CATCHMENT_FIELDS + cat_native_fields)
    riv_schema = pa.schema(RIVER_FIELDS + riv_native_fields)

    geod = Geod(ellps="WGS84")
    omitted_features: list[dict[str, Any]] = []
    crs_records: dict[str, str] = {}
    overall_bbox = [float("inf"), float("inf"), -float("inf"), -float("inf")]

    cat_count = 0
    riv_count = 0
    hillslope_count = 0
    channelized_count = 0
    omitted_cat_count = 0
    omitted_riv_count = 0

    part_stats: dict[str, dict[str, Any]] = {}

    def _pstat(pname: str) -> dict[str, Any]:
        return part_stats.setdefault(
            pname,
            {
                "total_source": 0,
                "already_valid": 0,
                "attempted": 0,
                "accepted": 0,
                "rejected": 0,
                "omitted_other": 0,
                "max_abs_area_diff_km2": 0.0,
                "max_rel_area_diff": 0.0,
            },
        )

    with _workspace_db(scratch_parent=output.parent) as (db, scratch_dir):
        # Pass 1a: Stream valid catchments/hillslopes to disk and collect eligible ring-normalization candidates.
        valid_cat_path = scratch_dir / "catchments_valid.arrow"
        repaired_cat_path = scratch_dir / "catchments_repaired.arrow"
        cat_records: list[tuple[bool, str, int | None]] = []
        seen_source_ids: set[tuple[str, str]] = set()
        seen_global_channelized_comids: set[int] = set()
        candidates: list[_RepairCandidate] = []

        with ipc.new_file(
            str(valid_cat_path), raw_cat_schema, options=ipc.IpcWriteOptions(compression="zstd")
        ) as valid_cat_writer:
            for pname, entry in supplied_cat_sources:
                is_hillslope = entry in hillslopes
                kind = "coastal_hillslope" if is_hillslope else "channelized"
                ns_domain = "coastal_hillslopes" if is_hillslope else "catchments"
                ns = f"{release}/{pname}/{ns_domain}"
                rename_map = cat_renames.get(pname, {})
                interned_pname = sys.intern(pname)
                pstat = _pstat(pname)
                with fiona.open(entry["shp"]) as src:
                    crs_records[f"{pname}/{ns_domain}"] = _check_crs(src, entry)
                    rows: list[dict[str, Any]] = []
                    for feat in src:
                        pstat["total_source"] += 1
                        props = dict(feat["properties"])
                        comid = int(props["COMID"]) if "COMID" in props and props["COMID"] is not None else None
                        raw_sid = str(comid) if comid is not None else str(props.get("FID", feat.id))
                        if (ns, raw_sid) in seen_source_ids:
                            raise ValueError(f"Duplicate source identity {(ns, raw_sid)}")
                        seen_source_ids.add((ns, raw_sid))
                        if not is_hillslope:
                            if comid is None:
                                raise ValueError(f"Missing COMID in channelized catchment {pname}")
                            if comid in seen_global_channelized_comids:
                                raise ValueError(f"Duplicate global channelized catchment COMID {comid}")
                            seen_global_channelized_comids.add(comid)

                        geom_dict = feat["geometry"]
                        area_km2 = 0.0
                        if geom_dict is None:
                            reason = "null geometry"
                            geom = None
                        else:
                            geom = shapely_shape(geom_dict)
                            if geom.is_empty:
                                reason = "empty geometry"
                            elif geom.geom_type not in {"Polygon", "MultiPolygon"}:
                                reason = f"unsupported geometry type {geom.geom_type}"
                            elif geom.has_z or geom.has_m:
                                reason = "unsupported Z/M coordinates"
                            else:
                                try:
                                    _validate_coordinates(geom)
                                except ValueError as exc:
                                    reason = str(exc)
                                else:
                                    if not geom.is_valid:
                                        reason = str(explain_validity(geom))
                                    else:
                                        area_km2 = _geodesic_area_km2(geom, geod)
                                        if not np.isfinite(area_km2) or area_km2 <= 0:
                                            reason = f"non-positive geodesic area {area_km2}"
                                        else:
                                            reason = ""

                        if reason:
                            if normalize_rings and geom is not None and reason.startswith("Ring Self-intersection["):
                                pstat["attempted"] += 1
                                rep_geom, metrics, fail_code = _attempt_ring_normalization(geom, reason, geod)
                                if rep_geom is not None and metrics is not None:
                                    native_vals: dict[str, Any] = {}
                                    for orig_name, val in props.items():
                                        mapped = rename_map.get(orig_name, orig_name)
                                        dtype = raw_cat_schema.field(mapped).type
                                        native_vals[mapped] = _coerce_value(val, dtype)
                                    candidates.append(
                                        _RepairCandidate(
                                            partition=pname,
                                            ns_domain=ns_domain,
                                            ns=ns,
                                            raw_sid=raw_sid,
                                            comid=comid,
                                            is_hillslope=is_hillslope,
                                            raw_reason=reason,
                                            repaired_geom=rep_geom,
                                            original_wkb=geom.wkb,
                                            metrics=metrics,
                                            native_row_values=native_vals,
                                        )
                                    )
                                    continue
                                pstat["rejected"] += 1
                                omitted_cat_count += 1
                                omitted_features.append(
                                    {
                                        "domain": ns_domain,
                                        "partition": pname,
                                        "source_namespace": ns,
                                        "source_id": raw_sid,
                                        "COMID": comid,
                                        "reason": reason,
                                        "repair_rejection": fail_code or "repair_check_failed",
                                    }
                                )
                                if not is_hillslope and comid is not None:
                                    db.execute(
                                        "INSERT INTO omitted_catchments VALUES (?, ?, ?)",
                                        (pname, comid, reason),
                                    )
                                continue

                            pstat["omitted_other"] += 1
                            omitted_cat_count += 1
                            omitted_features.append(
                                {
                                    "domain": ns_domain,
                                    "partition": pname,
                                    "source_namespace": ns,
                                    "source_id": raw_sid,
                                    "COMID": comid,
                                    "reason": reason,
                                }
                            )
                            if not is_hillslope and comid is not None:
                                db.execute(
                                    "INSERT INTO omitted_catchments VALUES (?, ?, ?)",
                                    (pname, comid, reason),
                                )
                            continue

                        assert geom is not None
                        pstat["already_valid"] += 1
                        wkb = geom.wkb
                        bbox = _narrowed_bbox(geom.bounds)
                        for idx, key in enumerate(("xmin", "ymin", "xmax", "ymax")):
                            overall_bbox[idx] = (
                                min(overall_bbox[idx], bbox[key]) if idx < 2 else max(overall_bbox[idx], bbox[key])
                            )
                        cat_count += 1
                        if is_hillslope:
                            hillslope_count += 1
                        else:
                            channelized_count += 1
                            assert comid is not None
                            db.execute("INSERT INTO catchments VALUES (?, ?, ?)", (pname, comid, cat_count))
                        cat_records.append((is_hillslope, interned_pname, None if is_hillslope else comid))

                        out_row: dict[str, Any] = {
                            "id": cat_count,
                            "level": 0,
                            "parent_id": None,
                            "source_namespace": ns,
                            "source_id": raw_sid,
                            "processing_partition": pname,
                            "catchment_kind": kind,
                            "geometry_normalization": "original",
                            "bbox": bbox,
                            "geometry": wkb,
                            "area_km2": float(np.float32(area_km2)),
                            "up_area_km2": None,
                            "outlet_lon": None,
                            "outlet_lat": None,
                        }
                        for field in cat_native_fields:
                            out_row[field.name] = None
                        for orig_name, val in props.items():
                            mapped = rename_map.get(orig_name, orig_name)
                            dtype = raw_cat_schema.field(mapped).type
                            out_row[mapped] = _coerce_value(val, dtype)
                        rows.append(out_row)
                        if len(rows) >= batch_rows:
                            valid_cat_writer.write_batch(pa.RecordBatch.from_pylist(rows, schema=raw_cat_schema))
                            rows = []
                    if rows:
                        valid_cat_writer.write_batch(pa.RecordBatch.from_pylist(rows, schema=raw_cat_schema))
                db.commit()

            if inherited_hillslopes:
                assert base_spatial_dir is not None and base_manifest is not None and base_meta is not None
                for pname in sorted(inherited_hillslopes):
                    crs_key = f"{pname}/coastal_hillslopes"
                    if crs_key in base_meta.get("crs_evidence", {}):
                        crs_records[crs_key] = base_meta["crs_evidence"][crs_key]
                    _pstat(pname)["source"] = "inherited_from_verified_bundle"
                base_cat_pf = pq.ParquetFile(base_spatial_dir / str(base_manifest["families"]["catchments"]["data"]))
                inherited_native_names = [f.name for f in cat_native_fields if f.name in base_cat_pf.schema_arrow.names]
                for batch in base_cat_pf.iter_batches(batch_size=batch_rows):
                    kinds = batch.column("catchment_kind").to_pylist()
                    pnames = batch.column("processing_partition").to_pylist()
                    indices = [
                        i
                        for i, (k_val, p_val) in enumerate(zip(kinds, pnames, strict=True))
                        if k_val == "coastal_hillslope" and p_val in inherited_hillslopes
                    ]
                    if not indices:
                        continue
                    matched = batch.take(indices)
                    rows = []
                    for row in matched.to_pylist():
                        if row.get("geometry_normalization", "original") != "original":
                            raise ValueError(
                                f"Unexpected geometry_normalization on inherited hillslope: "
                                f"{(row['source_namespace'], row['source_id'])}"
                            )
                        pname = str(row["processing_partition"])
                        ns = str(row["source_namespace"])
                        raw_sid = str(row["source_id"])
                        if (ns, raw_sid) in seen_source_ids:
                            raise ValueError(f"Duplicate source identity {(ns, raw_sid)}")
                        seen_source_ids.add((ns, raw_sid))
                        pstat = _pstat(pname)
                        pstat["total_source"] += 1
                        pstat["already_valid"] += 1
                        bbox = row["bbox"]
                        for idx, key in enumerate(("xmin", "ymin", "xmax", "ymax")):
                            overall_bbox[idx] = (
                                min(overall_bbox[idx], bbox[key]) if idx < 2 else max(overall_bbox[idx], bbox[key])
                            )
                        cat_count += 1
                        hillslope_count += 1
                        cat_records.append((True, sys.intern(pname), None))
                        out_row = {
                            "id": cat_count,
                            "level": 0,
                            "parent_id": None,
                            "source_namespace": ns,
                            "source_id": raw_sid,
                            "processing_partition": pname,
                            "catchment_kind": "coastal_hillslope",
                            "geometry_normalization": "original",
                            "bbox": bbox,
                            "geometry": row["geometry"],
                            "area_km2": float(row["area_km2"]),
                            "up_area_km2": None,
                            "outlet_lon": None,
                            "outlet_lat": None,
                        }
                        for field in cat_native_fields:
                            out_row[field.name] = row[field.name] if field.name in inherited_native_names else None
                        rows.append(out_row)
                    if rows:
                        valid_cat_writer.write_batch(pa.RecordBatch.from_pylist(rows, schema=raw_cat_schema))
                base_omitted_path = base_spatial_dir / "omitted-features.json"
                if base_omitted_path.is_file():
                    for item in json.loads(base_omitted_path.read_bytes()):
                        pname = str(item.get("partition", ""))
                        if item.get("domain") != "rivers" and pname in inherited_hillslopes:
                            omitted_features.append(item)
                            omitted_cat_count += 1
                            pstat = _pstat(pname)
                            pstat["total_source"] += 1
                            pstat["omitted_other"] += 1
        del seen_source_ids, seen_global_channelized_comids

        # Pass 1b: Global sibling interior non-overlap validation across all partitions and hillslopes.
        reject_reasons, neighbors_checked, overlap_stats = _validate_candidate_overlaps(candidates, valid_cat_path)
        has_repaired_rows = False
        repaired_gz_path = output / REPAIRED_FEATURES_ARTIFACT
        with (
            gzip.GzipFile(repaired_gz_path, mode="wb", mtime=0) as gz_stream,
            ipc.new_file(
                str(repaired_cat_path), raw_cat_schema, options=ipc.IpcWriteOptions(compression="zstd")
            ) as rep_cat_writer,
        ):
            repaired_rows: list[dict[str, Any]] = []
            for idx, cand in enumerate(candidates):
                pstat = _pstat(cand.partition)
                rej = reject_reasons.get(idx)
                if rej is not None:
                    pstat["rejected"] += 1
                    omitted_cat_count += 1
                    omitted_features.append(
                        {
                            "domain": cand.ns_domain,
                            "partition": cand.partition,
                            "source_namespace": cand.ns,
                            "source_id": cand.raw_sid,
                            "COMID": cand.comid,
                            "reason": cand.raw_reason,
                            "repair_rejection": rej,
                        }
                    )
                    if not cand.is_hillslope and cand.comid is not None:
                        db.execute(
                            "INSERT INTO omitted_catchments VALUES (?, ?, ?)",
                            (cand.partition, cand.comid, cand.raw_reason),
                        )
                    continue

                rep_geom = cand.repaired_geom
                bbox = _narrowed_bbox(rep_geom.bounds)
                for b_idx, key in enumerate(("xmin", "ymin", "xmax", "ymax")):
                    overall_bbox[b_idx] = (
                        min(overall_bbox[b_idx], bbox[key]) if b_idx < 2 else max(overall_bbox[b_idx], bbox[key])
                    )
                cat_count += 1
                if cand.is_hillslope:
                    hillslope_count += 1
                else:
                    channelized_count += 1
                    assert cand.comid is not None
                    db.execute("INSERT INTO catchments VALUES (?, ?, ?)", (cand.partition, cand.comid, cat_count))
                cat_records.append(
                    (cand.is_hillslope, sys.intern(cand.partition), None if cand.is_hillslope else cand.comid)
                )

                pstat["accepted"] += 1
                pstat["max_abs_area_diff_km2"] = max(
                    pstat["max_abs_area_diff_km2"], cand.metrics["geodesic_area_abs_diff_km2"]
                )
                pstat["max_rel_area_diff"] = max(pstat["max_rel_area_diff"], cand.metrics["geodesic_area_rel_diff"])

                audit_entry = {
                    "hfx_id": cat_count,
                    "domain": cand.ns_domain,
                    "partition": cand.partition,
                    "source_namespace": cand.ns,
                    "source_id": cand.raw_sid,
                    "COMID": cand.comid,
                    "rule": GEOMETRY_NORMALIZATION_RULE,
                    "algorithm": GEOMETRY_NORMALIZATION_ALGORITHM,
                    "raw_validity_reason": cand.raw_reason,
                    **cand.metrics,
                    "overlap_neighbors_checked": neighbors_checked[idx],
                    "checks": {
                        "valid_polygonal": True,
                        "vertices_exact": True,
                        "edges_exact": True,
                        "geodesic_area_within_tolerance": True,
                        "zero_sibling_interior_overlap": True,
                    },
                    "original_wkb_base64": base64.b64encode(cand.original_wkb).decode("ascii"),
                }
                gz_stream.write((json.dumps(audit_entry, sort_keys=True) + "\n").encode("utf-8"))

                out_row = {
                    "id": cat_count,
                    "level": 0,
                    "parent_id": None,
                    "source_namespace": cand.ns,
                    "source_id": cand.raw_sid,
                    "processing_partition": cand.partition,
                    "catchment_kind": "coastal_hillslope" if cand.is_hillslope else "channelized",
                    "geometry_normalization": "ring_self_intersection_split",
                    "bbox": bbox,
                    "geometry": rep_geom.wkb,
                    "area_km2": float(np.float32(cand.metrics["repaired_geodesic_area_km2"])),
                    "up_area_km2": None,
                    "outlet_lon": None,
                    "outlet_lat": None,
                }
                for field in cat_native_fields:
                    out_row[field.name] = cand.native_row_values.get(field.name)
                repaired_rows.append(out_row)
                has_repaired_rows = True
                if len(repaired_rows) >= batch_rows:
                    rep_cat_writer.write_batch(pa.RecordBatch.from_pylist(repaired_rows, schema=raw_cat_schema))
                    repaired_rows = []

            if repaired_rows:
                rep_cat_writer.write_batch(pa.RecordBatch.from_pylist(repaired_rows, schema=raw_cat_schema))
                repaired_rows = []
        if not has_repaired_rows:
            repaired_cat_path.unlink(missing_ok=True)
        db.commit()
        del candidates, reject_reasons, neighbors_checked

        # Pass 2a: Validate river geometries and load connectivity into SQLite.
        omitted_river_comids: dict[int, str] = {}
        next_riv_id = 0
        for pname, entry in supplied_riv_sources:
            ns = f"{release}/{pname}/rivers"
            with fiona.open(entry["shp"]) as src:
                crs_records[f"{pname}/rivers"] = _check_crs(src, entry)
                for feat in src:
                    props = dict(feat["properties"])
                    comid = int(props["COMID"])
                    if comid <= 0:
                        raise ValueError(f"River COMID must be positive, got {comid}")
                    next_down = int(props["NextDownID"])
                    maxup = int(props["maxup"])
                    ups = [int(props[f"up{k}"]) for k in range(1, 5)]
                    geom_dict = feat["geometry"]
                    if geom_dict is None:
                        reason = "null geometry"
                    else:
                        geom = shapely_shape(geom_dict)
                        if geom.is_empty:
                            reason = "empty geometry"
                        elif geom.geom_type not in {"LineString", "MultiLineString"}:
                            reason = f"unsupported river geometry type {geom.geom_type}"
                        elif geom.has_z or geom.has_m:
                            reason = "unsupported Z/M coordinates"
                        else:
                            try:
                                _validate_coordinates(geom)
                            except ValueError as exc:
                                reason = str(exc)
                            else:
                                if not geom.is_valid:
                                    reason = f"invalid river geometry: {explain_validity(geom)}"
                                    if next_down == comid and maxup == 0 and ups == [comid, 0, 0, 0]:
                                        reason += "; degenerate self-reference NextDownID=up1=COMID with maxup=0"
                                else:
                                    reason = ""
                    if reason:
                        cat_row = db.execute(
                            "SELECT hfx_id FROM catchments WHERE partition = ? AND comid = ?",
                            (pname, comid),
                        ).fetchone()
                        if cat_row is not None:
                            paired_status = "delivered"
                            paired_hfx_id: int | None = int(cat_row[0])
                        elif db.execute(
                            "SELECT 1 FROM omitted_catchments WHERE partition = ? AND comid = ?",
                            (pname, comid),
                        ).fetchone():
                            paired_status = "omitted_invalid_polygon"
                            paired_hfx_id = None
                        elif pname in unavailable_catchments:
                            paired_status = "source_unavailable"
                            paired_hfx_id = None
                        else:
                            paired_status = "not_in_source"
                            paired_hfx_id = None
                        omitted_riv_count += 1
                        omitted_river_comids[comid] = reason
                        omitted_features.append(
                            {
                                "domain": "rivers",
                                "partition": pname,
                                "source_namespace": ns,
                                "source_id": str(comid),
                                "COMID": comid,
                                "reason": reason,
                                "paired_catchment_status": paired_status,
                                "paired_catchment_hfx_id": paired_hfx_id,
                                "native_connectivity": {
                                    "NextDownID": next_down,
                                    "maxup": maxup,
                                    "up1": ups[0],
                                    "up2": ups[1],
                                    "up3": ups[2],
                                    "up4": ups[3],
                                },
                            }
                        )
                        continue
                    non_zero_ups = [u for u in ups if u != 0]
                    if maxup != len(non_zero_ups) or len(set(non_zero_ups)) != len(non_zero_ups):
                        raise ValueError(f"River {comid}: invalid maxup={maxup} vs up1..up4={ups}")
                    if next_down == comid or comid in non_zero_ups:
                        raise ValueError(f"River {comid}: self-loop in supplied graph")
                    next_riv_id += 1
                    try:
                        db.execute(
                            "INSERT INTO rivers VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                            (comid, next_riv_id, pname, next_down, maxup, *ups),
                        )
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(f"Duplicate river COMID {comid}") from exc
        db.commit()
        riv_count = next_riv_id

        topo = _analyze_topology(
            db,
            cat_records,
            unavailable_catchments=unavailable_catchments,
            unavailable_rivers=unavailable_rivers,
            omitted_river_comids=omitted_river_comids,
        )
        del cat_records
        cat_conn = topo.cat_conn

        # Pass 1c: Stream valid_cat_path (and repaired_cat_path if present) -> catchments.arrow, lookup.arrow, graph.parquet.
        insert_idx = len(CATCHMENT_RAW_FIELDS)
        cat_arrow_inputs = [valid_cat_path, *([repaired_cat_path] if has_repaired_rows else [])]
        with pq.ParquetWriter(output / "graph.parquet", GRAPH_SCHEMA) as cat_graph_writer:

            def final_catchment_batches() -> Iterator[pa.RecordBatch]:
                for cat_input_path in cat_arrow_inputs:
                    with pa.memory_map(str(cat_input_path), "r") as raw_stream:
                        raw_reader = ipc.open_file(raw_stream)
                        for batch_idx in range(raw_reader.num_record_batches):
                            batch = raw_reader.get_batch(batch_idx)
                            ids = cast(list[int], batch.column("id").to_pylist())
                            bboxes = cast(list[dict[str, float]], batch.column("bbox").to_pylist())
                            down_arr = pa.array([cat_conn.downstream_status[cid - 1] for cid in ids], type=pa.string())
                            up_arr = pa.array([cat_conn.upstream_cuts[cid - 1] for cid in ids], type=pa.int32())
                            cat_graph_writer.write_table(
                                pa.Table.from_pylist(
                                    [
                                        {
                                            "id": cid,
                                            "level": 0,
                                            "upstream_ids": cat_conn.upstream_ids[cid - 1],
                                            "bbox_minx": bbox["xmin"],
                                            "bbox_miny": bbox["ymin"],
                                            "bbox_maxx": bbox["xmax"],
                                            "bbox_maxy": bbox["ymax"],
                                        }
                                        for cid, bbox in zip(ids, bboxes, strict=True)
                                    ],
                                    schema=GRAPH_SCHEMA,
                                )
                            )
                            cols = [*batch.columns[:insert_idx], down_arr, up_arr, *batch.columns[insert_idx:]]
                            yield pa.RecordBatch.from_arrays(cols, schema=cat_schema)

            cat_directory = write_vectors(
                cat_schema,
                final_catchment_batches(),
                output / "catchments.arrow",
                output / "lookup.arrow",
                lookup_binding="hfx:catchments_sha256",
            )
            (output / "lookup-directory.json").write_text(json.dumps(cat_directory) + "\n")
        valid_cat_path.unlink(missing_ok=True)
        repaired_cat_path.unlink(missing_ok=True)
        cat_conn.upstream_ids.clear()
        cat_conn.downstream_status.clear()
        cat_conn.upstream_cuts.clear()

        # Pass 2c: Stream delivered rivers -> rivers.arrow, river-graph.parquet, local-associations.parquet
        with (
            pq.ParquetWriter(output / "river-graph.parquet", GRAPH_SCHEMA) as graph_writer,
            pq.ParquetWriter(output / "local-associations.parquet", ASSOC_SCHEMA) as assoc_writer,
        ):

            def river_batches() -> Iterator[pa.RecordBatch]:
                for pname, entry in supplied_riv_sources:
                    ns = f"{release}/{pname}/rivers"
                    rename_map = riv_renames.get(pname, {})
                    cat_unavailable = pname in unavailable_catchments
                    cat_lookup = dict(
                        db.execute("SELECT comid, hfx_id FROM catchments WHERE partition = ?", (pname,)).fetchall()
                    )
                    omitted_lookup = {
                        r[0]
                        for r in db.execute(
                            "SELECT comid FROM omitted_catchments WHERE partition = ?", (pname,)
                        ).fetchall()
                    }
                    with fiona.open(entry["shp"]) as src:
                        rows: list[dict[str, Any]] = []
                        graph_rows: list[dict[str, Any]] = []
                        assoc_rows: list[dict[str, int]] = []
                        for feat in src:
                            props = dict(feat["properties"])
                            comid = int(props["COMID"])
                            if comid in omitted_river_comids:
                                continue
                            hfx_id = topo.comid_to_hfx[comid]
                            geom = shapely_shape(feat["geometry"])
                            bbox = _narrowed_bbox(geom.bounds)
                            for idx, key in enumerate(("xmin", "ymin", "xmax", "ymax")):
                                overall_bbox[idx] = (
                                    min(overall_bbox[idx], bbox[key]) if idx < 2 else max(overall_bbox[idx], bbox[key])
                                )

                            if comid in cat_lookup:
                                local_status = "supplied"
                                assoc_rows.append({"river_id": hfx_id, "catchment_id": cat_lookup[comid]})
                            elif comid in omitted_lookup:
                                local_status = "omitted_invalid_polygon"
                            elif cat_unavailable:
                                local_status = "source_unavailable"
                            else:
                                local_status = "not_in_source"

                            out_row: dict[str, Any] = {
                                "id": hfx_id,
                                "level": 0,
                                "parent_id": None,
                                "source_namespace": ns,
                                "source_id": str(comid),
                                "processing_partition": pname,
                                "bbox": bbox,
                                "geometry": geom.wkb,
                                "local_catchment_status": local_status,
                                "downstream_boundary_status": topo.downstream_status_by_id[hfx_id],
                                "upstream_boundary_cuts": topo.upstream_cuts_by_id[hfx_id],
                                "component_catchment_id": None,
                                "upstream_node_id": None,
                                "downstream_node_id": None,
                                "geometry_orientation": 0,
                            }
                            for field in riv_native_fields:
                                out_row[field.name] = None
                            for orig_name, val in props.items():
                                mapped = rename_map.get(orig_name, orig_name)
                                dtype = riv_schema.field(mapped).type
                                out_row[mapped] = _coerce_value(val, dtype)
                            rows.append(out_row)

                            ups_hfx = sorted(topo.comid_to_hfx[u] for u in topo.incoming_by_comid[comid])
                            graph_rows.append(
                                {
                                    "id": hfx_id,
                                    "level": 0,
                                    "upstream_ids": ups_hfx,
                                    "bbox_minx": bbox["xmin"],
                                    "bbox_miny": bbox["ymin"],
                                    "bbox_maxx": bbox["xmax"],
                                    "bbox_maxy": bbox["ymax"],
                                }
                            )
                            if len(rows) >= batch_rows:
                                yield pa.RecordBatch.from_pylist(rows, schema=riv_schema)
                                graph_writer.write_table(pa.Table.from_pylist(graph_rows, schema=GRAPH_SCHEMA))
                                if assoc_rows:
                                    assoc_writer.write_table(pa.Table.from_pylist(assoc_rows, schema=ASSOC_SCHEMA))
                                rows, graph_rows, assoc_rows = [], [], []
                        if rows:
                            yield pa.RecordBatch.from_pylist(rows, schema=riv_schema)
                            graph_writer.write_table(pa.Table.from_pylist(graph_rows, schema=GRAPH_SCHEMA))
                            if assoc_rows:
                                assoc_writer.write_table(pa.Table.from_pylist(assoc_rows, schema=ASSOC_SCHEMA))

            riv_directory = write_vectors(
                riv_schema,
                river_batches(),
                output / "rivers.arrow",
                output / "river-lookup.arrow",
                lookup_binding="hfx:data_sha256",
            )
            (output / "river-lookup-directory.json").write_text(json.dumps(riv_directory) + "\n")

    for path_str, expected_meta in before_hashes.items():
        if _sha256(Path(path_str)) != expected_meta["sha256"]:
            raise ValueError(f"Source file changed during normalization: {path_str}")

    geom_norm_meta = _build_geometry_normalization_meta(
        enabled=normalize_rings,
        part_stats=part_stats,
        overlap_stats=overlap_stats,
    )
    return _write_staging_outputs(
        output,
        release=release,
        organization=inventory.get("organization"),
        region=inventory.get("region", inventory.get("organization", "global")),
        original_inventory_bytes=inventory_path.read_bytes(),
        omitted_features=omitted_features,
        cut_edges=topo.cut_edges,
        crs_records=crs_records,
        cat_orig_schemas=cat_orig_schemas,
        riv_orig_schemas=riv_orig_schemas,
        cat_renames=cat_renames,
        riv_renames=riv_renames,
        unavailable_catchments=unavailable_catchments,
        unavailable_rivers=unavailable_rivers,
        cat_count=cat_count,
        channelized_count=channelized_count,
        hillslope_count=hillslope_count,
        riv_count=riv_count,
        omitted_cat_count=omitted_cat_count,
        omitted_riv_count=omitted_riv_count,
        cat_conn_counts=cat_conn.counts,
        not_in_source_rivers=cat_conn.not_in_source_rivers,
        not_in_source_rivers_truncated=cat_conn.not_in_source_rivers_truncated,
        overall_bbox=overall_bbox,
        overflow_confluences=topo.overflow_confluences,
        aux_declarations=aux_declarations,
        partial_reason=partial_reason,
        geometry_normalization_meta=geom_norm_meta,
        derived_from_verified_bundle_manifest_sha256=base_manifest_sha,
        expected_raster_sha256=raster_before,
    )


def produce_from_bundle(
    base_spatial_dir: Path,
    output: Path,
    *,
    auxiliary_spec_path: Path | None = None,
    batch_rows: int = 4096,
    partial_reason: str | None = None,
) -> dict[str, Any]:
    """Derive catchment connectivity directly from a verified HFX 0.4.0 MERIT spatial bundle."""
    base_manifest, base_manifest_sha = _open_verified_base_bundle(base_spatial_dir)

    base_meta = json.loads((base_spatial_dir / "source-metadata.json").read_bytes())
    omitted_features: list[dict[str, Any]] = json.loads((base_spatial_dir / "omitted-features.json").read_bytes())
    orig_inv_bytes = (base_spatial_dir / "original-inventory.json").read_bytes()

    release = str(base_meta["release"])
    unavailable_catchments = {str(k): str(v) for k, v in base_meta["unavailable_partitions"]["catchments"].items()}
    unavailable_rivers = {str(k): str(v) for k, v in base_meta["unavailable_partitions"]["rivers"].items()}

    output.mkdir(parents=True, exist_ok=False)
    if auxiliary_spec_path is not None:
        aux_declarations, raster_before = _check_and_stage_rasters(json.loads(auxiliary_spec_path.read_bytes()), output)
    else:
        aux_declarations, raster_before = _stage_base_rasters(base_spatial_dir, base_manifest, output)

    base_cat_pf = pq.ParquetFile(base_spatial_dir / str(base_manifest["families"]["catchments"]["data"]))
    base_riv_pf = pq.ParquetFile(base_spatial_dir / str(base_manifest["families"]["rivers"]["data"]))

    base_cat_fields = [
        f
        for f in base_cat_pf.schema_arrow.with_metadata({})
        if f.name not in {"downstream_boundary_status", "upstream_boundary_cuts"}
    ]
    base_raw_prefix = [
        f
        for f in CATCHMENT_RAW_FIELDS
        if f.name != "geometry_normalization" or any(bf.name == "geometry_normalization" for bf in base_cat_fields)
    ]
    if [f.name for f in base_cat_fields[: len(base_raw_prefix)]] != [f.name for f in base_raw_prefix]:
        raise ValueError("Unexpected base catchments.parquet prefix layout")
    cat_native_fields = base_cat_fields[len(base_raw_prefix) :]
    cat_schema = pa.schema(CATCHMENT_FIELDS + cat_native_fields)

    riv_structural_names = {f.name for f in RIVER_FIELDS}
    base_riv_fields_map = {f.name: f for f in base_riv_pf.schema_arrow.with_metadata({})}
    riv_native_fields = [f for f in base_riv_pf.schema_arrow.with_metadata({}) if f.name not in riv_structural_names]
    if not all(f.name in base_riv_fields_map for f in RIVER_FIELDS):
        raise ValueError("Unexpected base rivers.parquet structural columns")
    riv_schema = pa.schema(RIVER_FIELDS + riv_native_fields)

    cat_count = base_cat_pf.metadata.num_rows
    riv_count = base_riv_pf.metadata.num_rows
    overall_bbox = [float("inf"), float("inf"), -float("inf"), -float("inf")]
    delivered_part_counts: dict[str, int] = {}

    with _workspace_db() as (db, _scratch_dir):
        omitted_cat_count = 0
        omitted_riv_count = 0
        omitted_river_comids: dict[int, str] = {}
        for item in omitted_features:
            domain = item["domain"]
            if domain == "catchments":
                omitted_cat_count += 1
                if item.get("COMID") is not None:
                    db.execute(
                        "INSERT INTO omitted_catchments VALUES (?, ?, ?)",
                        (str(item["partition"]), int(item["COMID"]), str(item["reason"])),
                    )
            elif domain == "coastal_hillslopes":
                omitted_cat_count += 1
            elif domain == "rivers":
                omitted_riv_count += 1
                omitted_river_comids[int(item["COMID"])] = str(item["reason"])

        cat_records_opt: list[tuple[bool, str, int | None] | None] = [None] * cat_count
        channelized_count = 0
        hillslope_count = 0
        seen_global_channelized_comids: set[int] = set()
        for batch in base_cat_pf.iter_batches(
            batch_size=batch_rows * 4,
            columns=["id", "processing_partition", "catchment_kind", "source_id", "bbox"],
        ):
            for cid, pname, kind, sid, bbox in zip(
                cast(list[int], batch.column("id").to_pylist()),
                cast(list[str], batch.column("processing_partition").to_pylist()),
                cast(list[str], batch.column("catchment_kind").to_pylist()),
                cast(list[str], batch.column("source_id").to_pylist()),
                cast(list[dict[str, float]], batch.column("bbox").to_pylist()),
                strict=True,
            ):
                if not (1 <= cid <= cat_count) or cat_records_opt[cid - 1] is not None:
                    raise ValueError(f"Invalid or duplicate catchment id in base bundle: {cid}")
                for idx, key in enumerate(("xmin", "ymin", "xmax", "ymax")):
                    overall_bbox[idx] = (
                        min(overall_bbox[idx], bbox[key]) if idx < 2 else max(overall_bbox[idx], bbox[key])
                    )
                is_hillslope = kind == "coastal_hillslope"
                interned_pname = sys.intern(str(pname))
                delivered_part_counts[interned_pname] = delivered_part_counts.get(interned_pname, 0) + 1
                if is_hillslope:
                    hillslope_count += 1
                    comid: int | None = None
                else:
                    channelized_count += 1
                    comid = int(sid)
                    if comid in seen_global_channelized_comids:
                        raise ValueError(f"Duplicate global channelized catchment COMID {comid}")
                    seen_global_channelized_comids.add(comid)
                    db.execute("INSERT INTO catchments VALUES (?, ?, ?)", (interned_pname, comid, cid))
                cat_records_opt[cid - 1] = (is_hillslope, interned_pname, comid)
        del seen_global_channelized_comids
        if any(rec is None for rec in cat_records_opt):
            raise ValueError("Gap in catchment id sequence in base bundle")
        cat_records = cast(list[tuple[bool, str, int | None]], cat_records_opt)

        seen_riv_ids = 0
        for batch in base_riv_pf.iter_batches(
            batch_size=batch_rows * 4,
            columns=[
                "id",
                "processing_partition",
                "source_id",
                "NextDownID",
                "maxup",
                "up1",
                "up2",
                "up3",
                "up4",
                "bbox",
            ],
        ):
            for rid, pname, sid, next_down, maxup, up1, up2, up3, up4, bbox in zip(
                cast(list[int], batch.column("id").to_pylist()),
                cast(list[str], batch.column("processing_partition").to_pylist()),
                cast(list[str], batch.column("source_id").to_pylist()),
                cast(list[int], batch.column("NextDownID").to_pylist()),
                cast(list[int], batch.column("maxup").to_pylist()),
                cast(list[int], batch.column("up1").to_pylist()),
                cast(list[int], batch.column("up2").to_pylist()),
                cast(list[int], batch.column("up3").to_pylist()),
                cast(list[int], batch.column("up4").to_pylist()),
                cast(list[dict[str, float]], batch.column("bbox").to_pylist()),
                strict=True,
            ):
                if not (1 <= rid <= riv_count):
                    raise ValueError(f"Invalid river id in base bundle: {rid}")
                seen_riv_ids += 1
                for idx, key in enumerate(("xmin", "ymin", "xmax", "ymax")):
                    overall_bbox[idx] = (
                        min(overall_bbox[idx], bbox[key]) if idx < 2 else max(overall_bbox[idx], bbox[key])
                    )
                comid_r = int(sid)
                ups = [int(up1), int(up2), int(up3), int(up4)]
                non_zero_ups = [u for u in ups if u != 0]
                if int(maxup) != len(non_zero_ups) or len(set(non_zero_ups)) != len(non_zero_ups):
                    raise ValueError(f"River {comid_r}: invalid maxup={maxup} vs up1..up4={ups}")
                if int(next_down) == comid_r or comid_r in non_zero_ups:
                    raise ValueError(f"River {comid_r}: self-loop in base bundle")
                db.execute(
                    "INSERT INTO rivers VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (comid_r, int(rid), str(pname), int(next_down), int(maxup), *ups),
                )
        if seen_riv_ids != riv_count:
            raise ValueError("River count mismatch in base bundle")
        db.commit()

        topo = _analyze_topology(
            db,
            cat_records,
            unavailable_catchments=unavailable_catchments,
            unavailable_rivers=unavailable_rivers,
            omitted_river_comids=omitted_river_comids,
        )
        del cat_records, cat_records_opt
        cat_conn = topo.cat_conn

        base_cut_edges_path = base_spatial_dir / "cut-edges.json"
        if base_cut_edges_path.exists() and json.loads(base_cut_edges_path.read_bytes()) != topo.cut_edges:
            raise ValueError("Reconstructed river cut-edges.json differs from base bundle")
        if base_meta.get("up1_to_up4_overflow_confluences") != topo.overflow_confluences:
            raise ValueError("Reconstructed up1..up4 overflow confluences differ from base bundle")

        # Stream catchments.parquet -> catchments.arrow, lookup.arrow, graph.parquet.
        base_cat_names = list(base_cat_pf.schema_arrow.names)
        has_geom_norm = "geometry_normalization" in base_cat_names
        native_indices = [base_cat_names.index(f.name) for f in cat_native_fields]
        id_idx = base_cat_names.index("id")
        bbox_idx = base_cat_names.index("bbox")
        with pq.ParquetWriter(output / "graph.parquet", GRAPH_SCHEMA) as cat_graph_writer:

            def base_catchment_batches() -> Iterator[pa.RecordBatch]:
                for batch in base_cat_pf.iter_batches(batch_size=batch_rows):
                    ids = cast(list[int], batch.column(id_idx).to_pylist())
                    bboxes = cast(list[dict[str, float]], batch.column(bbox_idx).to_pylist())
                    down_arr = pa.array([cat_conn.downstream_status[cid - 1] for cid in ids], type=pa.string())
                    up_arr = pa.array([cat_conn.upstream_cuts[cid - 1] for cid in ids], type=pa.int32())
                    cat_graph_writer.write_table(
                        pa.Table.from_pylist(
                            [
                                {
                                    "id": cid,
                                    "level": 0,
                                    "upstream_ids": cat_conn.upstream_ids[cid - 1],
                                    **dict(
                                        zip(
                                            ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                            bbox.values(),
                                            strict=True,
                                        )
                                    ),
                                }
                                for cid, bbox in zip(ids, bboxes, strict=True)
                            ],
                            schema=GRAPH_SCHEMA,
                        )
                    )
                    prefix_cols: list[pa.Array] = []
                    for f in CATCHMENT_RAW_FIELDS:
                        if f.name == "geometry_normalization" and not has_geom_norm:
                            prefix_cols.append(pa.array(["original"] * batch.num_rows, type=pa.string()))
                        else:
                            prefix_cols.append(batch.column(base_cat_names.index(f.name)))
                    cols = [
                        *prefix_cols,
                        down_arr,
                        up_arr,
                        *[batch.column(i) for i in native_indices],
                    ]
                    yield pa.RecordBatch.from_arrays(cols, schema=cat_schema)

            cat_directory = write_vectors(
                cat_schema,
                base_catchment_batches(),
                output / "catchments.arrow",
                output / "lookup.arrow",
                lookup_binding="hfx:catchments_sha256",
            )
            (output / "lookup-directory.json").write_text(json.dumps(cat_directory) + "\n")

        cat_conn.upstream_ids.clear()
        cat_conn.downstream_status.clear()
        cat_conn.upstream_cuts.clear()

        # Stream rivers.parquet -> rivers.arrow, river-lookup.arrow, river-graph.parquet, local-associations.parquet.
        chan_lookup = dict(db.execute("SELECT comid, hfx_id FROM catchments").fetchall())
        base_riv_names = list(base_riv_pf.schema_arrow.names)
        riv_col_indices = [base_riv_names.index(f.name) for f in riv_schema]
        riv_id_idx = base_riv_names.index("id")
        riv_sid_idx = base_riv_names.index("source_id")
        riv_bbox_idx = base_riv_names.index("bbox")
        riv_ds_idx = base_riv_names.index("downstream_boundary_status")
        riv_uc_idx = base_riv_names.index("upstream_boundary_cuts")

        with (
            pq.ParquetWriter(output / "river-graph.parquet", GRAPH_SCHEMA) as graph_writer,
            pq.ParquetWriter(output / "local-associations.parquet", ASSOC_SCHEMA) as assoc_writer,
        ):

            def base_river_batches() -> Iterator[pa.RecordBatch]:
                for batch in base_riv_pf.iter_batches(batch_size=batch_rows):
                    ids = cast(list[int], batch.column(riv_id_idx).to_pylist())
                    sids = cast(list[str], batch.column(riv_sid_idx).to_pylist())
                    bboxes = cast(list[dict[str, float]], batch.column(riv_bbox_idx).to_pylist())
                    obs_ds = cast(list[str], batch.column(riv_ds_idx).to_pylist())
                    obs_uc = cast(list[int], batch.column(riv_uc_idx).to_pylist())
                    graph_rows: list[dict[str, Any]] = []
                    assoc_rows: list[dict[str, int]] = []
                    for hfx_id, sid, bbox, ds_val, uc_val in zip(ids, sids, bboxes, obs_ds, obs_uc, strict=True):
                        if topo.downstream_status_by_id[hfx_id] != ds_val or topo.upstream_cuts_by_id[hfx_id] != uc_val:
                            raise ValueError(f"River {hfx_id} topology disagrees with base bundle")
                        comid_r = int(sid)
                        cat_id = chan_lookup.get(comid_r)
                        if cat_id is not None:
                            assoc_rows.append({"river_id": hfx_id, "catchment_id": cat_id})
                        ups_hfx = sorted(topo.comid_to_hfx[u] for u in topo.incoming_by_comid[comid_r])
                        graph_rows.append(
                            {
                                "id": hfx_id,
                                "level": 0,
                                "upstream_ids": ups_hfx,
                                **dict(
                                    zip(
                                        ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy"),
                                        bbox.values(),
                                        strict=True,
                                    )
                                ),
                            }
                        )
                    graph_writer.write_table(pa.Table.from_pylist(graph_rows, schema=GRAPH_SCHEMA))
                    if assoc_rows:
                        assoc_writer.write_table(pa.Table.from_pylist(assoc_rows, schema=ASSOC_SCHEMA))
                    yield pa.RecordBatch.from_arrays([batch.column(i) for i in riv_col_indices], schema=riv_schema)

            riv_directory = write_vectors(
                riv_schema,
                base_river_batches(),
                output / "rivers.arrow",
                output / "river-lookup.arrow",
                lookup_binding="hfx:data_sha256",
            )
            (output / "river-lookup-directory.json").write_text(json.dumps(riv_directory) + "\n")

    base_repaired_src = base_spatial_dir / REPAIRED_FEATURES_ARTIFACT
    if base_repaired_src.is_file():
        shutil.copyfile(base_repaired_src, output / REPAIRED_FEATURES_ARTIFACT)
    else:
        with gzip.GzipFile(output / REPAIRED_FEATURES_ARTIFACT, mode="wb", mtime=0):
            pass
    geom_norm_meta = base_meta.get("geometry_normalization")
    if geom_norm_meta is None:
        fallback_part_stats: dict[str, dict[str, Any]] = {}

        def _fstat(pname: str) -> dict[str, Any]:
            return fallback_part_stats.setdefault(
                pname,
                {
                    "total_source": 0,
                    "already_valid": 0,
                    "attempted": 0,
                    "accepted": 0,
                    "rejected": 0,
                    "omitted_other": 0,
                    "max_abs_area_diff_km2": 0.0,
                    "max_rel_area_diff": 0.0,
                },
            )

        for pname, cnt in delivered_part_counts.items():
            fstat = _fstat(pname)
            fstat["total_source"] += cnt
            fstat["already_valid"] += cnt
        for item in omitted_features:
            if item.get("domain") == "rivers":
                continue
            pname = str(item["partition"])
            fstat = _fstat(pname)
            fstat["total_source"] += 1
            if "repair_rejection" in item:
                fstat["attempted"] += 1
                fstat["rejected"] += 1
            else:
                fstat["omitted_other"] += 1

        geom_norm_meta = _build_geometry_normalization_meta(
            enabled=False,
            part_stats=fallback_part_stats,
            overlap_stats={
                "same_partition_pairs": 0,
                "cross_partition_pairs": 0,
                "coastal_hillslope_pairs": 0,
                "candidate_candidate_pairs": 0,
                "candidate_interior_overlaps_detected": 0,
                "accepted_interior_overlaps": 0,
            },
        )

    return _write_staging_outputs(
        output,
        release=release,
        organization=base_meta.get("organization"),
        region=base_manifest.get("region", base_meta.get("organization", "global")),
        original_inventory_bytes=orig_inv_bytes,
        omitted_features=omitted_features,
        cut_edges=topo.cut_edges,
        crs_records=base_meta["crs_evidence"],
        cat_orig_schemas=base_meta["original_schemas"]["catchments"],
        riv_orig_schemas=base_meta["original_schemas"]["rivers"],
        cat_renames=base_meta["column_renames"]["catchments"],
        riv_renames=base_meta["column_renames"]["rivers"],
        unavailable_catchments=unavailable_catchments,
        unavailable_rivers=unavailable_rivers,
        cat_count=cat_count,
        channelized_count=channelized_count,
        hillslope_count=hillslope_count,
        riv_count=riv_count,
        omitted_cat_count=omitted_cat_count,
        omitted_riv_count=omitted_riv_count,
        cat_conn_counts=cat_conn.counts,
        not_in_source_rivers=cat_conn.not_in_source_rivers,
        not_in_source_rivers_truncated=cat_conn.not_in_source_rivers_truncated,
        overall_bbox=overall_bbox,
        overflow_confluences=topo.overflow_confluences,
        aux_declarations=aux_declarations,
        partial_reason=partial_reason,
        geometry_normalization_meta=geom_norm_meta,
        derived_from_verified_bundle_manifest_sha256=base_manifest_sha,
        expected_raster_sha256=raster_before,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_or_output", type=Path, help="Inventory JSON path, or output directory with --from-spatial-bundle"
    )
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--from-spatial-bundle", type=Path)
    parser.add_argument("--base-spatial-dir", type=Path)
    parser.add_argument("--auxiliary-json", type=Path)
    parser.add_argument("--batch-rows", type=int, default=4096)
    parser.add_argument("--no-normalize-rings", action="store_true")
    parser.add_argument("--partial-reason")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.from_spatial_bundle is not None:
        if args.output is not None:
            parser.error("pass only output directory when --from-spatial-bundle is specified")
        result = produce_from_bundle(
            args.from_spatial_bundle,
            args.source_or_output,
            auxiliary_spec_path=args.auxiliary_json,
            batch_rows=args.batch_rows,
            partial_reason=args.partial_reason,
        )
    else:
        if args.output is None:
            parser.error("output directory is required")
        result = produce(
            args.source_or_output,
            args.output,
            auxiliary_spec_path=args.auxiliary_json,
            base_spatial_dir=args.base_spatial_dir,
            normalize_rings=not args.no_normalize_rings,
            batch_rows=args.batch_rows,
            partial_reason=args.partial_reason,
        )
    if args.report:
        with args.report.open("x") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
