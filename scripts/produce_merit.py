"""Normalize MERIT Basins vectors (+ optional MERIT Hydro rasters) without source repair."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import fiona
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pyproj import Geod
from shapely import get_coordinates
from shapely.geometry import shape as shapely_shape
from shapely.validation import explain_validity

from hfx.candidates import fingerprint
from hfx.vector_storage import write_vectors

BOUNDS = pa.struct([pa.field(name, pa.float32(), nullable=False) for name in ("xmin", "ymin", "xmax", "ymax")])
MIN_WGS84_LON = -180.0
MAX_UNWRAPPED_ANTIMERIDIAN_LON = 195.0

STRUCTURAL_FIELDS = {
    "id",
    "level",
    "parent_id",
    "source_namespace",
    "source_id",
    "processing_partition",
    "catchment_kind",
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

CATCHMENT_FIELDS: list[pa.Field] = [
    pa.field("id", pa.int64(), nullable=False),
    pa.field("level", pa.int16(), nullable=False),
    pa.field("parent_id", pa.int64()),
    pa.field("source_namespace", pa.string(), nullable=False),
    pa.field("source_id", pa.string(), nullable=False),
    pa.field("processing_partition", pa.string(), nullable=False),
    pa.field("catchment_kind", pa.string(), nullable=False),
    pa.field("bbox", BOUNDS, nullable=False),
    pa.field("geometry", pa.binary(), nullable=False),
    pa.field("area_km2", pa.float32(), nullable=False),
    pa.field("up_area_km2", pa.float32()),
    pa.field("outlet_lon", pa.float64()),
    pa.field("outlet_lat", pa.float64()),
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


@contextmanager
def _workspace_db() -> Iterator[sqlite3.Connection]:
    with TemporaryDirectory(prefix="hfx-merit-work-") as scratch:
        conn = sqlite3.connect(Path(scratch) / "merit.sqlite")
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
            yield conn
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


def _check_and_stage_rasters(
    auxiliary_spec: dict[str, Any] | None, output: Path
) -> tuple[list[dict[str, Any]], dict[Path, str]]:
    if auxiliary_spec is None:
        return [], {}
    entries = auxiliary_spec.get("auxiliary", [])
    source_paths = cast(dict[str, str], auxiliary_spec.get("files", {}))
    raster_hashes: dict[Path, str] = {}
    for rel_path, expected_src in source_paths.items():
        src_path = Path(expected_src)
        digest = _sha256(src_path)
        raster_hashes[src_path] = digest
        dest = output / rel_path
        if Path(rel_path).is_absolute() or not dest.resolve().is_relative_to(output.resolve()):
            raise ValueError(f"Raster path escapes bundle: {rel_path}")
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
        "Vectors keep original coordinates, part order, and native attributes; invalid source polygons and degenerate "
        "zero-length self-referencing river records (NextDownID=up1=COMID, maxup=0) are omitted with explicit provenance "
        "rather than repaired. Upstream IDs in river-graph.parquet derive from verified incoming NextDownID links so "
        "confluences exceeding the shapefile's four up1..up4 slots retain all incoming edges while native maxup and up1..up4 "
        "columns stay intact. Eastern Chukotka vectors (pfaf_35 and hillslope_35) preserve native unwrapped antimeridian "
        "longitudes (>180 deg E, up to ~191.003 deg E) without clipping or splitting, while MERIT Hydro rasters span wrapped "
        "longitudes (~-180.0004 to +179.9996 deg), requiring -360 deg longitude wrapping and split antimeridian raster reads. "
        "Derived geodesic area_km2 uses the WGS84 ellipsoid and is distinct from native unitarea and uparea. Catchment "
        f"connectivity is not_supplied; outlets, nodes, and upstream polygon area are null/absent.{caveat}\n"
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


def produce(
    inventory_path: Path,
    output: Path,
    *,
    auxiliary_spec_path: Path | None = None,
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
        if h_entry.get("status") == "supplied":
            before_hashes.update(_verify_sidecars(h_entry))
            supplied_cat_sources.append((pname, h_entry))

    if not supplied_cat_sources or not supplied_riv_sources:
        raise ValueError("Both supplied catchments and rivers are required")

    auxiliary_spec = json.loads(auxiliary_spec_path.read_bytes()) if auxiliary_spec_path else None

    output.mkdir(parents=True, exist_ok=False)
    aux_declarations, raster_before = _check_and_stage_rasters(auxiliary_spec, output)

    cat_native_fields, cat_renames, cat_orig_schemas = _unify_native_schema(supplied_cat_sources)
    riv_native_fields, riv_renames, riv_orig_schemas = _unify_native_schema(supplied_riv_sources)

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

    with _workspace_db() as db:
        # Pass 1: Catchments + coastal hillslopes -> catchments.arrow
        def catchment_batches() -> Iterator[pa.RecordBatch]:
            nonlocal cat_count, hillslope_count, channelized_count, omitted_cat_count
            seen_source_ids: set[tuple[str, str]] = set()
            seen_global_channelized_comids: set[int] = set()
            for pname, entry in supplied_cat_sources:
                is_hillslope = entry in hillslopes
                kind = "coastal_hillslope" if is_hillslope else "channelized"
                ns_domain = "coastal_hillslopes" if is_hillslope else "catchments"
                ns = f"{release}/{pname}/{ns_domain}"
                rename_map = cat_renames.get(pname, {})
                with fiona.open(entry["shp"]) as src:
                    crs_records[f"{pname}/{ns_domain}"] = _check_crs(src, entry)
                    rows: list[dict[str, Any]] = []
                    for feat in src:
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

                        out_row: dict[str, Any] = {
                            "id": cat_count,
                            "level": 0,
                            "parent_id": None,
                            "source_namespace": ns,
                            "source_id": raw_sid,
                            "processing_partition": pname,
                            "catchment_kind": kind,
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
                            dtype = cat_schema.field(mapped).type
                            out_row[mapped] = _coerce_value(val, dtype)
                        rows.append(out_row)
                        if len(rows) >= batch_rows:
                            yield pa.RecordBatch.from_pylist(rows, schema=cat_schema)
                            rows = []
                    if rows:
                        yield pa.RecordBatch.from_pylist(rows, schema=cat_schema)
                db.commit()

        cat_directory = write_vectors(
            cat_schema,
            catchment_batches(),
            output / "catchments.arrow",
            output / "lookup.arrow",
            lookup_binding="hfx:catchments_sha256",
        )
        (output / "lookup-directory.json").write_text(json.dumps(cat_directory) + "\n")

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

        # Validate internal reciprocity, cuts, and global acyclicity.
        cut_edges: list[dict[str, Any]] = []
        downstream_status_by_id: dict[int, str] = {}
        upstream_cuts_by_id: dict[int, int] = {}

        # Load compact river graph arrays (3M * ~56 bytes ~ 170 MB in Python dicts or compact NumPy).
        all_rows = db.execute("SELECT comid, hfx_id, partition, next_down, up1, up2, up3, up4 FROM rivers").fetchall()
        comid_to_hfx: dict[int, int] = {r[0]: r[1] for r in all_rows}
        comid_to_ups: dict[int, tuple[int, ...]] = {r[0]: tuple(u for u in r[4:8] if u != 0) for r in all_rows}
        next_down_by_comid: dict[int, int] = {r[0]: r[3] for r in all_rows}

        # Check reverse reciprocity first: every delivered upstream reference u must have u.next_down == comid.
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

        # Global Kahn topological cycle check.
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
        del all_rows, successors, in_degree

        # Pass 2b: Stream delivered rivers -> rivers.arrow, river-graph.parquet, local-associations.parquet
        graph_schema = pa.schema(
            [
                pa.field("id", pa.int64(), nullable=False),
                pa.field("level", pa.int16(), nullable=False),
                pa.field("upstream_ids", pa.list_(pa.int64()), nullable=False),
            ]
            + [pa.field("bbox_" + name, pa.float32(), nullable=False) for name in ("minx", "miny", "maxx", "maxy")]
        )
        assoc_schema = pa.schema(
            [pa.field("river_id", pa.int64(), nullable=False), pa.field("catchment_id", pa.int64(), nullable=False)]
        )

        with (
            pq.ParquetWriter(output / "river-graph.parquet", graph_schema) as graph_writer,
            pq.ParquetWriter(output / "local-associations.parquet", assoc_schema) as assoc_writer,
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
                            hfx_id = comid_to_hfx[comid]
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
                                "downstream_boundary_status": downstream_status_by_id[hfx_id],
                                "upstream_boundary_cuts": upstream_cuts_by_id[hfx_id],
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

                            ups_hfx = sorted(comid_to_hfx[u] for u in incoming_by_comid[comid])
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
                            if len(rows) >= batch_rows:
                                yield pa.RecordBatch.from_pylist(rows, schema=riv_schema)
                                graph_writer.write_table(pa.Table.from_pylist(graph_rows, schema=graph_schema))
                                if assoc_rows:
                                    assoc_writer.write_table(pa.Table.from_pylist(assoc_rows, schema=assoc_schema))
                                rows, graph_rows, assoc_rows = [], [], []
                        if rows:
                            yield pa.RecordBatch.from_pylist(rows, schema=riv_schema)
                            graph_writer.write_table(pa.Table.from_pylist(graph_rows, schema=graph_schema))
                            if assoc_rows:
                                assoc_writer.write_table(pa.Table.from_pylist(assoc_rows, schema=assoc_schema))

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
    for path_obj, expected_digest in raster_before.items():
        if _sha256(path_obj) != expected_digest:
            raise ValueError(f"Raster source changed during normalization: {path_obj}")

    (output / "original-inventory.json").write_bytes(inventory_path.read_bytes())
    (output / "omitted-features.json").write_text(json.dumps(omitted_features, indent=2) + "\n")
    (output / "cut-edges.json").write_text(json.dumps(cut_edges, indent=2) + "\n")
    source_metadata = {
        "release": release,
        "organization": inventory.get("organization"),
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
        "unwrapped_antimeridian_max_lon": float(overall_bbox[2]) if overall_bbox[2] > 180.0 else None,
        "up1_to_up4_overflow_confluences": overflow_confluences,
        "transformations": (
            "Zero geometry repair, clipping, simplification, reversal, or dissolving. Invalid polygons and degenerate "
            "zero-length self-referencing river records (NextDownID=up1=COMID, maxup=0) omitted with explicit provenance. "
            "River graph upstream_ids derived from verified NextDownID links so >4-tributary confluences keep all incoming "
            "edges while native maxup and up1..up4 attributes remain intact. Eastern Chukotka features (pfaf_35 and hillslope_35) "
            "retain native unwrapped antimeridian longitudes (>180 deg E) without clipping or splitting, while MERIT Hydro "
            "rasters use wrapped ~[-180, +180) longitudes (consumers querying rasters across >180 deg E must wrap longitudes "
            "by -360 deg and split reads across the antimeridian). Geodesic WGS84 ellipsoid area_km2 computed over native "
            "parts. Native COMIDs preserved; dataset-local HFX IDs assigned independently per domain."
        ),
    }
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
        "adapter_version": "merit-original-1",
        "region": inventory.get("region", inventory.get("organization", "global")),
        "catchments": "catchments.arrow",
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
        "catchment_connectivity": [{"level": 0, "status": "not_supplied"}],
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
    return {
        "status": "passed",
        "delivered_catchments": cat_count,
        "delivered_channelized_catchments": channelized_count,
        "delivered_coastal_hillslopes": hillslope_count,
        "delivered_rivers": riv_count,
        "omitted_catchments": omitted_cat_count,
        "omitted_rivers": omitted_riv_count,
        "cut_edges": len(cut_edges),
        "manifest_sha256": _sha256(output / "manifest.json"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--auxiliary-json", type=Path)
    parser.add_argument("--batch-rows", type=int, default=4096)
    parser.add_argument("--partial-reason")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = produce(
        args.inventory,
        args.output,
        auxiliary_spec_path=args.auxiliary_json,
        batch_rows=args.batch_rows,
        partial_reason=args.partial_reason,
    )
    if args.report:
        with args.report.open("x") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
