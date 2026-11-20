"""Independent binary shapefile (.shp/.dbf) verification against final MERIT HFX spatial bundle."""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import math
import sqlite3
import struct
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq
import shapely
from pyproj import Geod
from shapely import STRtree, from_wkb, get_coordinates, prepare
from shapely import box as shapely_box
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import orient
from shapely.validation import explain_validity

from hfx.candidates import fingerprint

GEOMETRY_NORMALIZATION_RULE = "merit_ring_self_intersection_lossless_split_v1"
GEODESIC_AREA_ABS_TOL_KM2 = 1e-9
GEODESIC_AREA_REL_TOL = 1e-12
REPAIRED_FEATURES_ARTIFACT = "repaired-features.jsonl.gz"
MIN_WGS84_LON = -180.0
MAX_UNWRAPPED_ANTIMERIDIAN_LON = 195.0
ANTIMERIDIAN_QUERY_MIN_X = -165.0
ANTIMERIDIAN_MAX_X = 180.0

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


def _sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _read_dbf_fields(dbf_path: Path) -> tuple[list[tuple[str, str, int, int]], int, int, int]:
    with dbf_path.open("rb") as stream:
        header = stream.read(32)
        if len(header) != 32:
            raise ValueError(f"Truncated DBF header: {dbf_path}")
        _, _, _, _, num_records, header_bytes, record_bytes = struct.unpack("<BBBBLHH20x", header)
        fields: list[tuple[str, str, int, int]] = []
        pos = 32
        while pos + 1 < header_bytes:
            marker = stream.read(1)
            if marker == b"\r":
                break
            desc = marker + stream.read(31)
            if len(desc) != 32:
                raise ValueError(f"Truncated DBF field descriptor: {dbf_path}")
            name = desc[:11].split(b"\x00", 1)[0].decode("ascii")
            ftype = chr(desc[11])
            flen = desc[16]
            fdec = desc[17]
            fields.append((name, ftype, flen, fdec))
            pos += 32
    return fields, num_records, header_bytes, record_bytes


def _read_dbf(dbf_path: Path) -> Iterator[dict[str, Any]]:
    """Parse dBASE III/IV records independently without Fiona or GDAL."""
    fields, num_records, header_bytes, record_bytes = _read_dbf_fields(dbf_path)
    with dbf_path.open("rb") as stream:
        stream.seek(header_bytes)
        for _ in range(num_records):
            raw = stream.read(record_bytes)
            if len(raw) != record_bytes:
                raise ValueError(f"Truncated DBF record in {dbf_path}")
            if raw[0] == 0x2A:
                raise ValueError(f"Deleted DBF record present in {dbf_path}")
            row: dict[str, Any] = {}
            offset = 1
            for name, ftype, flen, fdec in fields:
                chunk = raw[offset : offset + flen].decode("utf-8", errors="strict").strip()
                offset += flen
                if not chunk or (ftype in ("N", "F") and set(chunk) == {"*"}):
                    row[name] = None
                elif ftype in ("N", "F"):
                    if fdec == 0 and "." not in chunk:
                        row[name] = int(chunk)
                    else:
                        row[name] = float(chunk)
                else:
                    row[name] = chunk
            yield row


def _read_shp_parts(shp_path: Path) -> Iterator[tuple[int, list[tuple[tuple[float, float], ...]]]]:
    """Parse ESRI Shapefile PolyLine/Polygon parts directly via struct."""
    with shp_path.open("rb") as stream:
        header = stream.read(100)
        if len(header) != 100 or struct.unpack_from(">I", header, 0)[0] != 9994:
            raise ValueError(f"Invalid SHP header: {shp_path}")
        while True:
            rec_hdr = stream.read(8)
            if not rec_hdr:
                break
            if len(rec_hdr) != 8:
                raise ValueError(f"Truncated SHP record header: {shp_path}")
            _, content_words = struct.unpack(">2I", rec_hdr)
            payload = stream.read(content_words * 2)
            if len(payload) != content_words * 2:
                raise ValueError(f"Truncated SHP record payload: {shp_path}")
            shape_type = struct.unpack_from("<I", payload, 0)[0]
            if shape_type == 0:
                yield shape_type, []
                continue
            if shape_type not in (3, 5):
                raise ValueError(f"Unsupported SHP shape type {shape_type} in {shp_path}")
            num_parts, num_points = struct.unpack_from("<2I", payload, 36)
            parts = struct.unpack_from(f"<{num_parts}i", payload, 44)
            pts_offset = 44 + 4 * num_parts
            raw_pts = struct.unpack_from(f"<{num_points * 2}d", payload, pts_offset)
            rings: list[tuple[tuple[float, float], ...]] = []
            for idx, start in enumerate(parts):
                end = parts[idx + 1] if idx + 1 < num_parts else num_points
                rings.append(tuple((raw_pts[2 * k], raw_pts[2 * k + 1]) for k in range(start, end)))
            yield shape_type, rings


def _wkb_rings(wkb: bytes) -> list[tuple[tuple[float, float], ...]]:
    """Decode WKB LineString/MultiLineString/Polygon/MultiPolygon into ordered coordinate tuples."""
    buf = memoryview(wkb)

    def read_geom(offset: int) -> tuple[int, list[tuple[tuple[float, float], ...]]]:
        endian = "<" if buf[offset] == 1 else ">"
        gtype = struct.unpack_from(endian + "I", buf, offset + 1)[0]
        pos = offset + 5
        if gtype == 2:  # LineString
            count = struct.unpack_from(endian + "I", buf, pos)[0]
            pos += 4
            coords = struct.unpack_from(endian + f"{count * 2}d", buf, pos)
            pos += count * 16
            return pos, [tuple((coords[2 * i], coords[2 * i + 1]) for i in range(count))]
        if gtype == 3:  # Polygon
            n_rings = struct.unpack_from(endian + "I", buf, pos)[0]
            pos += 4
            out: list[tuple[tuple[float, float], ...]] = []
            for _ in range(n_rings):
                n_pts = struct.unpack_from(endian + "I", buf, pos)[0]
                pos += 4
                coords = struct.unpack_from(endian + f"{n_pts * 2}d", buf, pos)
                pos += n_pts * 16
                out.append(tuple((coords[2 * i], coords[2 * i + 1]) for i in range(n_pts)))
            return pos, out
        if gtype in (5, 6):  # MultiLineString / MultiPolygon
            n_geoms = struct.unpack_from(endian + "I", buf, pos)[0]
            pos += 4
            out = []
            for _ in range(n_geoms):
                pos, sub = read_geom(pos)
                out.extend(sub)
            return pos, out
        raise ValueError(f"Unsupported WKB geometry type {gtype}")

    end, rings = read_geom(0)
    if end != len(wkb):
        raise ValueError("Trailing WKB bytes")
    return rings


def _canon_ring_signature(rings: Sequence[Sequence[tuple[float, float]]], *, is_polygon: bool) -> str:
    """Hash vertex coordinates preserving exact part/ring coordinate precision."""
    h = hashlib.sha256()
    if is_polygon:
        # Polygon ring orientation or exterior/hole grouping may normalize clockwise->CCW;
        # normalize cyclic start + direction per ring so exact vertex set and ring sequence match.
        sigs: list[bytes] = []
        for ring in rings:
            pts = list(ring[:-1] if len(ring) > 1 and ring[0] == ring[-1] else ring)
            fwd_idx = pts.index(min(pts))
            fwd_rot = pts[fwd_idx:] + pts[:fwd_idx]
            rev = list(reversed(pts))
            rev_idx = rev.index(min(rev))
            rev_rot = rev[rev_idx:] + rev[:rev_idx]
            best = min(fwd_rot, rev_rot)
            sigs.append(struct.pack(f"<{len(best) * 2}d", *[c for pt in best for c in pt]))
        for sig in sorted(sigs):
            h.update(struct.pack("<I", len(sig)))
            h.update(sig)
    else:
        for part in rings:
            h.update(struct.pack("<I", len(part)))
            h.update(struct.pack(f"<{len(part) * 2}d", *[c for pt in part for c in pt]))
    return h.hexdigest()


def _canon_attrs_digest(attrs: dict[str, Any]) -> str:
    norm = {}
    for k, v in sorted(attrs.items()):
        if isinstance(v, float):
            norm[k] = struct.pack("!d", v).hex()
        else:
            norm[k] = v
    return hashlib.sha256(json.dumps(norm, sort_keys=True).encode()).hexdigest()


def _raw_river_is_invalid_or_out_of_bounds(rings: list[tuple[tuple[float, float], ...]]) -> bool:
    if not rings:
        return True
    distinct_pts: set[tuple[float, float]] = set()
    for ring in rings:
        if len(ring) < 2:
            return True
        for x, y in ring:
            if not (math.isfinite(x) and math.isfinite(y)) or x < -180.0 or x > 195.0 or abs(y) > 90.0:
                return True
            distinct_pts.add((x, y))
    return len(distinct_pts) < 2


def _signed_ring_area(ring: tuple[tuple[float, float], ...]) -> float:
    total = 0.0
    for (x1, y1), (x2, y2) in zip(ring[:-1], ring[1:], strict=True):
        total += x1 * y2 - x2 * y1
    return 0.5 * total


def _rings_to_raw_polygon(rings: list[tuple[tuple[float, float], ...]]) -> Polygon | MultiPolygon | None:
    if not rings:
        return None
    parts: list[tuple[tuple[tuple[float, float], ...], list[tuple[tuple[float, float], ...]]]] = []
    for ring in rings:
        if len(ring) < 4:
            return None
        if _signed_ring_area(ring) <= 0 or not parts:
            parts.append((ring, []))
        else:
            parts[-1][1].append(ring)
    if len(parts) == 1:
        return Polygon(parts[0][0], parts[0][1])
    return MultiPolygon([Polygon(ext, holes) for ext, holes in parts])


def _geodesic_area_km2(geom: Any, geod: Geod) -> float:
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    area = 0.0
    for poly in parts:
        ext = abs(geod.polygon_area_perimeter(*poly.exterior.xy)[0])
        holes = sum(abs(geod.polygon_area_perimeter(*ring.xy)[0]) for ring in poly.interiors)
        area += (ext - holes) / 1e6
    return area


def _canon_edge_vertex_digest(rings: Sequence[Sequence[tuple[float, float]]]) -> tuple[str, int]:
    """Hash exact IEEE 754 float64 vertex set and undirected edge multiset across rings."""
    verts: set[tuple[float, float]] = set()
    edges: list[tuple[tuple[float, float], tuple[float, float]]] = []
    zero_len = 0
    for ring in rings:
        pts = [(0.0 if x == 0.0 else float(x), 0.0 if y == 0.0 else float(y)) for x, y in ring]
        verts.update(pts)
        for a, b in zip(pts[:-1], pts[1:], strict=True):
            if a == b:
                zero_len += 1
            edges.append((a, b) if a <= b else (b, a))
    h = hashlib.sha256()
    sorted_verts = sorted(verts)
    h.update(struct.pack("<I", len(sorted_verts)))
    for x, y in sorted_verts:
        h.update(struct.pack("<2d", x, y))
    edges.sort()
    h.update(struct.pack("<I", len(edges)))
    for (x1, y1), (x2, y2) in edges:
        h.update(struct.pack("<4d", x1, y1, x2, y2))
    return h.hexdigest(), zero_len


def _raw_catchment_is_invalid_or_out_of_bounds(rings: list[tuple[tuple[float, float], ...]], geod: Geod) -> bool:
    if not rings:
        return True
    for ring in rings:
        for x, y in ring:
            if (
                not (math.isfinite(x) and math.isfinite(y))
                or x < MIN_WGS84_LON
                or x > MAX_UNWRAPPED_ANTIMERIDIAN_LON
                or abs(y) > 90.0
            ):
                return True
    geom = _rings_to_raw_polygon(rings)
    if geom is None or geom.is_empty or not geom.is_valid:
        return True
    area_km2 = _geodesic_area_km2(geom, geod)
    return not math.isfinite(area_km2) or area_km2 <= 0


def _part_and_hole_counts(geom: Any) -> tuple[int, int]:
    parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    return len(parts), sum(len(poly.interiors) for poly in parts)


def _shift_x(geom: Any, dx: float) -> Any:
    return shapely.transform(geom, lambda coords: coords + np.array([dx, 0.0], dtype=np.float64))


def _check_acyclic_downstream(edges: set[tuple[int, int]], node_ids: Iterable[int], label: str) -> None:
    """Verify out-degree <= 1 and global acyclicity via downstream path traversal."""
    successor: dict[int, int] = {}
    for src_id, dst_id in edges:
        if src_id == dst_id:
            raise ValueError(f"{label} graph contains a self-loop at {src_id}")
        if src_id in successor:
            raise ValueError(f"{label} graph is not a tree: multiple downstream edges from {src_id}")
        successor[src_id] = dst_id
    completed: set[int] = set()
    for node_id in node_ids:
        path_ids: set[int] = set()
        cur = node_id
        while cur in successor and cur not in completed:
            if cur in path_ids:
                raise ValueError(f"{label} graph contains a cycle")
            path_ids.add(cur)
            cur = successor[cur]
        completed.update(path_ids)


@contextmanager
def _verify_db(scratch_parent: Path | None = None) -> Iterator[sqlite3.Connection]:
    with TemporaryDirectory(prefix="hfx-merit-verify-", dir=scratch_parent) as scratch:
        conn = sqlite3.connect(Path(scratch) / "verify.sqlite")
        try:
            conn.execute("PRAGMA cache_size = -32768")
            conn.execute("PRAGMA temp_store = FILE")
            conn.execute("PRAGMA synchronous = OFF")
            conn.execute("PRAGMA journal_mode = MEMORY")
            conn.execute(
                "CREATE TABLE expected_features ("
                "domain TEXT NOT NULL, source_namespace TEXT NOT NULL, source_id TEXT NOT NULL, "
                "attr_digest TEXT NOT NULL, geom_sig TEXT NOT NULL, is_repaired INTEGER NOT NULL DEFAULT 0, "
                "raw_area_km2 REAL NOT NULL DEFAULT 0.0, seen INTEGER NOT NULL DEFAULT 0, "
                "PRIMARY KEY(domain, source_namespace, source_id))"
            )
            conn.execute(
                "CREATE TABLE id_map ("
                "domain TEXT NOT NULL, source_namespace TEXT NOT NULL, source_id TEXT NOT NULL, "
                "hfx_id INTEGER NOT NULL, partition TEXT NOT NULL, comid INTEGER, "
                "PRIMARY KEY(domain, source_namespace, source_id))"
            )
            conn.execute("CREATE INDEX idx_id_hfx ON id_map(domain, hfx_id)")
            conn.execute(
                "CREATE TABLE observed_catchments ("
                "hfx_id INTEGER PRIMARY KEY, catchment_kind TEXT NOT NULL, partition TEXT NOT NULL, "
                "comid INTEGER, downstream_boundary_status TEXT NOT NULL, upstream_boundary_cuts INTEGER)"
            )
            conn.execute(
                "CREATE TABLE expected_rivers ("
                "comid INTEGER PRIMARY KEY, partition TEXT NOT NULL, next_down INTEGER NOT NULL, "
                "up1 INTEGER NOT NULL, up2 INTEGER NOT NULL, up3 INTEGER NOT NULL, up4 INTEGER NOT NULL)"
            )
            yield conn
        finally:
            conn.close()


def _validate_base_bundle_inventory(
    inventory: dict[str, Any],
    inventory_path_bytes: bytes,
    base_spatial_dir: Path,
    spatial_dir: Path,
    *,
    inherited_hillslopes: set[str],
    supplied_hillslopes: set[str],
) -> tuple[dict[str, Any], list[str]]:
    """Verify inherited fallback eligibility and provenance against base_spatial_dir."""
    base_manifest = json.loads((base_spatial_dir / "manifest.json").read_bytes())
    if (
        base_manifest.get("format_version") != "0.4.0"
        or base_manifest.get("fabric_name") != "merit"
        or base_manifest.get("representation") != "spatial-geoparquet-adjacency"
        or base_manifest.get("representation_version") != 1
    ):
        raise ValueError("Invalid base_spatial_dir MERIT HFX 0.4.0 bundle")

    for rel, identity in base_manifest["artifacts"].items():
        path = base_spatial_dir / rel
        if (
            Path(rel).is_absolute()
            or not path.resolve().is_relative_to(base_spatial_dir.resolve())
            or not path.is_file()
            or fingerprint(path) != identity
        ):
            raise ValueError(f"Base bundle artifact differs or is missing: {rel}")

    partitions = sorted(inventory["partitions"], key=lambda item: item["partition"])
    for part in partitions:
        pname = str(part["partition"])
        if part["catchments"].get("status") != "supplied" or part["rivers"].get("status") != "supplied":
            raise ValueError(f"All pfaf_* partitions must be supplied when verifying with base_spatial_dir: {pname}")

    base_orig_path = base_spatial_dir / "original-inventory.json"
    spatial_orig_bytes = (spatial_dir / "original-inventory.json").read_bytes()
    if spatial_orig_bytes not in {base_orig_path.read_bytes(), inventory_path_bytes}:
        raise ValueError("original-inventory.json differs from inventory and base_spatial_dir")

    base_orig = json.loads(base_orig_path.read_bytes())
    if base_orig.get("release") != inventory.get("release"):
        raise ValueError("Release differs from base original-inventory.json")

    base_parts = {str(p["partition"]): p for p in base_orig["partitions"]}
    if {str(p["partition"]) for p in partitions} != set(base_parts):
        raise ValueError("pfaf_* partition set differs from base original-inventory.json")
    for part in partitions:
        pname = str(part["partition"])
        base_part = base_parts.get(pname)
        if base_part is None:
            raise ValueError(f"Partition absent from base original-inventory.json: {pname}")
        for domain in ("catchments", "rivers"):
            if base_part[domain].get("status") != "supplied":
                raise ValueError(f"Base original-inventory.json lacks supplied {pname}/{domain}")
            inv_sides = {
                name: (int(meta["size_bytes"]), str(meta["sha256"])) for name, meta in part[domain]["sidecars"].items()
            }
            base_sides = {
                name: (int(meta["size_bytes"]), str(meta["sha256"]))
                for name, meta in base_part[domain]["sidecars"].items()
            }
            if inv_sides != base_sides:
                raise ValueError(f"Supplied sidecars differ from base original-inventory.json for {pname}/{domain}")

    base_hills = {str(h["partition"]): h for h in base_orig.get("coastal_hillslopes", [])}
    inv_hills = {str(h["partition"]): h for h in inventory.get("coastal_hillslopes", [])}
    if set(inv_hills) != set(base_hills):
        raise ValueError("Coastal hillslope partition set differs from base original-inventory.json")
    for pname in supplied_hillslopes:
        base_h = base_hills[pname]
        inv_h = inv_hills[pname]
        if base_h.get("status") != "supplied":
            raise ValueError(f"Base original-inventory.json lacks supplied {pname}")
        inv_sides = {name: (int(meta["size_bytes"]), str(meta["sha256"])) for name, meta in inv_h["sidecars"].items()}
        base_sides = {name: (int(meta["size_bytes"]), str(meta["sha256"])) for name, meta in base_h["sidecars"].items()}
        if inv_sides != base_sides:
            raise ValueError(f"Supplied hillslope sidecars differ from base original-inventory.json for {pname}")
    for pname in inherited_hillslopes:
        if base_hills[pname].get("status") != "supplied":
            raise ValueError(f"Inherited hillslope {pname} was not supplied in base original-inventory.json")

    return base_manifest, sorted(base_hills)


def verify_spatial_bundle(
    inventory_path: Path,
    spatial_dir: Path,
    *,
    base_spatial_dir: Path | None = None,
) -> dict[str, Any]:
    """Independently verify original shapefiles against the final GeoParquet/IPC spatial bundle."""
    manifest = json.loads((spatial_dir / "manifest.json").read_bytes())
    if manifest["representation"] != "spatial-geoparquet-adjacency" or manifest["representation_version"] != 1:
        raise ValueError("Unexpected representation")
    expected_connectivity = [
        {
            "level": 0,
            "status": "supplied",
            "derivation": "river_local_associations",
            "completeness": "partial",
            "upstream_cuts_column": "upstream_boundary_cuts",
        }
    ]
    if (
        manifest["catchment_connectivity"] != expected_connectivity
        or manifest["topology"] != "tree"
        or manifest["has_up_area"]
    ):
        raise ValueError("Unexpected catchment connectivity, topology, or has_up_area")

    inventory_bytes = inventory_path.read_bytes()
    inventory = json.loads(inventory_bytes)
    release = str(inventory["release"])
    partitions = sorted(inventory["partitions"], key=lambda item: item["partition"])
    hillslopes = sorted(inventory.get("coastal_hillslopes", []), key=lambda item: item["partition"])

    raw_verified_hillslopes: list[str] = []
    inherited_hillslopes: set[str] = set()
    for h in hillslopes:
        pname = str(h["partition"])
        status = h.get("status")
        if status == "supplied":
            raw_verified_hillslopes.append(pname)
        elif status == "inherited_from_verified_bundle":
            inherited_hillslopes.add(pname)
        else:
            raise ValueError(f"Unsupported coastal_hillslopes status {status!r} for {pname}")

    if inherited_hillslopes and base_spatial_dir is None:
        raise ValueError("Inherited coastal_hillslopes require base_spatial_dir")

    base_manifest: dict[str, Any] | None = None
    if base_spatial_dir is not None:
        base_manifest, _ = _validate_base_bundle_inventory(
            inventory,
            inventory_bytes,
            base_spatial_dir,
            spatial_dir,
            inherited_hillslopes=inherited_hillslopes,
            supplied_hillslopes=set(raw_verified_hillslopes),
        )
    elif (spatial_dir / "original-inventory.json").read_bytes() != inventory_bytes:
        raise ValueError("original-inventory.json differs from inventory")

    for rel, identity in manifest["artifacts"].items():
        path = spatial_dir / rel
        if not path.resolve().is_relative_to(spatial_dir.resolve()) or fingerprint(path) != identity:
            raise ValueError(f"Final artifact fingerprint differs: {rel}")

    if REPAIRED_FEATURES_ARTIFACT not in manifest["artifacts"]:
        raise ValueError(f"Missing {REPAIRED_FEATURES_ARTIFACT} in manifest artifacts")

    omitted_list = json.loads((spatial_dir / "omitted-features.json").read_bytes())
    omitted_keys = {(item["source_namespace"], item["source_id"]): item for item in omitted_list}
    omitted_river_comids = {
        int(item["COMID"]) for item in omitted_list if item["domain"] == "rivers" and item.get("COMID") is not None
    }
    cut_edges = json.loads((spatial_dir / "cut-edges.json").read_bytes())

    repaired_records: dict[tuple[str, str], dict[str, Any]] = {}
    with gzip.open(spatial_dir / REPAIRED_FEATURES_ARTIFACT, "rt", encoding="utf-8") as rep_stream:
        for line in rep_stream:
            stripped = line.strip()
            if not stripped:
                continue
            rec = json.loads(stripped)
            rkey = (str(rec["source_namespace"]), str(rec["source_id"]))
            if rkey in repaired_records:
                raise ValueError(f"Duplicate entry in {REPAIRED_FEATURES_ARTIFACT}: {rkey}")
            if rkey in omitted_keys:
                raise ValueError(f"Feature in both {REPAIRED_FEATURES_ARTIFACT} and omitted-features.json: {rkey}")
            if rec.get("rule") != GEOMETRY_NORMALIZATION_RULE:
                raise ValueError(f"Unexpected repair rule in {REPAIRED_FEATURES_ARTIFACT}: {rec.get('rule')!r}")
            repaired_records[rkey] = rec

    source_hashes: dict[str, str] = {}
    geod = Geod(ellps="WGS84")
    obs_part_counts: dict[str, dict[str, int]] = {}

    def _opstat(pname: str) -> dict[str, int]:
        return obs_part_counts.setdefault(pname, {"already_valid": 0, "accepted": 0, "rejected": 0, "omitted_other": 0})

    with _verify_db(scratch_parent=spatial_dir.parent) as db:
        # 1. Verify original sidecar hashes and record expected delivered/omitted source features directly from .shp/.dbf.
        cat_entries: list[tuple[str, str, dict[str, Any]]] = []
        riv_entries: list[tuple[str, dict[str, Any]]] = []
        unavailable_catchments: set[str] = set()
        unavailable_rivers: set[str] = set()
        for part in partitions:
            pname = str(part["partition"])
            if part["catchments"].get("status") == "supplied":
                cat_entries.append((pname, "catchments", part["catchments"]))
            else:
                unavailable_catchments.add(pname)
            if part["rivers"].get("status") == "supplied":
                riv_entries.append((pname, part["rivers"]))
            else:
                unavailable_rivers.add(pname)
        for h in hillslopes:
            if h.get("status") == "supplied":
                cat_entries.append((str(h["partition"]), "coastal_hillslopes", h))

        for _, _, entry in cat_entries:
            for side in entry["sidecars"].values():
                digest = _sha256(Path(side["path"]))
                if digest != side["sha256"]:
                    raise ValueError(f"Original sidecar hash mismatch: {side['path']}")
                source_hashes[side["path"]] = digest
        for _, entry in riv_entries:
            for side in entry["sidecars"].values():
                digest = _sha256(Path(side["path"]))
                if digest != side["sha256"]:
                    raise ValueError(f"Original sidecar hash mismatch: {side['path']}")
                source_hashes[side["path"]] = digest

        delivered_catchments_by_part_comid: set[tuple[str, int]] = set()
        omitted_catchments_by_part_comid: set[tuple[str, int]] = set()

        for pname, ns_domain, entry in cat_entries:
            ns = f"{release}/{pname}/{ns_domain}"
            shp_path = Path(entry["shp"])
            dbf_path = shp_path.with_suffix(".dbf")
            dbf_rows = _read_dbf(dbf_path)
            shp_rows = _read_shp_parts(shp_path)
            pstat = _opstat(pname)
            for attrs, (_, rings) in zip(dbf_rows, shp_rows, strict=True):
                comid = int(attrs["COMID"]) if "COMID" in attrs and attrs["COMID"] is not None else None
                raw_sid = str(comid) if comid is not None else str(attrs.get("FID"))
                key = (ns, raw_sid)
                if key in omitted_keys:
                    if not _raw_catchment_is_invalid_or_out_of_bounds(rings, geod):
                        raise ValueError(f"Valid catchment unexpectedly omitted: {key}")
                    if ns_domain == "catchments" and comid is not None:
                        omitted_catchments_by_part_comid.add((pname, comid))
                    omitted_entry = omitted_keys[key]
                    if "repair_rejection" in omitted_entry:
                        raw_poly = _rings_to_raw_polygon(rings)
                        if raw_poly is None or not str(explain_validity(raw_poly)).startswith(
                            "Ring Self-intersection["
                        ):
                            raise ValueError(f"Rejected repair feature lacked Ring Self-intersection reason: {key}")
                        if omitted_entry["repair_rejection"] == "geodesic_area_mismatch":
                            rep_poly = shapely.make_valid(raw_poly, method="structure", keep_collapsed=False)
                            raw_a = _geodesic_area_km2(raw_poly, geod)
                            rep_a = _geodesic_area_km2(rep_poly, geod)
                            if abs(rep_a - raw_a) <= max(GEODESIC_AREA_ABS_TOL_KM2, GEODESIC_AREA_REL_TOL * raw_a):
                                raise ValueError(f"Rejected feature marked geodesic_area_mismatch within tol: {key}")
                        pstat["rejected"] += 1
                    else:
                        pstat["omitted_other"] += 1
                    continue
                if ns_domain == "catchments" and comid is not None:
                    delivered_catchments_by_part_comid.add((pname, comid))
                renamed = {("native_" + k if k in STRUCTURAL_FIELDS else k): v for k, v in attrs.items()}
                if key in repaired_records:
                    rec = repaired_records[key]
                    shp_ev_digest, shp_zero = _canon_edge_vertex_digest(rings)
                    orig_wkb = base64.b64decode(rec["original_wkb_base64"])
                    orig_ev_digest, orig_zero = _canon_edge_vertex_digest(_wkb_rings(orig_wkb))
                    if shp_zero > 0 or orig_zero > 0 or shp_ev_digest != orig_ev_digest:
                        raise ValueError(f"repaired-features.jsonl.gz original_wkb differs from raw shapefile: {key}")
                    raw_geom = from_wkb(orig_wkb)
                    if raw_geom.is_valid:
                        raise ValueError(f"Repaired feature had valid original geometry: {key}")
                    raw_reason = str(explain_validity(raw_geom))
                    if not raw_reason.startswith("Ring Self-intersection["):
                        raise ValueError(f"Repaired feature had ineligible raw validity reason {raw_reason!r}: {key}")
                    raw_area = _geodesic_area_km2(raw_geom, geod)
                    if not math.isfinite(raw_area) or raw_area <= 0:
                        raise ValueError(f"Repaired feature had non-positive raw geodesic area: {key}")
                    raw_parts, raw_holes = _part_and_hole_counts(raw_geom)
                    if raw_parts != rec.get("raw_part_count") or raw_holes != rec.get("raw_hole_count"):
                        raise ValueError(f"Repaired feature raw part/hole counts disagree with provenance: {key}")
                    if rec.get("raw_validity_reason") != raw_reason:
                        raise ValueError(f"Repaired feature raw_validity_reason disagrees with raw geometry: {key}")
                    expected_checks = {
                        "valid_polygonal": True,
                        "vertices_exact": True,
                        "edges_exact": True,
                        "geodesic_area_within_tolerance": True,
                        "zero_sibling_interior_overlap": True,
                    }
                    if rec.get("checks") != expected_checks:
                        raise ValueError(f"Repaired feature checks disagree with contract: {key}")
                    rec_raw_area = float(rec.get("raw_geodesic_area_km2", -1.0))
                    rec_rep_area = float(rec.get("repaired_geodesic_area_km2", -1.0))
                    if abs(rec_raw_area - raw_area) > max(1e-12, 1e-12 * raw_area):
                        raise ValueError(f"Repaired feature raw_geodesic_area_km2 disagrees with raw shapefile: {key}")
                    if abs(rec_rep_area - rec_raw_area) > max(
                        GEODESIC_AREA_ABS_TOL_KM2, GEODESIC_AREA_REL_TOL * rec_raw_area
                    ):
                        raise ValueError(f"Repaired feature provenance area diff exceeds tolerance: {key}")
                    rec.pop("original_wkb_base64", None)
                    pstat["accepted"] += 1
                    db.execute(
                        "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 1, ?, 0)",
                        ("catchments", ns, raw_sid, _canon_attrs_digest(renamed), shp_ev_digest, raw_area),
                    )
                else:
                    pstat["already_valid"] += 1
                    db.execute(
                        "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 0, 0.0, 0)",
                        (
                            "catchments",
                            ns,
                            raw_sid,
                            _canon_attrs_digest(renamed),
                            _canon_ring_signature(rings, is_polygon=True),
                        ),
                    )

        for pname, entry in riv_entries:
            ns = f"{release}/{pname}/rivers"
            shp_path = Path(entry["shp"])
            dbf_path = shp_path.with_suffix(".dbf")
            dbf_rows = _read_dbf(dbf_path)
            shp_rows = _read_shp_parts(shp_path)
            for attrs, (_, rings) in zip(dbf_rows, shp_rows, strict=True):
                comid = int(attrs["COMID"])
                raw_sid = str(comid)
                if (ns, raw_sid) in omitted_keys:
                    if not _raw_river_is_invalid_or_out_of_bounds(rings):
                        raise ValueError(f"Valid river unexpectedly omitted: {(ns, raw_sid)}")
                    continue
                renamed = {("native_" + k if k in STRUCTURAL_FIELDS else k): v for k, v in attrs.items()}
                db.execute(
                    "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 0, 0.0, 0)",
                    (
                        "rivers",
                        ns,
                        raw_sid,
                        _canon_attrs_digest(renamed),
                        _canon_ring_signature(rings, is_polygon=False),
                    ),
                )
                db.execute(
                    "INSERT INTO expected_rivers VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        comid,
                        pname,
                        int(attrs["NextDownID"]),
                        int(attrs["up1"]),
                        int(attrs["up2"]),
                        int(attrs["up3"]),
                        int(attrs["up4"]),
                    ),
                )
        db.commit()

        # 2. Stream catchments.parquet & rivers.parquet against expected features.
        cat_native_names: dict[tuple[str, str], set[str]] = {}
        for pname, ns_domain, entry in cat_entries:
            fields, _, _, _ = _read_dbf_fields(Path(entry["shp"]).with_suffix(".dbf"))
            cat_native_names[(pname, ns_domain)] = {
                ("native_" + n if n in STRUCTURAL_FIELDS else n) for n, _, _, _ in fields
            }

        cat_parquet = pq.ParquetFile(spatial_dir / manifest["families"]["catchments"]["data"])
        cat_structural = {
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
            "downstream_boundary_status",
            "upstream_boundary_cuts",
        }
        cat_all_native = [n for n in cat_parquet.schema_arrow.names if n not in cat_structural]

        norm_dict = next(
            (
                entry
                for entry in manifest.get("categorical_columns", [])
                if entry.get("artifact") == manifest["families"]["catchments"]["data"]
                and entry.get("column") == "geometry_normalization"
            ),
            None,
        )
        if (
            norm_dict is None
            or norm_dict.get("physical_type") != "string"
            or norm_dict.get("null_meaning") != "not permitted"
            or {item["code"] for item in norm_dict.get("codes", [])} != {"original", "ring_self_intersection_split"}
        ):
            raise ValueError("Missing or invalid categorical_columns dictionary for geometry_normalization")

        if inherited_hillslopes:
            assert base_spatial_dir is not None and base_manifest is not None
            base_cat_parquet = pq.ParquetFile(base_spatial_dir / base_manifest["families"]["catchments"]["data"])
            base_all_native = [n for n in base_cat_parquet.schema_arrow.names if n not in cat_structural]
            if set(base_all_native) != set(cat_all_native):
                raise ValueError("Catchment native columns differ from base_spatial_dir")
            for pname in inherited_hillslopes:
                cat_native_names[(pname, "coastal_hillslopes")] = set(base_all_native)
            for batch in base_cat_parquet.iter_batches(batch_size=8192):
                kinds = batch.column("catchment_kind").to_pylist()
                pnames = batch.column("processing_partition").to_pylist()
                match_indices = [
                    idx
                    for idx, (kind, pname) in enumerate(zip(kinds, pnames, strict=True))
                    if kind == "coastal_hillslope" and pname in inherited_hillslopes
                ]
                if not match_indices:
                    continue
                matched = batch.take(match_indices)
                for row in matched.to_pylist():
                    ns, sid = row["source_namespace"], row["source_id"]
                    if (ns, sid) in omitted_keys:
                        raise ValueError(f"Unexpected inherited hillslope in omitted-features.json: {(ns, sid)}")
                    if (ns, sid) in repaired_records:
                        raise ValueError(f"Inherited hillslope unexpectedly marked repaired: {(ns, sid)}")
                    if row.get("geometry_normalization", "original") != "original":
                        raise ValueError(f"Unexpected geometry_normalization on inherited hillslope: {(ns, sid)}")
                    _opstat(str(row["processing_partition"]))["already_valid"] += 1
                    row_attrs = {k: row[k] for k in base_all_native}
                    db.execute(
                        "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 0, 0.0, 0)",
                        (
                            "catchments",
                            ns,
                            sid,
                            _canon_attrs_digest(row_attrs),
                            _canon_ring_signature(_wkb_rings(row["geometry"]), is_polygon=True),
                        ),
                    )
            db.commit()

        delivered_catchments = 0
        verified_repaired: list[tuple[int, str, str, Any]] = []
        for rg_idx in range(cat_parquet.num_row_groups):
            for row in cat_parquet.read_row_group(rg_idx).to_pylist():
                delivered_catchments += 1
                ns, sid, pname = row["source_namespace"], row["source_id"], row["processing_partition"]
                if (ns, sid) in omitted_keys:
                    raise ValueError(f"Omitted feature leaked into catchments.parquet: {(ns, sid)}")
                exp = db.execute(
                    "SELECT attr_digest, geom_sig, is_repaired, raw_area_km2, seen "
                    "FROM expected_features WHERE domain='catchments' AND source_namespace=? AND source_id=?",
                    (ns, sid),
                ).fetchone()
                if exp is None or exp[4] != 0:
                    raise ValueError(f"Unexpected or duplicate catchment row {(ns, sid)}")
                ns_domain = "coastal_hillslopes" if row["catchment_kind"] == "coastal_hillslope" else "catchments"
                present_cols = cat_native_names[(pname, ns_domain)]
                row_attrs = {k: row[k] for k in cat_all_native if k in present_cols}
                if _canon_attrs_digest(row_attrs) != exp[0]:
                    raise ValueError(f"Catchment native attributes differ: {(ns, sid)}")
                wkb = row["geometry"]
                geom = from_wkb(wkb)
                if (
                    geom.is_empty
                    or not geom.is_valid
                    or geom.geom_type not in {"Polygon", "MultiPolygon"}
                    or geom.has_z
                    or geom.has_m
                ):
                    raise ValueError(f"Invalid geometry in catchments.parquet: {(ns, sid)}")
                coords = get_coordinates(geom)
                if (
                    not len(coords)
                    or not np.all(np.isfinite(coords))
                    or np.any(coords[:, 0] < MIN_WGS84_LON)
                    or np.any(coords[:, 0] > MAX_UNWRAPPED_ANTIMERIDIAN_LON)
                    or np.any(np.abs(coords[:, 1]) > 90.0)
                ):
                    raise ValueError(f"Catchment coordinates out of bounds: {(ns, sid)}")

                is_repaired = int(exp[2]) == 1
                obs_norm = row.get("geometry_normalization")
                expected_norm = "ring_self_intersection_split" if is_repaired else "original"
                if obs_norm != expected_norm:
                    raise ValueError(
                        f"Unexpected geometry_normalization {obs_norm!r} != {expected_norm!r} for {(ns, sid)}"
                    )
                if is_repaired:
                    rec = repaired_records[(ns, sid)]
                    if rec.get("hfx_id") != row["id"]:
                        raise ValueError(f"repaired-features.jsonl.gz hfx_id mismatch for {(ns, sid)}")
                    rep_ev_digest, rep_zero = _canon_edge_vertex_digest(_wkb_rings(wkb))
                    if rep_zero > 0 or rep_ev_digest != exp[1]:
                        raise ValueError(f"Repaired catchment geometry edge/vertex mismatch: {(ns, sid)}")
                    raw_area = float(exp[3])
                    rep_area = _geodesic_area_km2(geom, geod)
                    if abs(rep_area - raw_area) > max(GEODESIC_AREA_ABS_TOL_KM2, GEODESIC_AREA_REL_TOL * raw_area):
                        raise ValueError(f"Repaired catchment geodesic area mismatch: {(ns, sid)}")
                    rep_parts, rep_holes = _part_and_hole_counts(geom)
                    if rep_parts != rec.get("repaired_part_count") or rep_holes != rec.get("repaired_hole_count"):
                        raise ValueError(f"Repaired catchment part/hole count mismatch: {(ns, sid)}")
                    expected_area = rep_area
                    verified_repaired.append((int(row["id"]), pname, sid, geom))
                else:
                    if _canon_ring_signature(_wkb_rings(wkb), is_polygon=True) != exp[1]:
                        raise ValueError(f"Catchment geometry vertices differ: {(ns, sid)}")
                    expected_area = abs(geod.geometry_area_perimeter(orient(geom, sign=1.0))[0]) / 1e6

                if abs(row["area_km2"] - expected_area) > max(1e-5, expected_area * 1e-6):
                    raise ValueError(f"Geodesic area differs for {(ns, sid)}")
                b = row["bbox"]
                xmin, ymin, xmax, ymax = geom.bounds
                if not (b["xmin"] <= xmin and b["ymin"] <= ymin and b["xmax"] >= xmax and b["ymax"] >= ymax):
                    raise ValueError(f"Catchment bbox does not cover geometry: {(ns, sid)}")
                if row["level"] != 0 or row["parent_id"] is not None or row["up_area_km2"] is not None:
                    raise ValueError("Invented hierarchy or upstream area")
                if row["outlet_lon"] is not None or row["outlet_lat"] is not None:
                    raise ValueError("Invented outlet coordinates")
                db.execute(
                    "UPDATE expected_features SET seen=1 WHERE domain='catchments' AND source_namespace=? AND source_id=?",
                    (ns, sid),
                )
                db.execute(
                    "INSERT INTO id_map VALUES ('catchments', ?, ?, ?, ?, ?)",
                    (ns, sid, row["id"], pname, row.get("COMID")),
                )
                db.execute(
                    "INSERT INTO observed_catchments VALUES (?, ?, ?, ?, ?, ?)",
                    (
                        row["id"],
                        row["catchment_kind"],
                        pname,
                        row.get("COMID"),
                        row["downstream_boundary_status"],
                        row["upstream_boundary_cuts"],
                    ),
                )

        riv_parquet = pq.ParquetFile(spatial_dir / manifest["families"]["rivers"]["data"])
        if "local_catchment_id" in riv_parquet.schema_arrow.names:
            raise ValueError("Scalar local_catchment_id must be omitted")
        riv_structural = {
            "id",
            "level",
            "parent_id",
            "source_namespace",
            "source_id",
            "processing_partition",
            "bbox",
            "geometry",
            "local_catchment_status",
            "downstream_boundary_status",
            "upstream_boundary_cuts",
            "component_catchment_id",
            "upstream_node_id",
            "downstream_node_id",
            "geometry_orientation",
        }
        riv_all_native = [n for n in riv_parquet.schema_arrow.names if n not in riv_structural]
        delivered_rivers = 0
        for rg_idx in range(riv_parquet.num_row_groups):
            for row in riv_parquet.read_row_group(rg_idx).to_pylist():
                delivered_rivers += 1
                ns, sid, pname = row["source_namespace"], row["source_id"], row["processing_partition"]
                if (ns, sid) in omitted_keys:
                    raise ValueError(f"Omitted river leaked into rivers.parquet: {(ns, sid)}")
                exp = db.execute(
                    "SELECT attr_digest, geom_sig, seen FROM expected_features WHERE domain='rivers' AND source_namespace=? AND source_id=?",
                    (ns, sid),
                ).fetchone()
                if exp is None or exp[2] != 0:
                    raise ValueError(f"Unexpected or duplicate river row {(ns, sid)}")
                row_attrs = {k: row[k] for k in riv_all_native}
                if _canon_attrs_digest(row_attrs) != exp[0]:
                    raise ValueError(f"River native attributes differ: {(ns, sid)}")
                wkb = row["geometry"]
                if _canon_ring_signature(_wkb_rings(wkb), is_polygon=False) != exp[1]:
                    raise ValueError(f"River geometry vertices differ: {(ns, sid)}")
                geom = from_wkb(wkb)
                b = row["bbox"]
                xmin, ymin, xmax, ymax = geom.bounds
                if not (b["xmin"] <= xmin and b["ymin"] <= ymin and b["xmax"] >= xmax and b["ymax"] >= ymax):
                    raise ValueError(f"River bbox does not cover geometry: {(ns, sid)}")
                comid = int(row["COMID"])
                if (pname, comid) in delivered_catchments_by_part_comid:
                    expected_loc = "supplied"
                elif (pname, comid) in omitted_catchments_by_part_comid:
                    expected_loc = "omitted_invalid_polygon"
                elif pname in unavailable_catchments:
                    expected_loc = "source_unavailable"
                else:
                    expected_loc = "not_in_source"
                if row["local_catchment_status"] != expected_loc:
                    raise ValueError(
                        f"River {comid}: local_catchment_status {row['local_catchment_status']} != {expected_loc}"
                    )
                db.execute(
                    "UPDATE expected_features SET seen=1 WHERE domain='rivers' AND source_namespace=? AND source_id=?",
                    (ns, sid),
                )
                db.execute(
                    "INSERT INTO id_map VALUES ('rivers', ?, ?, ?, ?, ?)",
                    (ns, sid, row["id"], pname, comid),
                )
        db.commit()

        unseen = db.execute("SELECT COUNT(*) FROM expected_features WHERE seen=0").fetchone()[0]
        if unseen:
            raise ValueError(f"{unseen} expected source features missing from final bundle")

        # 3. Verify complete river flow graph, acyclicity, cuts, derived catchment connectivity, and relationships IPC.
        riv_comid_to_id = dict(db.execute("SELECT comid, hfx_id FROM id_map WHERE domain='rivers'").fetchall())
        cat_part_comid_to_id = {
            (r[0], r[1]): r[2]
            for r in db.execute(
                "SELECT partition, comid, hfx_id FROM id_map WHERE domain='catchments' AND comid IS NOT NULL "
                "AND source_namespace NOT LIKE '%/coastal_hillslopes'"
            ).fetchall()
        }

        cat_part_by_comid: dict[int, str] = {}
        for part, comid in cat_part_comid_to_id:
            if comid in cat_part_by_comid:
                raise ValueError(f"Duplicate global channelized catchment COMID {comid}")
            cat_part_by_comid[comid] = part

        expected_edges: set[tuple[int, int]] = set()
        expected_downstream_cuts = 0
        expected_upstream_cuts = 0
        all_rivs = db.execute("SELECT comid, partition, next_down, up1, up2, up3, up4 FROM expected_rivers").fetchall()
        riv_part_map = {comid: pname for comid, pname, *_ in all_rivs}
        next_down_map = {comid: next_down for comid, _, next_down, *_ in all_rivs}
        ups_map = {comid: (up1, up2, up3, up4) for comid, _, _, up1, up2, up3, up4 in all_rivs}
        riv_incoming_by_comid: dict[int, list[int]] = {comid: [] for comid in riv_comid_to_id}
        riv_upstream_cuts_by_comid: dict[int, int] = {}

        for comid, pname, next_down, up1, up2, up3, up4 in all_rivs:
            if comid in cat_part_by_comid and cat_part_by_comid[comid] != pname:
                raise ValueError(f"Cross-partition COMID collision for river {pname}/{comid}")
            rid = riv_comid_to_id[comid]
            if next_down != 0:
                if next_down in riv_comid_to_id:
                    if comid not in ups_map[next_down] and not all(v != 0 for v in ups_map[next_down]):
                        raise ValueError(
                            f"NextDownID overflow target {next_down} lacks full up1..up4 slots for {comid}"
                        )
                    expected_edges.add((rid, riv_comid_to_id[next_down]))
                    riv_incoming_by_comid[next_down].append(comid)
                else:
                    expected_downstream_cuts += 1
            up_cuts = 0
            for u in (up1, up2, up3, up4):
                if u == 0:
                    continue
                if u not in riv_comid_to_id:
                    expected_upstream_cuts += 1
                    up_cuts += 1
                elif next_down_map[u] != comid:
                    raise ValueError(f"Upstream {u} disagrees with target {comid} NextDownID")
            riv_upstream_cuts_by_comid[comid] = up_cuts

        if len([c for c in cut_edges if c["direction"] == "downstream"]) != expected_downstream_cuts:
            raise ValueError("Cut downstream edge count differs")
        if len([c for c in cut_edges if c["direction"] == "upstream"]) != expected_upstream_cuts:
            raise ValueError("Cut upstream edge count differs")

        _check_acyclic_downstream(expected_edges, riv_comid_to_id.values(), "River")

        expected_local_pairs = {
            (riv_comid_to_id[comid], cat_part_comid_to_id[(pname, comid)])
            for comid, pname, *_ in all_rivs
            if (pname, comid) in cat_part_comid_to_id
        }

        actual_down_edges: set[tuple[int, int]] = set()
        actual_up_edges: set[tuple[int, int]] = set()
        actual_riv_local: set[tuple[int, int]] = set()
        with ipc.open_file(str(spatial_dir / manifest["relationships"]["rivers"]["data"])) as reader:
            for b_idx in range(reader.num_record_batches):
                for row in reader.get_batch(b_idx).to_pylist():
                    rid = row["id"]
                    for u in row["upstream_ids"]:
                        actual_up_edges.add((u, rid))
                    for d in row["downstream_ids"]:
                        actual_down_edges.add((rid, d))
                    for cid in row["local_catchment_ids"]:
                        actual_riv_local.add((rid, cid))
        if actual_down_edges != expected_edges or actual_up_edges != expected_edges:
            raise ValueError("Final river flow relationships differ from original DBFs")
        if actual_riv_local != expected_local_pairs:
            raise ValueError("Final river local associations differ from original DBFs")

        # Independently derive catchment connectivity and completeness expectations from original .shp/.dbf inputs.
        expected_cat_edges: set[tuple[int, int]] = set()
        for comid, pname, next_down, *_ in all_rivs:
            if next_down != 0 and next_down in riv_comid_to_id:
                src_cid = cat_part_comid_to_id.get((pname, comid))
                dst_cid = cat_part_comid_to_id.get((riv_part_map[next_down], next_down))
                if src_cid is not None and dst_cid is not None:
                    expected_cat_edges.add((src_cid, dst_cid))

        obs_cats = db.execute(
            "SELECT hfx_id, catchment_kind, partition, comid, downstream_boundary_status, upstream_boundary_cuts "
            "FROM observed_catchments"
        ).fetchall()
        _check_acyclic_downstream(expected_cat_edges, [r[0] for r in obs_cats], "Catchment")

        paired_count = 0
        unpaired_omitted_riv = 0
        unpaired_riv_not_in_source = 0
        unpaired_riv_unavailable = 0
        hillslope_count = 0
        cuts_cat_count = 0
        total_up_cut_edges = 0
        downstream_cuts = 0
        exp_status_counts: Counter[str] = Counter()
        exp_cat_up_cuts_by_id: dict[int, int | None] = {}

        for cid, ckind, pname, comid, obs_down_status, obs_up_cuts in obs_cats:
            if ckind == "coastal_hillslope":
                hillslope_count += 1
                exp_down_status = "unassociated_coastal_hillslope"
                exp_up_cuts: int | None = None
            elif comid is None:
                raise ValueError(f"Channelized catchment {cid} missing COMID")
            elif comid not in riv_comid_to_id:
                exp_up_cuts = None
                if comid in omitted_river_comids:
                    exp_down_status = "omitted_degenerate_river"
                    unpaired_omitted_riv += 1
                elif pname in unavailable_rivers:
                    exp_down_status = "river_source_unavailable"
                    unpaired_riv_unavailable += 1
                else:
                    exp_down_status = "river_not_in_source"
                    unpaired_riv_not_in_source += 1
            else:
                paired_count += 1
                exp_up_cuts = riv_upstream_cuts_by_comid[comid] + sum(
                    1 for u in riv_incoming_by_comid[comid] if (riv_part_map[u], u) not in cat_part_comid_to_id
                )
                if exp_up_cuts > 0:
                    cuts_cat_count += 1
                    total_up_cut_edges += exp_up_cuts
                nd = next_down_map[comid]
                if nd == 0:
                    exp_down_status = "source_zero_sentinel"
                elif nd in riv_comid_to_id and (riv_part_map[nd], nd) in cat_part_comid_to_id:
                    exp_down_status = "internal"
                else:
                    exp_down_status = "cut_missing_downstream_catchment"
                    downstream_cuts += 1

            if obs_down_status != exp_down_status:
                raise ValueError(
                    f"Catchment {cid} (COMID {comid}): downstream_boundary_status {obs_down_status!r} != {exp_down_status!r}"
                )
            if obs_up_cuts != exp_up_cuts:
                raise ValueError(
                    f"Catchment {cid} (COMID {comid}): upstream_boundary_cuts {obs_up_cuts!r} != {exp_up_cuts!r}"
                )
            exp_cat_up_cuts_by_id[cid] = exp_up_cuts
            exp_status_counts[exp_down_status] += 1

        cat_down_dict = next(
            (
                entry
                for entry in manifest.get("categorical_columns", [])
                if entry.get("artifact") == manifest["families"]["catchments"]["data"]
                and entry.get("column") == "downstream_boundary_status"
            ),
            None,
        )
        if (
            cat_down_dict is None
            or cat_down_dict.get("physical_type") != "string"
            or cat_down_dict.get("null_meaning") != "not permitted"
        ):
            raise ValueError("Missing or invalid categorical_columns dictionary for downstream_boundary_status")
        declared_codes = {item["code"] for item in cat_down_dict["codes"]}
        if not set(exp_status_counts).issubset(declared_codes):
            raise ValueError("Undeclared downstream_boundary_status code in catchments.parquet")

        actual_cat_up_edges: set[tuple[int, int]] = set()
        actual_cat_down_edges: set[tuple[int, int]] = set()
        actual_cat_local: set[tuple[int, int]] = set()
        with ipc.open_file(str(spatial_dir / manifest["relationships"]["catchments"]["data"])) as reader:
            if not {"upstream_ids", "downstream_ids", "local_river_ids"}.issubset(reader.schema.names):
                raise ValueError("Catchment relationships missing expected flow/local columns")
            for b_idx in range(reader.num_record_batches):
                for row in reader.get_batch(b_idx).to_pylist():
                    cid = row["id"]
                    for u in row["upstream_ids"]:
                        actual_cat_up_edges.add((u, cid))
                    for d in row["downstream_ids"]:
                        actual_cat_down_edges.add((cid, d))
                    for rid in row["local_river_ids"]:
                        actual_cat_local.add((rid, cid))
        if actual_cat_up_edges != expected_cat_edges or actual_cat_down_edges != expected_cat_edges:
            raise ValueError("Final catchment flow relationships differ from river-derived edges")
        if actual_cat_local != expected_local_pairs:
            raise ValueError("Final catchment reverse local associations differ from original DBFs")

        # Verify source-metadata.json catchment connectivity counts via independent closure check.
        riv_in_deg = {comid: len(riv_incoming_by_comid[comid]) for comid in riv_comid_to_id}
        anc_has_gap = {
            comid: ((riv_part_map[comid], comid) not in cat_part_comid_to_id) or (riv_upstream_cuts_by_comid[comid] > 0)
            for comid in riv_comid_to_id
        }
        anc_has_delivered = {comid: ((riv_part_map[comid], comid) in cat_part_comid_to_id) for comid in riv_comid_to_id}
        anc_has_stranded = dict.fromkeys(riv_comid_to_id, False)
        riv_queue = [comid for comid, deg in riv_in_deg.items() if deg == 0]
        while riv_queue:
            cur = riv_queue.pop()
            nxt = next_down_map[cur]
            if nxt not in riv_comid_to_id:
                continue
            if anc_has_gap[cur]:
                anc_has_gap[nxt] = True
            if anc_has_delivered[cur]:
                anc_has_delivered[nxt] = True
            if ((riv_part_map[nxt], nxt) not in cat_part_comid_to_id) and anc_has_delivered[cur]:
                anc_has_stranded[nxt] = True
            if anc_has_stranded[cur]:
                anc_has_stranded[nxt] = True
            riv_in_deg[nxt] -= 1
            if riv_in_deg[nxt] == 0:
                riv_queue.append(nxt)

        complete_count = 0
        incomplete_paired_count = 0
        stranded_count = 0
        unverified_count = 0
        for cid, _ckind, _pname, comid, _obs_ds, _obs_uc in obs_cats:
            if exp_cat_up_cuts_by_id[cid] is None:
                unverified_count += 1
                continue
            assert comid is not None
            if not anc_has_gap[comid]:
                complete_count += 1
            else:
                incomplete_paired_count += 1
                if anc_has_stranded[comid]:
                    stranded_count += 1

        rivers_omitted_poly = 0
        rivers_cat_unavail = 0
        rivers_cat_not_in_src = 0
        for comid, pname, *_ in all_rivs:
            if (pname, comid) in cat_part_comid_to_id:
                continue
            if (pname, comid) in omitted_catchments_by_part_comid:
                rivers_omitted_poly += 1
            elif pname in unavailable_catchments:
                rivers_cat_unavail += 1
            else:
                rivers_cat_not_in_src += 1

        expected_meta_counts = {
            "derived_flow_edges": len(expected_cat_edges),
            "paired_channelized_catchments": paired_count,
            "unpaired_channelized_catchments_omitted_river": unpaired_omitted_riv,
            "unpaired_channelized_catchments_river_not_in_source": unpaired_riv_not_in_source,
            "unpaired_channelized_catchments_river_unavailable": unpaired_riv_unavailable,
            "unassociated_coastal_hillslopes": hillslope_count,
            "rivers_omitted_invalid_polygon": rivers_omitted_poly,
            "rivers_catchment_not_in_source": rivers_cat_not_in_src,
            "rivers_catchment_unavailable": rivers_cat_unavail,
            "downstream_boundary_status_counts": dict(sorted(exp_status_counts.items())),
            "catchments_with_upstream_cuts": cuts_cat_count,
            "total_upstream_catchment_cut_edges": total_up_cut_edges,
            "catchments_with_downstream_cuts": downstream_cuts,
            "upstream_complete_catchments": complete_count,
            "upstream_incomplete_paired_catchments": incomplete_paired_count,
            "upstream_incomplete_leaf_only_catchments": incomplete_paired_count - stranded_count,
            "upstream_incomplete_stranded_delivered_catchments": stranded_count,
            "upstream_unverified_catchments": unverified_count,
        }
        source_meta = json.loads((spatial_dir / "source-metadata.json").read_bytes())
        if source_meta.get("catchment_connectivity", {}).get("counts") != expected_meta_counts:
            raise ValueError("source-metadata.json catchment_connectivity.counts differs from independent verification")

        flow_edges_count = len(expected_edges)
        catchment_flow_edges_count = len(expected_cat_edges)
        local_associations_count = len(expected_local_pairs)

        # Release Step 3 river/graph memory before allocating global Step 4 spatial indexes.
        del (
            all_rivs,
            expected_edges,
            actual_down_edges,
            actual_up_edges,
            actual_riv_local,
            expected_local_pairs,
            expected_cat_edges,
            obs_cats,
            riv_incoming_by_comid,
            riv_part_map,
            next_down_map,
            ups_map,
            riv_comid_to_id,
            cat_part_comid_to_id,
            anc_has_gap,
            anc_has_delivered,
            anc_has_stranded,
            delivered_catchments_by_part_comid,
            omitted_catchments_by_part_comid,
        )

        # 4. Independent global sibling interior non-overlap verification for all repaired catchments
        #    and reconciliation of source-metadata.json["geometry_normalization"] counts.
        repaired_overlap_pairs_checked = 0
        if verified_repaired:
            rep_geoms = [item[3] for item in verified_repaired]
            for g in rep_geoms:
                prepare(g)
            rep_bounds = [g.bounds for g in rep_geoms]
            has_antimeridian = any(b[0] < ANTIMERIDIAN_QUERY_MIN_X or b[2] > ANTIMERIDIAN_MAX_X for b in rep_bounds)
            rep_tree = STRtree(rep_geoms)

            for i, (hid_i, _pname_i, _sid_i, g_i) in enumerate(verified_repaired):
                b_i = rep_bounds[i]
                queries: list[tuple[float, Any]] = [(0.0, g_i)]
                if has_antimeridian:
                    if b_i[2] > ANTIMERIDIAN_MAX_X:
                        queries.append((-360.0, _shift_x(g_i, -360.0)))
                    if b_i[0] < ANTIMERIDIAN_QUERY_MIN_X:
                        queries.append((360.0, _shift_x(g_i, 360.0)))
                for shift_dx, q_geom in queries:
                    for j in rep_tree.query(q_geom, predicate="intersects"):
                        j_idx = int(j)
                        if j_idx <= i:
                            continue
                        g_j = _shift_x(rep_geoms[j_idx], -shift_dx) if shift_dx else rep_geoms[j_idx]
                        if shift_dx != 0.0 and not g_i.intersects(g_j):
                            continue
                        repaired_overlap_pairs_checked += 1
                        if g_i.relate_pattern(g_j, "2********"):
                            raise ValueError(
                                f"Repaired catchments have sibling interior overlap: {hid_i} vs {verified_repaired[j_idx][0]}"
                            )

            for batch in cat_parquet.iter_batches(
                batch_size=8192,
                columns=["id", "geometry_normalization", "bbox", "geometry"],
            ):
                norms = cast(list[str], batch.column("geometry_normalization").to_pylist())
                orig_mask = np.array([v == "original" for v in norms], dtype=bool)
                if not np.any(orig_mask):
                    continue
                orig_indices = np.flatnonzero(orig_mask)
                bounds_arr = cast(pa.StructArray, batch.column("bbox"))
                xmins = bounds_arr.field("xmin").to_numpy().astype(np.float64)[orig_indices]
                ymins = bounds_arr.field("ymin").to_numpy().astype(np.float64)[orig_indices]
                xmaxs = bounds_arr.field("xmax").to_numpy().astype(np.float64)[orig_indices]
                ymaxs = bounds_arr.field("ymax").to_numpy().astype(np.float64)[orig_indices]
                ids_col = cast(list[int], batch.column("id").to_pylist())
                wkbs = batch.column("geometry")

                hit_batches: list[tuple[float, np.ndarray, np.ndarray]] = []
                r0, c0 = rep_tree.query(shapely_box(xmins, ymins, xmaxs, ymaxs))
                if len(r0):
                    hit_batches.append((0.0, orig_indices[r0], c0))
                if has_antimeridian:
                    west_sub = np.flatnonzero(xmins < ANTIMERIDIAN_QUERY_MIN_X)
                    if len(west_sub):
                        wr, wc = rep_tree.query(
                            shapely_box(
                                xmins[west_sub] + 360.0, ymins[west_sub], xmaxs[west_sub] + 360.0, ymaxs[west_sub]
                            )
                        )
                        if len(wr):
                            hit_batches.append((360.0, orig_indices[west_sub[wr]], wc))
                    east_sub = np.flatnonzero(xmaxs > ANTIMERIDIAN_MAX_X)
                    if len(east_sub):
                        er, ec = rep_tree.query(
                            shapely_box(
                                xmins[east_sub] - 360.0, ymins[east_sub], xmaxs[east_sub] - 360.0, ymaxs[east_sub]
                            )
                        )
                        if len(er):
                            hit_batches.append((-360.0, orig_indices[east_sub[er]], ec))

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
                        g_rep = rep_geoms[cand_idx]
                        if not g_rep.intersects(v_geom):
                            continue
                        repaired_overlap_pairs_checked += 1
                        if g_rep.relate_pattern(v_geom, "2********"):
                            raise ValueError(
                                f"Repaired catchment {verified_repaired[cand_idx][0]} has sibling interior overlap "
                                f"with original catchment {ids_col[row_pos]}"
                            )

        if len(verified_repaired) != len(repaired_records):
            raise ValueError("repaired-features.jsonl.gz contains unmatched records absent from catchments.parquet")

        gnorm_meta = source_meta.get("geometry_normalization")
        if not isinstance(gnorm_meta, dict):
            raise ValueError("Missing source-metadata.json geometry_normalization section")
        gcounts = gnorm_meta.get("global_counts", {})
        exp_accepted = sum(s["accepted"] for s in obs_part_counts.values())
        exp_valid = sum(s["already_valid"] for s in obs_part_counts.values())
        exp_rejected = sum(s["rejected"] for s in obs_part_counts.values())
        exp_other = sum(s["omitted_other"] for s in obs_part_counts.values())
        if (
            gcounts.get("accepted") != exp_accepted
            or gcounts.get("already_valid") != exp_valid
            or gcounts.get("rejected") != exp_rejected
            or gcounts.get("omitted_other") != exp_other
            or gcounts.get("attempted") != exp_accepted + exp_rejected
            or gcounts.get("total_source") != exp_valid + exp_accepted + exp_rejected + exp_other
            or gnorm_meta.get("overlap_validation", {}).get("accepted_interior_overlaps") != 0
        ):
            raise ValueError("source-metadata.json geometry_normalization counts differ from independent verification")
        meta_parts = gnorm_meta.get("by_partition", {})
        for pname, counts in obs_part_counts.items():
            mp = meta_parts.get(pname)
            if (
                mp is None
                or any(mp.get(k) != v for k, v in counts.items())
                or mp.get("attempted") != counts["accepted"] + counts["rejected"]
                or mp.get("total_source")
                != counts["already_valid"] + counts["accepted"] + counts["rejected"] + counts["omitted_other"]
            ):
                raise ValueError(f"source-metadata.json geometry_normalization partition {pname} counts differ")

        if (
            base_spatial_dir is not None
            and "derived_from_verified_bundle_manifest_sha256" in source_meta
            and source_meta["derived_from_verified_bundle_manifest_sha256"]
            != _sha256(base_spatial_dir / "manifest.json")
        ):
            raise ValueError("source-metadata.json derived_from_verified_bundle_manifest_sha256 differs from base")

    for path_str, expected_digest in source_hashes.items():
        if _sha256(Path(path_str)) != expected_digest:
            raise ValueError(f"Source file changed during verification: {path_str}")

    return {
        "status": "passed",
        "delivered_catchments": delivered_catchments,
        "repaired_catchments": len(verified_repaired),
        "repaired_overlap_pairs_checked": repaired_overlap_pairs_checked,
        "delivered_rivers": delivered_rivers,
        "omitted_catchments": len([x for x in omitted_list if x["domain"] != "rivers"]),
        "omitted_rivers": len([x for x in omitted_list if x["domain"] == "rivers"]),
        "flow_edges": flow_edges_count,
        "catchment_flow_edges": catchment_flow_edges_count,
        "local_associations": local_associations_count,
        "cut_edges": len(cut_edges),
        "raw_shapefile_verified_hillslope_partitions": sorted(raw_verified_hillslopes),
        "base_bundle_verified_hillslope_partitions": sorted(inherited_hillslopes),
        "manifest_sha256": _sha256(spatial_dir / "manifest.json"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("spatial_dir", type=Path)
    parser.add_argument("--base-spatial-dir", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    result = verify_spatial_bundle(args.inventory, args.spatial_dir, base_spatial_dir=args.base_spatial_dir)
    args.report.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
