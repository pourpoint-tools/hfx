"""Independent binary shapefile (.shp/.dbf) verification against final MERIT HFX spatial bundle."""

from __future__ import annotations

import argparse
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
from typing import Any

import pyarrow.ipc as ipc
import pyarrow.parquet as pq
from pyproj import Geod
from shapely import from_wkb
from shapely.ops import orient

from hfx.candidates import fingerprint

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
def _verify_db() -> Iterator[sqlite3.Connection]:
    with TemporaryDirectory(prefix="hfx-merit-verify-") as scratch:
        conn = sqlite3.connect(Path(scratch) / "verify.sqlite")
        try:
            conn.execute("PRAGMA cache_size = -32768")
            conn.execute("PRAGMA temp_store = FILE")
            conn.execute("PRAGMA synchronous = OFF")
            conn.execute("PRAGMA journal_mode = MEMORY")
            conn.execute(
                "CREATE TABLE expected_features ("
                "domain TEXT NOT NULL, source_namespace TEXT NOT NULL, source_id TEXT NOT NULL, "
                "attr_digest TEXT NOT NULL, geom_sig TEXT NOT NULL, seen INTEGER NOT NULL DEFAULT 0, "
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
    if (spatial_dir / "original-inventory.json").read_bytes() != base_orig_path.read_bytes():
        raise ValueError("original-inventory.json differs from base_spatial_dir")

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

    inventory = json.loads(inventory_path.read_bytes())
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
            base_spatial_dir,
            spatial_dir,
            inherited_hillslopes=inherited_hillslopes,
            supplied_hillslopes=set(raw_verified_hillslopes),
        )
    elif (spatial_dir / "original-inventory.json").read_bytes() != inventory_path.read_bytes():
        raise ValueError("original-inventory.json differs from inventory")

    for rel, identity in manifest["artifacts"].items():
        path = spatial_dir / rel
        if not path.resolve().is_relative_to(spatial_dir.resolve()) or fingerprint(path) != identity:
            raise ValueError(f"Final artifact fingerprint differs: {rel}")

    omitted_list = json.loads((spatial_dir / "omitted-features.json").read_bytes())
    omitted_keys = {(item["source_namespace"], item["source_id"]): item for item in omitted_list}
    omitted_river_comids = {
        int(item["COMID"]) for item in omitted_list if item["domain"] == "rivers" and item.get("COMID") is not None
    }
    cut_edges = json.loads((spatial_dir / "cut-edges.json").read_bytes())

    source_hashes: dict[str, str] = {}
    geod = Geod(ellps="WGS84")

    with _verify_db() as db:
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
            for attrs, (_, rings) in zip(dbf_rows, shp_rows, strict=True):
                comid = int(attrs["COMID"]) if "COMID" in attrs and attrs["COMID"] is not None else None
                raw_sid = str(comid) if comid is not None else str(attrs.get("FID"))
                if (ns, raw_sid) in omitted_keys:
                    if ns_domain == "catchments" and comid is not None:
                        omitted_catchments_by_part_comid.add((pname, comid))
                    continue
                if ns_domain == "catchments" and comid is not None:
                    delivered_catchments_by_part_comid.add((pname, comid))
                renamed = {("native_" + k if k in STRUCTURAL_FIELDS else k): v for k, v in attrs.items()}
                db.execute(
                    "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 0)",
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
                    "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 0)",
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
                    row_attrs = {k: row[k] for k in base_all_native}
                    db.execute(
                        "INSERT INTO expected_features VALUES (?, ?, ?, ?, ?, 0)",
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
        for rg_idx in range(cat_parquet.num_row_groups):
            for row in cat_parquet.read_row_group(rg_idx).to_pylist():
                delivered_catchments += 1
                ns, sid, pname = row["source_namespace"], row["source_id"], row["processing_partition"]
                if (ns, sid) in omitted_keys:
                    raise ValueError(f"Omitted feature leaked into catchments.parquet: {(ns, sid)}")
                exp = db.execute(
                    "SELECT attr_digest, geom_sig, seen FROM expected_features WHERE domain='catchments' AND source_namespace=? AND source_id=?",
                    (ns, sid),
                ).fetchone()
                if exp is None or exp[2] != 0:
                    raise ValueError(f"Unexpected or duplicate catchment row {(ns, sid)}")
                ns_domain = "coastal_hillslopes" if row["catchment_kind"] == "coastal_hillslope" else "catchments"
                present_cols = cat_native_names[(pname, ns_domain)]
                row_attrs = {k: row[k] for k in cat_all_native if k in present_cols}
                if _canon_attrs_digest(row_attrs) != exp[0]:
                    raise ValueError(f"Catchment native attributes differ: {(ns, sid)}")
                wkb = row["geometry"]
                if _canon_ring_signature(_wkb_rings(wkb), is_polygon=True) != exp[1]:
                    raise ValueError(f"Catchment geometry vertices differ: {(ns, sid)}")
                geom = from_wkb(wkb)
                if not geom.is_valid or geom.is_empty:
                    raise ValueError(f"Invalid geometry in catchments.parquet: {(ns, sid)}")
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
        "delivered_rivers": delivered_rivers,
        "omitted_catchments": len([x for x in omitted_list if x["domain"] != "rivers"]),
        "omitted_rivers": len([x for x in omitted_list if x["domain"] == "rivers"]),
        "flow_edges": len(expected_edges),
        "catchment_flow_edges": len(expected_cat_edges),
        "local_associations": len(expected_local_pairs),
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
