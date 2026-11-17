"""Independent binary shapefile (.shp/.dbf) verification against final MERIT HFX spatial bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
import struct
from collections.abc import Iterator
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


def _canon_ring_signature(rings: list[tuple[tuple[float, float], ...]], *, is_polygon: bool) -> str:
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
                "CREATE TABLE expected_rivers ("
                "comid INTEGER PRIMARY KEY, partition TEXT NOT NULL, next_down INTEGER NOT NULL, "
                "up1 INTEGER NOT NULL, up2 INTEGER NOT NULL, up3 INTEGER NOT NULL, up4 INTEGER NOT NULL)"
            )
            yield conn
        finally:
            conn.close()


def verify_spatial_bundle(inventory_path: Path, spatial_dir: Path) -> dict[str, Any]:
    """Independently verify original shapefiles against the final GeoParquet/IPC spatial bundle."""
    manifest = json.loads((spatial_dir / "manifest.json").read_bytes())
    if manifest["representation"] != "spatial-geoparquet-adjacency" or manifest["representation_version"] != 1:
        raise ValueError("Unexpected representation")
    if manifest["catchment_connectivity"] != [{"level": 0, "status": "not_supplied"}] or manifest["has_up_area"]:
        raise ValueError("Unexpected catchment connectivity or has_up_area")
    if (spatial_dir / "original-inventory.json").read_bytes() != inventory_path.read_bytes():
        raise ValueError("original-inventory.json differs from inventory")

    for rel, identity in manifest["artifacts"].items():
        path = spatial_dir / rel
        if not path.resolve().is_relative_to(spatial_dir.resolve()) or fingerprint(path) != identity:
            raise ValueError(f"Final artifact fingerprint differs: {rel}")

    inventory = json.loads(inventory_path.read_bytes())
    release = str(inventory["release"])
    partitions = sorted(inventory["partitions"], key=lambda item: item["partition"])
    hillslopes = sorted(inventory.get("coastal_hillslopes", []), key=lambda item: item["partition"])
    omitted_list = json.loads((spatial_dir / "omitted-features.json").read_bytes())
    omitted_keys = {(item["source_namespace"], item["source_id"]): item for item in omitted_list}
    cut_edges = json.loads((spatial_dir / "cut-edges.json").read_bytes())

    source_hashes: dict[str, str] = {}
    geod = Geod(ellps="WGS84")

    with _verify_db() as db:
        # 1. Verify original sidecar hashes and record expected delivered/omitted source features directly from .shp/.dbf.
        cat_entries: list[tuple[str, str, dict[str, Any]]] = []
        riv_entries: list[tuple[str, dict[str, Any]]] = []
        unavailable_catchments: set[str] = set()
        for part in partitions:
            pname = str(part["partition"])
            if part["catchments"].get("status") == "supplied":
                cat_entries.append((pname, "catchments", part["catchments"]))
            else:
                unavailable_catchments.add(pname)
            if part["rivers"].get("status") == "supplied":
                riv_entries.append((pname, part["rivers"]))
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
        }
        cat_all_native = [n for n in cat_parquet.schema_arrow.names if n not in cat_structural]

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

        # 3. Verify complete river flow graph, acyclicity, cuts, and relationships IPC against expected_rivers.
        riv_comid_to_id = dict(db.execute("SELECT comid, hfx_id FROM id_map WHERE domain='rivers'").fetchall())
        cat_part_comid_to_id = {
            (r[0], r[1]): r[2]
            for r in db.execute(
                "SELECT partition, comid, hfx_id FROM id_map WHERE domain='catchments' AND comid IS NOT NULL "
                "AND source_namespace NOT LIKE '%/coastal_hillslopes'"
            ).fetchall()
        }

        expected_edges: set[tuple[int, int]] = set()
        expected_downstream_cuts = 0
        expected_upstream_cuts = 0
        all_rivs = db.execute("SELECT comid, partition, next_down, up1, up2, up3, up4 FROM expected_rivers").fetchall()
        next_down_map = {comid: next_down for comid, _, next_down, *_ in all_rivs}
        ups_map = {comid: (up1, up2, up3, up4) for comid, _, _, up1, up2, up3, up4 in all_rivs}
        for comid, _, next_down, up1, up2, up3, up4 in all_rivs:
            rid = riv_comid_to_id[comid]
            if next_down != 0:
                if next_down in riv_comid_to_id:
                    if comid not in ups_map[next_down] and not all(v != 0 for v in ups_map[next_down]):
                        raise ValueError(
                            f"NextDownID overflow target {next_down} lacks full up1..up4 slots for {comid}"
                        )
                    expected_edges.add((rid, riv_comid_to_id[next_down]))
                else:
                    expected_downstream_cuts += 1
            for u in (up1, up2, up3, up4):
                if u == 0:
                    continue
                if u not in riv_comid_to_id:
                    expected_upstream_cuts += 1
                elif next_down_map[u] != comid or (riv_comid_to_id[u], rid) not in expected_edges:
                    if u in next_down_map and next_down_map[u] != comid:
                        raise ValueError(f"Upstream {u} disagrees with target {comid} NextDownID")
                    if (riv_comid_to_id[u], rid) not in expected_edges:
                        expected_edges.add((riv_comid_to_id[u], rid))

        if len([c for c in cut_edges if c["direction"] == "downstream"]) != expected_downstream_cuts:
            raise ValueError("Cut downstream edge count differs")
        if len([c for c in cut_edges if c["direction"] == "upstream"]) != expected_upstream_cuts:
            raise ValueError("Cut upstream edge count differs")

        # Independent cycle check via downstream path traversal.
        successor = dict(expected_edges)
        completed: set[int] = set()
        for rid in riv_comid_to_id.values():
            path_ids: set[int] = set()
            cur = rid
            while cur in successor and cur not in completed:
                if cur in path_ids:
                    raise ValueError("River graph contains a cycle")
                path_ids.add(cur)
                cur = successor[cur]
            completed.update(path_ids)

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

        actual_cat_local: set[tuple[int, int]] = set()
        with ipc.open_file(str(spatial_dir / manifest["relationships"]["catchments"]["data"])) as reader:
            for b_idx in range(reader.num_record_batches):
                for row in reader.get_batch(b_idx).to_pylist():
                    cid = row["id"]
                    if row.get("upstream_ids") or row.get("downstream_ids"):
                        raise ValueError("Catchment flow graph must be empty / not_supplied")
                    for rid in row["local_river_ids"]:
                        actual_cat_local.add((rid, cid))
        if actual_cat_local != expected_local_pairs:
            raise ValueError("Final catchment reverse local associations differ from original DBFs")

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
        "local_associations": len(expected_local_pairs),
        "cut_edges": len(cut_edges),
        "manifest_sha256": _sha256(spatial_dir / "manifest.json"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("spatial_dir", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    result = verify_spatial_bundle(args.inventory, args.spatial_dir)
    args.report.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
