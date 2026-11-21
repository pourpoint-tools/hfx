"""Independently verify normalized TDX staging and derived catchment-connectivity bundles."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
import struct
from collections import Counter, defaultdict
from collections.abc import Generator, Iterator
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq
from pyproj import Geod
from shapely import MultiPolygon, from_wkb, get_coordinates, to_wkb

from hfx.candidates import fingerprint

CATCHMENT_DERIVATION = "river_local_associations_with_zero_area_virtual_stub_contraction"
SOURCE_FIDS_TYPE = pa.list_(pa.field("item", pa.int64(), nullable=False))


def _rows(path: Path) -> Generator[dict[str, Any]]:
    with ipc.open_file(str(path)) as reader:
        for index in range(reader.num_record_batches):
            yield from reader.get_batch(index).to_pylist()


def _digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _geodesic_area_km2(geod: Geod, shape: Any) -> float:
    parts = shape.geoms if shape.geom_type == "MultiPolygon" else [shape]
    area = 0.0
    for polygon in parts:
        x, y = polygon.exterior.xy
        part_area = abs(geod.polygon_area_perimeter(x, y)[0])
        for ring in polygon.interiors:
            x, y = ring.xy
            part_area -= abs(geod.polygon_area_perimeter(x, y)[0])
        area += part_area / 1e6
    return area


def _independent_is_virtual_stub(
    linkno: int,
    length: float | None,
    geom_length: float | None,
    ds_area: float | None,
    us_area: float | None,
    up_deg: int,
    down_deg: int,
) -> bool:
    if linkno == 0 or up_deg < 1 or down_deg != 1:
        return False
    for val in (length, geom_length, ds_area, us_area):
        if val is None or isinstance(val, bool) or not isinstance(val, (int, float)) or not math.isfinite(val):
            return False
    return length == 0.0 and geom_length == 0.0 and ds_area == us_area


def _independent_derive_catchment_graph(
    river_ids: list[int],
    river_upstream: dict[int, list[int]],
    river_downstream: dict[int, list[int]],
    river_to_cat: dict[int, int],
    polygonless_attrs: dict[int, tuple[int, float | None, float | None, float | None, float | None]],
) -> tuple[dict[int, list[int]], dict[int, int], dict[str, int]]:
    """Independent recursive/memoized upstream walk (not Kahn's algorithm) to verify contraction."""
    contracted: set[int] = set()
    linkno0 = 0
    other = 0
    for rid in river_ids:
        if rid in river_to_cat:
            continue
        linkno, length, geom_len, ds_area, us_area = polygonless_attrs[rid]
        if _independent_is_virtual_stub(
            linkno, length, geom_len, ds_area, us_area, len(river_upstream[rid]), len(river_downstream[rid])
        ):
            contracted.add(rid)
        elif linkno == 0:
            linkno0 += 1
        else:
            other += 1

    # Iterative upstream resolution through contracted stubs only (no recursion limit).
    memo: dict[int, tuple[list[int], int]] = {}

    def resolve_incoming(start_rid: int) -> tuple[list[int], int]:
        if start_rid in memo:
            return memo[start_rid]
        stack: list[tuple[int, bool]] = [(start_rid, False)]
        while stack:
            rid, expanded = stack.pop()
            if rid in memo:
                continue
            if not expanded:
                stack.append((rid, True))
                for up in river_upstream[rid]:
                    if up in contracted and up not in memo:
                        stack.append((up, False))
            else:
                cats: list[int] = []
                cuts = 0
                for up in river_upstream[rid]:
                    if up in river_to_cat:
                        cats.append(river_to_cat[up])
                    elif up in contracted:
                        sub_cats, sub_cuts = memo[up]
                        cats.extend(sub_cats)
                        cuts += sub_cuts
                    else:
                        cuts += 1
                memo[rid] = (sorted(cats), cuts)
        return memo[start_rid]

    cat_upstream: dict[int, list[int]] = {}
    cat_cuts: dict[int, int] = {}
    for rid, cid in river_to_cat.items():
        ups, cuts = resolve_incoming(rid)
        if cid in ups or len(set(ups)) != len(ups):
            raise ValueError(f"Duplicate or self-linked catchment upstream at {cid}")
        cat_upstream[cid] = ups
        cat_cuts[cid] = cuts

    # Independent downstream-chain cycle check on catchment graph.
    cat_successor: dict[int, int] = {}
    for target, origins in cat_upstream.items():
        for origin in origins:
            if origin in cat_successor:
                raise ValueError(f"Non-tree catchment downstream split at {origin}")
            cat_successor[origin] = target
    completed: set[int] = set()
    for cid in cat_upstream:
        seen: set[int] = set()
        curr = cid
        while curr in cat_successor and curr not in completed:
            if curr in seen:
                raise ValueError("Derived catchment graph contains a cycle")
            seen.add(curr)
            curr = cat_successor[curr]
        completed.update(seen)

    # Independent iterative downstream propagation from cut catchments (immune to deep chains).
    affected: set[int] = set()
    for cid, cuts in cat_cuts.items():
        if cuts > 0:
            curr_cid: int | None = cid
            while curr_cid is not None and curr_cid not in affected:
                affected.add(curr_cid)
                curr_cid = cat_successor.get(curr_cid)
    affected_count = len(affected)
    stats = {
        "contracted_stub_count": len(contracted),
        "non_contracted_linkno0_count": linkno0,
        "non_contracted_other_count": other,
        "catchments_with_cuts_count": sum(1 for c in cat_cuts.values() if c > 0),
        "affected_downstream_catchments_count": affected_count,
    }
    return cat_upstream, cat_cuts, stats


@contextmanager
def _actual_relations(output: Path) -> Iterator[sqlite3.Connection]:
    """Disk-backed actual evidence preserves arbitrary physical relation order."""
    with TemporaryDirectory(prefix="hfx-verify-") as scratch:
        connection = sqlite3.connect(Path(scratch) / "actual.sqlite")
        try:
            connection.execute("PRAGMA cache_size = -8192")
            connection.execute("PRAGMA temp_store = FILE")
            connection.execute(
                "CREATE TABLE graph (id INTEGER PRIMARY KEY, level INTEGER NOT NULL, xmin REAL NOT NULL, ymin REAL NOT NULL, xmax REAL NOT NULL, ymax REAL NOT NULL)"
            )
            connection.execute(
                "CREATE TABLE edge (origin INTEGER NOT NULL, target INTEGER NOT NULL, PRIMARY KEY(origin, target))"
            )
            connection.execute("CREATE INDEX edge_target ON edge(target)")
            connection.execute(
                "CREATE TABLE cat_graph (id INTEGER PRIMARY KEY, level INTEGER NOT NULL, xmin REAL NOT NULL, ymin REAL NOT NULL, xmax REAL NOT NULL, ymax REAL NOT NULL)"
            )
            connection.execute(
                "CREATE TABLE cat_edge (origin INTEGER NOT NULL, target INTEGER NOT NULL, PRIMARY KEY(origin, target))"
            )
            connection.execute("CREATE INDEX cat_edge_target ON cat_edge(target)")
            connection.execute(
                "CREATE TABLE association (river INTEGER NOT NULL, catchment INTEGER NOT NULL, PRIMARY KEY(river, catchment))"
            )
            connection.execute("CREATE INDEX association_catchment ON association(catchment)")
            for table_prefix, filename in (("", "river-graph.parquet"), ("cat_", "graph.parquet")):
                for batch in pq.ParquetFile(output / filename).iter_batches(batch_size=4096):
                    with connection:
                        for row in batch.to_pylist():
                            identity, level, upstream = row["id"], row["level"], row["upstream_ids"]
                            bounds = [row[name] for name in ("bbox_minx", "bbox_miny", "bbox_maxx", "bbox_maxy")]
                            if (
                                type(identity) is not int
                                or type(level) is not int
                                or not isinstance(upstream, list)
                                or any(type(value) is not int for value in upstream)
                            ):
                                raise ValueError("Invalid graph identity or adjacency")
                            if any(type(value) is not float or not np.isfinite(value) for value in bounds):
                                raise ValueError("Invalid graph bounds")
                            connection.execute(
                                f"INSERT INTO {table_prefix}graph VALUES (?, ?, ?, ?, ?, ?)",
                                (identity, level, *bounds),
                            )
                            connection.executemany(
                                f"INSERT INTO {table_prefix}edge VALUES (?, ?)",
                                ((origin, identity) for origin in upstream),
                            )
            for batch in pq.ParquetFile(output / "local-associations.parquet").iter_batches(batch_size=4096):
                with connection:
                    for row in batch.to_pylist():
                        pair = row["river_id"], row["catchment_id"]
                        if any(type(value) is not int for value in pair):
                            raise ValueError("Invalid local association")
                        connection.execute("INSERT INTO association VALUES (?, ?)", pair)
            yield connection
        except sqlite3.IntegrityError as error:
            raise ValueError("Duplicate or invalid normalized relation") from error
        finally:
            connection.close()


def _verify_basin(
    item: dict[str, Any],
    schemas: dict[str, pa.Schema],
    iterators: dict[str, Generator[dict[str, Any]]],
    counts: dict[str, int],
    source_hashes: dict[str, str],
    relations: sqlite3.Connection,
) -> tuple[int, int, int, int]:
    """Independently prove one basin, releasing its topology before the next."""
    rivers: dict[tuple[str, int], int] = {}
    polygons: list[tuple[str, int, int]] = []
    downstream: set[tuple[str, int, int]] = set()
    upstream: set[tuple[str, int, int]] = set()
    river_bounds: dict[int, dict[str, float]] = {}
    river_attrs: dict[int, tuple[int, float | None, float | None, float | None, float | None]] = {}
    cat_bounds: dict[int, dict[str, float]] = {}
    actual_cat_cuts: dict[int, int | None] = {}
    degenerate = 0
    unavailable = item["catchments"].get("status") == "unavailable"
    geod = Geod(ellps="WGS84")
    basin = item["basin"]

    # Verify catchments (grouped by streamID in ascending min(fid) order).
    if not unavailable:
        schema = schemas["catchments"]
        actual_rows = iterators["catchments"]
        count = counts["catchments"]
        source = Path(item["catchments"]["path"])
        digest = _digest(source)
        if digest != item["catchments"]["sha256"]:
            raise ValueError("Original source digest differs")
        source_hashes[str(source)] = digest
        connection = sqlite3.connect(source.resolve().as_uri() + "?mode=ro", uri=True)
        try:
            (layer,) = connection.execute("SELECT table_name FROM gpkg_contents WHERE data_type='features'").fetchone()
            quoted = '"' + layer.replace('"', '""') + '"'
            fields = connection.execute(f"PRAGMA table_info({quoted})").fetchall()
            geometry, srs, z, m = connection.execute(
                "SELECT column_name, srs_id, z, m FROM gpkg_geometry_columns WHERE table_name=?", (layer,)
            ).fetchone()
            if (srs, z, m) != (4326, 0, 0):
                raise ValueError("Original geometry declaration differs")
            expected_types = {"INTEGER": pa.int64(), "MEDIUMINT": pa.int32(), "REAL": pa.float64()}
            for _, name, declared, nonnull, _, _ in fields:
                if name == geometry:
                    continue
                field = schema.field(name)
                mapped = name == "streamID"
                if mapped and declared not in ("INTEGER", "MEDIUMINT"):
                    raise ValueError("Unsupported original streamID declaration")
                expected_type = pa.int64() if mapped else expected_types[declared]
                metadata = (
                    {b"hfx:source_mapping": b"streamID INTEGER/MEDIUMINT to int64, lossless"}
                    if mapped
                    else {b"gpkg:declared_type": declared.encode()}
                )
                if field.type != expected_type or field.nullable != (not nonnull) or field.metadata != metadata:
                    raise ValueError(f"Native schema differs: {name}")

            by_stream: dict[int, list[tuple[int, bytes]]] = defaultdict(list)
            for values in connection.execute(f"SELECT * FROM {quoted} ORDER BY fid"):
                original = dict(zip((field[1] for field in fields), values, strict=True))
                blob = original[geometry]
                if blob[:3] != b"GP\x00" or blob[3] & 0xF0:
                    raise ValueError("Unsupported original geometry header")
                envelope = (blob[3] >> 1) & 7
                if envelope not in (0, 1):
                    raise ValueError("Unsupported source dimensions")
                wkb = blob[8 + (32 if envelope == 1 else 0) :]
                by_stream[int(original["streamID"])].append((int(original["fid"]), wkb))
        finally:
            connection.close()

        for stream_id, fid_wkbs in sorted(by_stream.items(), key=lambda kv: min(f for f, _ in kv[1])):
            actual = next(actual_rows, None)
            if actual is None:
                raise ValueError("Missing source feature")
            count += 1
            ordered = sorted(fid_wkbs, key=lambda pair: pair[0])
            expected_fids = [fid for fid, _ in ordered]
            if (
                actual["id"] != count
                or actual["source_id"] != str(expected_fids[0])
                or actual["fid"] != expected_fids[0]
                or actual["source_fids"] != expected_fids
                or actual["streamID"] != stream_id
                or actual["processing_basin"] != basin
                or actual["source_namespace"] != f"NGA-TDX-Hydro-20230126/{basin}/catchments"
            ):
                raise ValueError("Native value or source feature identity differs")
            if len(ordered) == 1 and actual["geometry"] != ordered[0][1]:
                raise ValueError("Original geometry bytes/order differ")
            shape = from_wkb(actual["geometry"])
            if len(ordered) > 1:
                expected_parts = []
                for _, part_wkb in ordered:
                    p = from_wkb(part_wkb)
                    if isinstance(p, MultiPolygon):
                        expected_parts.extend(p.geoms)
                    else:
                        expected_parts.append(p)
                if not isinstance(shape, MultiPolygon) or len(shape.geoms) != len(expected_parts):
                    raise ValueError("Grouped MultiPolygon part count differs")
                for actual_part, exp_part in zip(shape.geoms, expected_parts, strict=True):
                    if not actual_part.equals_exact(exp_part, tolerance=0.0):
                        raise ValueError("Grouped MultiPolygon part coordinates differ")
            coordinates = get_coordinates(shape)
            if (
                shape.geom_type not in {"Polygon", "MultiPolygon"}
                or shape.is_empty
                or not np.isfinite(coordinates).all()
                or (np.abs(coordinates[:, 0]) > 180).any()
                or (np.abs(coordinates[:, 1]) > 90).any()
            ):
                raise ValueError("Original geometry coordinates/type do not conform")
            bounds = actual["bbox"]
            xmin, ymin, xmax, ymax = shape.bounds
            if not all(np.isfinite(list(bounds.values()))) or not (
                bounds["xmin"] <= xmin <= xmax <= bounds["xmax"] and bounds["ymin"] <= ymin <= ymax <= bounds["ymax"]
            ):
                raise ValueError("Geometry bounds do not cover source")
            if actual["level"] != 0 or actual["parent_id"] is not None:
                raise ValueError("Invented hierarchy")
            area = _geodesic_area_km2(geod, shape)
            if not shape.is_valid or actual["area_km2"] != float(np.float32(area)):
                raise ValueError("Original polygon/derived geodesic area differs")
            if any(actual[name] is not None for name in ("up_area_km2", "outlet_lon", "outlet_lat")):
                raise ValueError("Invented upstream area or outlet")
            polygons.append((basin, stream_id, count))
            cat_bounds[count] = bounds
            actual_cat_cuts[count] = actual["upstream_boundary_cuts"]
        counts["catchments"] = count

    # Verify rivers.
    schema = schemas["rivers"]
    actual_rows = iterators["rivers"]
    count = counts["rivers"]
    source = Path(item["rivers"]["path"])
    digest = _digest(source)
    if digest != item["rivers"]["sha256"]:
        raise ValueError("Original source digest differs")
    source_hashes[str(source)] = digest
    connection = sqlite3.connect(source.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        (layer,) = connection.execute("SELECT table_name FROM gpkg_contents WHERE data_type='features'").fetchone()
        quoted = '"' + layer.replace('"', '""') + '"'
        fields = connection.execute(f"PRAGMA table_info({quoted})").fetchall()
        geometry, srs, z, m = connection.execute(
            "SELECT column_name, srs_id, z, m FROM gpkg_geometry_columns WHERE table_name=?", (layer,)
        ).fetchone()
        if (srs, z, m) != (4326, 0, 0):
            raise ValueError("Original geometry declaration differs")
        expected_types = {"INTEGER": pa.int64(), "MEDIUMINT": pa.int32(), "REAL": pa.float64()}
        for _, name, declared, nonnull, _, _ in fields:
            if name == geometry:
                continue
            field = schema.field(name)
            expected_type = expected_types[declared]
            metadata = {b"gpkg:declared_type": declared.encode()}
            if field.type != expected_type or field.nullable != (not nonnull) or field.metadata != metadata:
                raise ValueError(f"Native schema differs: {name}")
        for values in connection.execute(f"SELECT * FROM {quoted} ORDER BY fid"):
            original = dict(zip((field[1] for field in fields), values, strict=True))
            actual = next(actual_rows, None)
            if actual is None:
                raise ValueError("Missing source feature")
            count += 1
            if (
                actual["id"] != count
                or actual["source_id"] != str(original["fid"])
                or actual["processing_basin"] != basin
                or actual["source_namespace"] != f"NGA-TDX-Hydro-20230126/{basin}/rivers"
            ):
                raise ValueError("Source feature identity differs")
            for name, value in original.items():
                if name == geometry:
                    continue
                stored = actual[name]
                equal = (
                    struct.pack("!d", value) == struct.pack("!d", stored)
                    if isinstance(value, float) and isinstance(stored, float)
                    else value == stored and type(value) is type(stored)
                )
                if not equal:
                    raise ValueError(f"Native value or null differs: {name}")
            blob = original[geometry]
            if blob[:3] != b"GP\x00" or blob[3] & 0xF0:
                raise ValueError("Unsupported original geometry header")
            envelope = (blob[3] >> 1) & 7
            if envelope not in (0, 1):
                raise ValueError("Unsupported source dimensions")
            wkb = blob[8 + (32 if envelope == 1 else 0) :]
            if actual["geometry"] != wkb:
                raise ValueError("Original geometry bytes/order differ")
            shape = from_wkb(wkb)
            coordinates = get_coordinates(shape)
            if (
                shape.geom_type not in {"LineString", "MultiLineString"}
                or shape.is_empty
                or not np.isfinite(coordinates).all()
                or (np.abs(coordinates[:, 0]) > 180).any()
                or (np.abs(coordinates[:, 1]) > 90).any()
            ):
                raise ValueError("Original geometry coordinates/type do not conform")
            bounds = actual["bbox"]
            xmin, ymin, xmax, ymax = shape.bounds
            if not all(np.isfinite(list(bounds.values()))) or not (
                bounds["xmin"] <= xmin <= xmax <= bounds["xmax"] and bounds["ymin"] <= ymin <= ymax <= bounds["ymax"]
            ):
                raise ValueError("Geometry bounds do not cover source")
            if actual["level"] != 0 or actual["parent_id"] is not None:
                raise ValueError("Invented hierarchy")
            if actual["local_catchment_coverage"] != ("unavailable" if unavailable else "supplied"):
                raise ValueError("Local catchment coverage differs from original inventory")
            if actual["geometry_orientation"] != 0 or any(
                actual[name] is not None
                for name in ("component_catchment_id", "upstream_node_id", "downstream_node_id")
            ):
                raise ValueError("Invented river orientation or incidence")
            if shape.length == 0:
                degenerate += 1
            link = original["LINKNO"]
            if (basin, link) in rivers:
                raise ValueError("Duplicate source LINKNO")
            rivers[basin, link] = count
            river_bounds[count] = bounds
            river_attrs[count] = (
                link,
                original.get("Length"),
                float(shape.length),
                original.get("DSContArea"),
                original.get("USContArea"),
            )
            if original["DSLINKNO"] != -1:
                downstream.add((basin, link, original["DSLINKNO"]))
            for name in ("USLINKNO1", "USLINKNO2"):
                if original[name] != -1:
                    edge = (basin, original[name], link)
                    if edge in upstream:
                        raise ValueError("Duplicate original upstream edge")
                    upstream.add(edge)
    finally:
        connection.close()
    counts["rivers"] = count

    if not rivers or (not unavailable and not polygons):
        raise ValueError("Empty original basin domain")
    if downstream != upstream:
        raise ValueError("Original supplied graph directions disagree")
    if any((b, origin) not in rivers or (b, target) not in rivers for b, origin, target in downstream):
        raise ValueError("Unresolved original river relationship")
    expected_edges = {(rivers[b, origin], rivers[b, target]) for b, origin, target in downstream}
    successor = dict(expected_edges)
    completed: set[int] = set()
    for identity in rivers.values():
        path_ids: set[int] = set()
        current = identity
        while current in successor and current not in completed:
            if current in path_ids:
                raise ValueError("Original graph contains a cycle")
            path_ids.add(current)
            current = successor[current]
        completed.update(path_ids)

    low, high = min(river_bounds), max(river_bounds)
    graph_count = 0
    for identity, level, xmin, ymin, xmax, ymax in relations.execute(
        "SELECT * FROM graph WHERE id BETWEEN ? AND ? ORDER BY id", (low, high)
    ):
        if identity not in river_bounds or level != 0:
            raise ValueError("Graph feature coverage differs")
        if list(river_bounds[identity].values()) != [xmin, ymin, xmax, ymax]:
            raise ValueError("Graph bounds differ")
        graph_count += 1
    if graph_count != len(river_bounds):
        raise ValueError("Graph feature coverage differs")
    actual_edges = set(relations.execute("SELECT origin, target FROM edge WHERE target BETWEEN ? AND ?", (low, high)))
    if actual_edges != expected_edges:
        raise ValueError("Complete supplied flow differs")

    if any((b, stream) not in rivers for b, stream, _ in polygons):
        raise ValueError("Unmatched original polygon streamID")
    expected_associations = {(rivers[b, stream], identity) for b, stream, identity in polygons}
    if unavailable:
        if relations.execute("SELECT 1 FROM association WHERE river BETWEEN ? AND ? LIMIT 1", (low, high)).fetchone():
            raise ValueError("Invented associations for unavailable catchment source")
        return len(expected_edges), 0, 0, degenerate

    polygon_ids = [identity for _, _, identity in polygons]
    cat_low, cat_high = min(polygon_ids), max(polygon_ids)
    actual_associations = set(
        relations.execute(
            "SELECT river, catchment FROM association WHERE catchment BETWEEN ? AND ?",
            (cat_low, cat_high),
        )
    )
    if actual_associations != expected_associations:
        raise ValueError("Complete local associations differ")

    # Independent catchment connectivity verification.
    riv_up: dict[int, list[int]] = {rid: [] for rid in rivers.values()}
    riv_down: dict[int, list[int]] = {rid: [] for rid in rivers.values()}
    for origin, target in expected_edges:
        riv_up[target].append(origin)
        riv_down[origin].append(target)
    riv_to_cat = dict(expected_associations)
    exp_cat_up, exp_cat_cuts, _ = _independent_derive_catchment_graph(
        list(rivers.values()), riv_up, riv_down, riv_to_cat, river_attrs
    )
    for cid, expected_cut in exp_cat_cuts.items():
        if actual_cat_cuts[cid] != expected_cut:
            raise ValueError(f"upstream_boundary_cuts differs at catchment {cid}")
    cat_graph_count = 0
    for identity, level, xmin, ymin, xmax, ymax in relations.execute(
        "SELECT * FROM cat_graph WHERE id BETWEEN ? AND ? ORDER BY id", (cat_low, cat_high)
    ):
        if identity not in cat_bounds or level != 0:
            raise ValueError("Catchment graph feature coverage differs")
        if list(cat_bounds[identity].values()) != [xmin, ymin, xmax, ymax]:
            raise ValueError("Catchment graph bounds differ")
        cat_graph_count += 1
    if cat_graph_count != len(cat_bounds):
        raise ValueError("Catchment graph feature coverage differs")
    exp_cat_edges = {(origin, target) for target, origins in exp_cat_up.items() for origin in origins}
    actual_cat_edges = set(
        relations.execute("SELECT origin, target FROM cat_edge WHERE target BETWEEN ? AND ?", (cat_low, cat_high))
    )
    if actual_cat_edges != exp_cat_edges:
        raise ValueError("Derived catchment flow differs")
    return len(expected_edges), len(expected_associations), len(exp_cat_edges), degenerate


def _verify_source_metadata(inventory: dict[str, Any], output: Path) -> None:
    """Compare provenance declarations directly with each original SQLite file."""
    recorded_raw = json.loads((output / "source-metadata.json").read_bytes())
    recorded = recorded_raw["sources"] if isinstance(recorded_raw, dict) else recorded_raw
    expected: list[dict[str, Any]] = []
    sources = inventory["sources"]
    basins = [item["basin"] for item in sources]
    if not basins or len(set(basins)) != len(basins):
        raise ValueError("Invalid source basin coverage")
    supplied_polygons = 0
    for item in sorted(sources, key=lambda value: value["basin"]):
        for domain in ("catchments", "rivers"):
            entry = item.get(domain)
            if not isinstance(entry, dict):
                raise ValueError("Missing source coverage declaration")
            base = {"basin": item["basin"], "domain": domain}
            if "status" in entry:
                if (
                    domain != "catchments"
                    or entry.get("status") != "unavailable"
                    or set(entry) != {"status", "reason"}
                    or not isinstance(entry["reason"], str)
                    or not entry["reason"].strip()
                ):
                    raise ValueError("Invalid unavailable source coverage declaration")
                expected.append({**base, **entry})
                continue
            if not {"path", "sha256"}.issubset(entry):
                raise ValueError("Missing supplied source coverage declaration")
            supplied_polygons += domain == "catchments"
            source = Path(entry["path"])
            digest = entry["sha256"]
            connection = sqlite3.connect(source.resolve().as_uri() + "?mode=ro", uri=True)
            try:
                connection.row_factory = sqlite3.Row
                contents = [dict(row) for row in connection.execute("SELECT * FROM gpkg_contents")]
                schemas = {}
                for layer in contents:
                    if layer["data_type"] == "features":
                        name = layer["table_name"]
                        quoted = '"' + name.replace('"', '""') + '"'
                        schemas[name] = [dict(row) for row in connection.execute(f"PRAGMA table_info({quoted})")]
                expected.append(
                    {
                        **base,
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
    if not supplied_polygons:
        raise ValueError("No supplied catchment source coverage")
    if recorded != expected:
        raise ValueError("Source metadata differs from originals or coverage inventory")


def verify_originals(inventory_path: Path, output: Path) -> dict[str, object]:
    """Full source/value/geometry/edge comparison; no polygon-overlap certification."""
    inventory = json.loads(inventory_path.read_bytes())
    manifest = json.loads((output / "manifest.json").read_bytes())
    if (output / "original-inventory.json").read_bytes() != inventory_path.read_bytes():
        raise ValueError("Original inventory bytes differ")
    for name, identity in manifest["artifacts"].items():
        path = output / name
        if not path.resolve().is_relative_to(output.resolve()) or fingerprint(path) != identity:
            raise ValueError("Normalized artifact identity differs")
    if manifest.get("river_nodes") or manifest.get("graph") != "graph.parquet":
        raise ValueError("Invalid graph or node domain declaration")
    cc = manifest.get("catchment_connectivity")
    if (
        not isinstance(cc, list)
        or len(cc) != 1
        or cc[0].get("level") != 0
        or cc[0].get("status") != "supplied"
        or cc[0].get("completeness") not in ("complete", "partial")
        or cc[0].get("upstream_cuts_column") != "upstream_boundary_cuts"
        or cc[0].get("derivation") != CATCHMENT_DERIVATION
    ):
        raise ValueError("Invalid catchment_connectivity declaration")
    expected_categories = [
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
    ]
    if manifest.get("categorical_columns") != expected_categories:
        raise ValueError("Manifest local catchment coverage semantics differ")
    _verify_source_metadata(inventory, output)
    counts = {"catchments": 0, "rivers": 0}
    source_hashes: dict[str, str] = {}
    schemas: dict[str, pa.Schema] = {}
    iterators: dict[str, Generator[dict[str, Any]]] = {}
    for domain in counts:
        with ipc.open_file(str(output / f"{domain}.arrow")) as reader:
            schemas[domain] = reader.schema
        iterators[domain] = _rows(output / f"{domain}.arrow")
    coverage_field = schemas["rivers"].field("local_catchment_coverage")
    if coverage_field.type != pa.string() or coverage_field.nullable:
        raise ValueError("Invalid local catchment coverage schema")
    flow_edges = associations = catchment_edges = degenerate = 0
    try:
        with _actual_relations(output) as actual:
            for item in sorted(inventory["sources"], key=lambda value: value["basin"]):
                edges, pairs, cedges, zeros = _verify_basin(item, schemas, iterators, counts, source_hashes, actual)
                flow_edges += edges
                associations += pairs
                catchment_edges += cedges
                degenerate += zeros
            for table, expected in (
                ("graph", counts["rivers"]),
                ("edge", flow_edges),
                ("cat_graph", counts["catchments"]),
                ("cat_edge", catchment_edges),
                ("association", associations),
            ):
                if actual.execute(f"SELECT count(*) FROM {table}").fetchone()[0] != expected:
                    raise ValueError("Extra normalized relations")
            if any(next(rows, None) is not None for rows in iterators.values()):
                raise ValueError("Extra normalized feature")
    finally:
        for rows in iterators.values():
            rows.close()
    if counts["catchments"] != manifest["unit_count"] or counts["rivers"] != manifest["rivers"]["count"]:
        raise ValueError("Manifest feature counts differ")
    for path, digest in source_hashes.items():
        if _digest(Path(path)) != digest:
            raise ValueError("Source changed during verification")
    return {
        "status": "passed",
        "features": counts,
        "flow_edges": flow_edges,
        "catchment_flow_edges": catchment_edges,
        "local_associations": associations,
        "zero_length_rivers": degenerate,
        "manifest_sha256": _digest(output / "manifest.json"),
        "inventory_sha256": _digest(inventory_path),
        "geometry_coverage": "all features, exact embedded WKB/parts and bounds",
        "polygon_overlap_coverage": "grouped MultiPolygons exhaustively OGC-validated",
        "scope": "Original-to-normalized staging only; final spatial artifacts and publication require separate verification.",
    }


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


def _row_hash(basin: str, stream_id: int, fid: int, wkb: bytes, area_km2: float, bbox: dict[str, float]) -> int:
    h = hashlib.blake2b(digest_size=16)
    h.update(basin.encode())
    h.update(b"\x00")
    h.update(struct.pack("<qqf4f", stream_id, fid, area_km2, bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]))
    h.update(wkb)
    return int.from_bytes(h.digest(), "little")


def verify_derived_bundle(
    base_spatial_dir: Path,
    spatial_dir: Path,
) -> dict[str, Any]:
    """Independently verify a derived TDX connectivity spatial bundle against its verified base bundle."""
    base_manifest = json.loads((base_spatial_dir / "manifest.json").read_bytes())
    base_manifest_sha = _digest(base_spatial_dir / "manifest.json")
    manifest = json.loads((spatial_dir / "manifest.json").read_bytes())
    manifest_sha = _digest(spatial_dir / "manifest.json")

    if (
        manifest.get("format_version") != "0.4.0"
        or manifest.get("fabric_name") != "tdx"
        or manifest.get("representation") != "spatial-geoparquet-adjacency"
        or manifest.get("representation_version") != 1
    ):
        raise ValueError("Invalid derived spatial manifest header")

    cc = manifest.get("catchment_connectivity")
    if (
        not isinstance(cc, list)
        or len(cc) != 1
        or cc[0].get("level") != 0
        or cc[0].get("status") != "supplied"
        or cc[0].get("completeness") not in ("complete", "partial")
        or cc[0].get("upstream_cuts_column") != "upstream_boundary_cuts"
        or cc[0].get("derivation") != CATCHMENT_DERIVATION
    ):
        raise ValueError("Invalid derived catchment_connectivity declaration")

    for rel, identity in manifest["artifacts"].items():
        path = spatial_dir / rel
        if (
            Path(rel).is_absolute()
            or not path.resolve().is_relative_to(spatial_dir.resolve())
            or not path.is_file()
            or fingerprint(path) != identity
        ):
            raise ValueError(f"Derived bundle artifact differs or is missing: {rel}")

    # Verify source-manifest.json binding.
    if manifest.get("source_manifest") != "source-manifest.json" or _digest(
        spatial_dir / "source-manifest.json"
    ) != manifest.get("input_manifest_sha256"):
        raise ValueError("source-manifest.json binding differs")

    # Verify reused river feature and lookup artifacts are byte-identical to base.
    for rel in ("rivers.parquet", "rivers-lookup.arrow", "rivers-lookup-directory.json"):
        if manifest["artifacts"].get(rel) != base_manifest["artifacts"].get(rel):
            raise ValueError(f"Reused river artifact differs from base bundle: {rel}")

    # Verify base provenance files (except updated source-metadata.json and NOTICE) are preserved byte-for-byte.
    for rel, identity in base_manifest["artifacts"].items():
        if rel in {"manifest.json", "source-manifest.json", "source-metadata.json", "NOTICE"}:
            continue
        suffix = Path(rel).suffix.lower()
        if suffix in {".arrow", ".parquet", ".tif", ".tiff"} or "directory" in rel:
            continue
        if manifest["artifacts"].get(rel) != identity:
            raise ValueError(f"Base provenance file not preserved byte-for-byte: {rel}")

    # Verify rivers-relationships.arrow upstream_ids/downstream_ids are completely untouched.
    with (
        ipc.open_file(str(base_spatial_dir / "rivers-relationships.arrow")) as base_rrel,
        ipc.open_file(str(spatial_dir / "rivers-relationships.arrow")) as new_rrel,
    ):
        if base_rrel.num_record_batches != new_rrel.num_record_batches:
            raise ValueError("rivers-relationships batch count differs")
        for idx in range(base_rrel.num_record_batches):
            b0 = base_rrel.get_batch(idx)
            b1 = new_rrel.get_batch(idx)
            for col in ("id", "upstream_ids", "downstream_ids"):
                if not b0.column(col).equals(b1.column(col)):
                    raise ValueError(f"River flow graph modified in rivers-relationships.{col}")

    # Pass 1: verify derived catchments.parquet geometry, validity, area, bbox, and build key maps.
    geod = Geod(ellps="WGS84")
    stream_to_cid: dict[str, dict[int, int]] = defaultdict(dict)
    derived_cuts: dict[int, int | None] = {}
    ungrouped_xor = 0
    ungrouped_count = 0
    grouped_expected: dict[tuple[str, int], tuple[list[int], list[bytes]]] = {}
    grouped_basins_counts: Counter[str] = Counter()

    with pq.ParquetFile(spatial_dir / "catchments.parquet") as cat_pf:
        if cat_pf.metadata.num_rows != manifest["unit_count"]:
            raise ValueError("Catchment row count differs from manifest unit_count")
        for batch in _iter_parquet_batches(cat_pf):
            for row in batch.to_pylist():
                cid = int(row["id"])
                basin = str(row["processing_basin"])
                sid = int(row["streamID"])
                fid = int(row["fid"])
                fids = cast(list[int], row["source_fids"])
                wkb = bytes(row["geometry"])
                bbox = cast(dict[str, float], row["bbox"])
                area_km2 = float(row["area_km2"])
                if sid in stream_to_cid[basin]:
                    raise ValueError(f"Duplicate (processing_basin, streamID) in derived catchments: {(basin, sid)}")
                stream_to_cid[basin][sid] = cid
                derived_cuts[cid] = row["upstream_boundary_cuts"]
                if not fids or fids != sorted(set(fids)) or fid != fids[0] or row["source_id"] != str(fids[0]):
                    raise ValueError(f"Invalid fid/source_fids provenance at catchment {cid}")
                if len(fids) == 1:
                    ungrouped_xor ^= _row_hash(basin, sid, fid, wkb, area_km2, bbox)
                    ungrouped_count += 1
                else:
                    grouped_basins_counts[basin] += 1
                    shape = from_wkb(wkb)
                    if not isinstance(shape, MultiPolygon) or not shape.is_valid:
                        raise ValueError(f"Grouped catchment {cid} is not a valid MultiPolygon")
                    xmin, ymin, xmax, ymax = shape.bounds
                    if not (
                        bbox["xmin"] <= xmin <= xmax <= bbox["xmax"] and bbox["ymin"] <= ymin <= ymax <= bbox["ymax"]
                    ):
                        raise ValueError(f"Grouped catchment {cid} bbox does not cover geometry")
                    exp_area = _geodesic_area_km2(geod, shape)
                    if area_km2 != float(np.float32(exp_area)):
                        raise ValueError(f"Grouped catchment {cid} geodesic area differs")
                    part_hashes = [
                        hashlib.blake2b(to_wkb(p, output_dimension=2, byte_order=1), digest_size=16).digest()
                        for p in shape.geoms
                    ]
                    grouped_expected[(basin, sid)] = (fids, part_hashes)

    # Pass 2: verify every base catchment row matches either ungrouped XOR or grouped MultiPolygon parts.
    base_ungrouped_xor = 0
    base_ungrouped_count = 0
    base_grouped_parts: dict[tuple[str, int], list[tuple[int, list[bytes]]]] = defaultdict(list)

    with pq.ParquetFile(base_spatial_dir / "catchments.parquet") as base_cat_pf:
        for batch in _iter_parquet_batches(
            base_cat_pf,
            columns=["processing_basin", "streamID", "fid", "geometry", "area_km2", "bbox"],
        ):
            for basin, sid, fid, wkb, area_km2, bbox in zip(
                cast(list[str], batch.column("processing_basin").to_pylist()),
                cast(list[int], batch.column("streamID").to_pylist()),
                cast(list[int], batch.column("fid").to_pylist()),
                cast(list[bytes], batch.column("geometry").to_pylist()),
                cast(list[float], batch.column("area_km2").to_pylist()),
                cast(list[dict[str, float]], batch.column("bbox").to_pylist()),
                strict=True,
            ):
                key = (basin, sid)
                if key in grouped_expected:
                    shape = from_wkb(wkb)
                    parts = list(shape.geoms) if isinstance(shape, MultiPolygon) else [shape]
                    hashes = [
                        hashlib.blake2b(to_wkb(p, output_dimension=2, byte_order=1), digest_size=16).digest()
                        for p in parts
                    ]
                    base_grouped_parts[key].append((fid, hashes))
                else:
                    if sid not in stream_to_cid.get(basin, {}):
                        raise ValueError(f"Base catchment missing in derived bundle: {key}")
                    base_ungrouped_xor ^= _row_hash(basin, sid, fid, wkb, area_km2, bbox)
                    base_ungrouped_count += 1

    if base_ungrouped_count != ungrouped_count or base_ungrouped_xor != ungrouped_xor:
        raise ValueError("Ungrouped catchment rows differ from base bundle")
    if set(base_grouped_parts) != set(grouped_expected):
        raise ValueError("Grouped catchment key set differs from base bundle")
    extra_fragments_collapsed = 0
    for key, (exp_fids, exp_part_hashes) in grouped_expected.items():
        ordered = sorted(base_grouped_parts[key], key=lambda item: item[0])
        actual_fids = [fid for fid, _ in ordered]
        actual_part_hashes = [h for _, hs in ordered for h in hs]
        if actual_fids != exp_fids or actual_part_hashes != exp_part_hashes:
            raise ValueError(f"Grouped MultiPolygon parts differ from base polygons for {key}")
        extra_fragments_collapsed += len(actual_fids) - 1
    del base_grouped_parts, grouped_expected

    # Pass 3: independently derive catchment connectivity and local associations from base rivers.
    rivers_by_basin: dict[str, list[tuple[int, int, int, int, int, float | None, float | None, float | None]]] = (
        defaultdict(list)
    )
    polygonless_ids: set[int] = set()
    with pq.ParquetFile(base_spatial_dir / "rivers.parquet") as riv_pf:
        for batch in _iter_parquet_batches(
            riv_pf,
            columns=[
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
            ],
        ):
            for rid, basin, linkno, ds, u1, u2, length, ds_area, us_area, cov in zip(
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
                strict=True,
            ):
                rivers_by_basin[basin].append((rid, linkno, ds, u1, u2, length, ds_area, us_area))
                if cov == "supplied" and linkno not in stream_to_cid.get(basin, {}):
                    polygonless_ids.add(rid)

    # Read geometry only for polygonless rivers via rivers-lookup.arrow.
    poly_geom_lengths: dict[int, float] = {}
    if polygonless_ids:
        by_block: dict[int, list[tuple[int, int]]] = defaultdict(list)
        min_pid, max_pid = min(polygonless_ids), max(polygonless_ids)
        with ipc.open_file(str(base_spatial_dir / "rivers-lookup.arrow")) as reader:
            for bidx in range(reader.num_record_batches):
                batch = reader.get_batch(bidx)
                ids = batch.column("id").to_numpy()
                if ids[-1] < min_pid or ids[0] > max_pid:
                    continue
                blocks = batch.column("block").to_numpy()
                rows = batch.column("row").to_numpy()
                for i, rid in enumerate(ids.tolist()):
                    if rid in polygonless_ids:
                        by_block[int(blocks[i])].append((int(rows[i]), rid))
        with pq.ParquetFile(base_spatial_dir / "rivers.parquet") as riv_pf:
            for block, entries in sorted(by_block.items()):
                table = riv_pf.read_row_group(block, columns=["id", "geometry"]).combine_chunks()
                id_col = table.column("id")
                geom_col = table.column("geometry")
                for rpos, exp_rid in entries:
                    if id_col[rpos].as_py() != exp_rid:
                        raise ValueError(f"rivers-lookup mismatch at {exp_rid}")
                    poly_geom_lengths[exp_rid] = float(from_wkb(geom_col[rpos].as_py()).length)

    expected_cat_up: dict[int, list[int]] = {}
    expected_cat_down: dict[int, list[int]] = defaultdict(list)
    expected_cat_local_riv: dict[int, list[int]] = {}
    expected_riv_local_cat: dict[int, list[int]] = {}
    totals = Counter[str]()

    for basin in sorted(rivers_by_basin):
        b_rivers = rivers_by_basin[basin]
        link_to_rid = {r[1]: r[0] for r in b_rivers}
        riv_up: dict[int, list[int]] = {}
        riv_down: dict[int, list[int]] = {r[0]: [] for r in b_rivers}
        riv_to_cat: dict[int, int] = {}
        poly_attrs: dict[int, tuple[int, float | None, float | None, float | None, float | None]] = {}
        b_s2c = stream_to_cid.get(basin, {})
        for rid, linkno, ds, u1, u2, length, ds_area, us_area in b_rivers:
            ups = sorted(link_to_rid[u] for u in (u1, u2) if u != -1)
            riv_up[rid] = ups
            if ds != -1:
                riv_down[rid].append(link_to_rid[ds])
            cid = b_s2c.get(linkno)
            if cid is not None:
                riv_to_cat[rid] = cid
                expected_riv_local_cat[rid] = [cid]
                expected_cat_local_riv[cid] = [rid]
            else:
                expected_riv_local_cat[rid] = []
                if rid in poly_geom_lengths:
                    poly_attrs[rid] = (linkno, length, poly_geom_lengths[rid], ds_area, us_area)
        if b_s2c:
            c_up, c_cuts, stats = _independent_derive_catchment_graph(
                [r[0] for r in b_rivers], riv_up, riv_down, riv_to_cat, poly_attrs
            )
            for k, v in stats.items():
                totals[k] += v
            for cid, ups in c_up.items():
                expected_cat_up[cid] = ups
                for u in ups:
                    expected_cat_down[u].append(cid)
                if derived_cuts[cid] != c_cuts[cid]:
                    raise ValueError(f"upstream_boundary_cuts differs at catchment {cid}")
    del rivers_by_basin, poly_geom_lengths, derived_cuts

    # Pass 4: verify catchments-relationships.arrow and rivers-relationships.arrow local_catchment_ids.
    seen_cats = 0
    with ipc.open_file(str(spatial_dir / "catchments-relationships.arrow")) as reader:
        for bidx in range(reader.num_record_batches):
            batch = reader.get_batch(bidx)
            for cid, ups, downs, lrivs in zip(
                cast(list[int], batch.column("id").to_pylist()),
                cast(list[list[int]], batch.column("upstream_ids").to_pylist()),
                cast(list[list[int]], batch.column("downstream_ids").to_pylist()),
                cast(list[list[int]], batch.column("local_river_ids").to_pylist()),
                strict=True,
            ):
                if ups != expected_cat_up[cid]:
                    raise ValueError(f"catchments-relationships.upstream_ids differs at {cid}")
                if downs != sorted(expected_cat_down.get(cid, [])):
                    raise ValueError(f"catchments-relationships.downstream_ids differs at {cid}")
                if lrivs != expected_cat_local_riv[cid]:
                    raise ValueError(f"catchments-relationships.local_river_ids differs at {cid}")
                seen_cats += 1
    if seen_cats != manifest["unit_count"]:
        raise ValueError("catchments-relationships row count differs")

    seen_rivs = 0
    with ipc.open_file(str(spatial_dir / "rivers-relationships.arrow")) as reader:
        for bidx in range(reader.num_record_batches):
            batch = reader.get_batch(bidx)
            for rid, lcats in zip(
                cast(list[int], batch.column("id").to_pylist()),
                cast(list[list[int]], batch.column("local_catchment_ids").to_pylist()),
                strict=True,
            ):
                if lcats != expected_riv_local_cat[rid]:
                    raise ValueError(f"rivers-relationships.local_catchment_ids differs at {rid}")
                seen_rivs += 1
    if seen_rivs != manifest["families"]["rivers"]["count"]:
        raise ValueError("rivers-relationships row count differs")

    # Pass 5: verify source-metadata.json audit counts match independent verification.
    smeta = json.loads((spatial_dir / "source-metadata.json").read_bytes())
    if smeta.get("derived_from_verified_bundle_manifest_sha256") != base_manifest_sha:
        raise ValueError("source-metadata.json derived_from_verified_bundle_manifest_sha256 differs")
    cc_audit = smeta.get("catchment_connectivity", {})
    for key in (
        "contracted_stub_count",
        "non_contracted_linkno0_count",
        "non_contracted_other_count",
        "catchments_with_cuts_count",
        "affected_downstream_catchments_count",
    ):
        if cc_audit.get(key) != totals[key]:
            raise ValueError(f"source-metadata.json catchment_connectivity.{key} differs")
    if cc_audit.get("grouped_multipolygon_count") != sum(grouped_basins_counts.values()):
        raise ValueError("source-metadata.json grouped_multipolygon_count differs")
    if cc_audit.get("extra_fragment_rows_collapsed") != extra_fragments_collapsed:
        raise ValueError("source-metadata.json extra_fragment_rows_collapsed differs")

    return {
        "status": "passed",
        "base_manifest_sha256": base_manifest_sha,
        "manifest_sha256": manifest_sha,
        "catchment_count": seen_cats,
        "river_count": seen_rivs,
        "grouped_multipolygon_count": sum(grouped_basins_counts.values()),
        "extra_fragment_rows_collapsed": extra_fragments_collapsed,
        "grouped_basins": dict(sorted(grouped_basins_counts.items())),
        "contracted_stub_count": totals["contracted_stub_count"],
        "non_contracted_linkno0_count": totals["non_contracted_linkno0_count"],
        "non_contracted_other_count": totals["non_contracted_other_count"],
        "catchments_with_cuts_count": totals["catchments_with_cuts_count"],
        "affected_downstream_catchments_count": totals["affected_downstream_catchments_count"],
        "scope": "Independent derived-bundle comparison against verified base spatial bundle.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Inventory JSON path, or base spatial dir with --from-bundle")
    parser.add_argument("output", type=Path, help="Normalized staging dir, or derived spatial dir with --from-bundle")
    parser.add_argument("--from-bundle", action="store_true", help="Verify derived spatial bundle against base bundle")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists() or args.report.resolve().is_relative_to(args.output.resolve()):
        raise ValueError("Report must be new and outside bundle")
    result = (
        verify_derived_bundle(args.source, args.output)
        if args.from_bundle
        else verify_originals(args.source, args.output)
    )
    with args.report.open("x") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
