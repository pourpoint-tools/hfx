"""Test-first coverage for TDX polygon fragment collection and catchment connectivity."""

import hashlib
import json
import math
import sqlite3
import struct
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.parquet as pq
import pytest
from shapely import MultiPolygon, Polygon, from_wkb, to_wkb

from hfx.spatial_bundle import build as build_spatial
from hfx.spatial_bundle import verify as verify_spatial


def _poly_wkb(xmin: float, ymin: float, xmax: float, ymax: float) -> bytes:
    return to_wkb(
        Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]),
        output_dimension=2,
        byte_order=1,
    )


def _line_wkb(x0: float, y0: float, x1: float, y1: float) -> bytes:
    return struct.pack("<BII4d", 1, 2, 2, x0, y0, x1, y1)


def _gpkg_geom(wkb: bytes) -> bytes:
    return b"GP" + bytes([0, 1]) + struct.pack("<i", 4326) + wkb


def _write_tdx_gpkg(
    path: Path,
    *,
    domain: str,
    rows: list[tuple[Any, ...]],
    geom_type: str | None = None,
) -> dict[str, str]:
    layer = "basins" if domain == "catchments" else "TDX_streamnet_test_01"
    kind = geom_type or ("POLYGON" if domain == "catchments" else "LINESTRING")
    with sqlite3.connect(path) as conn:
        conn.execute("CREATE TABLE gpkg_contents (table_name TEXT, data_type TEXT, srs_id INTEGER)")
        conn.execute("INSERT INTO gpkg_contents VALUES (?, 'features', 4326)", (layer,))
        conn.execute(
            "CREATE TABLE gpkg_geometry_columns "
            "(table_name TEXT, column_name TEXT, geometry_type_name TEXT, srs_id INTEGER, z INTEGER, m INTEGER)"
        )
        conn.execute("INSERT INTO gpkg_geometry_columns VALUES (?, 'geom', ?, 4326, 0, 0)", (layer, kind))
        if domain == "catchments":
            conn.execute("CREATE TABLE basins (fid INTEGER PRIMARY KEY NOT NULL, geom POLYGON, streamID INTEGER)")
            conn.executemany("INSERT INTO basins VALUES (?, ?, ?)", rows)
        else:
            conn.execute(
                f"CREATE TABLE {layer} ("
                "fid INTEGER PRIMARY KEY NOT NULL, geom LINESTRING, "
                "LINKNO MEDIUMINT, DSLINKNO MEDIUMINT, USLINKNO1 MEDIUMINT, USLINKNO2 MEDIUMINT, "
                "Length REAL, DSContArea REAL, USContArea REAL)"
            )
            conn.executemany(f"INSERT INTO {layer} VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def test_collect_polygon_fragments_corner_touching_valid() -> None:
    from scripts.produce_tdx import collect_polygon_fragments

    p1 = _poly_wkb(0.0, 0.0, 1.0, 1.0)
    p2 = _poly_wkb(1.0, 1.0, 2.0, 2.0)  # corner-touching at (1, 1) -> OGC valid MultiPolygon
    wkb, area, bbox, fids = collect_polygon_fragments("1020018110", 42, [(17, p2), (5, p1)])
    assert fids == [5, 17]
    mp = from_wkb(wkb)
    assert isinstance(mp, MultiPolygon)
    assert mp.is_valid
    assert len(mp.geoms) == 2
    assert mp.geoms[0].equals_exact(from_wkb(p1), tolerance=0.0)
    assert mp.geoms[1].equals_exact(from_wkb(p2), tolerance=0.0)
    assert bbox["xmin"] <= 0.0 and bbox["ymin"] <= 0.0 and bbox["xmax"] >= 2.0 and bbox["ymax"] >= 2.0
    # WGS84 geodesic area of [0,1]x[0,1] (12308.778 km2) + [1,2]x[1,2] (12305.129 km2).
    assert area == pytest.approx(24613.907, abs=0.01)


def test_collect_polygon_fragments_rejects_overlapping_or_edge_sharing() -> None:
    from scripts.produce_tdx import collect_polygon_fragments

    p1 = _poly_wkb(0.0, 0.0, 1.0, 1.0)
    p_edge = _poly_wkb(1.0, 0.0, 2.0, 1.0)  # shares edge x=1 -> invalid MultiPolygon under OGC
    p_overlap = _poly_wkb(0.5, 0.5, 1.5, 1.5)
    with pytest.raises(ValueError, match="Invalid grouped MultiPolygon"):
        collect_polygon_fragments("1020018110", 42, [(1, p1), (2, p_edge)])
    with pytest.raises(ValueError, match="Invalid grouped MultiPolygon"):
        collect_polygon_fragments("1020018110", 42, [(1, p1), (2, p_overlap)])


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": 0.0,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            True,
        ),
        # LINKNO == 0 is never contracted even if attributes are zero
        (
            {
                "linkno": 0,
                "length": 0.0,
                "geom_length": 0.0,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
        # Any single check failing -> False
        (
            {
                "linkno": 10,
                "length": 0.001,
                "geom_length": 0.0,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": 1e-6,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": 0.0,
                "ds_cont_area": 5.1,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": 0.0,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 0,
                "downstream_degree": 1,
            },
            False,
        ),
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": 0.0,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 0,
            },
            False,
        ),
        # Fail closed on None / NaN / inf
        (
            {
                "linkno": 10,
                "length": None,
                "geom_length": 0.0,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": math.nan,
                "ds_cont_area": 5.0,
                "us_cont_area": 5.0,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
        (
            {
                "linkno": 10,
                "length": 0.0,
                "geom_length": 0.0,
                "ds_cont_area": math.inf,
                "us_cont_area": math.inf,
                "upstream_degree": 2,
                "downstream_degree": 1,
            },
            False,
        ),
    ],
)
def test_is_contractible_virtual_stub_fail_closed(kwargs: dict[str, Any], expected: bool) -> None:
    from scripts.produce_tdx import is_contractible_virtual_stub

    assert is_contractible_virtual_stub(**kwargs) is expected


def test_derive_basin_catchment_connectivity_contraction_chains_and_cuts() -> None:
    """Test 3-way stub, chained stubs, cut propagation through stubs, and non-bridging across real gaps."""
    from scripts.produce_tdx import derive_basin_catchment_connectivity

    # Network:
    # 1 (cat 101) \
    #              -> 3 (virtual stub) \
    # 2 (cat 102) /                     -> 5 (virtual stub) -> 6 (cat 106)
    # 4 (cat 104) ----------------------/
    #
    # 7 (LINKNO=0, no cat) \
    #                       -> 9 (virtual stub) -> 10 (cat 110) -> 11 (cat 111)
    # 8 (cat 108) ---------/
    #
    # 12 (cat 112) -> 13 (nonzero polygonless, NOT contracted) -> 14 (cat 114)
    river_flow = {
        1: [],
        2: [],
        3: [1, 2],
        4: [],
        5: [3, 4],
        6: [5],
        7: [],
        8: [],
        9: [7, 8],
        10: [9],
        11: [10],
        12: [],
        13: [12],
        14: [13],
    }
    river_to_catchment = {1: 101, 2: 102, 4: 104, 6: 106, 8: 108, 10: 110, 11: 111, 12: 112, 14: 114}
    polygonless_attrs = {
        3: (30, 0.0, 0.0, 10.0, 10.0),
        5: (50, 0.0, 0.0, 15.0, 15.0),
        7: (0, 2.5, 0.02, 5.0, 0.0),  # LINKNO == 0 headwater
        9: (90, 0.0, 0.0, 8.0, 8.0),  # virtual stub below LINKNO == 0
        13: (130, 1.0, 0.01, 9.0, 5.0),  # nonzero polygonless internal river
    }
    result = derive_basin_catchment_connectivity(river_flow, river_to_catchment, polygonless_attrs)

    # Chained virtual stubs 3 -> 5 contract into 3-way confluence at cat 106 with 0 cuts
    assert result.catchment_upstream[106] == [101, 102, 104]
    assert result.catchment_cuts[106] == 0

    # Virtual stub 9 propagates cat 108 AND the 1 cut from LINKNO==0 river 7 into cat 110
    assert result.catchment_upstream[110] == [108]
    assert result.catchment_cuts[110] == 1
    assert result.catchment_upstream[111] == [110]
    assert result.catchment_cuts[111] == 0

    # Non-contracted polygonless river 13 MUST NOT bridge cat 112 into cat 114
    assert result.catchment_upstream[114] == []
    assert result.catchment_cuts[114] == 1

    # Audit counts
    assert result.contracted_stub_count == 3
    assert result.max_stub_chain_depth == 2
    assert result.non_contracted_linkno0_count == 1
    assert result.non_contracted_other_count == 1
    assert result.catchments_with_cuts_count == 2  # 110 and 114
    assert result.affected_downstream_catchments_count == 3  # 110, 111, 114


def test_produce_from_bundle_and_reuse_rivers_end_to_end(tmp_path: Path) -> None:
    """Build a base spatial bundle (not_supplied, multi-polygon), derive connectivity, reuse rivers, verify."""
    from scripts.produce_tdx import produce_from_bundle
    from scripts.verify_tdx import verify_derived_bundle

    # Build a synthetic base spatial bundle mimicking full62-20261027:
    # Basin 1020018110:
    #   rivers:
    #     LINKNO 10 (headwater, 2 polygon fragments fids 5 & 17)
    #     LINKNO 20 (headwater, 1 polygon fid 6)
    #     LINKNO 30 (virtual stub: Length=0, geom_len=0, DS=US=10.0, up=[10,20], down=40)
    #     LINKNO 40 (downstream of 30, 1 polygon fid 8)
    #     LINKNO 0  (real headwater, no polygon, Length=3.0, DS=5.0, US=0.0, down=50)
    #     LINKNO 50 (downstream of 0, 1 polygon fid 9)
    base_staging = tmp_path / "base-staging"
    base_staging.mkdir()

    p10a = _poly_wkb(0.0, 0.0, 0.5, 0.5)
    p10b = _poly_wkb(0.5, 0.5, 1.0, 1.0)
    p20 = _poly_wkb(2.0, 0.0, 2.5, 0.5)
    p40 = _poly_wkb(3.0, 0.0, 3.5, 0.5)
    p50 = _poly_wkb(4.0, 0.0, 4.5, 0.5)

    line_real = _line_wkb(0.0, 0.0, 0.1, 0.1)
    line_zero = _line_wkb(1.0, 1.0, 1.0, 1.0)

    basin = "1020018110"
    cat_gpkg = tmp_path / f"{basin}-cat.gpkg"
    riv_gpkg = tmp_path / f"{basin}-riv.gpkg"
    cat_inv = _write_tdx_gpkg(
        cat_gpkg,
        domain="catchments",
        rows=[
            (5, _gpkg_geom(p10a), 10),
            (6, _gpkg_geom(p20), 20),
            (8, _gpkg_geom(p40), 40),
            (9, _gpkg_geom(p50), 50),
            (17, _gpkg_geom(p10b), 10),
        ],
    )
    riv_inv = _write_tdx_gpkg(
        riv_gpkg,
        domain="rivers",
        rows=[
            (1, _gpkg_geom(line_real), 0, 50, -1, -1, 3.0, 5.0, 0.0),
            (2, _gpkg_geom(line_real), 10, 30, -1, -1, 2.0, 4.0, 0.0),
            (3, _gpkg_geom(line_real), 20, 30, -1, -1, 2.0, 6.0, 0.0),
            (4, _gpkg_geom(line_zero), 30, 40, 10, 20, 0.0, 10.0, 10.0),
            (5, _gpkg_geom(line_real), 40, -1, 30, -1, 4.0, 15.0, 10.0),
            (6, _gpkg_geom(line_real), 50, -1, 0, -1, 4.0, 9.0, 5.0),
        ],
    )
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(
        json.dumps(
            {
                "release": "NGA-TDX-Hydro-20230126",
                "sources": [{"basin": basin, "catchments": cat_inv, "rivers": riv_inv}],
            }
        )
    )

    # First produce a legacy-style base staging (un-grouped, not_supplied) via helper or directly,
    # then convert to base spatial bundle.
    from hfx.candidates import fingerprint
    from hfx.tdx_source import read_layer
    from hfx.vector_storage import write_vectors
    from scripts.produce_tdx import normalized_batch, source_relationships

    cat_batch = normalized_batch(next(read_layer(cat_gpkg, "catchments")), "catchments", basin, first_id=1)
    riv_batch = normalized_batch(next(read_layer(riv_gpkg, "rivers")), "rivers", basin, first_id=1)
    cat_dir = write_vectors(
        cat_batch.schema,
        iter([cat_batch]),
        base_staging / "catchments.arrow",
        base_staging / "lookup.arrow",
        lookup_binding="hfx:catchments_sha256",
    )
    (base_staging / "lookup-directory.json").write_text(json.dumps(cat_dir) + "\n")
    riv_dir = write_vectors(
        riv_batch.schema,
        iter([riv_batch]),
        base_staging / "rivers.arrow",
        base_staging / "river-lookup.arrow",
        lookup_binding="hfx:data_sha256",
    )
    (base_staging / "river-lookup-directory.json").write_text(json.dumps(riv_dir) + "\n")

    riv_rows = riv_batch.to_pylist()
    cat_rows = cat_batch.to_pylist()
    rflow, assocs = source_relationships(
        [(r["id"], r["LINKNO"], r["DSLINKNO"], r["USLINKNO1"], r["USLINKNO2"]) for r in riv_rows],
        [(c["id"], c["streamID"]) for c in cat_rows],
    )
    graph_schema = pa.schema(
        [
            pa.field("id", pa.int64(), False),
            pa.field("level", pa.int16(), False),
            pa.field("upstream_ids", pa.list_(pa.int64()), False),
        ]
        + [pa.field("bbox_" + n, pa.float32(), False) for n in ("minx", "miny", "maxx", "maxy")]
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "id": r["id"],
                    "level": 0,
                    "upstream_ids": rflow[r["id"]],
                    "bbox_minx": r["bbox"]["xmin"],
                    "bbox_miny": r["bbox"]["ymin"],
                    "bbox_maxx": r["bbox"]["xmax"],
                    "bbox_maxy": r["bbox"]["ymax"],
                }
                for r in riv_rows
            ],
            schema=graph_schema,
        ),
        base_staging / "river-graph.parquet",
    )
    pq.write_table(
        pa.Table.from_pylist(
            [{"river_id": r, "catchment_id": c} for r, c in assocs],
            schema=pa.schema([pa.field("river_id", pa.int64(), False), pa.field("catchment_id", pa.int64(), False)]),
        ),
        base_staging / "local-associations.parquet",
    )
    (base_staging / "original-inventory.json").write_bytes(inv_path.read_bytes())
    (base_staging / "source-metadata.json").write_text(
        json.dumps([{"basin": basin, "domain": "catchments", "status": "supplied"}]) + "\n"
    )
    (base_staging / "NOTICE").write_text("Base NOTICE\n")
    (base_staging / "provenance").mkdir()
    (base_staging / "provenance" / "raw-NOTICE.txt").write_text("Raw notice\n")
    base_manifest_obj = {
        "format_version": "0.4.0",
        "fabric_name": "tdx",
        "fabric_version": "NGA-TDX-Hydro-20230126",
        "crs": "EPSG:4326",
        "has_up_area": False,
        "topology": "tree",
        "bbox": [0.0, 0.0, 4.5, 1.0],
        "unit_count": 5,
        "created_at": "2026-10-27T04:02:06+00:00",
        "adapter_version": "tdx-original-2",
        "region": basin,
        "catchments": "catchments.arrow",
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
        "catchment_connectivity": [{"level": 0, "status": "not_supplied"}],
        "rivers": {
            "data": "rivers.arrow",
            "graph": "river-graph.parquet",
            "lookup": "river-lookup.arrow",
            "directory": "river-lookup-directory.json",
            "count": 6,
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
            p.relative_to(base_staging).as_posix(): fingerprint(p)
            for p in sorted(base_staging.rglob("*"))
            if p.is_file()
        },
    }
    (base_staging / "manifest.json").write_text(json.dumps(base_manifest_obj, indent=2) + "\n")

    base_spatial = tmp_path / "base-spatial"
    build_spatial(
        base_staging,
        base_spatial,
        source_manifest_sha256=hashlib.sha256((base_staging / "manifest.json").read_bytes()).hexdigest(),
    )

    # Now derive new staging from base_spatial, build new spatial reusing rivers, and verify!
    derived_staging = tmp_path / "derived-staging"
    produce_from_bundle(base_spatial, derived_staging)

    derived_spatial = tmp_path / "derived-spatial"
    staging_sha = hashlib.sha256((derived_staging / "manifest.json").read_bytes()).hexdigest()
    build_spatial(
        derived_staging,
        derived_spatial,
        source_manifest_sha256=staging_sha,
        reuse_families={"rivers": base_spatial},
    )

    # 1. Full spatial_bundle.verify passes (including reused rivers!)
    sp_report = verify_spatial(derived_staging, derived_spatial, source_manifest_sha256=staging_sha)
    assert sp_report["status"] == "passed"

    # 2. Reused river artifacts are byte-identical to base_spatial
    for rel in ("rivers.parquet", "rivers-lookup.arrow", "rivers-lookup-directory.json"):
        assert (derived_spatial / rel).read_bytes() == (base_spatial / rel).read_bytes()

    # 3. Extra provenance preserved
    assert (derived_spatial / "provenance" / "raw-NOTICE.txt").read_text() == "Raw notice\n"

    # 4. Derived catchments count == 4 (5 base rows minus 1 extra fragment for streamID 10)
    manifest = json.loads((derived_spatial / "manifest.json").read_text())
    assert manifest["unit_count"] == 4
    assert manifest["catchment_connectivity"] == [
        {
            "level": 0,
            "status": "supplied",
            "completeness": "partial",
            "upstream_cuts_column": "upstream_boundary_cuts",
            "derivation": "river_local_associations_with_zero_area_virtual_stub_contraction",
        }
    ]
    cats = pq.read_table(derived_spatial / "catchments.parquet").sort_by("id")
    assert cats["streamID"].to_pylist() == [10, 20, 40, 50]
    assert cats["source_fids"].to_pylist() == [[5, 17], [6], [8], [9]]
    assert cats["upstream_boundary_cuts"].to_pylist() == [0, 0, 0, 1]

    # Catchment relationships: cat 3 (streamID 40) has upstream [1, 2] via contracted stub 30
    with ipc.open_file(str(derived_spatial / "catchments-relationships.arrow")) as reader:
        crel = reader.read_all()
    assert crel["upstream_ids"].to_pylist() == [[], [], [1, 2], []]

    # 5. Independent derived-bundle verification passes
    derived_report = verify_derived_bundle(base_spatial, derived_spatial)
    assert derived_report["status"] == "passed"
    assert derived_report["grouped_multipolygon_count"] == 1
    assert derived_report["contracted_stub_count"] == 1
    assert derived_report["non_contracted_linkno0_count"] == 1


def test_independent_derive_handles_deep_chain_without_recursion_error() -> None:
    """Catchment chains deeper than Python's default recursion limit (1000) must not fail."""
    from scripts.verify_tdx import _independent_derive_catchment_graph

    depth = 1500
    river_ids = list(range(1, depth + 2))
    river_upstream = {1: []}
    river_downstream: dict[int, list[int]] = {i: [] for i in river_ids}
    for i in range(2, depth + 2):
        river_upstream[i] = [i - 1]
        river_downstream[i - 1].append(i)
    # River 1 is LINKNO==0 headwater (no catchment); rivers 2..depth+1 have catchments.
    river_to_cat = {i: i - 1 for i in range(2, depth + 2)}
    poly_attrs = {1: (0, 2.0, 0.01, 5.0, 0.0)}
    cat_up, cat_cuts, stats = _independent_derive_catchment_graph(
        river_ids, river_upstream, river_downstream, river_to_cat, poly_attrs
    )
    assert len(cat_up) == depth
    assert cat_cuts[1] == 1
    assert stats["affected_downstream_catchments_count"] == depth
