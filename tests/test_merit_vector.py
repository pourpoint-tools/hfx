"""Test-first verification for MERIT Basins vector normalization and independent binary shapefile readback."""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, cast

import pytest

for dependency in ("fiona", "numpy", "pyarrow", "pyproj", "shapely"):
    pytest.importorskip(dependency)

import fiona  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_shapefile(
    path: Path,
    *,
    geom_type: str,
    properties_schema: dict[str, str],
    features: list[dict[str, Any]],
    write_prj: bool = True,
) -> dict[str, dict[str, Any]]:
    schema = {"geometry": geom_type, "properties": properties_schema}
    kwargs: dict[str, Any] = {"driver": "ESRI Shapefile", "schema": schema, "encoding": "UTF-8"}
    if write_prj:
        kwargs["crs"] = "EPSG:4326"
    with fiona.open(path, "w", **kwargs) as dst:
        for feat in features:
            dst.write(feat)
    if not write_prj:
        path.with_suffix(".prj").unlink(missing_ok=True)
    sidecars: dict[str, dict[str, Any]] = {}
    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg"):
        side = path.with_suffix(ext)
        if side.exists():
            sidecars[side.name] = {"path": str(side), "size_bytes": side.stat().st_size, "sha256": _sha256(side)}
    return sidecars


RIV_PROPS = {
    "COMID": "int:10",
    "lengthkm": "float:16.6",
    "lengthdir": "float:16.6",
    "sinuosity": "float:16.6",
    "slope": "float:16.6",
    "uparea": "float:16.6",
    "order": "int:4",
    "strmDrop_t": "float:16.6",
    "slope_taud": "float:16.6",
    "NextDownID": "int:10",
    "maxup": "int:4",
    "up1": "int:10",
    "up2": "int:10",
    "up3": "int:10",
    "up4": "int:10",
}

CAT_PROPS = {"COMID": "int:10", "unitarea": "float:16.6"}
HILL_PROPS = {"FID": "int:10", "area_km2": "float:16.6"}  # Intentional collision test on native column name!


def _river_feature(
    comid: int, coords: list[tuple[float, float]], next_down: int, up: tuple[int, ...]
) -> dict[str, Any]:
    padded = list(up) + [0] * (4 - len(up))
    return {
        "geometry": {"type": "LineString", "coordinates": coords},
        "properties": {
            "COMID": comid,
            "lengthkm": 1.5,
            "lengthdir": 1.2,
            "sinuosity": 1.25,
            "slope": 0.01,
            "uparea": 30.0,
            "order": 1 if not up else 2,
            "strmDrop_t": 5.0,
            "slope_taud": 0.008,
            "NextDownID": next_down,
            "maxup": len(up),
            "up1": padded[0],
            "up2": padded[1],
            "up3": padded[2],
            "up4": padded[3],
        },
    }


def _box_poly(x0: float, y0: float, x1: float, y1: float) -> dict[str, Any]:
    return {"type": "Polygon", "coordinates": [[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]]}


def _bowtie_poly(x0: float, y0: float) -> dict[str, Any]:
    return {
        "type": "Polygon",
        "coordinates": [[(x0, y0), (x0 + 0.1, y0 + 0.1), (x0 + 0.1, y0), (x0, y0 + 0.1), (x0, y0)]],
    }


def _pinched_hole_poly(x0: float, y0: float, cell: float = 0.05) -> dict[str, Any]:
    """3x3 grid minus bottom-left [0,1]x[0,1] and center [1,2]x[1,2] hole pinched at (1,1)."""

    def pt(ix: int, iy: int) -> tuple[float, float]:
        return (x0 + ix * cell, y0 + iy * cell)

    return {
        "type": "Polygon",
        "coordinates": [
            [
                pt(1, 0),
                pt(3, 0),
                pt(3, 3),
                pt(0, 3),
                pt(0, 1),
                pt(1, 1),
                pt(1, 2),
                pt(2, 2),
                pt(2, 1),
                pt(1, 1),
                pt(1, 0),
            ]
        ],
    }


def _opposite_winding_fig8_poly(x0: float, y0: float, cell: float = 0.05) -> dict[str, Any]:
    """Two diagonal squares touching at (1,1) traced with opposite winding so raw area cancels."""

    def pt(ix: int, iy: int) -> tuple[float, float]:
        return (x0 + ix * cell, y0 + iy * cell)

    return {
        "type": "Polygon",
        "coordinates": [
            [
                pt(0, 0),
                pt(1, 0),
                pt(1, 1),
                pt(1, 2),
                pt(2, 2),
                pt(2, 1),
                pt(1, 1),
                pt(0, 1),
                pt(0, 0),
            ]
        ],
    }


def _same_winding_fig8_poly(x0: float, y0: float, cell: float = 0.0625) -> dict[str, Any]:
    """Two diagonal squares touching at (1,1) traced with same winding -> splits into 2-part MultiPolygon."""

    def pt(ix: int, iy: int) -> tuple[float, float]:
        return (x0 + ix * cell, y0 + iy * cell)

    return {
        "type": "Polygon",
        "coordinates": [
            [
                pt(0, 0),
                pt(1, 0),
                pt(1, 1),
                pt(2, 1),
                pt(2, 2),
                pt(1, 2),
                pt(1, 1),
                pt(0, 1),
                pt(0, 0),
            ]
        ],
    }


def _build_fixture(tmp_path: Path, *, cycle: bool = False, mismatch: bool = False, cut_target: bool = False) -> Path:
    src = tmp_path / "source"
    src.mkdir(parents=True)

    # Main basin 27: valid 3-reach tree + 1 river with invalid bowtie catchment.
    rivers = [
        _river_feature(
            27000035,
            [(-16.0, 66.0), (-15.9, 65.9)],
            27000015 if not cycle else 27000058,
            (27000015,) if cycle else (),
        ),
        _river_feature(27000058, [(-15.8, 66.0), (-15.9, 65.9)], 27000015, (27000035,) if cycle else ()),
        _river_feature(
            27000015,
            [(-15.9, 65.9), (-15.9, 65.8)],
            0 if not (cycle or cut_target) else (27000035 if cycle else 28000001),
            (27000035,) if mismatch else ((27000058,) if cycle else (27000035, 27000058)),
        ),
        _river_feature(27000099, [(-15.6, 65.9), (-15.6, 65.8)], 0, ()),
    ]
    cats = [
        {"geometry": _box_poly(-16.05, 65.9, -15.95, 66.05), "properties": {"COMID": 27000035, "unitarea": 10.0}},
        {"geometry": _box_poly(-15.85, 65.9, -15.75, 66.05), "properties": {"COMID": 27000058, "unitarea": 10.0}},
        {"geometry": _box_poly(-15.95, 65.75, -15.85, 65.9), "properties": {"COMID": 27000015, "unitarea": 12.0}},
        {"geometry": _bowtie_poly(-15.65, 65.8), "properties": {"COMID": 27000099, "unitarea": 5.0}},
    ]
    hills = [
        {"geometry": _box_poly(-15.5, 65.8, -15.4, 65.9), "properties": {"FID": 9001, "area_km2": 3.5}},
    ]

    riv_side = _write_shapefile(
        src / "riv_pfaf_27_MERIT_Hydro_v07_Basins_v01_bugfix1.shp",
        geom_type="LineString",
        properties_schema=RIV_PROPS,
        features=rivers,
        write_prj=True,
    )
    cat_side = _write_shapefile(
        src / "cat_pfaf_27_MERIT_Hydro_v07_Basins_v01_bugfix1.shp",
        geom_type="Polygon",
        properties_schema=CAT_PROPS,
        features=cats,
        write_prj=False,
    )
    hill_side = _write_shapefile(
        src / "hillslope_pfaf_27.shp",
        geom_type="Polygon",
        properties_schema=HILL_PROPS,
        features=hills,
        write_prj=False,
    )

    inventory: dict[str, Any] = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_27",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_pfaf_27_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": riv_side,
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_pfaf_27_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": cat_side,
                },
            }
        ],
        "coastal_hillslopes": [
            {
                "partition": "hillslope_pfaf_27",
                "status": "supplied",
                "shp": str(src / "hillslope_pfaf_27.shp"),
                "crs_evidence": "paired_river_wgs84_degree_grid",
                "sidecars": hill_side,
            }
        ],
    }
    if cut_target:
        inventory["partitions"].append(
            {
                "partition": "pfaf_28",
                "rivers": {"status": "unavailable", "reason": "Test excluded partition"},
                "catchments": {"status": "unavailable", "reason": "Test excluded partition"},
            }
        )
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(json.dumps(inventory, indent=2))
    return inv_path


def test_merit_vector_production_omission_hillslope_and_independent_verify(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build, verify
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    inv = _build_fixture(tmp_path)
    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    produce_receipt = produce(inv, staging)
    assert produce_receipt["omitted_catchments"] == 1
    assert produce_receipt["delivered_catchments"] == 4
    assert produce_receipt["delivered_rivers"] == 4

    src_digest = hashlib.sha256((staging / "manifest.json").read_bytes()).hexdigest()
    build(staging, spatial, source_manifest_sha256=src_digest)
    assert verify(staging, spatial, source_manifest_sha256=src_digest)["status"] == "passed"

    spatial_manifest = json.loads((spatial / "manifest.json").read_text())
    assert spatial_manifest["catchment_connectivity"] == [
        {
            "level": 0,
            "status": "supplied",
            "derivation": "river_local_associations",
            "completeness": "partial",
            "upstream_cuts_column": "upstream_boundary_cuts",
        }
    ]
    assert not spatial_manifest["has_up_area"]

    cats = {r["COMID"] or r["native_FID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist()}
    assert set(cats) == {27000015, 27000035, 27000058, 9001}
    assert cats[27000015]["catchment_kind"] == "channelized"
    assert cats[27000015]["downstream_boundary_status"] == "source_zero_sentinel"
    assert cats[27000015]["upstream_boundary_cuts"] == 0
    assert cats[27000035]["downstream_boundary_status"] == "internal"
    assert cats[27000035]["upstream_boundary_cuts"] == 0
    assert cats[27000058]["downstream_boundary_status"] == "internal"
    assert cats[27000058]["upstream_boundary_cuts"] == 0
    assert cats[9001]["catchment_kind"] == "coastal_hillslope"
    assert cats[9001]["downstream_boundary_status"] == "unassociated_coastal_hillslope"
    assert cats[9001]["upstream_boundary_cuts"] is None
    assert cats[9001]["native_area_km2"] == pytest.approx(3.5)  # Collision-renamed column preserved!

    with ipc.open_file(str(spatial / "catchments-relationships.arrow")) as reader:
        cat_rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    assert cat_rels[cats[27000015]["id"]]["upstream_ids"] == sorted([cats[27000035]["id"], cats[27000058]["id"]])
    assert cat_rels[cats[27000015]["id"]]["downstream_ids"] == []
    assert cat_rels[cats[27000035]["id"]]["downstream_ids"] == [cats[27000015]["id"]]
    assert cat_rels[cats[27000058]["id"]]["downstream_ids"] == [cats[27000015]["id"]]
    assert cat_rels[cats[9001]["id"]]["upstream_ids"] == []
    assert cat_rels[cats[9001]["id"]]["downstream_ids"] == []

    rivs = {r["COMID"]: r for r in pq.read_table(spatial / "rivers.parquet").to_pylist()}
    assert set(rivs) == {27000015, 27000035, 27000058, 27000099}
    assert "local_catchment_id" not in pq.read_schema(spatial / "rivers.parquet").names
    assert rivs[27000099]["local_catchment_status"] == "omitted_invalid_polygon"
    assert rivs[27000015]["local_catchment_status"] == "supplied"
    assert rivs[27000015]["downstream_boundary_status"] == "source_zero_sentinel"
    assert rivs[27000035]["downstream_boundary_status"] == "internal"

    with ipc.open_file(str(spatial / "rivers-relationships.arrow")) as reader:
        riv_rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    assert riv_rels[rivs[27000099]["id"]]["local_catchment_ids"] == []
    assert riv_rels[rivs[27000015]["id"]]["local_catchment_ids"] == [cats[27000015]["id"]]

    omitted = json.loads((spatial / "omitted-features.json").read_text())
    assert len(omitted) == 1 and omitted[0]["COMID"] == 27000099
    assert "Self-intersection" in omitted[0]["reason"]

    verified = verify_spatial_bundle(inv, spatial)
    assert verified["status"] == "passed"
    assert verified["delivered_rivers"] == 4
    assert verified["delivered_catchments"] == 4
    assert verified["omitted_catchments"] == 1
    assert verified["catchment_flow_edges"] == 2


def test_cut_downstream_edge_recorded_explicitly_not_as_terminal(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    inv = _build_fixture(tmp_path, cut_target=True)
    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    produce(inv, staging)
    build(staging, spatial, source_manifest_sha256=hashlib.sha256((staging / "manifest.json").read_bytes()).hexdigest())

    rivs = {r["COMID"]: r for r in pq.read_table(spatial / "rivers.parquet").to_pylist()}
    assert rivs[27000015]["downstream_boundary_status"] == "cut_omitted_target"
    cats = {r["COMID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist() if r["COMID"] is not None}
    assert cats[27000015]["downstream_boundary_status"] == "cut_missing_downstream_catchment"
    assert cats[27000015]["upstream_boundary_cuts"] == 0
    cuts = json.loads((spatial / "cut-edges.json").read_text())
    assert cuts == [
        {
            "direction": "downstream",
            "origin_comid": 27000015,
            "target_comid": 28000001,
            "origin_partition": "pfaf_27",
            "reason": "target_not_in_delivered_partitions",
        }
    ]
    assert verify_spatial_bundle(inv, spatial)["status"] == "passed"


def test_no_bridging_across_omitted_middle_polygon_and_unmatched_associations(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)

    # Chain A (101) -> B (102, bowtie) -> C (103) -> E (105), unmatched D (104) -> C (103),
    # plus leaf-only headwater gap F (201, bowtie) -> G (202).
    rivers = [
        _river_feature(101, [(0.0, 1.0), (0.0, 0.5)], 102, ()),
        _river_feature(102, [(0.0, 0.5), (0.0, 0.2)], 103, (101,)),
        _river_feature(103, [(0.0, 0.2), (0.0, 0.1)], 105, (102, 104)),
        _river_feature(104, [(0.2, 0.4), (0.0, 0.2)], 103, ()),
        _river_feature(105, [(0.0, 0.1), (0.0, 0.0)], 0, (103,)),
        _river_feature(201, [(1.0, 1.0), (1.0, 0.5)], 202, ()),
        _river_feature(202, [(1.0, 0.5), (1.0, 0.0)], 0, (201,)),
    ]
    cats = [
        {"geometry": _box_poly(-0.1, 0.5, 0.1, 1.0), "properties": {"COMID": 101, "unitarea": 10.0}},
        {"geometry": _bowtie_poly(-0.1, 0.2), "properties": {"COMID": 102, "unitarea": 10.0}},
        {"geometry": _box_poly(-0.1, 0.1, 0.1, 0.2), "properties": {"COMID": 103, "unitarea": 10.0}},
        {"geometry": _box_poly(-0.1, 0.0, 0.1, 0.1), "properties": {"COMID": 105, "unitarea": 10.0}},
        {"geometry": _bowtie_poly(0.9, 0.5), "properties": {"COMID": 201, "unitarea": 5.0}},
        {"geometry": _box_poly(0.9, 0.0, 1.1, 0.5), "properties": {"COMID": 202, "unitarea": 5.0}},
    ]
    riv_side = _write_shapefile(
        src / "riv_pfaf_11.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=rivers, write_prj=True
    )
    cat_side = _write_shapefile(
        src / "cat_pfaf_11.shp", geom_type="Polygon", properties_schema=CAT_PROPS, features=cats, write_prj=False
    )
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(
        json.dumps(
            {
                "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
                "organization": "pfaf_level_02",
                "partitions": [
                    {
                        "partition": "pfaf_11",
                        "rivers": {
                            "status": "supplied",
                            "shp": str(src / "riv_pfaf_11.shp"),
                            "crs_evidence": "supplied_prj",
                            "sidecars": riv_side,
                        },
                        "catchments": {
                            "status": "supplied",
                            "shp": str(src / "cat_pfaf_11.shp"),
                            "crs_evidence": "paired_river_wgs84_degree_grid",
                            "sidecars": cat_side,
                        },
                    }
                ],
            },
            indent=2,
        )
    )

    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    receipt = produce(inv_path, staging)
    assert receipt["catchment_flow_edges"] == 1
    assert receipt["upstream_complete_catchments"] == 1
    build(staging, spatial, source_manifest_sha256=receipt["manifest_sha256"])

    cats_tbl = {r["COMID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist()}
    assert set(cats_tbl) == {101, 103, 105, 202}
    assert cats_tbl[101]["downstream_boundary_status"] == "cut_missing_downstream_catchment"
    assert cats_tbl[101]["upstream_boundary_cuts"] == 0
    assert cats_tbl[103]["downstream_boundary_status"] == "internal"
    assert cats_tbl[103]["upstream_boundary_cuts"] == 2
    assert cats_tbl[105]["downstream_boundary_status"] == "source_zero_sentinel"
    assert cats_tbl[105]["upstream_boundary_cuts"] == 0
    assert cats_tbl[202]["downstream_boundary_status"] == "source_zero_sentinel"
    assert cats_tbl[202]["upstream_boundary_cuts"] == 1

    with ipc.open_file(str(spatial / "catchments-relationships.arrow")) as reader:
        cat_rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    assert cat_rels[cats_tbl[101]["id"]]["downstream_ids"] == []
    assert cat_rels[cats_tbl[103]["id"]]["upstream_ids"] == []
    assert cat_rels[cats_tbl[103]["id"]]["downstream_ids"] == [cats_tbl[105]["id"]]
    assert cat_rels[cats_tbl[105]["id"]]["upstream_ids"] == [cats_tbl[103]["id"]]
    assert cat_rels[cats_tbl[202]["id"]]["upstream_ids"] == []

    rivs_tbl = {r["COMID"]: r for r in pq.read_table(spatial / "rivers.parquet").to_pylist()}
    assert rivs_tbl[102]["local_catchment_status"] == "omitted_invalid_polygon"
    assert rivs_tbl[201]["local_catchment_status"] == "omitted_invalid_polygon"
    assert rivs_tbl[104]["local_catchment_status"] == "not_in_source"

    smeta = json.loads((spatial / "source-metadata.json").read_text())
    cc_counts = smeta["catchment_connectivity"]["counts"]
    assert cc_counts["derived_flow_edges"] == 1
    assert cc_counts["upstream_complete_catchments"] == 1
    assert cc_counts["upstream_incomplete_paired_catchments"] == 3
    assert cc_counts["upstream_incomplete_stranded_delivered_catchments"] == 2
    assert cc_counts["upstream_incomplete_leaf_only_catchments"] == 1
    assert cc_counts["rivers_catchment_not_in_source"] == 1
    assert [r["COMID"] for r in smeta["catchment_connectivity"]["not_in_source_rivers"]] == [104]

    verified = verify_spatial_bundle(inv_path, spatial)
    assert verified["status"] == "passed"
    assert verified["catchment_flow_edges"] == 1


def test_rejects_mismatched_and_cyclic_river_graphs(tmp_path: Path) -> None:
    from scripts.produce_merit import produce

    with pytest.raises(ValueError, match="disagree"):
        produce(_build_fixture(tmp_path / "m", mismatch=True), tmp_path / "out_m")

    with pytest.raises(ValueError, match="cycle"):
        produce(_build_fixture(tmp_path / "c", cycle=True), tmp_path / "out_c")


def test_auxiliary_rasters_bound_through_spatial_bundle(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build, verify
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    inv = _build_fixture(tmp_path)
    raw_dir = tmp_path / "rasters"
    raw_dir.mkdir()
    d_file = raw_dir / "flow_dir.tif"
    a_file = raw_dir / "flow_acc.tif"
    d_file.write_bytes(b"fake-dir-tiff-bytes")
    a_file.write_bytes(b"fake-acc-tiff-bytes")

    aux_spec = {
        "files": {"aux/d8/flow_dir.tif": str(d_file), "aux/d8/flow_acc.tif": str(a_file)},
        "auxiliary": [
            {
                "schema": "hfx.aux.d8_raster.v2",
                "artifacts": {"flow_dir": "aux/d8/flow_dir.tif", "flow_acc": "aux/d8/flow_acc.tif"},
                "metadata": {"crs": "EPSG:4326", "flow_dir_encoding": "merit", "flow_acc_units": "km2"},
            }
        ],
    }
    aux_path = tmp_path / "aux.json"
    aux_path.write_text(json.dumps(aux_spec))

    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    produce(inv, staging, auxiliary_spec_path=aux_path)
    digest = hashlib.sha256((staging / "manifest.json").read_bytes()).hexdigest()
    build(staging, spatial, source_manifest_sha256=digest)
    assert verify(staging, spatial, source_manifest_sha256=digest)["status"] == "passed"
    manifest = json.loads((spatial / "manifest.json").read_text())
    assert manifest["auxiliary"] == aux_spec["auxiliary"]
    assert (spatial / "aux/d8/flow_dir.tif").read_bytes() == b"fake-dir-tiff-bytes"
    assert verify_spatial_bundle(inv, spatial)["status"] == "passed"


def test_preserves_unwrapped_chukotka_antimeridian_longitudes_and_rejects_out_of_range(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)

    riv_side = _write_shapefile(
        src / "riv_pfaf_35.shp",
        geom_type="LineString",
        properties_schema=RIV_PROPS,
        features=[_river_feature(35000001, [(185.05, 66.05), (185.15, 66.15)], 0, ())],
        write_prj=True,
    )
    cat_side = _write_shapefile(
        src / "cat_pfaf_35.shp",
        geom_type="Polygon",
        properties_schema=CAT_PROPS,
        features=[
            {"geometry": _box_poly(185.0, 66.0, 185.2, 66.2), "properties": {"COMID": 35000001, "unitarea": 100.0}},
            {"geometry": _box_poly(200.0, 66.0, 200.1, 66.1), "properties": {"COMID": 35000002, "unitarea": 50.0}},
        ],
        write_prj=False,
    )
    hill_side = _write_shapefile(
        src / "hillslope_35_clean.shp",
        geom_type="Polygon",
        properties_schema=HILL_PROPS,
        features=[
            {"geometry": _box_poly(190.5, 66.0, 191.0, 66.1), "properties": {"FID": 3501, "area_km2": 200.0}},
        ],
        write_prj=True,
    )

    inventory = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_35",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_pfaf_35.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": riv_side,
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_pfaf_35.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": cat_side,
                },
            }
        ],
        "coastal_hillslopes": [
            {
                "partition": "hillslope_35",
                "status": "supplied",
                "shp": str(src / "hillslope_35_clean.shp"),
                "crs_evidence": "supplied_prj",
                "sidecars": hill_side,
            }
        ],
    }
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(json.dumps(inventory, indent=2))

    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    receipt = produce(inv_path, staging)
    assert receipt["delivered_catchments"] == 2
    assert receipt["delivered_rivers"] == 1
    assert receipt["omitted_catchments"] == 1

    build(staging, spatial, source_manifest_sha256=receipt["manifest_sha256"])

    cats = pq.read_table(spatial / "catchments.parquet").to_pylist()
    assert len(cats) == 2
    assert all(row["area_km2"] > 0 for row in cats)
    assert max(row["bbox"]["xmax"] for row in cats) >= 191.0

    rivs = pq.read_table(spatial / "rivers.parquet").to_pylist()
    assert len(rivs) == 1
    assert rivs[0]["bbox"]["xmax"] >= 185.15

    omitted = json.loads((spatial / "omitted-features.json").read_text())
    assert omitted == [
        {
            "domain": "catchments",
            "partition": "pfaf_35",
            "source_namespace": "MERIT-Hydro-v07-Basins-v01-bugfix1/pfaf_35/catchments",
            "source_id": "35000002",
            "COMID": 35000002,
            "reason": "Coordinates outside finite WGS84 longitude/latitude bounds",
        }
    ]
    assert verify_spatial_bundle(inv_path, spatial)["status"] == "passed"


def test_omits_degenerate_self_loop_river_and_preserves_five_tributary_confluence(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)

    tribs = (64072580, 64072915, 64076420, 64076434)
    all_tribs = (*tribs, 64076441)
    rivers = [
        _river_feature(64072579, [(-52.90, -32.55), (-52.89, -32.56)], 0, tribs),
        *(
            _river_feature(c, [(-52.92 + 0.005 * i, -32.53), (-52.90, -32.55)], 64072579, ())
            for i, c in enumerate(all_tribs)
        ),
        {
            "geometry": {"type": "LineString", "coordinates": [(-15.0, 65.0), (-15.0, 65.0)]},
            "properties": {
                **_river_feature(11021766, [(-15.0, 65.0), (-15.0, 65.0)], 11021766, (11021766,))["properties"],
                "maxup": 0,
            },
        },
    ]
    cats = [
        {"geometry": _box_poly(-52.95, -32.60, -52.85, -32.54), "properties": {"COMID": 64072579, "unitarea": 40.0}},
        *(
            {
                "geometry": _box_poly(-52.95 + 0.01 * i, -32.54, -52.945 + 0.01 * i, -32.50),
                "properties": {"COMID": c, "unitarea": 10.0},
            }
            for i, c in enumerate(all_tribs)
        ),
        {"geometry": _box_poly(-15.1, 64.9, -14.9, 65.1), "properties": {"COMID": 11021766, "unitarea": 30.0}},
    ]

    riv_side = _write_shapefile(
        src / "riv_pfaf_64.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=rivers, write_prj=True
    )
    cat_side = _write_shapefile(
        src / "cat_pfaf_64.shp", geom_type="Polygon", properties_schema=CAT_PROPS, features=cats, write_prj=False
    )
    inventory = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_64",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_pfaf_64.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": riv_side,
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_pfaf_64.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": cat_side,
                },
            }
        ],
    }
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(json.dumps(inventory, indent=2))

    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    receipt = produce(inv_path, staging)
    assert receipt["delivered_rivers"] == 6
    assert receipt["omitted_rivers"] == 1
    assert receipt["delivered_catchments"] == 7
    assert receipt["catchment_flow_edges"] == 5
    assert receipt["upstream_complete_catchments"] == 6

    build(staging, spatial, source_manifest_sha256=receipt["manifest_sha256"])

    cats_tbl = {r["COMID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist()}
    assert set(cats_tbl) == {64072579, *all_tribs, 11021766}
    assert cats_tbl[11021766]["downstream_boundary_status"] == "omitted_degenerate_river"
    assert cats_tbl[11021766]["upstream_boundary_cuts"] is None
    assert cats_tbl[64072579]["downstream_boundary_status"] == "source_zero_sentinel"
    assert cats_tbl[64072579]["upstream_boundary_cuts"] == 0

    with ipc.open_file(str(spatial / "catchments-relationships.arrow")) as reader:
        cat_rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    expected_cat_up_ids = sorted(cats_tbl[c]["id"] for c in all_tribs)
    assert cat_rels[cats_tbl[64072579]["id"]]["upstream_ids"] == expected_cat_up_ids
    assert cat_rels[cats_tbl[64076441]["id"]]["downstream_ids"] == [cats_tbl[64072579]["id"]]
    assert cat_rels[cats_tbl[11021766]["id"]]["upstream_ids"] == []
    assert cat_rels[cats_tbl[11021766]["id"]]["downstream_ids"] == []

    omitted = json.loads((spatial / "omitted-features.json").read_text())
    assert len(omitted) == 1
    assert omitted[0]["domain"] == "rivers"
    assert omitted[0]["COMID"] == 11021766
    assert omitted[0]["paired_catchment_status"] == "delivered"
    assert omitted[0]["paired_catchment_hfx_id"] == cats_tbl[11021766]["id"]
    assert omitted[0]["native_connectivity"] == {
        "NextDownID": 11021766,
        "maxup": 0,
        "up1": 11021766,
        "up2": 0,
        "up3": 0,
        "up4": 0,
    }

    smeta = json.loads((spatial / "source-metadata.json").read_text())
    assert len(smeta["up1_to_up4_overflow_confluences"]) == 1
    assert smeta["up1_to_up4_overflow_confluences"][0]["upstream_comid"] == 64076441
    assert smeta["up1_to_up4_overflow_confluences"][0]["downstream_comid"] == 64072579

    rivs = {r["COMID"]: r for r in pq.read_table(spatial / "rivers.parquet").to_pylist()}
    assert rivs[64072579]["maxup"] == 4
    with ipc.open_file(str(spatial / "rivers-relationships.arrow")) as reader:
        rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    expected_up_ids = sorted(rivs[c]["id"] for c in all_tribs)
    assert rels[rivs[64072579]["id"]]["upstream_ids"] == expected_up_ids
    assert rels[rivs[64076441]["id"]]["downstream_ids"] == [rivs[64072579]["id"]]

    assert verify_spatial_bundle(inv_path, spatial)["status"] == "passed"


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ("catchment_edges", "Final catchment flow relationships differ"),
        ("downstream_status", "downstream_boundary_status"),
        ("upstream_cuts", "upstream_boundary_cuts"),
        ("metadata_counts", r"catchment_connectivity\.counts differs"),
    ],
)
def test_verifier_rejects_tampered_catchment_connectivity(tmp_path: Path, mutation: str, match: str) -> None:
    import pyarrow as pa

    from hfx.candidates import fingerprint
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    inv = _build_fixture(tmp_path)
    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    receipt = produce(inv, staging)
    build(staging, spatial, source_manifest_sha256=receipt["manifest_sha256"])

    manifest = json.loads((spatial / "manifest.json").read_text())
    if mutation == "catchment_edges":
        rel_path = spatial / "catchments-relationships.arrow"
        table = ipc.open_file(str(rel_path)).read_all()
        col_idx = table.schema.get_field_index("upstream_ids")
        values = cast(list[list[int]], table.column(col_idx).to_pylist())
        ids = cast(list[int], table.column("id").to_pylist())
        values[0] = [ids[1]] if not values[0] else []
        table = table.set_column(col_idx, table.schema.field(col_idx), pa.array(values, type=pa.list_(pa.int64())))
        with ipc.new_file(str(rel_path), table.schema) as writer:
            writer.write_table(table)
        manifest["artifacts"]["catchments-relationships.arrow"] = fingerprint(rel_path)
    elif mutation in ("downstream_status", "upstream_cuts"):
        cat_path = spatial / "catchments.parquet"
        table = pq.read_table(cat_path)
        if mutation == "downstream_status":
            col_idx = table.schema.get_field_index("downstream_boundary_status")
            vals = table.column(col_idx).to_pylist()
            vals[0] = "internal" if vals[0] != "internal" else "source_zero_sentinel"
            table = table.set_column(col_idx, table.schema.field(col_idx), pa.array(vals, type=pa.string()))
        else:
            col_idx = table.schema.get_field_index("upstream_boundary_cuts")
            vals = table.column(col_idx).to_pylist()
            vals[0] = 99 if vals[0] != 99 else 0
            table = table.set_column(col_idx, table.schema.field(col_idx), pa.array(vals, type=pa.int32()))
        pq.write_table(table, cat_path)
        digest = fingerprint(cat_path)
        manifest["artifacts"]["catchments.parquet"] = digest
        manifest["families"]["catchments"]["data_sha256"] = digest["sha256"]
    else:
        meta_path = spatial / "source-metadata.json"
        smeta = json.loads(meta_path.read_text())
        smeta["catchment_connectivity"]["counts"]["derived_flow_edges"] += 1
        meta_path.write_text(json.dumps(smeta, indent=2) + "\n")
        manifest["artifacts"]["source-metadata.json"] = fingerprint(meta_path)
    (spatial / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    with pytest.raises(ValueError, match=match):
        verify_spatial_bundle(inv, spatial)


def test_dbf_reader_treats_asterisk_fill_as_null(tmp_path: Path) -> None:
    import struct

    from scripts.verify_merit import _read_dbf

    dbf_path = tmp_path / "null_num.dbf"
    header = struct.pack("<BBBBLHH20x", 3, 26, 11, 17, 2, 65, 13)
    field_desc = b"slope".ljust(11, b"\x00") + b"F" + b"\x00" * 4 + bytes([12, 4]) + b"\x00" * 14
    rec_valid = b" " + b"     1.2500 "
    rec_null = b" " + b"************"
    dbf_path.write_bytes(header + field_desc + b"\r" + rec_valid + rec_null)

    assert list(_read_dbf(dbf_path)) == [{"slope": 1.25}, {"slope": None}]


def test_canon_ring_signature_invariant_to_rotation_and_reversal() -> None:
    from scripts.verify_merit import _canon_ring_signature

    ring_ccw = [((0.0, 0.0), (0.5, 1.0), (0.2, 2.0), (2.0, 0.0), (0.0, 0.0))]
    ring_cw_rotated = [((0.2, 2.0), (0.5, 1.0), (0.0, 0.0), (2.0, 0.0), (0.2, 2.0))]
    assert _canon_ring_signature(ring_ccw, is_polygon=True) == _canon_ring_signature(ring_cw_rotated, is_polygon=True)


def test_produce_from_bundle_and_inherited_hillslope_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import sys

    import pyarrow as pa

    from hfx.candidates import fingerprint
    from hfx.spatial_bundle import build
    from scripts.produce_merit import main as produce_main
    from scripts.produce_merit import produce
    from scripts.verify_merit import main as verify_main
    from scripts.verify_merit import verify_spatial_bundle

    inv_path = _build_fixture(tmp_path)
    staging_ref = tmp_path / "staging_ref"
    spatial_ref = tmp_path / "spatial_ref"
    ref_receipt = produce(inv_path, staging_ref)
    build(staging_ref, spatial_ref, source_manifest_sha256=ref_receipt["manifest_sha256"])

    # Build a #281-style base bundle (no derived catchment connectivity columns or artifact).
    base_spatial = tmp_path / "base_spatial"
    base_staging = tmp_path / "base_staging"
    base_receipt = produce(inv_path, base_staging)
    build(base_staging, base_spatial, source_manifest_sha256=base_receipt["manifest_sha256"])
    base_cat_path = base_spatial / "catchments.parquet"
    base_cat = pq.read_table(base_cat_path).drop(["downstream_boundary_status", "upstream_boundary_cuts"])
    pq.write_table(base_cat, base_cat_path)
    (base_spatial / "catchments-relationships.arrow").unlink()
    base_meta_path = base_spatial / "source-metadata.json"
    base_meta = json.loads(base_meta_path.read_text())
    base_meta.pop("catchment_connectivity", None)
    base_meta_path.write_text(json.dumps(base_meta, indent=2) + "\n")
    base_manifest = json.loads((base_spatial / "manifest.json").read_text())
    base_manifest["catchment_connectivity"] = [{"level": 0, "status": "not_supplied"}]
    base_manifest["relationships"].pop("catchments", None)
    base_manifest["artifacts"].pop("catchments-relationships.arrow", None)
    base_manifest["categorical_columns"] = [
        item for item in base_manifest.get("categorical_columns", []) if item.get("artifact") != "catchments.parquet"
    ]
    cat_digest = fingerprint(base_cat_path)
    base_manifest["artifacts"]["catchments.parquet"] = cat_digest
    base_manifest["families"]["catchments"]["data_sha256"] = cat_digest["sha256"]
    base_manifest["artifacts"]["source-metadata.json"] = fingerprint(base_meta_path)
    (base_spatial / "manifest.json").write_text(json.dumps(base_manifest, indent=2) + "\n")
    base_manifest_sha = hashlib.sha256((base_spatial / "manifest.json").read_bytes()).hexdigest()

    # Re-derive full connectivity strictly from the verified base bundle via CLI entrypoint.
    staging2 = tmp_path / "staging2"
    spatial2 = tmp_path / "spatial2"
    produce_report = tmp_path / "produce_report.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "produce_merit.py",
            "--from-spatial-bundle",
            str(base_spatial),
            "--report",
            str(produce_report),
            str(staging2),
        ],
    )
    produce_main()
    derived_receipt = json.loads(produce_report.read_text())
    assert derived_receipt["derived_from_verified_bundle_manifest_sha256"] == base_manifest_sha
    assert derived_receipt["delivered_catchments"] == ref_receipt["delivered_catchments"]
    assert derived_receipt["delivered_rivers"] == ref_receipt["delivered_rivers"]
    build(staging2, spatial2, source_manifest_sha256=derived_receipt["manifest_sha256"])

    # Compare derived bundle directly against fresh shapefile production.
    assert (
        pq.read_table(spatial2 / "catchments.parquet").to_pylist()
        == pq.read_table(spatial_ref / "catchments.parquet").to_pylist()
    )
    assert (
        pq.read_table(spatial2 / "rivers.parquet").to_pylist()
        == pq.read_table(spatial_ref / "rivers.parquet").to_pylist()
    )
    with (
        ipc.open_file(str(spatial2 / "catchments-relationships.arrow")) as r2,
        ipc.open_file(str(spatial_ref / "catchments-relationships.arrow")) as r1,
    ):
        assert r2.read_all().to_pylist() == r1.read_all().to_pylist()

    smeta2 = json.loads((spatial2 / "source-metadata.json").read_text())
    assert smeta2["derived_from_verified_bundle_manifest_sha256"] == base_manifest_sha
    assert (
        smeta2["catchment_connectivity"]
        == json.loads((spatial_ref / "source-metadata.json").read_text())["catchment_connectivity"]
    )

    # Mark hillslope_pfaf_27 as inherited and delete its raw .shp so verification cannot read it from disk.
    inv_data = json.loads(inv_path.read_text())
    for side in inv_data["coastal_hillslopes"][0]["sidecars"].values():
        Path(side["path"]).unlink()
    inv_data["coastal_hillslopes"] = [
        {
            "partition": "hillslope_pfaf_27",
            "status": "inherited_from_verified_bundle",
            "reason": "Simulated quota block on hillslope raw shapefile",
        }
    ]
    inv_inherited = tmp_path / "inventory_inherited.json"
    inv_inherited.write_text(json.dumps(inv_data, indent=2))

    with pytest.raises(ValueError, match="base_spatial_dir"):
        verify_spatial_bundle(inv_inherited, spatial2)

    bad_pfaf_inv = json.loads(inv_inherited.read_text())
    bad_pfaf_inv["partitions"][0]["catchments"]["sidecars"]["cat_pfaf_27_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"][
        "sha256"
    ] = "0" * 64
    bad_pfaf_path = tmp_path / "inventory_bad_pfaf.json"
    bad_pfaf_path.write_text(json.dumps(bad_pfaf_inv, indent=2))
    with pytest.raises(ValueError, match="Supplied sidecars differ"):
        verify_spatial_bundle(bad_pfaf_path, spatial2, base_spatial_dir=base_spatial)

    missing_pfaf_inv = json.loads(inv_inherited.read_text())
    missing_pfaf_inv["partitions"] = []
    missing_pfaf_path = tmp_path / "inventory_missing_pfaf.json"
    missing_pfaf_path.write_text(json.dumps(missing_pfaf_inv, indent=2))
    with pytest.raises(ValueError, match="partition set differs"):
        verify_spatial_bundle(missing_pfaf_path, spatial2, base_spatial_dir=base_spatial)

    verify_report = tmp_path / "verify_report.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "verify_merit.py",
            str(inv_inherited),
            str(spatial2),
            "--base-spatial-dir",
            str(base_spatial),
            "--report",
            str(verify_report),
        ],
    )
    verify_main()
    verified = json.loads(verify_report.read_text())
    assert verified["status"] == "passed"
    assert verified["raw_shapefile_verified_hillslope_partitions"] == []
    assert verified["base_bundle_verified_hillslope_partitions"] == ["hillslope_pfaf_27"]

    # Tampering with an inherited hillslope row in spatial2 must fail verification against base_spatial_dir.
    tampered_path = spatial2 / "catchments.parquet"
    tampered_tbl = pq.read_table(tampered_path)
    fid_idx = tampered_tbl.schema.get_field_index("native_FID")
    fids = tampered_tbl.column(fid_idx).to_pylist()
    kinds = tampered_tbl.column("catchment_kind").to_pylist()
    for i, kind in enumerate(kinds):
        if kind == "coastal_hillslope":
            fids[i] = 999999
    tampered_tbl = tampered_tbl.set_column(
        fid_idx, tampered_tbl.schema.field(fid_idx), pa.array(fids, type=tampered_tbl.schema.field(fid_idx).type)
    )
    pq.write_table(tampered_tbl, tampered_path)
    m2 = json.loads((spatial2 / "manifest.json").read_text())
    tdig = fingerprint(tampered_path)
    m2["artifacts"]["catchments.parquet"] = tdig
    m2["families"]["catchments"]["data_sha256"] = tdig["sha256"]
    (spatial2 / "manifest.json").write_text(json.dumps(m2, indent=2) + "\n")

    with pytest.raises(ValueError, match="Catchment native attributes differ"):
        verify_spatial_bundle(inv_inherited, spatial2, base_spatial_dir=base_spatial)


def _read_repaired_jsonl(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


def test_restores_ring_self_intersection_reconnects_upstream_and_writes_provenance(tmp_path: Path) -> None:
    from shapely import from_wkb

    from hfx.spatial_bundle import build, verify
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)

    # Chain A (101) -> B (102, pinched hole at (0.0625, 0.0625)) -> C (103), with valid neighbor N (104)
    # exactly filling B's 1-cell hole [0.0625, 0.125] x [0.0625, 0.125] (shared border only, 0 interior overlap).
    rivers = [
        _river_feature(101, [(0.02, 0.22), (0.0625, 0.15)], 102, ()),
        _river_feature(102, [(0.0625, 0.15), (0.15, 0.02)], 103, (101,)),
        _river_feature(103, [(0.15, 0.02), (0.25, 0.01)], 0, (102,)),
        _river_feature(104, [(0.09, 0.09), (0.10, 0.10)], 0, ()),
    ]
    cats = [
        {"geometry": _box_poly(0.0, 0.1875, 0.1875, 0.28), "properties": {"COMID": 101, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(0.0, 0.0, cell=0.0625), "properties": {"COMID": 102, "unitarea": 12.0}},
        {"geometry": _box_poly(0.1875, 0.0, 0.28, 0.125), "properties": {"COMID": 103, "unitarea": 10.0}},
        {"geometry": _box_poly(0.0625, 0.0625, 0.125, 0.125), "properties": {"COMID": 104, "unitarea": 3.0}},
    ]
    riv_35 = [_river_feature(3501, [(185.01, 66.01), (185.05, 66.05)], 0, ())]
    cat_35 = [
        {"geometry": _same_winding_fig8_poly(185.0, 66.0, cell=0.0625), "properties": {"COMID": 3501, "unitarea": 96.0}}
    ]

    riv_side = _write_shapefile(
        src / "riv_pfaf_11.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=rivers, write_prj=True
    )
    cat_side = _write_shapefile(
        src / "cat_pfaf_11.shp", geom_type="Polygon", properties_schema=CAT_PROPS, features=cats, write_prj=False
    )
    riv35_side = _write_shapefile(
        src / "riv_pfaf_35.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=riv_35, write_prj=True
    )
    cat35_side = _write_shapefile(
        src / "cat_pfaf_35.shp", geom_type="Polygon", properties_schema=CAT_PROPS, features=cat_35, write_prj=False
    )
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(
        json.dumps(
            {
                "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
                "organization": "pfaf_level_02",
                "partitions": [
                    {
                        "partition": "pfaf_11",
                        "rivers": {
                            "status": "supplied",
                            "shp": str(src / "riv_pfaf_11.shp"),
                            "crs_evidence": "supplied_prj",
                            "sidecars": riv_side,
                        },
                        "catchments": {
                            "status": "supplied",
                            "shp": str(src / "cat_pfaf_11.shp"),
                            "crs_evidence": "paired_river_wgs84_degree_grid",
                            "sidecars": cat_side,
                        },
                    },
                    {
                        "partition": "pfaf_35",
                        "rivers": {
                            "status": "supplied",
                            "shp": str(src / "riv_pfaf_35.shp"),
                            "crs_evidence": "supplied_prj",
                            "sidecars": riv35_side,
                        },
                        "catchments": {
                            "status": "supplied",
                            "shp": str(src / "cat_pfaf_35.shp"),
                            "crs_evidence": "paired_river_wgs84_degree_grid",
                            "sidecars": cat35_side,
                        },
                    },
                ],
            },
            indent=2,
        )
    )

    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    receipt = produce(inv_path, staging)
    assert receipt["delivered_catchments"] == 5
    assert receipt["omitted_catchments"] == 0
    assert receipt["catchment_flow_edges"] == 2
    assert receipt["upstream_complete_catchments"] == 5
    assert (staging / "repaired-features.jsonl.gz").is_file()

    src_digest = hashlib.sha256((staging / "manifest.json").read_bytes()).hexdigest()
    build(staging, spatial, source_manifest_sha256=src_digest)
    assert verify(staging, spatial, source_manifest_sha256=src_digest)["status"] == "passed"

    manifest = json.loads((spatial / "manifest.json").read_text())
    assert "repaired-features.jsonl.gz" in manifest["artifacts"]
    norm_cat_col = [
        item
        for item in manifest["categorical_columns"]
        if item.get("artifact") == "catchments.parquet" and item.get("column") == "geometry_normalization"
    ]
    assert len(norm_cat_col) == 1
    assert {c["code"] for c in norm_cat_col[0]["codes"]} == {"original", "ring_self_intersection_split"}

    cats_tbl = {r["COMID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist()}
    assert set(cats_tbl) == {101, 102, 103, 104, 3501}
    assert cats_tbl[101]["geometry_normalization"] == "original"
    assert cats_tbl[102]["geometry_normalization"] == "ring_self_intersection_split"
    assert cats_tbl[3501]["geometry_normalization"] == "ring_self_intersection_split"
    assert cats_tbl[103]["geometry_normalization"] == "original"
    assert cats_tbl[104]["geometry_normalization"] == "original"
    assert cats_tbl[103]["upstream_boundary_cuts"] == 0

    repaired_geom = from_wkb(cats_tbl[102]["geometry"])
    assert repaired_geom.is_valid
    assert repaired_geom.geom_type == "Polygon"
    assert len(repaired_geom.interiors) == 1

    repaired_mp = from_wkb(cats_tbl[3501]["geometry"])
    assert repaired_mp.is_valid
    assert repaired_mp.geom_type == "MultiPolygon"
    assert len(repaired_mp.geoms) == 2

    with ipc.open_file(str(spatial / "catchments-relationships.arrow")) as reader:
        cat_rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    assert cat_rels[cats_tbl[103]["id"]]["upstream_ids"] == [cats_tbl[102]["id"]]
    assert cat_rels[cats_tbl[102]["id"]]["upstream_ids"] == [cats_tbl[101]["id"]]

    rivs_tbl = {r["COMID"]: r for r in pq.read_table(spatial / "rivers.parquet").to_pylist()}
    assert rivs_tbl[102]["local_catchment_status"] == "supplied"
    with ipc.open_file(str(spatial / "rivers-relationships.arrow")) as reader:
        riv_rels = {r["id"]: r for r in reader.read_all().to_pylist()}
    assert riv_rels[rivs_tbl[102]["id"]]["local_catchment_ids"] == [cats_tbl[102]["id"]]

    repaired_records = {r["COMID"]: r for r in _read_repaired_jsonl(spatial / "repaired-features.jsonl.gz")}
    assert set(repaired_records) == {102, 3501}
    rec = repaired_records[102]
    assert rec["hfx_id"] == cats_tbl[102]["id"]
    assert rec["raw_validity_reason"].startswith("Ring Self-intersection[")
    assert rec["repaired_hole_count"] == 1
    assert rec["checks"]["vertices_exact"] is True
    assert rec["checks"]["edges_exact"] is True
    assert rec["checks"]["geodesic_area_within_tolerance"] is True
    assert rec["checks"]["zero_sibling_interior_overlap"] is True
    orig_geom = from_wkb(base64.b64decode(rec["original_wkb_base64"]))
    assert not orig_geom.is_valid
    assert repaired_records[3501]["repaired_geom_type"] == "MultiPolygon"
    assert repaired_records[3501]["repaired_part_count"] == 2

    smeta = json.loads((spatial / "source-metadata.json").read_text())
    gnorm = smeta["geometry_normalization"]
    assert gnorm["global_counts"]["attempted"] == 2
    assert gnorm["global_counts"]["accepted"] == 2
    assert gnorm["global_counts"]["rejected"] == 0
    assert gnorm["by_partition"]["pfaf_11"]["accepted"] == 1
    assert gnorm["by_partition"]["pfaf_35"]["accepted"] == 1

    verified = verify_spatial_bundle(inv_path, spatial)
    assert verified["status"] == "passed"
    assert verified["repaired_catchments"] == 2


def test_ring_normalization_rejects_area_mismatch_and_sibling_overlaps(tmp_path: Path) -> None:
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)

    # pfaf_11:
    # - 1101: valid neighbor at [0.05, 0.20] x [0.05, 0.20] (overlaps 1102 and cross-partition 1201!)
    # - 1102: pinched hole at (0, 0) overlapping 1101 -> REJECTED
    # - 1103: opposite-winding figure-8 at (1.0, 0.0) -> REJECTED for geodesic area mismatch
    # - 1104 & 1105: two pinched-hole candidates overlapping each other at (2.0, 0.0) / (2.02, 0.02) -> BOTH REJECTED
    # - 1106: pinched hole at (3.0, 0.0) overlapping coastal hillslope 9901 -> REJECTED
    riv_11 = [_river_feature(c, [(0.0, 0.0), (0.01, 0.01)], 0, ()) for c in (1101, 1102, 1103, 1104, 1105, 1106)]
    cat_11 = [
        {"geometry": _box_poly(0.05, 0.05, 0.20, 0.20), "properties": {"COMID": 1101, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(0.0, 0.0, cell=0.05), "properties": {"COMID": 1102, "unitarea": 10.0}},
        {"geometry": _opposite_winding_fig8_poly(1.0, 0.0, cell=0.05), "properties": {"COMID": 1103, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(2.0, 0.0, cell=0.05), "properties": {"COMID": 1104, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(2.02, 0.02, cell=0.05), "properties": {"COMID": 1105, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(3.0, 0.0, cell=0.05), "properties": {"COMID": 1106, "unitarea": 10.0}},
    ]

    # pfaf_12:
    # - 1201: pinched hole at (0.10, 0.10) overlapping cross-partition valid neighbor 1101 -> REJECTED
    riv_12 = [_river_feature(1201, [(0.1, 0.1), (0.12, 0.12)], 0, ())]
    cat_12 = [{"geometry": _pinched_hole_poly(0.10, 0.10, cell=0.05), "properties": {"COMID": 1201, "unitarea": 10.0}}]

    # pfaf_35 (unwrapped antimeridian >180 E) & pfaf_36 (wrapped -180 W):
    # - 3501: pinched hole around 185.00..185.15 overlapping wrapped neighbor 3601 at -174.98..-174.90 -> REJECTED
    # - 3502: unwrapped valid at [186.02, 186.12] overlapping wrapped candidate 3602 at -174.0..-173.85 -> 3602 REJECTED
    # - 3503 & 3603: two pinched-hole candidates overlapping across antimeridian (187.0 & -172.98) -> BOTH REJECTED
    riv_35 = [
        _river_feature(3501, [(185.01, 66.01), (185.02, 66.02)], 0, ()),
        _river_feature(3502, [(186.05, 66.05), (186.08, 66.08)], 0, ()),
        _river_feature(3503, [(187.01, 66.01), (187.05, 66.05)], 0, ()),
    ]
    cat_35 = [
        {"geometry": _pinched_hole_poly(185.0, 66.0, cell=0.05), "properties": {"COMID": 3501, "unitarea": 10.0}},
        {"geometry": _box_poly(186.02, 66.02, 186.12, 66.12), "properties": {"COMID": 3502, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(187.0, 66.0, cell=0.05), "properties": {"COMID": 3503, "unitarea": 10.0}},
    ]
    riv_36 = [
        _river_feature(3601, [(-174.95, 66.05), (-174.94, 66.06)], 0, ()),
        _river_feature(3602, [(-173.95, 66.05), (-173.90, 66.08)], 0, ()),
        _river_feature(3603, [(-172.95, 66.05), (-172.90, 66.08)], 0, ()),
    ]
    cat_36 = [
        {"geometry": _box_poly(-174.98, 66.02, -174.90, 66.12), "properties": {"COMID": 3601, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(-174.0, 66.0, cell=0.05), "properties": {"COMID": 3602, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(-172.98, 66.02, cell=0.05), "properties": {"COMID": 3603, "unitarea": 10.0}},
    ]

    # Coastal hillslope overlapping 1106
    hills = [{"geometry": _box_poly(3.01, 0.01, 3.04, 0.12), "properties": {"FID": 9901, "area_km2": 3.0}}]

    inv = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_11",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_11.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        src / "riv_11.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=riv_11
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_11.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        src / "cat_11.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cat_11,
                        write_prj=False,
                    ),
                },
            },
            {
                "partition": "pfaf_12",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_12.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        src / "riv_12.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=riv_12
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_12.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        src / "cat_12.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cat_12,
                        write_prj=False,
                    ),
                },
            },
            {
                "partition": "pfaf_35",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_35.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        src / "riv_35.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=riv_35
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_35.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        src / "cat_35.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cat_35,
                        write_prj=False,
                    ),
                },
            },
            {
                "partition": "pfaf_36",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_36.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        src / "riv_36.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=riv_36
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_36.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        src / "cat_36.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cat_36,
                        write_prj=False,
                    ),
                },
            },
        ],
        "coastal_hillslopes": [
            {
                "partition": "hillslope_11",
                "status": "supplied",
                "shp": str(src / "hill_11.shp"),
                "crs_evidence": "supplied_prj",
                "sidecars": _write_shapefile(
                    src / "hill_11.shp", geom_type="Polygon", properties_schema=HILL_PROPS, features=hills
                ),
            }
        ],
    }
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(json.dumps(inv, indent=2))

    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    receipt = produce(inv_path, staging)
    assert receipt["delivered_catchments"] == 4  # 1101, 3502, 3601, 9901 only
    assert receipt["omitted_catchments"] == 10

    from hfx.spatial_bundle import build

    build(staging, spatial, source_manifest_sha256=receipt["manifest_sha256"])

    omitted = {item["COMID"]: item for item in json.loads((spatial / "omitted-features.json").read_text())}
    assert set(omitted) == {1102, 1103, 1104, 1105, 1106, 1201, 3501, 3503, 3602, 3603}
    assert omitted[1103]["repair_rejection"] == "geodesic_area_mismatch"
    for comid in (1102, 1104, 1105, 1106, 1201, 3501, 3503, 3602, 3603):
        assert "overlap" in omitted[comid]["repair_rejection"]

    smeta = json.loads((spatial / "source-metadata.json").read_text())
    gnorm = smeta["geometry_normalization"]
    assert gnorm["global_counts"]["attempted"] == 10
    assert gnorm["global_counts"]["accepted"] == 0
    assert gnorm["global_counts"]["rejected"] == 10
    assert _read_repaired_jsonl(spatial / "repaired-features.jsonl.gz") == []
    assert verify_spatial_bundle(inv_path, spatial)["status"] == "passed"


def test_produce_hybrid_with_inherited_hillslopes_and_restored_catchments(tmp_path: Path) -> None:
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    # 1. Build base bundle containing a coastal hillslope at [3.01, 0.01, 3.04, 0.12].
    base_src = tmp_path / "base_src"
    base_src.mkdir(parents=True)
    riv_base = [
        _river_feature(1101, [(0.02, 0.20), (0.05, 0.12)], 0, ()),
        _river_feature(1102, [(3.1, 0.1), (3.15, 0.15)], 0, ()),
    ]
    cat_base = [
        {"geometry": _pinched_hole_poly(0.0, 0.0, cell=0.05), "properties": {"COMID": 1101, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(3.0, 0.0, cell=0.05), "properties": {"COMID": 1102, "unitarea": 10.0}},
    ]
    hill_base = [{"geometry": _box_poly(3.01, 0.01, 3.04, 0.12), "properties": {"FID": 9901, "area_km2": 3.0}}]

    base_inv = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_11",
                "rivers": {
                    "status": "supplied",
                    "shp": str(base_src / "riv_11.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        base_src / "riv_11.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=riv_base
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(base_src / "cat_11.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        base_src / "cat_11.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cat_base,
                        write_prj=False,
                    ),
                },
            }
        ],
        "coastal_hillslopes": [
            {
                "partition": "hillslope_11",
                "status": "supplied",
                "shp": str(base_src / "hill_11.shp"),
                "crs_evidence": "supplied_prj",
                "sidecars": _write_shapefile(
                    base_src / "hill_11.shp", geom_type="Polygon", properties_schema=HILL_PROPS, features=hill_base
                ),
            }
        ],
    }
    base_inv_path = tmp_path / "base_inv.json"
    base_inv_path.write_text(json.dumps(base_inv, indent=2))
    base_staging = tmp_path / "base_staging"
    base_spatial = tmp_path / "base_spatial"
    base_rec = produce(base_inv_path, base_staging, normalize_rings=False)
    build(base_staging, base_spatial, source_manifest_sha256=base_rec["manifest_sha256"])

    # Convert base_spatial into an exact #285-schema bundle (without geometry_normalization or
    # repaired-features.jsonl.gz) so hybrid inheritance exercises the actual #285 schema.
    from hfx.candidates import fingerprint as hfx_fingerprint

    base_cat_path = base_spatial / "catchments.parquet"
    pq.write_table(pq.read_table(base_cat_path).drop(["geometry_normalization"]), base_cat_path)
    (base_spatial / "repaired-features.jsonl.gz").unlink()
    bm = json.loads((base_spatial / "manifest.json").read_text())
    bm["artifacts"].pop("repaired-features.jsonl.gz", None)
    bm["categorical_columns"] = [
        item for item in bm.get("categorical_columns", []) if item.get("column") != "geometry_normalization"
    ]
    base_cat_dig = hfx_fingerprint(base_cat_path)
    bm["artifacts"]["catchments.parquet"] = base_cat_dig
    bm["families"]["catchments"]["data_sha256"] = base_cat_dig["sha256"]
    (base_spatial / "manifest.json").write_text(json.dumps(bm, indent=2) + "\n")

    # 2. Relocate pfaf_11 sidecars to a new directory and delete raw hillslope files to prove
    # path-independent sidecar matching and inherited hillslope overlap checking.
    new_src = tmp_path / "relocated_src"
    new_src.mkdir(parents=True)
    for part_dom in ("rivers", "catchments"):
        entry = base_inv["partitions"][0][part_dom]
        new_sidecars = {}
        for sname, meta in entry["sidecars"].items():
            dst = new_src / sname
            dst.write_bytes(Path(meta["path"]).read_bytes())
            new_sidecars[sname] = {**meta, "path": str(dst)}
        entry["shp"] = str(new_src / Path(entry["shp"]).name)
        entry["sidecars"] = new_sidecars
    for meta in base_inv["coastal_hillslopes"][0]["sidecars"].values():
        Path(meta["path"]).unlink()
    base_inv["coastal_hillslopes"] = [
        {"partition": "hillslope_11", "status": "inherited_from_verified_bundle", "reason": "Quota test"}
    ]
    hybrid_inv_path = tmp_path / "hybrid_inv.json"
    hybrid_inv_path.write_text(json.dumps(base_inv, indent=2))

    staging = tmp_path / "hybrid_staging"
    spatial = tmp_path / "hybrid_spatial"
    rec = produce(hybrid_inv_path, staging, base_spatial_dir=base_spatial)
    assert rec["delivered_catchments"] == 2  # restored 1101 + inherited hillslope 9901 (1102 overlaps 9901!)
    assert rec["omitted_catchments"] == 1
    build(staging, spatial, source_manifest_sha256=rec["manifest_sha256"])

    verified = verify_spatial_bundle(hybrid_inv_path, spatial, base_spatial_dir=base_spatial)
    assert verified["status"] == "passed"
    assert verified["repaired_catchments"] == 1
    assert verified["base_bundle_verified_hillslope_partitions"] == ["hillslope_11"]


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ("tamper_repaired_geom", "Catchment geometry|edge|vertex|hole"),
        ("tamper_repaired_jsonl_wkb", "repaired-features|original_wkb|raw"),
        ("tamper_geom_norm_col", "geometry_normalization"),
        ("silent_valid_omission", "Valid catchment unexpectedly omitted"),
        ("tamper_norm_metadata", "geometry_normalization"),
    ],
)
def test_verifier_rejects_tampered_repair_provenance_and_silent_omissions(
    tmp_path: Path, mutation: str, match: str
) -> None:
    import pyarrow as pa
    from shapely import from_wkb, to_wkb
    from shapely.geometry import Polygon

    from hfx.candidates import fingerprint
    from hfx.spatial_bundle import build
    from scripts.produce_merit import produce
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)
    rivers = [
        _river_feature(101, [(0.02, 0.20), (0.05, 0.12)], 102, ()),
        _river_feature(102, [(0.05, 0.12), (0.12, 0.02)], 0, (101,)),
    ]
    cats = [
        {"geometry": _box_poly(0.0, 0.16, 0.15, 0.25), "properties": {"COMID": 101, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(0.0, 0.0, cell=0.05), "properties": {"COMID": 102, "unitarea": 12.0}},
    ]
    inv = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_11",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_11.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        src / "riv_11.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=rivers
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_11.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        src / "cat_11.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cats,
                        write_prj=False,
                    ),
                },
            }
        ],
    }
    inv_path = tmp_path / "inv.json"
    inv_path.write_text(json.dumps(inv, indent=2))
    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    rec = produce(inv_path, staging)
    build(staging, spatial, source_manifest_sha256=rec["manifest_sha256"])

    manifest = json.loads((spatial / "manifest.json").read_text())
    if mutation == "tamper_repaired_geom":
        cat_path = spatial / "catchments.parquet"
        tbl = pq.read_table(cat_path)
        rows = tbl.to_pylist()
        for r in rows:
            if r["COMID"] == 102:
                # Fill the hole (dropping interior ring!) while remaining valid.
                g = from_wkb(r["geometry"])
                r["geometry"] = to_wkb(Polygon(g.exterior))
        new_tbl = pa.Table.from_pylist(rows, schema=tbl.schema)
        pq.write_table(new_tbl, cat_path)
        dig = fingerprint(cat_path)
        manifest["artifacts"]["catchments.parquet"] = dig
        manifest["families"]["catchments"]["data_sha256"] = dig["sha256"]
    elif mutation == "tamper_repaired_jsonl_wkb":
        rep_path = spatial / "repaired-features.jsonl.gz"
        recs = _read_repaired_jsonl(rep_path)
        recs[0]["original_wkb_base64"] = base64.b64encode(
            to_wkb(Polygon([(0.0, 0.0), (0.1, 0.0), (0.1, 0.1), (0.0, 0.1), (0.0, 0.0)]))
        ).decode("ascii")
        with gzip.GzipFile(rep_path, mode="wb", mtime=0) as gz:
            gz.write((json.dumps(recs[0]) + "\n").encode("utf-8"))
        manifest["artifacts"]["repaired-features.jsonl.gz"] = fingerprint(rep_path)
    elif mutation == "tamper_geom_norm_col":
        cat_path = spatial / "catchments.parquet"
        tbl = pq.read_table(cat_path)
        col_idx = tbl.schema.get_field_index("geometry_normalization")
        vals = ["original" for _ in range(tbl.num_rows)]
        tbl = tbl.set_column(col_idx, tbl.schema.field(col_idx), pa.array(vals, type=pa.string()))
        pq.write_table(tbl, cat_path)
        dig = fingerprint(cat_path)
        manifest["artifacts"]["catchments.parquet"] = dig
        manifest["families"]["catchments"]["data_sha256"] = dig["sha256"]
    elif mutation == "silent_valid_omission":
        omit_path = spatial / "omitted-features.json"
        omitted = json.loads(omit_path.read_text())
        omitted.append(
            {
                "domain": "catchments",
                "partition": "pfaf_11",
                "source_namespace": "MERIT-Hydro-v07-Basins-v01-bugfix1/pfaf_11/catchments",
                "source_id": "101",
                "COMID": 101,
                "reason": "Simulated silent omission",
            }
        )
        omit_path.write_text(json.dumps(omitted, indent=2) + "\n")
        manifest["artifacts"]["omitted-features.json"] = fingerprint(omit_path)
    else:
        meta_path = spatial / "source-metadata.json"
        smeta = json.loads(meta_path.read_text())
        smeta["geometry_normalization"]["global_counts"]["accepted"] += 1
        meta_path.write_text(json.dumps(smeta, indent=2) + "\n")
        manifest["artifacts"]["source-metadata.json"] = fingerprint(meta_path)

    (spatial / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with pytest.raises(ValueError, match=match):
        verify_spatial_bundle(inv_path, spatial)


def test_verifier_step4_rejects_overlapping_repaired_catchment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import scripts.produce_merit as pm
    from hfx.spatial_bundle import build
    from scripts.verify_merit import verify_spatial_bundle

    src = tmp_path / "source"
    src.mkdir(parents=True)
    rivers = [
        _river_feature(101, [(0.06, 0.06), (0.10, 0.10)], 0, ()),
        _river_feature(102, [(0.15, 0.05), (0.18, 0.08)], 0, ()),
    ]
    cats = [
        {"geometry": _box_poly(0.02, 0.02, 0.14, 0.14), "properties": {"COMID": 101, "unitarea": 10.0}},
        {"geometry": _pinched_hole_poly(0.0, 0.0, cell=0.0625), "properties": {"COMID": 102, "unitarea": 12.0}},
    ]
    inv = {
        "release": "MERIT-Hydro-v07-Basins-v01-bugfix1",
        "organization": "pfaf_level_02",
        "partitions": [
            {
                "partition": "pfaf_11",
                "rivers": {
                    "status": "supplied",
                    "shp": str(src / "riv_11.shp"),
                    "crs_evidence": "supplied_prj",
                    "sidecars": _write_shapefile(
                        src / "riv_11.shp", geom_type="LineString", properties_schema=RIV_PROPS, features=rivers
                    ),
                },
                "catchments": {
                    "status": "supplied",
                    "shp": str(src / "cat_11.shp"),
                    "crs_evidence": "paired_river_wgs84_degree_grid",
                    "sidecars": _write_shapefile(
                        src / "cat_11.shp",
                        geom_type="Polygon",
                        properties_schema=CAT_PROPS,
                        features=cats,
                        write_prj=False,
                    ),
                },
            }
        ],
    }
    inv_path = tmp_path / "inv.json"
    inv_path.write_text(json.dumps(inv, indent=2))

    monkeypatch.setattr(
        pm,
        "_validate_candidate_overlaps",
        lambda cands, _valid: (
            {},
            dict.fromkeys(range(len(cands)), 0),
            {
                "same_partition_pairs": 0,
                "cross_partition_pairs": 0,
                "coastal_hillslope_pairs": 0,
                "candidate_candidate_pairs": 0,
                "candidate_interior_overlaps_detected": 0,
                "accepted_interior_overlaps": 0,
            },
        ),
    )
    staging = tmp_path / "staging"
    spatial = tmp_path / "spatial"
    rec = pm.produce(inv_path, staging)
    build(staging, spatial, source_manifest_sha256=rec["manifest_sha256"])

    with pytest.raises(ValueError, match="sibling interior overlap"):
        verify_spatial_bundle(inv_path, spatial)
