"""Test-first verification for MERIT Basins vector normalization and independent binary shapefile readback."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

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
    assert spatial_manifest["catchment_connectivity"] == [{"level": 0, "status": "not_supplied"}]
    assert not spatial_manifest["has_up_area"]

    cats = {r["COMID"] or r["native_FID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist()}
    assert set(cats) == {27000015, 27000035, 27000058, 9001}
    assert cats[27000015]["catchment_kind"] == "channelized"
    assert cats[9001]["catchment_kind"] == "coastal_hillslope"
    assert cats[9001]["native_area_km2"] == pytest.approx(3.5)  # Collision-renamed column preserved!

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
    rivers = [
        _river_feature(64072579, [(-52.90, -32.55), (-52.89, -32.56)], 0, tribs),
        *(
            _river_feature(c, [(-52.92 + 0.005 * i, -32.53), (-52.90, -32.55)], 64072579, ())
            for i, c in enumerate((*tribs, 64076441))
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
    assert receipt["delivered_catchments"] == 1

    build(staging, spatial, source_manifest_sha256=receipt["manifest_sha256"])

    cats_tbl = {r["COMID"]: r for r in pq.read_table(spatial / "catchments.parquet").to_pylist()}
    assert set(cats_tbl) == {11021766}

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
    expected_up_ids = sorted(rivs[c]["id"] for c in (*tribs, 64076441))
    assert rels[rivs[64072579]["id"]]["upstream_ids"] == expected_up_ids
    assert rels[rivs[64076441]["id"]]["downstream_ids"] == [rivs[64072579]["id"]]

    assert verify_spatial_bundle(inv_path, spatial)["status"] == "passed"


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
