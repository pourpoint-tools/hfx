"""Test-first verification for MERIT Hydro raster normalization and layout proof."""

from __future__ import annotations

from pathlib import Path

import pytest

for dependency in ("numpy", "rasterio", "rio_cogeo", "tifffile"):
    pytest.importorskip(dependency)

import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from rasterio.enums import MaskFlags  # noqa: E402
from rasterio.transform import Affine  # noqa: E402
from rio_cogeo.cogeo import cog_validate  # noqa: E402


def _write_merit_source_tile(
    dir_path: Path,
    upa_path: Path,
    *,
    transform: Affine,
    dir_data: np.ndarray,
    upa_data: np.ndarray,
    dir_nodata: float | None = None,
    upa_nodata: float | None = -9999.0,
    upa_transform: Affine | None = None,
) -> None:
    common = {
        "driver": "GTiff",
        "height": dir_data.shape[0],
        "width": dir_data.shape[1],
        "count": 1,
        "crs": "EPSG:4326",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
        "compress": "deflate",
    }
    with rasterio.open(dir_path, "w", dtype="uint8", nodata=dir_nodata, transform=transform, **common) as dst:
        dst.write(dir_data, 1)
        dst.update_tags(AREA_OR_POINT="Area")
    with rasterio.open(
        upa_path, "w", dtype="float32", nodata=upa_nodata, transform=upa_transform or transform, **common
    ) as dst:
        dst.write(upa_data, 1)
        dst.update_tags(AREA_OR_POINT="Area")


def _sample_arrays(shape: tuple[int, int] = (800, 800)) -> tuple[np.ndarray, np.ndarray]:
    direction = np.full(shape, 247, dtype=np.uint8)
    direction[10:20, 10:20] = np.array(
        [[1, 2, 4, 8, 16, 32, 64, 128, 0, 255]] * 10,
        dtype=np.uint8,
    )
    bits = np.full(shape, np.float32(-9999.0).view(np.uint32), dtype=np.uint32)
    bits[10:20, 10:20] = np.float32(12.5).view(np.uint32)
    # Preserve negative zero float32 bit pattern cleanly.
    bits[10, 10] = np.uint32(0x80000000)
    return direction, bits.view(np.float32)


def test_normalizes_two_tiles_preserves_terminals_float_bits_and_verifies(tmp_path: Path) -> None:
    from hfx.merit_raster import TilePair, build_merit_rasters
    from scripts.verify_merit_rasters import verify_merit_raster_pair

    dx = 1.0 / 1200.0
    t0 = Affine(dx, 0.0, -20.0 - dx * 0.5, 0.0, -dx, 70.0 - dx * 0.5)
    t1 = Affine(dx, 0.0, -20.0 - dx * 0.5 + 800 * dx, 0.0, -dx, 70.0 - dx * 0.5)

    d0, u0 = _sample_arrays()
    d1, u1 = _sample_arrays()
    d1[30, 30] = 255
    u1[30, 30] = np.float32(0.0081)

    p_d0, p_u0 = tmp_path / "tile0_dir.tif", tmp_path / "tile0_upa.tif"
    p_d1, p_u1 = tmp_path / "tile1_dir.tif", tmp_path / "tile1_upa.tif"
    _write_merit_source_tile(p_d0, p_u0, transform=t0, dir_data=d0, upa_data=u0)
    _write_merit_source_tile(p_d1, p_u1, transform=t1, dir_data=d1, upa_data=u1)

    out_dir = tmp_path / "out"
    receipt = build_merit_rasters(
        [
            TilePair(name="tile0", dir_path=p_d0, upa_path=p_u0),
            TilePair(name="tile1", dir_path=p_d1, upa_path=p_u1),
        ],
        out_dir,
        complete_global=True,
        tile_size=400,
    )

    flow_dir = out_dir / "flow_dir.tif"
    flow_acc = out_dir / "flow_acc.tif"
    assert flow_dir.is_file() and flow_acc.is_file()

    with rasterio.open(flow_dir) as dir_ds, rasterio.open(flow_acc) as acc_ds:
        assert str(dir_ds.crs) == "EPSG:4326"
        assert str(acc_ds.crs) == "EPSG:4326"
        assert dir_ds.nodata == 247.0
        assert acc_ds.nodata == -9999.0
        assert dir_ds.mask_flag_enums == ([MaskFlags.nodata],)
        assert acc_ds.mask_flag_enums == ([MaskFlags.nodata],)
        assert dir_ds.width == 1600 and dir_ds.height == 800
        assert dir_ds.transform == t0 and acc_ds.transform == t0
        assert dir_ds.block_shapes == [(400, 400)]
        assert acc_ds.block_shapes == [(400, 400)]

        out_dir_arr = dir_ds.read(1)
        out_dir_mask = dir_ds.read_masks(1)
        assert out_dir_arr[10, 18] == 0 and out_dir_mask[10, 18] == 255  # river mouth
        assert out_dir_arr[10, 19] == 255 and out_dir_mask[10, 19] == 255  # depression
        assert out_dir_arr[0, 0] == 247 and out_dir_mask[0, 0] == 0  # ocean nodata

        out_acc_arr = acc_ds.read(1)
        # Signed negative zero bit pattern preserved bit-for-bit.
        assert out_acc_arr.view(np.uint32)[10, 10] == np.uint32(0x80000000)

    for path in (flow_dir, flow_acc):
        valid, errors, warnings = cog_validate(str(path), strict=False, quiet=True)
        assert valid and not errors
        assert warnings == ["The file is greater than 512xH or 512xW, it is recommended to include internal overviews"]

    assert receipt["counts"]["255"] >= 11
    assert receipt["counts"]["0"] >= 20
    assert len(receipt["depression_witnesses"]) >= 1

    verified = verify_merit_raster_pair(
        flow_dir,
        flow_acc,
        [
            TilePair(name="tile0", dir_path=p_d0, upa_path=p_u0),
            TilePair(name="tile1", dir_path=p_d1, upa_path=p_u1),
        ],
        tile_size=400,
    )
    assert verified["status"] == "passed"
    assert verified["dir_bits_equal"] is True
    assert verified["upa_bits_equal"] is True
    assert verified["masks_equal"] is True


def test_rejects_unexpected_direction_code(tmp_path: Path) -> None:
    from hfx.merit_raster import TilePair, build_merit_rasters

    dx = 1.0 / 1200.0
    t = Affine(dx, 0.0, -20.0 - dx * 0.5, 0.0, -dx, 70.0 - dx * 0.5)
    d, u = _sample_arrays((400, 400))
    d[5, 5] = 3  # Invalid D8 code
    u[5, 5] = np.float32(1.0)
    p_d, p_u = tmp_path / "d.tif", tmp_path / "u.tif"
    _write_merit_source_tile(p_d, p_u, transform=t, dir_data=d, upa_data=u)

    with pytest.raises(ValueError, match="direction code"):
        build_merit_rasters([TilePair(name="bad", dir_path=p_d, upa_path=p_u)], tmp_path / "out", tile_size=400)


def test_rejects_ocean_mask_disagreement(tmp_path: Path) -> None:
    from hfx.merit_raster import TilePair, build_merit_rasters

    dx = 1.0 / 1200.0
    t = Affine(dx, 0.0, -20.0 - dx * 0.5, 0.0, -dx, 70.0 - dx * 0.5)
    d, u = _sample_arrays((400, 400))
    d[5, 5] = 247
    u[5, 5] = np.float32(1.0)  # valid upa where dir says ocean
    p_d, p_u = tmp_path / "d.tif", tmp_path / "u.tif"
    _write_merit_source_tile(p_d, p_u, transform=t, dir_data=d, upa_data=u)

    with pytest.raises(ValueError, match="247|mask"):
        build_merit_rasters([TilePair(name="mismatch", dir_path=p_d, upa_path=p_u)], tmp_path / "out", tile_size=400)


def test_rejects_misaligned_tile_grid(tmp_path: Path) -> None:
    from hfx.merit_raster import TilePair, build_merit_rasters

    dx = 1.0 / 1200.0
    t0 = Affine(dx, 0.0, -20.0 - dx * 0.5, 0.0, -dx, 70.0 - dx * 0.5)
    t1 = Affine(dx, 0.0, -20.0 - dx * 0.25 + 400 * dx, 0.0, -dx, 70.0 - dx * 0.5)
    d, u = _sample_arrays((400, 400))
    p_d0, p_u0 = tmp_path / "d0.tif", tmp_path / "u0.tif"
    p_d1, p_u1 = tmp_path / "d1.tif", tmp_path / "u1.tif"
    _write_merit_source_tile(p_d0, p_u0, transform=t0, dir_data=d, upa_data=u)
    _write_merit_source_tile(p_d1, p_u1, transform=t1, dir_data=d, upa_data=u)

    with pytest.raises(ValueError, match="(?i)grid|align"):
        build_merit_rasters(
            [TilePair(name="a", dir_path=p_d0, upa_path=p_u0), TilePair(name="b", dir_path=p_d1, upa_path=p_u1)],
            tmp_path / "out",
            complete_global=True,
            tile_size=400,
        )


def test_partial_coverage_emits_only_acquired_tiles_without_ocean_fill(tmp_path: Path) -> None:
    from hfx.merit_raster import TilePair, build_merit_rasters
    from scripts.verify_merit_rasters import verify_merit_raster_pair

    dx = 1.0 / 1200.0
    # Leave a gap between tile0 and tile2; partial mode must NOT fabricate ocean in between.
    t0 = Affine(dx, 0.0, -20.0 - dx * 0.5, 0.0, -dx, 70.0 - dx * 0.5)
    t2 = Affine(dx, 0.0, -20.0 - dx * 0.5 + 800 * dx, 0.0, -dx, 70.0 - dx * 0.5)
    d, u = _sample_arrays((400, 400))
    p_d0, p_u0 = tmp_path / "d0.tif", tmp_path / "u0.tif"
    p_d2, p_u2 = tmp_path / "d2.tif", tmp_path / "u2.tif"
    _write_merit_source_tile(p_d0, p_u0, transform=t0, dir_data=d, upa_data=u)
    _write_merit_source_tile(p_d2, p_u2, transform=t2, dir_data=d, upa_data=u)

    receipt = build_merit_rasters(
        [TilePair(name="t0", dir_path=p_d0, upa_path=p_u0), TilePair(name="t2", dir_path=p_d2, upa_path=p_u2)],
        tmp_path / "out",
        complete_global=False,
        tile_size=400,
    )
    assert len(receipt["pairs"]) == 2
    for entry in receipt["pairs"]:
        verified = verify_merit_raster_pair(
            tmp_path / "out" / entry["flow_dir"],
            tmp_path / "out" / entry["flow_acc"],
            [
                TilePair(
                    name=entry["name"],
                    dir_path=tmp_path / f"{'d0' if entry['name'] == 't0' else 'd2'}.tif",
                    upa_path=tmp_path / f"{'u0' if entry['name'] == 't0' else 'u2'}.tif",
                )
            ],
            tile_size=400,
        )
        assert verified["status"] == "passed"


def test_tar_archive_roundtrip_and_verification(tmp_path: Path) -> None:
    import tarfile

    from hfx.merit_raster import TarArchivePair, build_merit_rasters
    from scripts.verify_merit_rasters import verify_merit_raster_pair

    dx = 1.0 / 1200.0
    t = Affine(dx, 0.0, -20.0 - dx * 0.5, 0.0, -dx, 70.0 - dx * 0.5)
    d, u = _sample_arrays((400, 400))
    p_d, p_u = tmp_path / "n65w020_dir.tif", tmp_path / "n65w020_upa.tif"
    _write_merit_source_tile(p_d, p_u, transform=t, dir_data=d, upa_data=u)

    dir_tar, upa_tar = tmp_path / "dir_n60w030.tar", tmp_path / "upa_n60w030.tar"
    with tarfile.open(dir_tar, "w") as tf:
        tf.add(p_d, arcname="dir_n60w030/n65w020_dir.tif")
    with tarfile.open(upa_tar, "w") as tf:
        tf.add(p_u, arcname="upa_n60w030/n65w020_upa.tif")

    sources = [TarArchivePair(block="n60w030", dir_tar=dir_tar, upa_tar=upa_tar)]
    out_dir = tmp_path / "out"
    receipt = build_merit_rasters(sources, out_dir, complete_global=True, tile_size=400)
    assert receipt["counts"]["255"] == 10
    verified = verify_merit_raster_pair(out_dir / "flow_dir.tif", out_dir / "flow_acc.tif", sources, tile_size=400)
    assert verified["status"] == "passed"
    assert verified["verified_source_tiles"] == 1


def test_accepts_float64_ulp_transform_rounding_and_rejects_subpixel_shift(tmp_path: Path) -> None:
    from hfx.merit_raster import TilePair, build_merit_rasters
    from scripts.verify_merit_rasters import verify_merit_raster_pair

    dx = 1.0 / 1200.0
    t_dir = Affine(dx, 0.0, 4.999583333333334, 0.0, -dx, 4.999583333333333)
    t_upa_ulp = Affine(dx, 0.0, 4.999583333333334, 0.0, -dx, 4.999583333333334)
    d, u = _sample_arrays((400, 400))

    p_d, p_u = tmp_path / "ulp_dir.tif", tmp_path / "ulp_upa.tif"
    _write_merit_source_tile(p_d, p_u, transform=t_dir, upa_transform=t_upa_ulp, dir_data=d, upa_data=u)

    out_dir = tmp_path / "out_ulp"
    sources = [TilePair(name="n00e005", dir_path=p_d, upa_path=p_u)]
    build_merit_rasters(sources, out_dir, complete_global=True, tile_size=400)
    assert (
        verify_merit_raster_pair(out_dir / "flow_dir.tif", out_dir / "flow_acc.tif", sources, tile_size=400)["status"]
        == "passed"
    )

    t_upa_shifted = Affine(dx, 0.0, 4.999583433333334, 0.0, -dx, 4.999583333333334)
    p_bd, p_bu = tmp_path / "bad_dir.tif", tmp_path / "bad_upa.tif"
    _write_merit_source_tile(p_bd, p_bu, transform=t_dir, upa_transform=t_upa_shifted, dir_data=d, upa_data=u)
    with pytest.raises(ValueError, match="grid mismatch"):
        build_merit_rasters([TilePair(name="bad", dir_path=p_bd, upa_path=p_bu)], tmp_path / "out_bad", tile_size=400)
