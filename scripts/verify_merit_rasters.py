"""Independent verification and remote range checks for MERIT Hydro HFX rasters."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, BinaryIO, cast

import numpy as np
import numpy.typing as npt
import rasterio
import tifffile
from rasterio.enums import MaskFlags
from rasterio.windows import Window
from rio_cogeo.cogeo import cog_validate

from hfx.candidates import fingerprint
from hfx.merit_raster import (
    DEFAULT_TILE_SIZE,
    DIR_NODATA,
    TRANSFORM_TOLERANCE,
    UPA_NODATA_BITS,
    GridPlan,
    TarArchivePair,
    TileHeader,
    TilePair,
    inspect_tar_archives,
    inspect_tile_pair,
    plan_grid,
)

OVERVIEW_WARNING = "The file is greater than 512xH or 512xW, it is recommended to include internal overviews"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def check_tiff_layout(
    path: Path, *, expected_tile_size: int, expected_nodata: float, expected_dtype: str
) -> dict[str, Any]:
    """Verify BigTIFF structure, IFD-before-data ordering, zero trailing gap, and COG validator."""
    with tifffile.TiffFile(path) as tif:
        _require(tif.is_bigtiff, f"{path.name}: expected BigTIFF")
        pages = list(tif.pages)
        _require(len(pages) == 1, f"{path.name}: expected single IFD without overviews")
        page = pages[0]
        _require(isinstance(page, tifffile.TiffPage), f"{path.name}: missing native TiffPage")
        assert isinstance(page, tifffile.TiffPage)
        _require(
            page.is_tiled and page.tilewidth == expected_tile_size and page.tilelength == expected_tile_size,
            f"{path.name}: expected {expected_tile_size}x{expected_tile_size} tiles",
        )
        _require(
            page.compression in (8, 32946) and page.predictor == 1, f"{path.name}: expected DEFLATE without predictor"
        )
        _require(page.samplesperpixel == 1 and not page.subifds, f"{path.name}: unexpected bands or SubIFDs")
        _require(page.offset < 300, f"{path.name}: IFD offset {page.offset} >= 300")
        offsets, lengths = page.dataoffsets, page.databytecounts
        expected_tiles = ((page.imagewidth + expected_tile_size - 1) // expected_tile_size) * (
            (page.imagelength + expected_tile_size - 1) // expected_tile_size
        )
        _require(
            len(offsets) == len(lengths) == expected_tiles and expected_tiles > 0, f"{path.name}: tile count mismatch"
        )
        _require(page.offset < offsets[0], f"{path.name}: IFD follows tile data")
        for tag in page.tags.values():
            _require(
                tag.valueoffset + tag.valuebytecount <= offsets[0], f"{path.name}: tag {tag.code} follows tile data"
            )
        end = 0
        size = path.stat().st_size
        for offset, length in zip(offsets, lengths, strict=True):
            _require(length > 0 and offset >= end and offset + length <= size, f"{path.name}: invalid tile range")
            end = offset + length
        _require(end == size, f"{path.name}: unexpected trailing bytes ({size - end})")

    with rasterio.open(path) as ds:
        _require(
            ds.files == [str(path.resolve())] or ds.files == [str(path)], f"{path.name}: external sidecars unsupported"
        )
        _require(str(ds.crs) == "EPSG:4326", f"{path.name}: CRS is {ds.crs}")
        _require(
            ds.count == 1 and ds.dtypes[0] == expected_dtype, f"{path.name}: dtype {ds.dtypes[0]} != {expected_dtype}"
        )
        _require(ds.nodata == expected_nodata, f"{path.name}: nodata {ds.nodata} != {expected_nodata}")
        _require(
            ds.mask_flag_enums == ([MaskFlags.nodata],), f"{path.name}: unexpected mask flags {ds.mask_flag_enums}"
        )
        _require(not ds.overviews(1), f"{path.name}: unexpected overviews")

    valid, errors, warnings = cog_validate(str(path), strict=False, quiet=True)
    _require(
        valid and not errors and all(w == OVERVIEW_WARNING for w in warnings),
        f"{path.name}: cog_validate failed: errors={errors}, warnings={warnings}",
    )
    return {"tiles": expected_tiles, "cog_warnings": warnings, "identity": fingerprint(path)}


def _verify_band_strips(
    hdr: TileHeader,
    src_dir: rasterio.DatasetReader,
    src_upa: rasterio.DatasetReader,
    dst_dir: rasterio.DatasetReader,
    dst_upa: rasterio.DatasetReader,
    plan: GridPlan,
    row_idx: int,
    col_idx: int,
) -> tuple[int, int]:
    tile_h, tile_w = plan.tile_source_shape
    ts = plan.tile_size
    base_y = row_idx * tile_h
    base_x = col_idx * tile_w
    checked_cells = 0
    for sr in range(tile_h // ts):
        y0 = sr * ts
        src_win = Window.from_slices((y0, y0 + ts), (0, tile_w))
        dst_win = Window.from_slices((base_y + y0, base_y + y0 + ts), (base_x, base_x + tile_w))
        sd = cast(npt.NDArray[np.uint8], src_dir.read(1, window=src_win))
        su = cast(npt.NDArray[np.float32], src_upa.read(1, window=src_win))
        su_mask = cast(npt.NDArray[np.uint8], src_upa.read_masks(1, window=src_win))

        dd = cast(npt.NDArray[np.uint8], dst_dir.read(1, window=dst_win))
        dd_mask = cast(npt.NDArray[np.uint8], dst_dir.read_masks(1, window=dst_win))
        du = cast(npt.NDArray[np.float32], dst_upa.read(1, window=dst_win))
        du_mask = cast(npt.NDArray[np.uint8], dst_upa.read_masks(1, window=dst_win))

        _require(np.array_equal(sd, dd), f"{hdr.name}: direction bytes differ at strip {sr}")
        expected_dir_mask = np.where(sd == DIR_NODATA, np.uint8(0), np.uint8(255))
        _require(np.array_equal(expected_dir_mask, dd_mask), f"{hdr.name}: direction mask differs at strip {sr}")
        _require(
            np.array_equal(su.view(np.uint32), du.view(np.uint32)),
            f"{hdr.name}: accumulation float32 bits differ at strip {sr}",
        )
        _require(np.array_equal(su_mask, du_mask), f"{hdr.name}: accumulation mask differs at strip {sr}")
        _require(
            np.array_equal(sd == DIR_NODATA, su_mask == 0),
            f"{hdr.name}: ocean-247 disagrees with upa mask at strip {sr}",
        )
        checked_cells += sd.size
    return checked_cells, 0


def verify_merit_raster_pair(
    flow_dir: Path,
    flow_acc: Path,
    sources: Sequence[TilePair] | Sequence[TarArchivePair],
    *,
    tile_size: int = DEFAULT_TILE_SIZE,
) -> dict[str, Any]:
    """Independently verify output BigTIFF layout, exact values, float bits, and nodata masks."""
    dir_layout = check_tiff_layout(
        flow_dir, expected_tile_size=tile_size, expected_nodata=float(DIR_NODATA), expected_dtype="uint8"
    )
    acc_layout = check_tiff_layout(
        flow_acc, expected_tile_size=tile_size, expected_nodata=-9999.0, expected_dtype="float32"
    )
    if not sources:
        raise ValueError("At least one source is required")
    headers = (
        [inspect_tile_pair(cast(TilePair, item)) for item in sources]
        if isinstance(sources[0], TilePair)
        else inspect_tar_archives(cast(Sequence[TarArchivePair], sources))
    )
    plan = plan_grid(headers, tile_size=tile_size)

    tile_h, tile_w = plan.tile_source_shape
    ocean_cells = 0
    land_cells = 0
    with rasterio.open(flow_dir) as dst_dir, rasterio.open(flow_acc) as dst_upa:
        _require(
            dst_dir.width == plan.width and dst_dir.height == plan.height, "Output dir dimensions differ from grid"
        )
        _require(
            dst_upa.width == plan.width and dst_upa.height == plan.height, "Output acc dimensions differ from grid"
        )
        _require(
            dst_dir.transform == dst_upa.transform
            and dst_dir.transform.almost_equals(plan.transform, precision=TRANSFORM_TOLERANCE),
            "Output transform differs",
        )

        for row_idx in range(plan.tile_rows):
            with TemporaryDirectory(prefix="hfx-verify-raster-") as scratch:
                scratch_dir = Path(scratch)
                band_paths: dict[int, tuple[TileHeader, Path, Path]] = {}
                by_archive: dict[tuple[Path, Path], list[tuple[int, TileHeader]]] = {}
                for col_idx in range(plan.tile_cols):
                    hdr = plan.placed_tiles.get((row_idx, col_idx))
                    if hdr is None:
                        continue
                    if hdr.dir_path is not None and hdr.upa_path is not None:
                        band_paths[col_idx] = (hdr, hdr.dir_path, hdr.upa_path)
                    else:
                        assert hdr.dir_archive and hdr.upa_archive
                        by_archive.setdefault((hdr.dir_archive, hdr.upa_archive), []).append((col_idx, hdr))

                for (d_tar_path, u_tar_path), items in by_archive.items():
                    with tarfile.open(d_tar_path, "r") as d_tar, tarfile.open(u_tar_path, "r") as u_tar:
                        for col_idx, hdr in items:
                            assert hdr.dir_member and hdr.upa_member
                            d_path = scratch_dir / f"{hdr.name}_dir.tif"
                            u_path = scratch_dir / f"{hdr.name}_upa.tif"
                            d_in, u_in = d_tar.extractfile(hdr.dir_member), u_tar.extractfile(hdr.upa_member)
                            if d_in is None or u_in is None:
                                raise ValueError(f"{hdr.name}: missing tar member")
                            with d_in, d_path.open("wb") as out:
                                shutil.copyfileobj(d_in, out)
                            with u_in, u_path.open("wb") as out:
                                shutil.copyfileobj(u_in, out)
                            band_paths[col_idx] = (hdr, d_path, u_path)

                for col_idx in range(plan.tile_cols):
                    if col_idx not in band_paths:
                        base_y = row_idx * tile_h
                        base_x = col_idx * tile_w
                        for sr in range(tile_h // tile_size):
                            y0 = base_y + sr * tile_size
                            win = Window.from_slices((y0, y0 + tile_size), (base_x, base_x + tile_w))
                            dd = dst_dir.read(1, window=win)
                            dm = dst_dir.read_masks(1, window=win)
                            du = dst_upa.read(1, window=win)
                            um = dst_upa.read_masks(1, window=win)
                            _require(bool(np.all(dd == DIR_NODATA)), f"Ocean fill dir mismatch at {(row_idx, col_idx)}")
                            _require(bool(np.all(dm == 0)), f"Ocean fill dir mask mismatch at {(row_idx, col_idx)}")
                            _require(
                                bool(np.all(du.view(np.uint32) == UPA_NODATA_BITS)),
                                f"Ocean fill upa bits mismatch at {(row_idx, col_idx)}",
                            )
                            _require(bool(np.all(um == 0)), f"Ocean fill upa mask mismatch at {(row_idx, col_idx)}")
                            ocean_cells += dd.size
                        continue

                    hdr, d_path, u_path = band_paths[col_idx]
                    with rasterio.open(d_path) as s_dir, rasterio.open(u_path) as s_upa:
                        c_land, _ = _verify_band_strips(hdr, s_dir, s_upa, dst_dir, dst_upa, plan, row_idx, col_idx)
                        land_cells += c_land
                    if hdr.dir_path is None:
                        d_path.unlink(missing_ok=True)
                        u_path.unlink(missing_ok=True)

    return {
        "status": "passed",
        "dir_bits_equal": True,
        "upa_bits_equal": True,
        "masks_equal": True,
        "dir_layout": dir_layout,
        "acc_layout": acc_layout,
        "verified_source_tiles": len(plan.placed_tiles),
        "land_tile_cells": land_cells,
        "ocean_fill_cells": ocean_cells,
    }


def probe_remote_raster_windows(
    bundle_url: str,
    local_bundle: Path | None,
    windows: Sequence[dict[str, Any]],
    *,
    aws_profile: str | None = None,
    s3_endpoint: str | None = None,
    aux_index: int = 0,
) -> dict[str, Any]:
    """Perform authenticated HTTP 206 range reads over paired HFX raster windows."""
    from urllib.parse import urljoin

    from scripts.check_vector_access import HTTP, S3, Fingerprint, RangeFile

    http: HTTP = HTTP()
    if aws_profile:
        import boto3
        from botocore.config import Config

        if not s3_endpoint:
            raise ValueError("--s3-endpoint is required with --aws-profile")
        http = S3(
            boto3.Session(profile_name=aws_profile).client(
                "s3", endpoint_url=s3_endpoint, config=Config(s3={"addressing_style": "path"})
            ),
            s3_endpoint,
        )

    base = bundle_url.rstrip("/") + "/"
    manifest_bytes = http.get(urljoin(base, "manifest.json"))
    manifest = json.loads(manifest_bytes)
    entry = manifest["auxiliary"][aux_index]
    _require(entry["schema"] == "hfx.aux.d8_raster.v2", "unexpected raster schema")
    _require(entry["metadata"]["flow_dir_encoding"] == "merit", "expected flow_dir_encoding=merit")

    def open_remote(rel_path: str) -> rasterio.DatasetReader:
        fp = cast(Fingerprint, manifest["artifacts"][rel_path])
        url = urljoin(base, rel_path)
        rf = RangeFile(http, url, fp)

        def opener(name: str, mode: str = "rb") -> BinaryIO:
            if name != rel_path or mode not in ("r", "rb"):
                raise FileNotFoundError(name)
            return cast(BinaryIO, rf)

        return rasterio.open(rel_path, opener=opener)

    dir_rel = entry["artifacts"]["flow_dir"]
    acc_rel = entry["artifacts"]["flow_acc"]
    results: list[dict[str, Any]] = []
    with (
        rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"),
        open_remote(dir_rel) as r_dir,
        open_remote(acc_rel) as r_acc,
    ):
        _require(r_dir.nodata == 247.0 and r_acc.nodata == -9999.0, "remote nodata mismatch")
        _require(str(r_dir.crs) == "EPSG:4326" and str(r_acc.crs) == "EPSG:4326", "remote CRS mismatch")
        local_dir = rasterio.open(local_bundle / dir_rel) if local_bundle else None
        local_acc = rasterio.open(local_bundle / acc_rel) if local_bundle else None
        try:
            for item in windows:
                col, row, width, height = item["window"]
                win = Window.from_slices((row, row + height), (col, col + width))
                req_before = len(http.requests)
                d_arr = r_dir.read(1, window=win)
                d_mask = r_dir.read_masks(1, window=win)
                a_arr = r_acc.read(1, window=win)
                a_mask = r_acc.read_masks(1, window=win)
                reqs = http.requests[req_before:]
                _require(np.array_equal(d_arr == DIR_NODATA, a_mask == 0), "remote ocean mask mismatch")
                if local_dir is not None and local_acc is not None:
                    _require(np.array_equal(d_arr, local_dir.read(1, window=win)), "remote dir differs from local")
                    _require(
                        np.array_equal(d_mask, local_dir.read_masks(1, window=win)),
                        "remote dir mask differs from local",
                    )
                    _require(
                        np.array_equal(a_arr.view(np.uint32), local_acc.read(1, window=win).view(np.uint32)),
                        "remote acc bits differ from local",
                    )
                    _require(
                        np.array_equal(a_mask, local_acc.read_masks(1, window=win)),
                        "remote acc mask differs from local",
                    )
                vals, cnts = np.unique(d_arr, return_counts=True)
                results.append(
                    {
                        "name": item["name"],
                        "window": item["window"],
                        "http_requests": len(reqs),
                        "http_body_bytes": sum(r["bytes"] for r in reqs),
                        "dir_sha256": hashlib.sha256(d_arr.tobytes()).hexdigest(),
                        "dir_mask_sha256": hashlib.sha256(d_mask.tobytes()).hexdigest(),
                        "acc_sha256": hashlib.sha256(a_arr.tobytes()).hexdigest(),
                        "acc_mask_sha256": hashlib.sha256(a_mask.tobytes()).hexdigest(),
                        "dir_value_counts": {str(int(k)): int(v) for k, v in zip(vals, cnts, strict=True)},
                        "matches_local": local_dir is not None,
                    }
                )
        finally:
            if local_dir is not None:
                local_dir.close()
            if local_acc is not None:
                local_acc.close()

    return {
        "status": "passed",
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "total_http_requests": len(http.requests),
        "total_http_body_bytes": sum(r["bytes"] for r in http.requests),
        "windows": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    verify_cmd = sub.add_parser("verify-local")
    verify_cmd.add_argument("flow_dir", type=Path)
    verify_cmd.add_argument("flow_acc", type=Path)
    verify_cmd.add_argument("tar_manifest", type=Path, help="JSON list of {block, dir_tar, upa_tar}")
    verify_cmd.add_argument("--tile-size", type=int, default=DEFAULT_TILE_SIZE)
    verify_cmd.add_argument("--report", type=Path, required=True)

    probe_cmd = sub.add_parser("probe-remote")
    probe_cmd.add_argument("bundle_url")
    probe_cmd.add_argument("windows_json", type=Path)
    probe_cmd.add_argument("--local-bundle", type=Path)
    probe_cmd.add_argument("--aws-profile")
    probe_cmd.add_argument("--s3-endpoint")
    probe_cmd.add_argument("--aux-index", type=int, default=0)
    probe_cmd.add_argument("--report", type=Path, required=True)

    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    if args.command == "verify-local":
        raw_pairs = json.loads(args.tar_manifest.read_text())
        if raw_pairs and "dir_path" in raw_pairs[0]:
            sources: list[TilePair] | list[TarArchivePair] = [
                TilePair(p["name"], Path(p["dir_path"]), Path(p["upa_path"])) for p in raw_pairs
            ]
        else:
            sources = [TarArchivePair(p["block"], Path(p["dir_tar"]), Path(p["upa_tar"])) for p in raw_pairs]
        result = verify_merit_raster_pair(args.flow_dir, args.flow_acc, sources, tile_size=args.tile_size)
    else:
        windows = json.loads(args.windows_json.read_text())
        result = probe_remote_raster_windows(
            args.bundle_url,
            args.local_bundle,
            windows,
            aws_profile=args.aws_profile,
            s3_endpoint=args.s3_endpoint,
            aux_index=args.aux_index,
        )
    args.report.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
