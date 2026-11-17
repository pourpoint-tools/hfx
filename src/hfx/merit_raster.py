"""MERIT Hydro native-grid D8 raster validation and BigTIFF normalization."""

from __future__ import annotations

import argparse
import json
import shutil
import tarfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from queue import Queue
from tempfile import TemporaryDirectory
from threading import Thread
from time import perf_counter
from typing import Any, Literal, cast

import numpy as np
import numpy.typing as npt
import rasterio
import tifffile
from rasterio.transform import Affine
from rasterio.windows import Window

from hfx.candidates import fingerprint

ALLOWED_DIR_CODES = frozenset({0, 1, 2, 4, 8, 16, 32, 64, 128, 247, 255})
DIR_NODATA = 247
UPA_NODATA = np.float32(-9999.0)
UPA_NODATA_BITS = np.uint32(np.float32(-9999.0).view(np.uint32))
DEFAULT_TILE_SIZE = 400
TRANSFORM_TOLERANCE = 1e-12


@dataclass(frozen=True)
class TilePair:
    """One paired MERIT Hydro flow-direction and upstream-area source GeoTIFF."""

    name: str
    dir_path: Path
    upa_path: Path


@dataclass(frozen=True)
class TarArchivePair:
    """One paired 30-degree MERIT Hydro dir and upa tar archive."""

    block: str
    dir_tar: Path
    upa_tar: Path


@dataclass(frozen=True)
class TileHeader:
    """Inspected GeoTIFF header metadata read directly from source files."""

    name: str
    width: int
    height: int
    transform: tuple[float, float, float, float, float, float]
    crs: str
    dir_dtype: str
    dir_nodata: float | None
    upa_dtype: str
    upa_nodata: float | None
    dir_archive: Path | None = None
    dir_member: str | None = None
    upa_archive: Path | None = None
    upa_member: str | None = None
    dir_path: Path | None = None
    upa_path: Path | None = None

    @property
    def affine(self) -> Affine:
        return Affine(*self.transform)


@dataclass(frozen=True)
class GridPlan:
    """Exact grid geometry and tile-index layout derived from source headers."""

    west: float
    north: float
    dx: float
    dy: float
    width: int
    height: int
    tile_size: int
    tile_cols: int
    tile_rows: int
    tile_source_shape: tuple[int, int]
    placed_tiles: dict[tuple[int, int], TileHeader]

    @property
    def transform(self) -> Affine:
        return Affine(self.dx, 0.0, self.west, 0.0, self.dy, self.north)


def _inspect_open_pair(
    name: str,
    dir_ds: rasterio.DatasetReader,
    upa_ds: rasterio.DatasetReader,
    *,
    dir_archive: Path | None = None,
    dir_member: str | None = None,
    upa_archive: Path | None = None,
    upa_member: str | None = None,
    dir_path: Path | None = None,
    upa_path: Path | None = None,
) -> TileHeader:
    if dir_ds.count != 1 or upa_ds.count != 1:
        raise ValueError(f"{name}: expected single-band rasters")
    if str(dir_ds.crs) != "EPSG:4326" or str(upa_ds.crs) != "EPSG:4326":
        raise ValueError(f"{name}: expected EPSG:4326 CRS, got {dir_ds.crs} and {upa_ds.crs}")
    if dir_ds.dtypes[0] != "uint8":
        raise ValueError(f"{name}: expected uint8 direction, got {dir_ds.dtypes[0]}")
    if dir_ds.nodata not in (None, float(DIR_NODATA)):
        raise ValueError(f"{name}: unexpected direction nodata {dir_ds.nodata}")
    if upa_ds.dtypes[0] != "float32" or upa_ds.nodata != float(UPA_NODATA):
        raise ValueError(f"{name}: expected float32 upa with nodata -9999, got {upa_ds.dtypes[0]} / {upa_ds.nodata}")
    if dir_ds.shape != upa_ds.shape or not dir_ds.transform.almost_equals(
        upa_ds.transform, precision=TRANSFORM_TOLERANCE
    ):
        raise ValueError(f"{name}: dir/upa grid mismatch")
    t = dir_ds.transform
    if t.b != 0.0 or t.d != 0.0 or t.a <= 0.0 or t.e >= 0.0:
        raise ValueError(f"{name}: expected north-up axis-aligned affine transform, got {t}")
    return TileHeader(
        name=name,
        width=dir_ds.width,
        height=dir_ds.height,
        transform=(t.a, t.b, t.c, t.d, t.e, t.f),
        crs="EPSG:4326",
        dir_dtype=dir_ds.dtypes[0],
        dir_nodata=dir_ds.nodata,
        upa_dtype=upa_ds.dtypes[0],
        upa_nodata=upa_ds.nodata,
        dir_archive=dir_archive,
        dir_member=dir_member,
        upa_archive=upa_archive,
        upa_member=upa_member,
        dir_path=dir_path,
        upa_path=upa_path,
    )


def inspect_tile_pair(pair: TilePair) -> TileHeader:
    with rasterio.open(pair.dir_path) as dir_ds, rasterio.open(pair.upa_path) as upa_ds:
        return _inspect_open_pair(
            pair.name,
            dir_ds,
            upa_ds,
            dir_path=pair.dir_path,
            upa_path=pair.upa_path,
        )


def _archive_tile_key(member_name: str, suffix: str) -> str | None:
    base = Path(member_name).name
    if not base.endswith(suffix) or base.startswith("."):
        return None
    return base[: -len(suffix)]


def inspect_tar_archives(pairs: Sequence[TarArchivePair]) -> list[TileHeader]:
    """Read authoritative GeoTIFF headers for every paired tile across tar archives."""
    headers: list[TileHeader] = []
    seen_names: set[str] = set()
    with TemporaryDirectory(prefix="hfx-merit-hdr-") as scratch:
        scratch_dir = Path(scratch)
        for archive_pair in pairs:
            with tarfile.open(archive_pair.dir_tar, "r") as d_tar, tarfile.open(archive_pair.upa_tar, "r") as u_tar:
                dir_members: dict[str, tarfile.TarInfo] = {}
                for member in d_tar.getmembers():
                    if not member.isfile():
                        continue
                    key = _archive_tile_key(member.name, "_dir.tif")
                    if key is not None:
                        if key in dir_members:
                            raise ValueError(f"Duplicate dir tile {key} in {archive_pair.dir_tar}")
                        dir_members[key] = member
                upa_members: dict[str, tarfile.TarInfo] = {}
                for member in u_tar.getmembers():
                    if not member.isfile():
                        continue
                    key = _archive_tile_key(member.name, "_upa.tif")
                    if key is not None:
                        if key in upa_members:
                            raise ValueError(f"Duplicate upa tile {key} in {archive_pair.upa_tar}")
                        upa_members[key] = member
                if set(dir_members) != set(upa_members):
                    raise ValueError(
                        f"{archive_pair.block}: mismatched dir/upa members: "
                        f"dir-only={sorted(set(dir_members) - set(upa_members))}, "
                        f"upa-only={sorted(set(upa_members) - set(dir_members))}"
                    )
                for key in sorted(dir_members):
                    if key in seen_names:
                        raise ValueError(f"Duplicate tile {key} across archives")
                    seen_names.add(key)
                    d_tmp = scratch_dir / f"{key}_dir.tif"
                    u_tmp = scratch_dir / f"{key}_upa.tif"
                    d_src = d_tar.extractfile(dir_members[key])
                    u_src = u_tar.extractfile(upa_members[key])
                    if d_src is None or u_src is None:
                        raise ValueError(f"Unable to extract member for {key}")
                    with d_src, d_tmp.open("wb") as out:
                        shutil.copyfileobj(d_src, out)
                    with u_src, u_tmp.open("wb") as out:
                        shutil.copyfileobj(u_src, out)
                    try:
                        with rasterio.open(d_tmp) as dir_ds, rasterio.open(u_tmp) as upa_ds:
                            headers.append(
                                _inspect_open_pair(
                                    key,
                                    dir_ds,
                                    upa_ds,
                                    dir_archive=archive_pair.dir_tar,
                                    dir_member=dir_members[key].name,
                                    upa_archive=archive_pair.upa_tar,
                                    upa_member=upa_members[key].name,
                                )
                            )
                    finally:
                        d_tmp.unlink(missing_ok=True)
                        u_tmp.unlink(missing_ok=True)
    if not headers:
        raise ValueError("No paired tiles found in archives")
    return headers


def plan_grid(headers: Sequence[TileHeader], *, tile_size: int = DEFAULT_TILE_SIZE) -> GridPlan:
    """Validate grid alignment across tiles and compute exact tile placement."""
    if not headers:
        raise ValueError("At least one tile header is required")
    if tile_size <= 0 or tile_size % 16 != 0:
        raise ValueError("TIFF tile_size must be a positive multiple of 16")
    ref = headers[0]
    dx, _, _, _, dy, _ = ref.transform
    tile_w, tile_h = ref.width, ref.height
    if tile_w % tile_size != 0 or tile_h % tile_size != 0:
        raise ValueError(f"Source shape {(tile_h, tile_w)} must be divisible by tile_size={tile_size}")
    for h in headers:
        if h.width != tile_w or h.height != tile_h:
            raise ValueError(f"{h.name}: inconsistent tile dimensions {(h.height, h.width)} vs {(tile_h, tile_w)}")
        if abs(h.transform[0] - dx) > TRANSFORM_TOLERANCE or abs(h.transform[4] - dy) > TRANSFORM_TOLERANCE:
            raise ValueError(f"{h.name}: inconsistent pixel size {h.transform} vs {ref.transform}")

    west = min(h.transform[2] for h in headers)
    north = max(h.transform[5] for h in headers)
    east = max(h.transform[2] + h.width * dx for h in headers)
    south = min(h.transform[5] + h.height * dy for h in headers)

    total_w_exact = (east - west) / dx
    total_h_exact = (north - south) / (-dy)
    total_w = round(total_w_exact)
    total_h = round(total_h_exact)
    if abs(total_w_exact - total_w) > 1e-6 or abs(total_h_exact - total_h) > 1e-6:
        raise ValueError("Grid extent is not an integer multiple of pixel dimensions")
    if total_w % tile_w != 0 or total_h % tile_h != 0:
        raise ValueError("Total extent is not aligned to source tile dimensions")

    placed: dict[tuple[int, int], TileHeader] = {}
    for h in headers:
        col_exact = (h.transform[2] - west) / dx
        row_exact = (north - h.transform[5]) / (-dy)
        col_off = round(col_exact)
        row_off = round(row_exact)
        if abs(col_exact - col_off) > 1e-6 or abs(row_exact - row_off) > 1e-6:
            raise ValueError(f"{h.name}: grid misalignment at col={col_exact}, row={row_exact}")
        if col_off % tile_w != 0 or row_off % tile_h != 0:
            raise ValueError(f"{h.name}: tile offset ({row_off}, {col_off}) not aligned to source tile grid")
        key = (row_off // tile_h, col_off // tile_w)
        if key in placed:
            raise ValueError(f"{h.name}: overlaps with {placed[key].name} at tile block {key}")
        placed[key] = h

    return GridPlan(
        west=west,
        north=north,
        dx=dx,
        dy=dy,
        width=total_w,
        height=total_h,
        tile_size=tile_size,
        tile_cols=total_w // tile_w,
        tile_rows=total_h // tile_h,
        tile_source_shape=(tile_h, tile_w),
        placed_tiles=placed,
    )


def validate_chunk_semantics(
    dir_arr: npt.NDArray[np.uint8],
    upa_arr: npt.NDArray[np.float32],
    upa_valid_mask: npt.NDArray[np.bool_],
    *,
    tile_name: str,
    row_offset: int,
    transform: Affine,
    counts: dict[int, int],
    depression_witnesses: list[dict[str, object]],
    mouth_witnesses: list[dict[str, object]],
    max_witnesses: int = 8,
) -> tuple[int, int, float | None, float | None]:
    """Validate MERIT direction codes, mask equality with upa, and upa value sanity."""
    unique_vals, unique_counts = np.unique(dir_arr, return_counts=True)
    for val_np, cnt_np in zip(unique_vals, unique_counts, strict=True):
        val = int(val_np)
        if val not in ALLOWED_DIR_CODES:
            r, c = np.argwhere(dir_arr == val_np)[0]
            raise ValueError(f"{tile_name}: unexpected direction code {val} at row={row_offset + int(r)}, col={int(c)}")
        counts[val] = counts.get(val, 0) + int(cnt_np)

    dir_ocean = dir_arr == DIR_NODATA
    upa_invalid = ~upa_valid_mask
    if not np.array_equal(dir_ocean, upa_invalid):
        mismatch = np.argwhere(dir_ocean != upa_invalid)
        r, c = int(mismatch[0, 0]), int(mismatch[0, 1])
        raise ValueError(
            f"{tile_name}: direction==247 disagrees with upa invalid mask ({len(mismatch)} cells, "
            f"first at row={row_offset + r}, col={c}: dir={int(dir_arr[r, c])}, upa={float(upa_arr[r, c])})"
        )

    valid_upa = upa_arr[upa_valid_mask]
    if valid_upa.size:
        if not np.all(np.isfinite(valid_upa)) or np.any(valid_upa < 0):
            bad = np.argwhere(upa_valid_mask & (~np.isfinite(upa_arr) | (upa_arr < 0)))
            r, c = int(bad[0, 0]), int(bad[0, 1])
            raise ValueError(f"{tile_name}: invalid valid upa value {upa_arr[r, c]} at row={row_offset + r}, col={c}")
        min_upa: float | None = float(valid_upa.min())
        max_upa: float | None = float(valid_upa.max())
    else:
        min_upa = None
        max_upa = None

    for code, witnesses in ((255, depression_witnesses), (0, mouth_witnesses)):
        if len(witnesses) < max_witnesses and np.any(dir_arr == code):
            pos = np.argwhere(dir_arr == code)[0]
            r, c = int(pos[0]), int(pos[1])
            lon, lat = transform @ (c + 0.5, row_offset + r + 0.5)
            witnesses.append(
                {
                    "tile": tile_name,
                    "row": row_offset + r,
                    "col": c,
                    "lon": float(lon),
                    "lat": float(lat),
                    "upa": float(upa_arr[r, c]),
                }
            )

    valid_count = int(upa_valid_mask.sum())
    invalid_count = int(upa_invalid.sum())
    return valid_count, invalid_count, min_upa, max_upa


def _geotiff_extratags(
    transform: Affine,
    dtype: Literal["uint8", "float32"],
    nodata: float,
) -> list[tuple[int, str, int, Any, bool]]:
    tag_dtype_map = {2: "s", 3: "H", 12: "d"}
    with TemporaryDirectory(prefix="hfx-geotiff-tpl-") as scratch:
        template = Path(scratch) / "template.tif"
        with rasterio.open(
            template,
            "w",
            driver="GTiff",
            width=16,
            height=16,
            count=1,
            dtype=dtype,
            crs="EPSG:4326",
            transform=transform,
            nodata=nodata,
        ) as dst:
            dst.write(np.full((16, 16), nodata, dtype=dtype), 1)
            dst.update_tags(AREA_OR_POINT="Area")
        extratags: list[tuple[int, str, int, Any, bool]] = []
        with tifffile.TiffFile(template) as tif:
            page = tif.pages[0]
            if not isinstance(page, tifffile.TiffPage):
                raise ValueError("Template GeoTIFF missing main page")
            for code in (33550, 33922, 34735, 34736, 34737, 42112, 42113):
                if code in page.tags:
                    tag = page.tags[code]
                    dt = tag_dtype_map[int(tag.dtype)]
                    val = str(tag.value) if dt == "s" else tag.value
                    count = 0 if dt == "s" else tag.count
                    extratags.append((tag.code, dt, count, val, True))
        return extratags


def _write_queued_bigtiff(
    output_path: Path,
    *,
    shape: tuple[int, int],
    dtype: Literal["uint8", "float32"],
    tile_size: int,
    transform: Affine,
    nodata: float,
    queue: Queue[npt.NDArray[np.generic] | None],
    maxworkers: int,
    errors: list[BaseException],
) -> None:
    extratags = _geotiff_extratags(transform, dtype, nodata)
    saw_sentinel = False

    def iterator() -> Iterator[npt.NDArray[np.generic]]:
        nonlocal saw_sentinel
        while True:
            item = queue.get()
            if item is None:
                saw_sentinel = True
                return
            yield item

    try:
        tifffile.imwrite(
            output_path,
            iterator(),
            shape=shape,
            dtype=dtype,
            bigtiff=True,
            tile=(tile_size, tile_size),
            compression="deflate",
            predictor=False,
            photometric="minisblack",
            metadata=None,
            extratags=extratags,
            maxworkers=maxworkers,
        )
    except BaseException as exc:  # pragma: no cover - surfaced to caller thread
        errors.append(exc)
        while not saw_sentinel:
            if queue.get() is None:
                saw_sentinel = True


def _build_single_grid(
    plan: GridPlan,
    dir_out: Path,
    acc_out: Path,
    *,
    maxworkers: int = 4,
) -> dict[str, Any]:
    tile_h, tile_w = plan.tile_source_shape
    ts = plan.tile_size
    sub_rows = tile_h // ts
    sub_cols = tile_w // ts

    ocean_dir_tile = np.full((ts, ts), DIR_NODATA, dtype=np.uint8)
    ocean_acc_tile = np.full((ts, ts), UPA_NODATA, dtype=np.float32)

    dir_queue: Queue[npt.NDArray[np.generic] | None] = Queue(maxsize=64)
    acc_queue: Queue[npt.NDArray[np.generic] | None] = Queue(maxsize=64)
    errors: list[BaseException] = []

    t_dir = Thread(
        target=_write_queued_bigtiff,
        kwargs={
            "output_path": dir_out,
            "shape": (plan.height, plan.width),
            "dtype": "uint8",
            "tile_size": ts,
            "transform": plan.transform,
            "nodata": float(DIR_NODATA),
            "queue": dir_queue,
            "maxworkers": maxworkers,
            "errors": errors,
        },
        daemon=True,
    )
    t_acc = Thread(
        target=_write_queued_bigtiff,
        kwargs={
            "output_path": acc_out,
            "shape": (plan.height, plan.width),
            "dtype": "float32",
            "tile_size": ts,
            "transform": plan.transform,
            "nodata": float(UPA_NODATA),
            "queue": acc_queue,
            "maxworkers": maxworkers,
            "errors": errors,
        },
        daemon=True,
    )
    t_dir.start()
    t_acc.start()

    counts: dict[int, int] = {}
    depression_witnesses: list[dict[str, object]] = []
    mouth_witnesses: list[dict[str, object]] = []
    valid_cells = 0
    invalid_cells = 0
    ocean_fill_tiles = 0
    upa_min: float | None = None
    upa_max: float | None = None

    try:
        for row_idx in range(plan.tile_rows):
            with TemporaryDirectory(prefix="hfx-merit-band-") as scratch:
                scratch_dir = Path(scratch)
                band_readers: dict[int, tuple[TileHeader, rasterio.DatasetReader, rasterio.DatasetReader]] = {}
                try:
                    # Extract and open only the current latitude strip of tiles.
                    by_archive: dict[tuple[Path, Path], list[tuple[int, TileHeader]]] = {}
                    direct_cols: list[tuple[int, TileHeader]] = []
                    for col_idx in range(plan.tile_cols):
                        hdr = plan.placed_tiles.get((row_idx, col_idx))
                        if hdr is None:
                            continue
                        if hdr.dir_path is not None and hdr.upa_path is not None:
                            direct_cols.append((col_idx, hdr))
                        elif (
                            hdr.dir_archive is not None
                            and hdr.upa_archive is not None
                            and hdr.dir_member is not None
                            and hdr.upa_member is not None
                        ):
                            by_archive.setdefault((hdr.dir_archive, hdr.upa_archive), []).append((col_idx, hdr))
                        else:  # pragma: no cover - invariant enforced by constructor
                            raise ValueError(f"{hdr.name}: incomplete tile source location")

                    for col_idx, hdr in direct_cols:
                        assert hdr.dir_path is not None and hdr.upa_path is not None
                        band_readers[col_idx] = (hdr, rasterio.open(hdr.dir_path), rasterio.open(hdr.upa_path))

                    for (d_tar_path, u_tar_path), items in by_archive.items():
                        with tarfile.open(d_tar_path, "r") as d_tar, tarfile.open(u_tar_path, "r") as u_tar:
                            for col_idx, hdr in items:
                                assert hdr.dir_member is not None and hdr.upa_member is not None
                                d_local = scratch_dir / f"{hdr.name}_dir.tif"
                                u_local = scratch_dir / f"{hdr.name}_upa.tif"
                                d_stream = d_tar.extractfile(hdr.dir_member)
                                u_stream = u_tar.extractfile(hdr.upa_member)
                                if d_stream is None or u_stream is None:
                                    raise ValueError(f"{hdr.name}: failed to extract from tar")
                                with d_stream, d_local.open("wb") as out:
                                    shutil.copyfileobj(d_stream, out)
                                with u_stream, u_local.open("wb") as out:
                                    shutil.copyfileobj(u_stream, out)
                                band_readers[col_idx] = (hdr, rasterio.open(d_local), rasterio.open(u_local))

                    for sr in range(sub_rows):
                        y0 = sr * ts
                        window = Window.from_slices((y0, y0 + ts), (0, tile_w))
                        strips: dict[int, tuple[npt.NDArray[np.uint8], npt.NDArray[np.float32]]] = {}
                        for col_idx, (hdr, d_ds, u_ds) in band_readers.items():
                            d_strip = cast(npt.NDArray[np.uint8], d_ds.read(1, window=window))
                            u_strip = cast(npt.NDArray[np.float32], u_ds.read(1, window=window))
                            u_mask = cast(npt.NDArray[np.uint8], u_ds.read_masks(1, window=window))
                            u_valid = (u_mask != 0) & (u_strip.view(np.uint32) != UPA_NODATA_BITS)
                            v_cnt, i_cnt, c_min, c_max = validate_chunk_semantics(
                                d_strip,
                                u_strip,
                                u_valid,
                                tile_name=hdr.name,
                                row_offset=y0,
                                transform=hdr.affine,
                                counts=counts,
                                depression_witnesses=depression_witnesses,
                                mouth_witnesses=mouth_witnesses,
                            )
                            valid_cells += v_cnt
                            invalid_cells += i_cnt
                            if c_min is not None:
                                upa_min = c_min if upa_min is None else min(upa_min, c_min)
                            if c_max is not None:
                                upa_max = c_max if upa_max is None else max(upa_max, c_max)
                            strips[col_idx] = (d_strip, u_strip)

                        for col_idx in range(plan.tile_cols):
                            if col_idx not in strips:
                                ocean_fill_tiles += sub_cols
                                counts[DIR_NODATA] = counts.get(DIR_NODATA, 0) + sub_cols * ts * ts
                                invalid_cells += sub_cols * ts * ts
                                for _ in range(sub_cols):
                                    if errors:
                                        raise errors[0]
                                    dir_queue.put(ocean_dir_tile)
                                    acc_queue.put(ocean_acc_tile)
                            else:
                                d_strip, u_strip = strips[col_idx]
                                for sc in range(sub_cols):
                                    if errors:
                                        raise errors[0]
                                    x0 = sc * ts
                                    dir_queue.put(np.ascontiguousarray(d_strip[:, x0 : x0 + ts]))
                                    acc_queue.put(np.ascontiguousarray(u_strip[:, x0 : x0 + ts]))
                finally:
                    for _, d_ds, u_ds in band_readers.values():
                        d_ds.close()
                        u_ds.close()
    except BaseException:
        dir_queue.put(None)
        acc_queue.put(None)
        t_dir.join()
        t_acc.join()
        dir_out.unlink(missing_ok=True)
        acc_out.unlink(missing_ok=True)
        raise

    dir_queue.put(None)
    acc_queue.put(None)
    t_dir.join()
    t_acc.join()
    if errors:
        raise errors[0]

    return {
        "flow_dir": dir_out.name,
        "flow_acc": acc_out.name,
        "flow_dir_identity": fingerprint(dir_out),
        "flow_acc_identity": fingerprint(acc_out),
        "width": plan.width,
        "height": plan.height,
        "transform": list(plan.transform)[:6],
        "crs": "EPSG:4326",
        "tile_size": ts,
        "source_tiles": len(plan.placed_tiles),
        "ocean_fill_output_tiles": ocean_fill_tiles,
        "valid_cells": valid_cells,
        "invalid_cells": invalid_cells,
        "upa_min": upa_min,
        "upa_max": upa_max,
        "counts": {str(k): v for k, v in sorted(counts.items())},
        "depression_witnesses": depression_witnesses,
        "mouth_witnesses": mouth_witnesses,
    }


def _headers_from_sources(sources: Sequence[TilePair] | Sequence[TarArchivePair]) -> list[TileHeader]:
    if not sources:
        raise ValueError("At least one raster source is required")
    first = sources[0]
    if isinstance(first, TilePair):
        return [inspect_tile_pair(cast(TilePair, item)) for item in sources]
    return inspect_tar_archives(cast(Sequence[TarArchivePair], sources))


def build_merit_rasters(
    sources: Sequence[TilePair] | Sequence[TarArchivePair],
    output_dir: Path,
    *,
    complete_global: bool = False,
    tile_size: int = DEFAULT_TILE_SIZE,
    maxworkers: int = 4,
) -> dict[str, Any]:
    """Validate and package MERIT Hydro dir/upa sources into HFX-conformant BigTIFF pairs."""
    headers = _headers_from_sources(sources)
    output_dir.mkdir(parents=True, exist_ok=True)
    if complete_global:
        plan = plan_grid(headers, tile_size=tile_size)
        pair_receipt = _build_single_grid(
            plan,
            output_dir / "flow_dir.tif",
            output_dir / "flow_acc.tif",
            maxworkers=maxworkers,
        )
        return {
            "mode": "complete_global",
            "pairs": [pair_receipt],
            "counts": pair_receipt["counts"],
            "depression_witnesses": pair_receipt["depression_witnesses"],
            "mouth_witnesses": pair_receipt["mouth_witnesses"],
        }

    # Partial mode: never fill missing land tiles with ocean nodata; validate cross-tile grid coherence
    # across all headers, then emit one verified pair per acquired source tile.
    _ = plan_grid(headers, tile_size=tile_size)
    pairs: list[dict[str, Any]] = []
    total_counts: dict[int, int] = {}
    depression_witnesses: list[dict[str, object]] = []
    mouth_witnesses: list[dict[str, object]] = []
    for hdr in headers:
        plan = plan_grid([hdr], tile_size=tile_size)
        pair_receipt = _build_single_grid(
            plan,
            output_dir / f"{hdr.name}_flow_dir.tif",
            output_dir / f"{hdr.name}_flow_acc.tif",
            maxworkers=maxworkers,
        )
        pair_receipt["name"] = hdr.name
        pairs.append(pair_receipt)
        for key, value in pair_receipt["counts"].items():
            code = int(key)
            total_counts[code] = total_counts.get(code, 0) + int(value)
        for w in pair_receipt["depression_witnesses"]:
            if len(depression_witnesses) < 8:
                depression_witnesses.append(w)
        for w in pair_receipt["mouth_witnesses"]:
            if len(mouth_witnesses) < 8:
                mouth_witnesses.append(w)

    return {
        "mode": "partial_per_tile",
        "pairs": pairs,
        "counts": {str(k): v for k, v in sorted(total_counts.items())},
        "depression_witnesses": depression_witnesses,
        "mouth_witnesses": mouth_witnesses,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tar_manifest", type=Path, help="JSON list of {block, dir_tar, upa_tar}")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--complete-global", action="store_true")
    parser.add_argument("--tile-size", type=int, default=DEFAULT_TILE_SIZE)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    started = perf_counter()
    raw_pairs = json.loads(args.tar_manifest.read_text())
    if raw_pairs and "dir_path" in raw_pairs[0]:
        sources: list[TilePair] | list[TarArchivePair] = [
            TilePair(p["name"], Path(p["dir_path"]), Path(p["upa_path"])) for p in raw_pairs
        ]
    else:
        sources = [TarArchivePair(p["block"], Path(p["dir_tar"]), Path(p["upa_tar"])) for p in raw_pairs]
    result = build_merit_rasters(
        sources,
        args.output_dir,
        complete_global=args.complete_global,
        tile_size=args.tile_size,
        maxworkers=args.workers,
    )
    result["seconds"] = perf_counter() - started
    args.report.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
