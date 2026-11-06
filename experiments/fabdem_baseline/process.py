"""Mosaic FABDEM, project elevations, then fill depressions and route D8."""

from __future__ import annotations

import argparse
import json
import math
from contextlib import ExitStack
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pyflwdir
import rasterio
from affine import Affine
from numpy.typing import NDArray
from pyflwdir.dem import fill_depressions
from pyproj import Transformer
from rasterio.merge import merge
from rasterio.warp import Resampling, reproject
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from shapely.ops import transform, unary_union

CRS = "EPSG:32616"
NODATA = -9999


def basin_geometry(path: Path) -> BaseGeometry:
    data = json.loads(path.read_text())
    return unary_union([shape(feature["geometry"]) for feature in data["features"]])


def route(elevation: NDArray[np.float32]) -> tuple[NDArray[np.float32], NDArray[np.uint8], NDArray[np.int32]]:
    filled, directions = fill_depressions(elevation, nodata=NODATA, outlets="edge", max_depth=-1)
    flow = pyflwdir.from_array(directions, ftype="d8")
    counts = flow.accuflux(np.ones(elevation.shape, dtype=np.int32))
    return filled, directions, counts


def write_raster(path: Path, values: NDArray, affine: Affine, nodata: int | float, crs: str = CRS) -> None:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=values.shape[1],
        height=values.shape[0],
        count=1,
        dtype=values.dtype,
        crs=crs,
        transform=affine,
        nodata=nodata,
        compress="deflate",
        tiled=True,
    ) as dst:
        dst.write(values, 1)


def project_elevation(
    source: NDArray[np.float32],
    native: Affine,
    source_crs: str,
    source_nodata: float,
    target: Affine,
    shape: tuple[int, int],
) -> NDArray[np.float32]:
    elevation = np.full(shape, NODATA, dtype=np.float32)
    reproject(
        source,
        elevation,
        src_transform=native,
        src_crs=source_crs,
        src_nodata=source_nodata,
        dst_transform=target,
        dst_crs=CRS,
        dst_nodata=NODATA,
        resampling=Resampling.bilinear,
        tolerance=0.0,
        XSCALE=1.0,
        YSCALE=1.0,
    )
    return elevation


def process(base: Path, basin: Path, buffer_km: float) -> None:
    output = base / "outputs" / f"buffer-{buffer_km:g}km"
    output.mkdir(parents=True, exist_ok=True)
    geometry = transform(Transformer.from_crs(4326, CRS, always_xy=True).transform, basin_geometry(basin))
    xmin, ymin, xmax, ymax = geometry.buffer(buffer_km * 1000).bounds
    xmin, ymin = math.floor(xmin / 30) * 30, math.floor(ymin / 30) * 30
    xmax, ymax = math.ceil(xmax / 30) * 30, math.ceil(ymax / 30) * 30
    affine = Affine(30, 0, xmin, 0, -30, ymax)
    elevation = np.full((int((ymax - ymin) / 30), int((xmax - xmin) / 30)), NODATA, dtype=np.float32)
    with ExitStack() as stack:
        sources = [stack.enter_context(rasterio.open(path)) for path in sorted((base / "inputs").glob("N*FABDEM*.tif"))]
        mosaic, native_transform = merge(sources)
        elevation = project_elevation(
            mosaic[0], native_transform, sources[0].crs, sources[0].nodata, affine, elevation.shape
        )
    if np.any(elevation == NODATA) or not np.isfinite(elevation).all():
        raise ValueError("Processing rectangle has missing elevations; acquire more tiles before routing")
    filled, directions, counts = route(elevation)
    for name, values, nodata in [
        ("elevation", elevation, NODATA),
        ("filled", filled, NODATA),
        ("direction", directions, 247),
        ("accumulation_cells", counts, -1),
    ]:
        write_raster(output / f"{name}.tif", values, affine, nodata)
    settings = {
        "buffer_km": buffer_km,
        "crs": CRS,
        "resolution_m": 30,
        "bounds": [xmin, ymin, xmax, ymax],
        "shape": list(elevation.shape),
        "resampling": "bilinear elevation before routing; tolerance=0; XSCALE=YSCALE=1 (fixed sampling kernel)",
        "method": "pyflwdir Wang & Liu priority flood; edge outlets; max_depth=-1; connectivity=8",
        "accumulation": "inclusive upstream cell count; multiply by 0.0009 for projected km2",
        "nodata_cells": 0,
        "filled_cells": int((filled > elevation).sum()),
        "max_fill_m": float((filled - elevation).max()),
        "gdal_version": rasterio.__gdal_version__,
        "versions": {name: version(name) for name in ["pyflwdir", "numpy", "numba", "rasterio", "pyproj", "shapely"]},
    }
    (output / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    print(json.dumps(settings), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("basin", type=Path)
    parser.add_argument("--buffer-km", type=float, required=True)
    args = parser.parse_args()
    process(args.base, args.basin, args.buffer_km)
