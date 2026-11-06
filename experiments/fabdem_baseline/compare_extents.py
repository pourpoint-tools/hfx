"""Compare exact aligned overlap across processed terrain extents."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import rasterio
from process import basin_geometry
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.windows import from_bounds
from shapely.ops import transform


def compare(base: Path, small: int, large: int) -> dict[str, object]:
    output = base / "outputs"
    result: dict[str, object] = {"small_km": small, "large_km": large}
    for name in ["elevation", "filled", "direction", "accumulation_cells", "reference_connected", "adjusted_connected"]:
        with (
            rasterio.open(output / f"buffer-{small}km" / f"{name}.tif") as a,
            rasterio.open(output / f"buffer-{large}km" / f"{name}.tif") as b,
        ):
            x = a.read(1)
            window = from_bounds(*a.bounds, b.transform)
            y = b.read(1, window=window)
            if a.transform != b.window_transform(window) or x.shape != y.shape:
                raise ValueError("Extent comparison requires exactly aligned grids")
            geometry = transform(
                Transformer.from_crs(4326, a.crs, always_xy=True).transform,
                basin_geometry(base / "inputs/coarse_union.geojson"),
            )
            mask = geometry_mask([geometry], x.shape, a.transform, invert=True)
            result[name] = {
                "overlap_changed_cells": int((x != y).sum()),
                "vector_changed_cells": int(((x != y) & mask).sum()),
                "max_absolute_difference": float(np.max(np.abs(x.astype(float) - y))),
            }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    args = parser.parse_args()
    result = [compare(args.base, 10, 20), compare(args.base, 20, 25)]
    (args.base / "report/extent-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
