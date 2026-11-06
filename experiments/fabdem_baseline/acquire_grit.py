"""Read a bounded native-grid window from the published HFX 0.4.0 GRIT bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import rasterio
from process import basin_geometry
from pyproj import Transformer
from rasterio.features import geometry_mask
from rasterio.windows import Window
from shapely.ops import transform

SOURCE = "https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/aux/d8/"


def acquire(base: Path, cached: Path | None) -> None:
    # Prior diagnostic window, enlarged beyond the basin; native GRIT cell indices.
    window = Window(233826, 145103, 2658, 3091)  # ty: ignore[too-many-positional-arguments] -- attrs-generated signature
    arrays = {}
    provenance = {"source": SOURCE, "window_col_row_width_height": list(window.flatten())}
    for name, key in [("flow_dir.tif", "raw_direction"), ("flow_acc.tif", "accumulation")]:
        with rasterio.open(SOURCE + name) as src:
            arrays[key] = src.read(1, window=window)
            arrays["direction_valid" if key == "raw_direction" else "accumulation_valid"] = (
                src.read_masks(1, window=window) > 0
            )
            affine, crs = src.window_transform(window), str(src.crs)
            provenance[name] = {
                "url": SOURCE + name,
                "dtype": src.dtypes[0],
                "nodata": str(src.nodata),
                "window_array_sha256": hashlib.sha256(arrays[key].tobytes()).hexdigest(),
            }
    geometry = transform(
        Transformer.from_crs(4326, crs, always_xy=True).transform, basin_geometry(base / "inputs/coarse_union.geojson")
    )
    arrays["vector_mask"] = geometry_mask([geometry], arrays["raw_direction"].shape, affine, invert=True)
    arrays["transform"] = np.array(affine)
    arrays["crs"] = np.array(crs)
    if cached is not None:
        with np.load(cached) as previous:
            provenance["previous_cache_comparison"] = {
                key: bool(np.array_equal(previous[key], arrays[key])) for key in arrays
            }
    np.savez_compressed(base / "inputs/grit-window.npz", allow_pickle=False, **arrays)
    (base / "inputs/grit-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("--cached", type=Path)
    args = parser.parse_args()
    acquire(args.base, args.cached)
