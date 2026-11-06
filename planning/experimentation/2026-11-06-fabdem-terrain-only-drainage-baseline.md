# FABDEM terrain-routing experiment decisions

[Published vision](../visions/2026-11-06-fabdem-terrain-only-drainage-baseline.md).

Question: can a terrain-only FABDEM V1-2 D8 baseline reproduce drainage coverage, without GRIT guidance?

- Selected pyflwdir 0.5.12 Wang–Liu priority flood (all edge outlets, unlimited fill) and inclusive D8 accumulation. Stock QGIS GRASS 7.8.6RC2 was available; pyflwdir keeps the bounded experiment in one locked uv environment and separate from the archived GRASS investigations. No claim that one method is best.
- Four official TIFF members suffice for 10/20/25 km metric buffered envelopes. Initial two-member input did not cover the eastern corner of the 20 km UTM envelope; processing refused nodata instead of inventing an internal outlet.
- Default GDAL bilinear reprojection varied with output extent even on aligned common cells: maximum difference 1.9403305 m (10 vs 20 km). `tolerance=0` reduced this to 0.0671234 m; fixed `XSCALE=YSCALE=1` removed it exactly. Commands: `process.py ... --buffer-km 10`, `process.py ... --buffer-km 20`, then aligned overlap subtraction. This is a fixed bilinear interpolation kernel, not an inferred hydrological tolerance. Native one-arcsecond spacing is approximately 26 m east-west /31 m north-south here; destination is 30 m UTM. No antialiasing/downsampling average is intended.
- The discarded default-warp 10 km result connected 99.8736% of vector centers at the reference cell, whereas the 20 km result connected one cell. Do not use that preliminary result. With fixed resampling both, and 25 km, connect only the reference cell; a separately reported radius-limited adjusted outlet connects 99.8519%.
- Final common raw elevations are bitwise equal for 10/20 and 20/25 km. Fill, direction and accumulation do change outside the vector basin; inside it they are identical. Both evaluated catchments are identical and contain no processing-boundary cells. This supports local extent stability, not global completeness.

Detailed metrics, commands, figures and preserved input provenance are local at `~/Desktop/hfx-269-investigation/fabdem-baseline/`. Results are in `experiments/fabdem_baseline/README.md`. No rasters, production changes or source-generation attribution are published.
