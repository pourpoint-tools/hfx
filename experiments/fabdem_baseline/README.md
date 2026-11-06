# FABDEM terrain-only drainage baseline

[Vision](../../planning/visions/2026-11-06-fabdem-terrain-only-drainage-baseline.md) · [HFX #270](https://github.com/pourpoint-tools/hfx/issues/270) · originating [#269](https://github.com/pourpoint-tools/hfx/issues/269)

**At the unchanged diagnostic coordinate, the terrain baseline selects only one cell, not the vector basin.** Its containing cell is centered at `(-89.3229324181, 31.6394414139)`, 11.073 m from the reference `(-89.3229956399, 31.6393574878)`. This result is identical with 10, 20 and 25 km buffers.

A **separate adjusted comparison** selects the greatest accumulation within 300 m of that coordinate. Its center is `(-89.3222329373, 31.6367484548)`, 298.260 m away. It connects 99.8519% of vector-mask centers. This radius-limited choice near the search boundary is neither a unique nor an established correct outlet. It must not replace the reference result silently.

## Results

All three extents give the same following metrics:

| Measurement | Reference cell | Adjusted comparison |
|---|---:|---:|
| Inclusive contributing cells | 1 | 2,587,918 |
| Projected contributing area, km² | 0.0009 | 2,329.1262 |
| Connected vector-mask centers / 2,243,052 | 1 | 2,239,730 |
| Vector coverage | 0.0000446% | 99.8518982% |
| Contributing cells outside vector mask | 0 | 348,188 |
| Outside-mask contributing area, km² | 0 | 313.3692 |
| Contributing processing-boundary cells | 0 | 0 |

All derived direction rasters have **zero cycle cells, zero invalid directions and zero accumulation-count mismatches** under an independent topological evaluator. Cell counts include the cell itself. The entire processing rectangle is valid elevation; missing values cause processing to stop, not become internal outlets. D8 coding is ESRI powers of two; 0 is a terminal, 247 is nodata. The evaluator does not alter directions.

Public GRIT HFX 0.4.0 data were fetched afresh from:
- `https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/aux/d8/flow_dir.tif`
- `https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/aux/d8/flow_acc.tif`

The native EPSG:8857 window starts at zero-based row 145103, column 233826, and is 3091×2658 cells. All 8,215,878 direction **and accumulation** values, masks, CRS and transform match the previous local cache exactly. This is a window verification, not a checksum of the global rasters. Local `inputs/grit-provenance.json` records array hashes and URLs. No original Zenodo tile is needed to rerun this experiment.

In that GRIT window, 2,021 cells are cycle members (whole rectangle, not just basin cycles). The unchanged GRIT outlet connects 616,288 / 2,241,901 native vector centers, about 27.49%, allowing paths outside the polygon. Its published accumulation is 2,017.65234375 km². Prior exhaustive downstream investigation covered continuations outside this window; this bounded evaluator does not repeat or supersede it.

Terrain accumulation is sampled at GRIT native cell-center locations. Different grids are not treated as interchangeable indices, and direction pointers are never warped. Terrain-minus-GRIT accumulation over valid vector centers has minimum −2,017.6514, median approximately 0, 95th percentile 0.0342 and maximum 2,329.0803 km². These pointwise differences are highly sensitive to channel location; the near-zero median mostly describes hillslopes. GRIT's width-partitioned accumulation at bifurcations is **not** ordinary inclusive D8 accumulation. These differences do not establish an accumulation error.

## Method and extent checks

FABDEM V1-2 native inputs retain one-arcsecond WGS84 horizontal georeferencing and EGM2008 vertical reference. Four members (`N31/N32W090` and `N31/N32W089`) were retrieved by HTTP ranges from the official `N30W090-N40W080_FABDEM_V1-2.zip`. ZIP CRC checks verify extraction; SHA-256, source URL, original README and version are retained locally. The second longitude pair covers the expanded envelope's eastern corner.

Elevation is mosaicked, then bilinearly projected to aligned 30 m EPSG:32616 cells. This is the local UTM zone, suitable for metric buffering and routing; km² here means projected grid area, not exact ellipsoidal polygon area. The vector union ending at GRIT source segment 190026628 defines only the comparison mask and buffered rectangular extent, never terrain routing constraints.

Pyflwdir 0.5.12 uses Wang & Liu (2006), DOI `10.1080/13658810500433453`, priority-flood depression filling and D8 routing through pits/flats. All valid edge cells can be outlets; fill depth is unlimited. No GRIT rivers, pointers or repairs enter the derivation. Filled terrain and fill statistics are retained. This is one established method, not an exhaustive conditioning comparison. FABDEM itself already includes corrections/post-processing.

Exact coordinate transformation (`tolerance=0`) and a fixed bilinear sampling kernel (`XSCALE=YSCALE=1`) avoid output-extent-dependent reprojection changes. Native spacing here is approximately 26×31 m, not square 30 m; the selected kernel interpolates elevations, without an area-average downsampling filter. These are explicit engineering choices, not tuned hydrological tolerances. [Experiment notes](../../planning/experimentation/2026-11-06-fabdem-terrain-only-drainage-baseline.md) preserve the rejected default-warp result and reason for changing the settings.

10→20 km and 20→25 km common raw elevations are bitwise identical. Filled elevations, pointers and counts change outside the vector basin (26,651 and 45,738 changed direction cells respectively), but **zero vector-mask cells change** in any of those products. Both outlet-connected masks remain identical; neither includes boundary cells. A preliminary 15 km directory remains locally as rejected pre-fix exploratory output; it is not part of these final results. The final 25 km run is a third check, not evidence of global completeness. Other watersheds near the processing edges remain extent-sensitive and should not be used as complete products.

Machine-readable summary metrics, extent checks and acquisition receipts are retained in [`results/`](results/); large rasters remain local.

Versions: Python 3.13.14, pyflwdir 0.5.12, NumPy 2.3.5, Numba 0.68.1, Rasterio 1.5.2, pyproj 3.8.0, Shapely 2.2.0; remaining packages are pinned in `uv.lock`. The isolated uv project avoids changing HFX production dependencies. Stock QGIS GRASS was available but not needed; no archived GRASS failure reproductions were rerun.

## Local evidence and figures

Durable evidence: `~/Desktop/hfx-269-investigation/fabdem-baseline/`:

- `inputs/`: original TIFFs, official README/license, checksums, basin GeoJSON and public GRIT native-window cache.
- `scripts/`: copies of the reviewed experiment scripts, project file and lock.
- `outputs/buffer-{10,20,25}km/`: elevation, filled elevation, D8, inclusive accumulation, reference/adjusted catchments, settings and metrics.
- `report/`: this report, command logs, extent comparison, provenance notes and figures.

Figures are local, not Git attachments: `accumulation-25km.png` compares drainage networks and pointwise area; `connectivity-25km.png` shows the large area connected only by the adjusted terrain result; `outlet-directions-25km.png` shows native-grid D8 arrows around the reference and adjusted cells, with vector boundary in orange. The terrain and GRIT channels do not occupy identical cells. Arrow coordinates, not pointer codes, are transformed for display. Figures and rasters carry the FABDEM non-commercial/ShareAlike attribution locally.

## Reproduce

From the HFX checkout, using the existing saved basin GeoJSON:

```sh
uv sync --project experiments/fabdem_baseline --frozen --python 3.13
B="$HOME/Desktop/hfx-269-investigation/fabdem-baseline"
uv run --project experiments/fabdem_baseline python experiments/fabdem_baseline/acquire.py "$B/inputs"
# Copy the saved coarse_union.geojson into $B/inputs before continuing.
uv run --project experiments/fabdem_baseline python experiments/fabdem_baseline/acquire_grit.py "$B"
uv run --project experiments/fabdem_baseline python experiments/fabdem_baseline/process.py "$B" "$B/inputs/coarse_union.geojson" --buffer-km 10
uv run --project experiments/fabdem_baseline python experiments/fabdem_baseline/evaluate.py "$B" "$B/inputs/grit-window.npz" --buffer-km 10
```

Repeat the last two commands for `--buffer-km 20` and `25`, then:

```sh
uv run --project experiments/fabdem_baseline python experiments/fabdem_baseline/compare_extents.py "$B"
uv run --project experiments/fabdem_baseline pytest experiments/fabdem_baseline -q
uv run --project experiments/fabdem_baseline ruff check experiments/fabdem_baseline
uv run --project experiments/fabdem_baseline ty check --project experiments/fabdem_baseline --config 'src.include=["experiments/fabdem_baseline"]' experiments/fabdem_baseline
```

The saved basin is a local prerequisite, not a published attachment; this script does not reconstruct the upstream vector union from cloud tables. Range acquisition preserves native files and verifies CRC each time. Raster generation overwrites named derived outputs; preserve prior settings/metrics first when exploring changes.

## Interpretation and limits

Loop-free terrain routing establishes a reproducible internally consistent baseline under this conditioning method. It does **not** identify which GRIT generation stage caused its published loops. The unchanged reference outlet fails to capture the terrain drainage basin, while a nearby, explicitly separate outlet captures nearly all vector centers and additional terrain outside the mask. This highlights grid/channel/outlet correspondence, not a license to snap production outlets automatically.

This is personal non-commercial testing. FABDEM inputs and derivatives are CC BY-NC-SA 4.0, with the required Copernicus attribution retained in `inputs/readme.txt` and `outputs/LICENSE.txt`. Do not publish these rasters or assume the repository MIT license applies to them. No production bundle, Pourpoint, GRIT repair, replacement adoption, upstream publication or cloud infrastructure change is made. Findings return to #269 for later discussion; the experiment does not settle that investigation.
