# FABDEM terrain-only drainage baseline

## Outcome

Derive an independent D8 flow-direction and accumulation pair from FABDEM for the
GRIT diagnostic basin. Determine whether the derived directions contain loops,
how much of the vector basin connects to the diagnostic outlet, and where drainage
paths and accumulation differ from GRIT.

This vision defines the current experiment in [HFX #270](https://github.com/pourpoint-tools/hfx/issues/270).
It supersedes that issue's later plan to reconstruct GRIT's processing: the
maintainer now wants a terrain-only baseline. Delivering the experiment and its
report completes #270, not the originating investigation
[HFX #269](https://github.com/pourpoint-tools/hfx/issues/269). Return findings to both.
Publication of this vision completes neither issue.

The checked GRIT source window already contains the published loops; HFX did not
introduce them. Almost all omitted vector-mask cell centers drain into cycles.
Those observations motivate the comparison, but do not establish the generation
cause. Standard terrain-routing methods can prevent cycles by construction, so a
loop-free baseline is not proof of which GRIT stage caused its loops.

## Inputs and acquisition

Use **FABDEM V1-2**, not the V1-0 input used by the archived GRIT workflow. This is
an independent baseline, not exact GRIT reproduction. Use the
[official dataset](https://doi.org/10.5523/bris.s5hqmjcdj8yo2ibzi9b4ew3sn) and retain
its [README](https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/readme.txt),
source URLs, version, download checksums, and attribution. V1-2 corrects earlier
tile-edge artifacts. Its native spacing is one arcsecond, horizontal CRS WGS84
(EPSG:4326), vertical reference EGM2008 (EPSG:3855); it is not a square 30 m grid.

The saved basin is `~/Desktop/pourpoint-223-investigation/coarse_union.geojson`.
Its approximate bounds are longitude -89.701773 to -89.204948 and latitude
31.628573 to 32.354305. Begin with a **10 km metric buffer** around the basin and
use its rectangular envelope for terrain processing, rather than masking the
DEM to the vector basin. This initially requires:

- `N31W090_FABDEM_V1-2.tif`
- `N32W090_FABDEM_V1-2.tif`

Both are in the official
[`N30W090-N40W080_FABDEM_V1-2.zip`](https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/N30W090-N40W080_FABDEM_V1-2.zip).
Acquisition research verified HTTP range support and the archive directory:
2,300,584,454 bytes for the archive, versus 46,969,818 bytes for the two TIFF
members combined. Retrieve only those members where practical; verify extraction
integrity. No third-party mirror is required. Mosaic and crop locally, preserving
native input files and their georeferencing.

This is **personal, non-commercial testing only**. FABDEM's README states
CC BY-NC-SA 4.0 and requires Copernicus attribution. Preserve those terms with the
inputs and derived data; do not publish rasters or assume unrestricted reuse.
FABDEM already includes corrections and post-processing. Terrain-only means no
added GRIT river guidance, not untouched elevation observations.

## Processing and comparison

Produce a reproducible GIS processing script and a separate evaluation script.
Use project-native tooling and `uv` for Python dependencies and execution. Choose
and document a suitable established terrain-routing method, pit/flat treatment,
software versions, grid, resampling, nodata handling, and accumulation convention.
These reversible engineering choices are left to implementation; no routing
library or exhaustive method comparison is prescribed.

- Derive single-direction D8 and its paired accumulation from the same processed
  terrain. Do not burn GRIT rivers, overlay GRIT directions, or repair GRIT data.
- If projection is needed, reproject the elevation first, then derive directions.
  Never warp derived direction pointers. Compare differing raster grids spatially
  rather than pretending their cell indices are interchangeable.
- Keep the existing diagnostic outlet as the reference: longitude/latitude
  `(-89.32299563991558, 31.639357487809487)`. This is not the original CAMELS gauge.
  Report the actual evaluated cell center and displacement. Any terrain-based
  outlet adjustment is a separate, explicitly reported comparison, not a silent
  replacement of the reference result.
- Use the upstream vector union ending at GRIT source segment `190026628` as the
  comparison mask, not as a constraint on routing through the terrain.
- Check for cycles, invalid directions, and consistency between derived directions
  and accumulation. Measure outlet-connected vector-mask coverage and contributing
  terrain outside that mask. Inspect drainage differences with maps.
- Treat 10 km as an initial extent, not proof of upstream completeness. Check
  boundary-connected contributing paths and repeat with a larger extent to assess
  stability. Expand further if needed; report unresolved edge effects rather than
  interpreting clipped accumulation as a data defect.
- Compare accumulation with explicit units and conventions. GRIT uses
  width-partitioned accumulation at bifurcations, so ordinary D8 accumulation is
  not an equivalent product and differences alone do not establish an error.

Existing GRIT data are at
`https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/`, with
`aux/d8/flow_dir.tif` and `aux/d8/flow_acc.tif`. Prior local evidence is under
`~/Desktop/hfx-269-investigation/` and `~/Desktop/pourpoint-223-investigation/`.
These are local artifacts, not published attachments. Some scripts retain old
`/tmp` paths; inspect before reuse. Do not rerun unrelated software-failure
investigations as a prerequisite to this terrain experiment.

## Deliverables and limits

Deliver the acquired DEM inputs with provenance, reproducible processing and
evaluation scripts, derived rasters, and a concise report with figures. Record
commands, versions, processing settings, extent checks, outlet handling, and
comparison metrics so a fresh agent can reproduce and interpret the results.
Keep large inputs and outputs outside Git; preserve useful research evidence in
a durable local location and link it in the report.

Success is a reproducible comparison, not a predetermined winner or agreement
with the vector basin. Explain what a loop-free result establishes, where the
terrain baseline disagrees with GRIT, and remaining limits. If acquisition or
processing cannot be completed, report the specific blocker without claiming
completion.

This authorizes bounded experimental work only. No production bundle changes,
Pourpoint changes, replacement-data adoption, GRIT repairs, upstream publication,
or cloud infrastructure spending are included. Findings inform a later decision
in #269; they do not resolve its source-generation question automatically.
