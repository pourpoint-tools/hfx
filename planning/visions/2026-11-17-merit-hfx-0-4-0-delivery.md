# MERIT HFX 0.4.0 delivery

Deliver MERIT Hydro native rasters and MERIT Basins rivers, catchments and supplied
relationships for private use through HFX 0.4.0. Aim for global source coverage,
but a verified partial delivery is acceptable if acquisition, source defects or
the budget prevent completeness. Account for every omission and its reason;
never present a sample or incomplete product as a complete global fabric.

This is the implementation handoff for [HFX #281](https://github.com/pourpoint-tools/hfx/issues/281),
following the [bounded trial and findings](https://github.com/pourpoint-tools/hfx/issues/282#issuecomment-6945979797).
Successful production delivery and its verified consumer handoff complete #281;
publication of this vision does not. Start implementation separately. No production
objects or contract files were changed during the trial.

## Scope and source fidelity

Use the current contract at [docs/site/index.md](../../docs/site/index.md):
`format_version=0.4.0`, representation `spatial-geoparquet-adjacency`, version 1.
Version number alone is not representation compatibility. Keep normalization
outside Pourpoint and preserve source provenance, attributes, types and meanings.

Target the currently documented MERIT Hydro v1.0.1 direction (`dir`) and upstream
area (`upa`, km²) products, plus MERIT Basins coastline-bugfix1 vectors. Pin actual
acquired release identities and source hashes rather than trusting filenames.
Confirm their relationship and document differences; the trial established
representation preservation, not scientific equivalence of vector and raster
hydrology. Other Hydro raster products are not required for this delivery.
Include useful supplied coastal hillslopes where practical, separately from the
main channelized catchments. Do not invent channels or derive missing catchments.

The Basins source uses a 25 km² channel threshold. `pfaf_level_01` and
`pfaf_level_02` are alternative geographic partitions of the same features,
not two hierarchy levels. Acquire one organization. Preserve source COMIDs and
namespaces separately from dataset-local HFX IDs; establish uniqueness before
reusing COMIDs globally. Preserve native stream order as an attribute, not a tier.
A single river/catchment tier with null parents is sufficient unless source
evidence establishes otherwise. Do not invent source nodes, feature outlets,
parents, flow edges or geometry orientation.

Check river/catchment associations by native identity, not row order or spatial
intersection. Validate supplied upstream and downstream relationships against
each other, including references, multiplicities, cycles and subset boundaries.
Preserve source area quantities separately from derived geodesic local polygon
area. Native river `uparea` is not automatically HFX upstream polygon-union area.
Keep unsupported upstream area and outlets null. Catchment connectivity may remain
`not_supplied`; any derivation from river associations needs explicit evidence and
provenance, not an implicit assumption.

Do not repair, clip, simplify, reverse or dissolve source geometry silently.
Invalid catchment polygons may be omitted under the accepted partial-delivery
policy. Record their source identities and failure reasons. Do not delete an
otherwise usable river merely because its polygon was omitted. Distinguish
omitted/unavailable polygons from naturally absent polygons, and preserve
association evidence in provenance when the target feature cannot be delivered.
If a region or feature is excluded, retain native references in provenance and
mark cut relationships as boundaries, never natural terminals. Keep the delivered
HFX relationship tables referentially valid. Partial delivery must still satisfy
the contract within its explicitly declared coverage.

## Agreed raster normalization

Preserve the native grid and raw direction codes. Explicitly set direction
**nodata to 247**, the documented unsigned ocean code. This is an intentional
metadata and validity-mask normalization, not byte-for-byte metadata preservation.
Preserve decoded direction values and accumulation float32 bit patterns, including
source invalid values; retain upstream-area nodata and units. Verify the declared
mask transformation independently against source semantics. Do not warp, resample,
recalculate flow, or reconstruct transforms from tile names.

Keep river-mouth code **0** and inland-depression code **255** distinct, along
with the ordinary ESRI D8 codes. Clarify their meaning in the HFX raster contract
and coordinate the exact interpretation with Pourpoint before claiming consumer
compatibility. The maintainer approved this source-preserving direction, not a
silent conversion to GRASS or a collapse of mouth/depression distinctions. Keep
this change narrow; do not redesign raster storage or implement the Pourpoint
engine in this effort. Record the consumer contract handoff and any outstanding
reader support explicitly; do not claim an untested reader can consume the result.
If an unexpected source convention contradicts this mapping, isolate the affected
coverage and resolve it rather than guessing.

Acquire a bounded original depression example if practical, since the Iceland
tile did not contain 255. Verify allowed codes and nodata throughout delivered
coverage, including coastline and tile seams. Document remaining untested cases.
Use lossless, internally tiled, bounded-range-readable GeoTIFFs with no required
overviews. Verify TIFF layout and actual selected remote reads; a successful local
open is not proof of range-read suitability. Tile dimensions, compression settings
and packaging boundaries are implementation choices, not inherited global rules
from the trial's 256×256 tiles. Respect HFX coverage rules: do not promise arbitrary
multi-tile reads if consumers require one entry to cover a request. Declare exact
raster coverage, coherent overlap and request-boundary limitations.

## Evidence to reuse and limits to resolve

The trial ran entirely on Hetzner using HFX revision
`0a4cde91bad65616e773e0b6704ca27936034f77`. Its temporary server and IPs were deleted;
no live server or source sample should be assumed to remain.

- Iceland Basins pfaf27 supplied 1,973 rivers and 1,973 uniquely matching
  catchments in 24,314,903 downloaded bytes. All 1,651 river edges were reciprocal
  and acyclic; zero was the absent-link sentinel in that region. There were
  **117 invalid catchment polygons**, all ring self-intersections. None repaired.
  Catchment `.prj` was absent; WGS84 was an explicit trial assumption, while the
  river `.prj` supplied it. Establish and document the catchment CRS from reliable
  source evidence before production; do not promote the trial assumption to fact.
- Complete valid basin `27000035 → 27000015 ← 27000058` passed normalization,
  current HFX spatial build/verification and independent Fiona-original-to-final
  readback. Six feature records, two river edges and three local associations
  were preserved. Native DBF integer/real fields exported as int64/float64;
  actual field names include `strmDrop_t` and `slope_taud`. Geometry checks compared
  reader-exported WKB coordinates and part order, not original shapefile bytes.
  All nine source sidecar hashes stayed unchanged. This does not establish
  global uniqueness, validity, disjointness, topology or performance.
- Hydro archives `dir_n60w030.tar` and `upa_n60w030.tar` totaled 598,753,280 bytes.
  SHA-256 values are respectively
  `7eb9a17044067216cdc721669db8fb46887ebb54c424b389d806cf3afce5a58f`
  and `6eabef181f3f5cfc644360faa2418ac9d7e24689d9b527975dcdaa0c72dd0779`.
  Tile `n65w020` has a half-cell-shifted native origin, not whole-degree edges.
  A 720×720 crop covering the vector basin passed exact cells, masks and grid
  checks after lossless repackaging. Across all 36,000,000 native tile cells,
  direction247 positions equaled upstream-area invalid cells exactly. Direction
  lacked a nodata tag and contained mouths but no depressions. The unchanged
  trial rasters were **not HFX-conformant**; tagging nodata and coordinated terminal
  semantics remain production work.
- The source listing contained 57 direction and 57 upstream-area archives,
  approximately **54.24 GB compressed** from rounded display sizes. This is not
  an exact transfer total, global output-size estimate or runtime prediction.
  Global vectors, hillslopes and extracted/temp storage remain to estimate.

Compact scripts and receipts were retained outside Git at
`~/.prime/agent/session-artifacts/01a1e8ec-b78c-73c1-8d60-8918df30db09/merit-trial/`.
They contain no source rasters/shapefiles or credentials. Use them if available;
the linked findings and identified originals remain the durable starting point.
Preserve useful production scripts and compact reproducible evidence in the
repository during implementation, not only in temporary server directories.

Reuse `hfx.vector_storage.write_vectors` for normalized staging and
`hfx.spatial_bundle.build` / `verify` for final storage where appropriate. The
latter verifies normalized inputs against final storage; it does not establish
original-source fidelity. Add a small MERIT normalizer and independent checks,
not a general adapter framework or supported reader library. Use strongly typed
Python, uv, Ruff, ty and focused test-first development.

## Acquisition, hosting and resources

All source acquisition, wrangling and verification must run on a **Hetzner server**,
not the laptop. Keep bulk data and credentials outside Git. Registration for
MERIT Hydro is complete under **CC BY-NC 4.0**; obtain the credential through the
maintainer's private channel if it is no longer available. Never place it in
issues, scripts, URLs in reports or logs. The official source pages are
[MERIT Hydro](https://global-hydrodynamics.github.io/MERIT_Hydro/) and
[MERIT Basins](https://www.reachhydro.org/home/params/merit-basins).

Use the existing R2 bucket only if the exact releases' terms permit the intended
redistribution. Otherwise deliver to the maintainer's private Hetzner storage
alongside TDX 0.4.0. The Hydro page requests permission for complete original-format
redistribution; do not assume format conversion removes all licensing concerns.
Basins has separate applicable terms. Confirm attribution, license and access
obligations for both products; private hosting is not a waiver. If public delivery
permission remains uncertain, use the approved private fallback, not a public
upload while waiting. Identify the exact existing bucket/prefix and access method
before writing; use a new disjoint MERIT prefix, preserve all TDX/GRIT and other
existing objects, and do not change bucket policy. Never publish signed URLs or
credentials in the handoff.

The **US$10 total acquisition/conversion cap includes the trial**, temporary
compute, IPv4, working storage and applicable transfer charges. The trial's estimated
compute/IP cost was about US$0.033 at observed gross rates and EUR/USD1.1014;
reconcile actual charges before allocating the remaining budget. This is not a
fresh US$10 allowance. No standing resources remain from the trial.

Recheck current prices, tax, currency conversion, capacity and existing storage
before provisioning. CPX22 in hel1 worked for the tiny trial (2 shared cores,
4 GB RAM, 80 GB disk; observed compute EUR0.0285384/hour gross plus
EUR0.0008648/hour IPv4), but this is not a mandated production machine or a
capacity/runtime guarantee. Estimate and monitor memory, scratch/output disk,
transfer, cumulative spend and completion time. Process incrementally and release
verified temporary inputs when safe instead of staging the entire globe by default.
Reserve enough budget to verify, save receipts and clean up. Track an early
stop point: the budget is not provider-enforced. Stop before exceeding it and
handoff verified partial coverage rather than silently expanding spend.

Existing hosting is the intended durable destination. Ongoing hosting charges
were not given an unlimited budget; report expected retention cost and obtain
approval for any new recurring obligation before committing to it. Delete this
job's temporary server, volumes and owned IPs after successful preservation and
verification; powering off is not cleanup. If blocked, save resumable evidence
and arrange cleanup within the cap, not indefinite billable retention.

## Completion and consumer handoff

Fully verify each delivered artifact against original sources and declared
transformations with bounded memory. Check identities, native attributes, geometry,
relationships, null/unknown meanings, raster values, masks and grids. Retain
source inventories, source/output hashes, commands, environment/code versions,
coverage and exclusion receipts. Selected read demonstrations supplement, not
replace, complete producer fidelity checks.

Upload verified immutable artifacts first and publish the root manifest last.
Verify hosted identities and focused range reads through the intended public or
private access path, including vector records/relationships and paired native
raster windows. Document exact coverage and limitations, rather than claiming
uniform global availability or source hydrological certification.

Return a handoff on #281 with the exact destination and manifest identity,
representation, release/license/provenance records, delivered and omitted products
and regions, verification receipts, outstanding consumer support, actual/estimated
spend and confirmed temporary-resource cleanup. A useful verified partial delivery
can complete the agreed outcome; a blocked run with no usable delivery cannot.
Do not close #281 for this vision or the successful tiny trial alone.
