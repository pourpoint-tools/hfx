# MERIT catchment connectivity

Add the best defensible catchment connectivity to the MERIT HFX dataset so
Pourpoint can select upstream catchments without MERIT-specific inference in the
engine. This is a known source conversion with exceptions to validate, not an
open-ended research project or a promise of complete global watersheds.

This vision addresses [HFX #285](https://github.com/pourpoint-tools/hfx/issues/285).
Publication leaves that issue open. Its intended completion includes the prepared
bundle and consumer verification described below; code changes alone do not
establish delivery. Production replacement and paid infrastructure require
separate authorization.

## Starting point

The [existing private delivery and verification receipts](https://github.com/pourpoint-tools/hfx/issues/281#issuecomment-6961363360)
contain MERIT Hydro v1.0.1, MERIT-Basins `v0.7/v1.0_bugfix1`, and Coastal Hillslopes.
The river graph and local river–catchment associations are present, but catchment
connectivity is declared `not_supplied`. The missing preparation step is to
translate source river topology into catchment topology, not reconstruct ordinary
river connectivity from scratch.

Delineator demonstrates the conversion. Its
[preparation code](https://github.com/mheberger/delineator/blob/dc147de5f831d7e8b2cc2d5668dac80939c196fc/data_mgmt/2_catchments/aggregate_basins.py#L449-L495)
joins catchment and river records on `COMID`, attaching the river's `NextDownID`
to the catchment. Its
[upstream query](https://github.com/mheberger/delineator/blob/dc147de5f831d7e8b2cc2d5668dac80939c196fc/src/delineator/queries.py#L14-L48)
then traverses those connections. This is concrete reference evidence, not a
requirement to copy its implementation or all its assumptions. In particular,
Delineator filters self-referencing and short rivers and fills unmatched
catchments' downstream values with zero; its coastal preparation also assigns
zero. Those choices do not independently prove natural termination.

Use the current HFX specification at `docs/site/index.md`, the existing production
and verification tools (`scripts/produce_merit.py`, `scripts/verify_merit.py`), and
`scripts/provenance/merit-delivery.json`. Preserve their production capabilities
and source provenance. The current representation is HFX 0.4.0
`spatial-geoparquet-adjacency`, version 1; format version alone does not establish
consumer compatibility.

## Required behavior

- Materialize catchment connections supported by shared source identities and
  river `NextDownID` relationships. Validate the mapping against this source
  release and the delivered local associations. Write immediate relationships
  through the existing HFX contract, with resolvable delivered IDs and consistent
  topology declarations.
- Try seriously to resolve exceptional cases using available source evidence.
  Do not abandon well-supported connections because some areas remain incomplete,
  and do not invent connections to claim global completeness. Leave reversible
  implementation and investigation methods to the implementing agent.
- Account explicitly for the delivery's 149,323 omitted invalid catchment
  polygons, 332 retained catchments whose degenerate rivers were omitted,
  277,371 coastal hillslopes, and any ambiguous or missing associations. These
  delivery counts are evidence to reconcile, not permanent constants or a
  substitute for inspecting the data.
- Preserve source relationships and geometry. Distinguish derived catchment
  relationships from supplied river relationships. Do not silently bridge omitted
  polygons as though the resulting edge were an immediate source connection,
  repair geometry, or turn missing information into a natural outlet. Preserve
  known external edge identities and document cuts and unresolved cases.
- Record derivation rules, source evidence, coverage, exceptions, and limitations
  in bundle provenance. Missing polygon area remains missing even if the network
  path is known; raster refinement within the input polygon footprint cannot
  restore it.

The quality target is the best supported connectivity obtainable from the source,
not the impossible guarantee that every watershed is complete. Keep unsupported
cases explicit and explain what prevents resolving them.

## Consumer behavior and verification

An incomplete watershed must not be presented as an ordinary complete result.
Partial results require explicit opt-in and clear missing-coverage limitations.
Check whether the existing contract and consumer behavior can carry and enforce
this distinction; coordinate any necessary generic consumer support rather than
introducing MERIT-specific inference into Pourpoint. Do not claim that provenance
alone enforces a runtime safeguard.

Demonstrate that source-supported river connections map correctly to catchment
connections and that upstream traversal selects the expected catchments. Cover
ordinary chains and confluences as well as the relevant omissions, unmatched
associations, and terminal/unknown distinctions. Validate the complete produced
graph structurally and report the scope of independent source checks; sampled
hydrologic checks do not establish worldwide completeness.

Verify a candidate bundle with Pourpoint upstream selection and raster refinement
when [MERIT raster reader support](https://github.com/pourpoint-tools/pourpoint/issues/224)
is available. Report supported cases, incomplete cases, and any consumer blockers
separately. Reader support alone is not proof of end-to-end compatibility.
Antimeridian handling is separate consumer work, not an implicit part of this
conversion. Do not block preparing valid connectivity merely because independent
consumer work is unfinished, but leave end-to-end verification explicitly pending.

## Delivery boundaries

Keep source normalization in HFX preparation. This work does not authorize a new
hydrology model, raster-grid changes, missing-polygon reconstruction, blanket
MERIT compatibility claims, or public redistribution of the private dataset.
Preserve the existing published bundle while preparing and verifying changes.
Before paid processing or replacing production, obtain explicit authorization;
the earlier delivery's budget is not renewed by this vision.

Provide reproducible preparation and verification changes, the resulting
connectivity coverage and exception report, and an exact candidate bundle identity
when produced. Coordinate publication and consumer verification to complete #285;
if authorization or consumer support blocks those steps, report the remaining
work rather than closing the issue as delivered.
