# TDX catchment connectivity

Deliver a fully HFX 0.4.0 conformant private TDX Hydro vector bundle with
evidence-backed catchment connectivity and explicit completeness metadata, so
Pourpoint can delineate upstream watersheds and enforce generic completeness
safeguards without TDX-specific engine changes.

This vision addresses [HFX #289](https://github.com/pourpoint-tools/hfx/issues/289).
Publication of this vision leaves #289 open. Completing #289 requires the
independently verified private bundle at a new immutable S3 prefix, updated
documentation and reproducible tooling, and end-to-end Pourpoint verification;
code changes alone do not establish delivery.

## Starting point and verified findings

The existing private full-coverage bundle is:

- `s3://pourpoint-hfx/tdx/hfx-0.4.0/NGA-TDX-Hydro-20230126-full62-20261027/`
  (endpoint `https://fsn1.your-objectstorage.com`, AWS CLI profile
  `hetzner-pourpoint`).
- HFX `0.4.0`, `representation: "spatial-geoparquet-adjacency"`,
  `representation_version: 1`, vector-only (`auxiliary: []`), all 62 processing
  basins supplied.
- `15,940,505` rivers (`rivers.parquet`, 15.9 GB) forming a verified acyclic
  binary tree (`upstream_degree` in `{0, 2}`, `downstream_degree` in `{0, 1}`,
  no unresolved cross-basin links), and `16,864,998` catchment rows
  (`catchments.parquet`, 60.0 GB) where every row has a nonnull `streamID`
  matching a river `LINKNO` in the same processing basin.
- Level 0 `catchment_connectivity` is currently `not_supplied`, causing native
  Pourpoint `delineate` at CAMELS gauge `01013500` (`-68.58264, 47.23739`) to
  fail with `selection.connectivity_unavailable`.

Prior investigation and unanswered NGA inquiry context live in
[HFX #268](https://github.com/pourpoint-tools/hfx/issues/268). Prior delivery
history is recorded in
`planning/visions/2026-10-22-tdx-hfx-0-4-0-private-vector-delivery.md`, and the
MERIT connectivity precedent is in
`planning/visions/2026-11-18-merit-catchment-connectivity.md`.

Global audit of `rivers-relationships.arrow`, `source-metadata.json`, and
targeted probes of `rivers.parquet`, `catchments.parquet`, and raw GeoPackages
(`5020055870`, `2020057170`) established two distinct source patterns:

1. **Multi-polygon `streamID` rows (`334,017` rivers, `928,570` extra rows)
   occur in only 2 of the 62 basins:**
   - 60 basins use GeoPackage geometry type `GEOMETRY` (`Polygon` and
     `MultiPolygon`) and have strictly at most 1 catchment row per river.
   - Only basins `1020018110` (`122,259` multi-polygon rivers) and `5020049720`
     (`211,758` multi-polygon rivers) were exported by NGA with GeoPackage type
     `POLYGON`, splitting multi-part catchments into separate `POLYGON` rows
     sharing the same `streamID`. Probed extra rows in both basins are
     corner-touching (`distance = 0.0`) 1-pixel (`~150 m²`) and 2-pixel
     (`~300 m²`) fragments of the same river's drainage unit.
2. **Polygonless rivers (`4,077` rivers, `0.026%`) split into two causes:**
   - **13 real rivers with `LINKNO == 0`** (`12` headwaters and `1` isolated):
     all 13 have real lengths (`1.9–8.9 km`) and real contributing areas
     (`5.7–17.5 km²`), likely dropped by NGA treating value `0` as NoData during
     raster-to-polygon export. The `12` headwaters sit upstream of `5,791`
     rivers (`0.036%`) and represent genuine missing catchment area.
   - **4,064 internal rivers** (all with `upstream_degree == 2` and
     `downstream_degree == 1`): in a 40-river global sample plus 2/2 in
     `2020057170`, 100% had `Length == 0.0`, `geom.length == 0.0`, and
     `DSContArea == USContArea`. Because TauDEM represents confluences as
     strictly binary junctions (`USLINKNO1`/`USLINKNO2`), it inserts a
     zero-length, zero-area virtual link wherever 3 or more tributaries meet at
     one pixel; that virtual link covers no cells and has no polygon. The
     40-river sample is **not** exhaustive proof and must be verified across all
     `4,064` links during production.

## Required behavior and derivation rules

- **Collect single-part `POLYGON` fragments in `1020018110` and `5020049720`
  into one `MultiPolygon` per `(processing_basin, streamID)`:**
  - Pack the original polygon parts into one `MultiPolygon` without dissolving,
    simplifying, repairing, or altering vertex coordinates.
  - Exhaustively verify OGC validity (`is_valid`) across all `334,017` grouped
    `MultiPolygon`s; if any edge-sharing or overlapping fragment pair produces
    an invalid `MultiPolygon`, stop and report for an explicit contract decision
    rather than silently repairing or dissolving.
  - Sum exact WGS84 geodesic `area_km2` across parts and compute the
    outward-covering `bbox`.
  - Preserve original `fid` provenance for grouped rows (exact column/metadata
    representation is left to implementation, provided every original `fid`
    remains traceable, `source_id` stays unique within `source_namespace`, and
    any schema change is documented in `docs/site/vector-data.md` and `NOTICE`).
  - Confirm that all 62 basins now have strictly at most 1 catchment per river
    (`15,936,428` catchments, `tree` topology).
  - Update the TDX-specific note in `docs/site/spatial-layout.md` so the
    specification accurately describes collecting NGA's single-part `POLYGON`
    fragments without implying arbitrary feature merging.
- **Contract only exhaustively verified zero-length, zero-area internal virtual
  confluence stubs:**
  - During production, exhaustively check every polygonless river.
  - Contract (bypass) a polygonless river in the **catchment** flow graph if and
    only if all five conditions hold (failing closed if any attribute is null or
    non-finite): `Length == 0`, `geom.length == 0`, `DSContArea == USContArea`,
    `upstream_degree >= 1`, and `downstream_degree == 1`. Handle chained virtual
    stubs transitively if any exist.
  - Wire the upstream catchments of a contracted stub directly into its
    downstream catchment (producing a 3-way or k-way confluence in the catchment
    tree). Leave the **river** flow graph (`rivers-relationships.arrow`
    `upstream_ids`/`downstream_ids`) completely untouched.
  - Never contract `LINKNO == 0` rivers or any polygonless river that fails any
    of the five checks. Transitive contraction must propagate both upstream
    catchments and upstream boundary cuts across contracted stubs onto the
    downstream catchment (and must never bridge across a downstream
    non-contracted polygonless river).
  - Set `upstream_boundary_cuts` (`int32`, nullable) to `0` (non-null) for every
    delivered catchment with no incoming cuts, and to the positive count of
    non-contracted incoming gaps otherwise.
- **Publish explicit completeness metadata and provenance:**
  - Declare level 0 `catchment_connectivity` as `status: "supplied"`,
    `completeness: "partial"`, `upstream_cuts_column: "upstream_boundary_cuts"`,
    and an explicit `derivation` string under `docs/site/bundle.md` and
    `docs/site/vector-data.md`.
  - Amend `docs/site/vector-data.md` (flow connectivity section) to explicitly
    distinguish exhaustively verified zero-length, zero-area virtual confluence
    stubs from omitted/missing polygons with drainage area, and document the
    exhaustive audit counts (grouped basins, contracted zero-area stubs,
    non-contracted `LINKNO == 0` omissions, any unexpected nonzero polygonless
    rivers, affected downstream closures) and remaining limitations in
    `source-metadata.json`, `NOTICE`, and a committed provenance receipt under
    `scripts/provenance/`.
  - Keep `scripts/produce_tdx.py`, `scripts/verify_tdx.py`, and tests updated so
    raw-source normalization remains reproducible alongside the derived-bundle
    workflow.

## Consumer boundary and verification

- **No Pourpoint changes are in scope.** The updated bundle must work with
  Pourpoint (`pourpoint-tools/pourpoint` `main`, commit
  `28b4ea153dd3ab7dfb300de3827e61dc1951c766` or newer carrying
  [Pourpoint #227](https://github.com/pourpoint-tools/pourpoint/issues/227))
  as-is. If verification uncovers a Pourpoint bug or missing generic extension,
  open a GitHub issue in `pourpoint-tools/pourpoint` and stop implementation.
- **Bundle integrity verification:** Independently verify exact feature/geometry
  preservation (including grouped `MultiPolygon` parts against base polygons),
  exact river-to-catchment `local_catchment_ids`/`local_river_ids` rebuild
  (required because catchment IDs change after grouping), catchment graph
  acyclicity and completeness counts, full spatial bundle verification, and full
  authenticated S3 readback before and after publishing `manifest.json` last.
- **Pourpoint end-to-end verification (`cache=False`):**
  - CAMELS gauge `01013500` (`-68.58264, 47.23739`), level 0, `refine=False`:
    succeeds with complete upstream delineation (`coverage.status == "complete"`).
  - Witness cases covering: an ordinary complete confluence, a contracted 3-way
    confluence (`upstream_boundary_cuts == 0`, returning all 3 upstream
    branches), grouped `MultiPolygon` catchments in both `1020018110` and
    `5020049720`, and a downstream closure below a missing `LINKNO == 0`
    headwater (rejected by default with `selection.incomplete_coverage`;
    returned with cut diagnostics under `allow_incomplete=True`).

## Execution, compute, and delivery boundaries

- **Derive from the verified live `full62-20261027` bundle:** Produce a fresh
  normalized staging manifest (`source-manifest.json` bound by
  `input_manifest_sha256`) whose semantic declarations (`unit_count`,
  `catchment_connectivity`, `categorical_columns`, etc.) match the final spatial
  `manifest.json`, while recording the base `full62-20261027` manifest SHA-256
  in provenance. Reuse unchanged river feature and lookup artifacts
  (`rivers.parquet`, `rivers-lookup.arrow`, `rivers-lookup-directory.json`) via
  verified S3 server-side copy where exact byte identity and provenance bindings
  hold, and rebuild `catchments.parquet`, catchment lookup, both relationship
  files (`catchments-relationships*` and `rivers-relationships*`), provenance,
  and manifests.
- **Compute:** Run the global build on the UZH VM (`ubuntu@172.23.211.16`,
  hostname `nik`, 15 vCPUs, ~58 GiB RAM, 4.9 TiB `/data` disk; follow the
  connection, workload check, and `/data` preservation rules in
  `../hydrologeez/README.md`). Paid cloud compute is not authorized.
- **Publication:** Publish to a **new unused immutable private prefix** under
  `s3://pourpoint-hfx/tdx/hfx-0.4.0/` using profile `hetzner-pourpoint`. Do not
  modify or delete `NGA-TDX-Hydro-20230126-full62-20261027/`, earlier bundles,
  or the raw archive. Keep private access policy intact; do not enable anonymous
  access.
- **Exclusions:** This is a dataset update under HFX 0.4.0, not a new HFX format
  version. Connectivity completeness does not certify vector–raster consistency
  ([HFX #264](https://github.com/pourpoint-tools/hfx/issues/264)).

## Risks and uncertainties

- Any of the `4,064` internal polygonless rivers outside the 40-river sample
  that fails the zero-length/zero-area checks will remain an explicit cut and
  increase the fraction of `partial` downstream closures; report the exact
  exhaustive counts honestly rather than relaxing the rule.
- Rewriting `catchments.parquet` (60 GB) and rebuilding global relationship
  indices requires bounded-memory streaming; verify disk and memory headroom on
  `/data` before running the full pass.
- Re-checking S3 bucket storage headroom before upload is required so the new
  ~76 GB bundle fits cleanly alongside existing deliveries.
