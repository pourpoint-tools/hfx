# MERIT source geometry restoration

Restore source-backed catchment coverage in a new global MERIT HFX bundle, verify
it with Pourpoint, and promote it as the current private MERIT delivery. Work on
the UZH VM. Keep existing immutable bundles untouched.

This vision addresses [HFX #292](https://github.com/pourpoint-tools/hfx/issues/292).
Publishing this document leaves #292 open. Completing the work requires the
verified global bundle, promotion, reproducible tooling and provenance, and a
findings handoff to the [originating delineation discussion](https://github.com/CooperBigFoot/hydrologeez/issues/203).
Code or a diagnostic repair alone does not complete delivery.

## Starting evidence

Use the [investigation findings](https://github.com/pourpoint-tools/hfx/issues/292#issuecomment-7010288963)
and their retained evidence, not assumptions that every invalid polygon is
repairable. The current promoted delivery from [#285](https://github.com/pourpoint-tools/hfx/issues/285) is:

- `s3://pourpoint-hfx/merit/hfx-0.4.0/merit-hydro-v1.0.1-basins-v1.0-bugfix1-20261118-connectivity-candidate/manifest.json`
- Manifest SHA-256: `9d3d48af7fe6f1dbd9661621cba307d662e6a35eee65f7de3911a6c31edecb2d`.
- Private Hetzner endpoint `https://fsn1.your-objectstorage.com`, profile
  `hetzner-pourpoint`; preserve private access.

Three of six CAMELS smoke-test gauges stopped on upstream coverage cuts. The
investigation traced these to four omitted MERIT source catchments:

| Gauge | Missing source COMIDs | Source partition |
| --- | --- | --- |
| 01013500 | 72042012, 72044278 | pfaf_72 |
| 02053200 | 73010886 | pfaf_73 |
| 05444000 | 74024994 | pfaf_74 |

These have delivered rivers and complete source river topology. Raw polygons
contain self-touching rings at raster-cell corners. Diagnostic ring repairs
preserved their exact vertices and undirected edges, changed geodesic area only
at floating-point scale, and introduced no neighbor interior overlap. Restoring
them also reconnects valid upstream catchments stranded behind the gaps.
The three restored diagnostic refinements matched unmasked native D8 traces
pixel-for-pixel. This did not test a rebuilt bundle through strict selection.

The same repair checks passed for all 9,926 invalid polygons in the three audited
partitions. All 149,323 globally omitted polygons have reported ring
self-intersections at raster-cell corners, but this does not prove the remaining
polygons have only that defect or pass the restoration checks. Global acceptance
must be established feature by feature, not extrapolated from the regional audit.

Evidence includes the smoke-test directory
`~/Desktop/uzh-results/2026-11-19-merit-camels-smoke`, mirrored on UZH at
`/data/2026-11-19-merit-camels-smoke`. The investigation comment identifies raw
source hashes, the diagnostic script, and `findings_receipt.json` under
`~/.gemini/antigravity-cli/brain/48c42b77-e76a-41c8-857f-a683f6bdb973/scratch/`.
Preserve durable reproduction evidence with the delivery rather than relying
solely on an agent scratch directory.

## Approved restoration policy

Explicitly document a narrow MERIT ring-normalization rule in the current HFX
specification at `docs/site/`, production tooling, and delivery provenance. This
is not permission for arbitrary geometry cleanup or a new HFX format version.

- Admit a repaired source polygon only when the documented checks establish
  valid Polygon/MultiPolygon output, exact source vertex and undirected edge
  preservation, geodesic area preservation within justified numerical tolerance,
  and no introduced interior overlaps with neighboring catchments. Retain holes
  and disconnected components; do not fill, simplify, move, or discard source land.
- Verify applicability across all global MERIT catchment partitions. Regional
  staging is allowed, but the deliverable is a global bundle with the same
  acceptance policy everywhere. Include cross-partition neighbors where relevant
  to overlap validation; do not mistake partition-local checks for global proof.
- Preserve original geometries, source identities, raw failure reasons, and
  per-feature repair evidence, including the chosen algorithm, versions, checks,
  and measured changes. Document changes in ring/part representation explicitly.
  The exact provenance layout is left to implementation; each repair must be
  auditable and originals recoverable under the applicable source license.
- Leave already valid source geometries unchanged. Features outside the approved
  repair rule or failing any acceptance check remain omitted with explicit reasons.
  Do not progressively relax checks to achieve a target restoration count.
- Choose and document the exact repair method and numerical tolerances using the
  evidence. The diagnostic agreement between Shapely structure and linework
  methods is useful evidence, not a mandate for a particular implementation.

Coordinate shared specification edits with
[TDX catchment connectivity](2026-11-19-tdx-catchment-connectivity.md), tracked in
[#289](https://github.com/pourpoint-tools/hfx/issues/289). That work contracts
exhaustively verified zero-length, zero-area virtual links in the catchment graph
while preserving the river graph. MERIT's missing units have real drainage area:
**restore their polygons, never bridge their gaps by graph contraction**. Preserve
TDX's distinct approved exception without broadening it or overwriting concurrent
contract changes. Source normalization stays in HFX, not in Pourpoint.

## Build and delivery

Update reproducible MERIT production and verification tools, tests, and current
documentation. Use the verified source releases and provenance from the existing
delivery. Re-acquire and hash-check originals where required. Preserve licensing,
source attributes, identities, and unchanged river geometry and topology.

Rebuild affected catchment features, local river–catchment associations,
connectivity, upstream completeness indicators, lookups, and artifact bindings.
Newly restored catchments must participate through their actual source-backed
relationships. Rebuild dependent river associations even where river geometry
is unchanged. Do not assume HFX IDs remain stable across bundles; preserve source
identity and ensure all references resolve within the new bundle.

Produce internally consistent staging and final manifests, retaining the base
manifest identity and transformation history in provenance. Reuse verified
unchanged artifacts only where byte identity and bindings remain valid. Copy the
unchanged verified native direction and accumulation rasters into the new bundle;
do not reprocess them. The current external-raster contract requires public HTTPS
availability and is not a private cross-prefix reference mechanism. Verified
server-side copying is appropriate where supported.

Publish at a new unused immutable private prefix under
`s3://pourpoint-hfx/merit/hfx-0.4.0/`. Verify artifacts before publishing the root
manifest last. Preserve all existing deliveries, their manifests, and access
policies. Check current bucket headroom, including concurrent TDX delivery, before
uploading; historical allowance estimates are not current authorization for extra
paid capacity.

## Verification and promotion

The implementing agent owns the necessary checks and promotion; a second human
promotion approval is not required after successful verification. At minimum:

- Independently verify source-to-output fidelity and every accepted repair against
  the approved invariants. Report exhaustive counts by partition: attempted,
  accepted, rejected, remaining omissions, and validation coverage. Do not label
  sampled evidence as exhaustive or assume all 149,323 polygons will be restored.
- Validate bundle schemas and artifact bindings, associations in both directions,
  graph integrity and acyclicity, and completeness/cut propagation. Confirm actual
  counts and distinguish restored polygons from unresolved source/topology limits.
- Rerun all six original smoke-test gauges through the actual rebuilt bundle:
  `01013500`, `02053200`, `05444000`, `08104900`, `09066000`, `11141280`.
  Use the retained USGS coordinates with explicit NAD83-to-EPSG:4326 handling,
  level 0, the original nearby-candidate search, 100 m maximum vector snap,
  3-cell raster search window, `refine=True`, and `allow_incomplete=False`.
  Use Pourpoint `28b4ea153dd3ab7dfb300de3827e61dc1951c766` or a documented newer
  compatible revision. Use fresh/disabled caches to avoid old-bundle results.
- Require strict complete-closure success for the three restored cases and no
  unexplained regressions for the other three. Compare boundaries, areas, outlet
  cells, and CAMELS agreement with retained evidence, including diagnostic
  restored results. Investigate differences rather than retuning parameters or
  replacing gauges. Existing diagnostic IoUs are evidence, not invented universal
  accuracy thresholds.
- Exercise genuine remaining incomplete or unverified closures: they must still
  reject by default and report missing coverage under explicit partial opt-in.
  Never use opt-in results as successful complete delivery checks.
- Verify uploaded identities and authenticated consumer range access against the
  exact published manifest. Retain independent verification and upload/readback
  receipts; successful upload alone does not establish integrity or usability.

If required checks fail, resolve the failure or report the blocker without
promotion. Expected, documented residual omissions do not by themselves prevent
promotion when the bundle correctly represents them and the agreed checks pass.

After successful verification, record the exact new URI and manifest digest as
the current private MERIT delivery, update relevant delivery documentation and
known consumer selectors where present, and verify those references. Do not move
objects, mutate old bundles, change access policy, or claim that an external
application configuration was updated without evidence.

Report the cause, global restoration results, remaining limitations, verification,
and promotion to #292 and the originating hydrologeez #203 discussion before its
full delineation vision is written. Complete #292 only when this delivery and
handoff are done; do not close the broader delineation discussion.

## Execution boundaries and remaining limits

Follow `../hydrologeez/README.md` for UZH connection, host-key checks, workload
inspection, and storage safety. Run on `ubuntu@172.23.211.16` (hostname `nik`).
Inspect live CPU, RAM, disk, and active workloads before work; README capacities
are snapshots. Use bounded-memory processing and put large data under `/data` or
`/data2`, not the system disk. Preserve existing datasets and workloads. Do not
format disks, reboot, expand infrastructure, or provision paid compute/storage.
Use `uv` and the repository's development guidelines.

This work does not authorize a 4,492-gauge national delineation run, native
raster-only Pourpoint implementation, arbitrary geometry repair, or suppression
of completeness safeguards. If consumer verification reveals a generic Pourpoint
blocker, report and coordinate it rather than adding MERIT-specific inference or
silently expanding this work into engine changes.

Unassociated coastal hillslopes, missing-source catchments, omitted degenerate
rivers, and antimeridian limitations remain explicit unless independently resolved
within the approved policy. Restoring polygons does not establish missing river
topology. Promotion certifies the verified delivery and its declared limits, not
universal watershed completeness, national success rates, or hydrological truth.
