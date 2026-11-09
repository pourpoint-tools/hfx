# GRIT cycle stream-mapping hypothesis

Investigate whether GRIT’s known closed D8 cycles are associated with the stream cells where vector directions were enforced. Return evidence that lets the maintainer decide how to proceed with Michel, not a raster repair.

This bounded experiment completes [#277](https://github.com/pourpoint-tools/hfx/issues/277). Return findings to the originating discussion [#269](https://github.com/pourpoint-tools/hfx/issues/269), which stays open. Publication of this vision completes neither issue.

## Question and context

Michel reports that, when preparing the drainage raster for external use, he forced directions inside raster streams to follow vector segment directions set using RivGraph, following the mainstem. He suspects imprecise raster-cell/vector mapping introduced the cycles. His immediate diagnostic question is whether all cycle-cell clusters touch stream cells. He also reports that published reach lines and catchment boundaries were slightly straightened, although segments and reach nodes are consistent. This is a reported processing explanation and causal hypothesis, not yet a demonstrated cause.

The investigated upstream footprint terminates at GRIT source segment `190026628`. Existing exhaustive tracing found 294 distinct cycles receiving 1,625,331 of 1,625,680 omitted vector-mask cell centers (99.98%). The selected terminal-channel diagnostic cell is global zero-based row 148032, column 235090, coordinate `(-89.32299563991558, 31.639357487809487)`. It is not the original CAMELS gauge or the declared vector outlet. Independent forward/reverse tracing agrees with Pourpoint. All 8,215,878 direction cells in the checked window match the original GRIT v1.0 source tile; HFX preparation did not introduce those loops.

Use #269 and its later comments as the evidence history, including the completed FABDEM comparison. That comparison is not a reproduction of GRIT processing and its nearby high-coverage outlet is not an established correct outlet. Do not repeat it as part of this experiment.

## Bounded investigation

Start with the 294 known cycles in this basin and the public river geometry. Recover and verify the existing cycle inventory rather than silently substituting a different population. Preserve native row/column identities and distinguish cycle-member cells from all cells draining into them. Report whether cycles have members outside the vector mask and any limits on their recovered extent.

Evaluate whether each cycle cluster intersects published river geometry and characterize distances for non-intersecting clusters. Distinguish at least one member touching a river from every member doing so. State the cell-footprint/intersection convention, coordinate system, distance units, river selection and coverage. Include spatially relevant river segments even where they cross the vector-mask boundary. Buffer or proximity comparisons may describe uncertainty; they do not establish membership in the original stream mask.

Inspect representative native-cell configurations, including non-overlaps and ambiguous junctions where present. Show raw directions, cycle members and river geometry together. Determine what the spatial evidence supports or challenges, without treating correlation as proof of processing history.

Inspect available provenance and published processing evidence proportionately. If original stream-cell mapping or before-enforcement directions are available, use them to distinguish original stream membership from published-line overlap and strengthen causal attribution. If unavailable, identify the specific missing evidence and the question it would resolve. Do not guess a mapping or widen a buffer until all cycles appear to match. A qualified or inconclusive answer is a valid outcome.

## Inputs and reproducibility

- Public bundle: `https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/`. Discover river artifacts through its manifest and current HFX specification, starting at `docs/site/index.md`.
- Direction raster: `aux/d8/flow_dir.tif` under that URL; native EPSG:8857, 30 m cells, declared GRASS direction encoding. Preserve original directions.
- Existing local evidence: `~/Desktop/hfx-269-investigation/` and `~/Desktop/pourpoint-223-investigation/`. These are local artifacts, not public attachments; some scripts retain old `/tmp` paths. Inspect rather than assume portability.
- Prior runnable experiment and report: `experiments/fabdem_baseline/`. Its GRIT provenance and cached-window comparison can help identify inputs; its derived terrain directions are not GRIT’s pre-enforcement raster.
- Standalone loop example in #269: global cells `(147888,234986)=5` and `(147889,234985)=1`, pointing back to each other. Original Zenodo source details and independent checks are linked there.

Use project-native tools (`uv`, typed Python, existing spatial dependencies) and a small reproducible experiment, not a new validation framework. Follow repository testing guidance for new analysis code. Keep large inputs and derived spatial artifacts outside Git; record their locations, identities, commands and preparation steps. Keep scripts and a concise report in a dedicated repository experiment directory. Report inaccessible inputs or evidence gaps explicitly rather than claiming a completed comparison.

## Evidence to return

Provide a per-cycle inventory with native cell identities, overlap classification and distance evidence, aggregate counts with explicit denominators, and representative cell-level figures. Record input identities, geometry/raster conventions, commands, and limitations so another agent can reproduce the result.

The report must answer separately:

1. How many known cycle clusters intersect the published rivers, and how close are the others?
2. What can be established about membership in the original stream cells that Michel modified?
3. Does available evidence establish direction enforcement as the cause, merely support it, challenge it, or leave it unresolved? What additional artifact would distinguish these possibilities?

Post a concise summary and durable report links to #277 and #269. Clearly distinguish observed results, Michel’s reported workflow, and remaining hypotheses. Leave the parent discussion open for the maintainer’s decision and contact with Michel.

## Exclusions

No direction correction, hole filling, production changes, global cycle survey, replacement raster adoption, Zenodo publication, or upstream contact. Do not implement Michel’s proposed fix. The prospective University of Zurich VM and possible future corrected release are outside this experiment. Eliminating cycles would not by itself prove vector–raster watershed consistency; outlet placement and straightened boundaries remain separate concerns.
