# GRASS Watershed Duplicate Seed and Invariant Investigation

Associated discussion: [Issue #273](https://github.com/pourpoint-tools/hfx/issues/273). Implementation of this vision will provide the complete evidentiary findings required to resolve the issue.

## Context and Purpose

During the investigation of published D8 drainage cycles in hydrofabric datasets ([#269](https://github.com/pourpoint-tools/hfx/issues/269)), source inspection of the GRASS GIS fork pinned by archived GRIT (`mwort/grass` commit [`6a32f679`](https://github.com/mwort/grass/commit/6a32f679c3b0defe092f69545bb2c7b4e51d2a84), parent commit [`e225508f`](https://github.com/mwort/grass/commit/e225508f5cedd15cc694fdf35b2850676ffa4ee1)) revealed two concrete anomalies in the segmented `r.watershed` engine (`-m`):

1. **Duplicate Watershed Seed Insertion (Upstream Defect):** In `raster/r.watershed/seg/init_vars.c`, valid cells adjacent to NULL or MASKed boundaries are assigned aspect (`af.asp`) and added to the seed min-heap via `add_pt()`, but the local variable `asp_value` remains initialized to `0`. If depression/sink input is provided, the subsequent check (`if (pit_flag && asp_value == 0)`) tests the stale local variable and inserts the same cell into the min-heap a second time. (In contrast, the RAM engine in `ram/init_vars.c` correctly checks `asp[seg_idx] == 0`). Upstream `main` (`grass-upstream` commit `aca52e2`) retains this exact discrepancy.
2. **Broken A* Count Invariant and Out-of-Bounds Memory Write:** Duplicate seeding causes the min-heap search in `seg/do_astar.c` to pop an extra point. The countdown variable `doer` decrements past 0 to `-1`. Crucially, before the count invariant check, `Segment_put(&astar_pts, (char *)&heap_p.pnt, 0, doer)` writes out-of-bounds at index `-1` into memory directly preceding `SEG->cache`.
3. **Fork Fail-Open Modification (`dc9ef282`):** Upstream GRASS treats the count mismatch as fatal (`G_fatal_error`). The archived fork downgraded this to a warning (`G_warning`), permitting execution to continue past memory corruption into SFD flow accumulation, causing downstream allocator aborts (`munmap_chunk(): invalid pointer`). Furthermore, the 1-based `search_heap` array is allocated for `do_points + 1` columns; duplicate insertion can overflow the segment cache during initialization itself (`free(): invalid size`).

Initial ARM64 testing revealed an architecture confounder: plain `char` is unsigned on ARM64 (`__CHAR_UNSIGNED__ 1`), converting negative drainage codes to unsigned 248–254. Running under an isolated Docker `linux/amd64` container (`hfx273-grass-build-amd64`) eliminated this artifact, reproducing clean signed controls and distinguishing the early initialization abort on NULL-boundary fixtures from the post-warning allocator abort on MASK-boundary fixtures.

This vision specifies the implementation required to instrument these failures, verify the upstream root fix, test whether D8 cycles can ever be produced by bypassed execution, and prepare the local findings report.

## Settlement and Scope Boundaries

* **No Production Changes:** No changes to production GRIT tools, workflows, or published hydrofabric datasets.
* **No Premature Publication:** Findings, reports, and patch diffs remain strictly local drafts. No issues or pull requests may be submitted to OSGeo/GRASS or GRIT repositories without explicit human authorization.
* **Separation of Causation:** Verified software defects in GRASS must be kept strictly separate from the question of whether this condition occurred in production GRIT runs, or whether it caused the D8 cycles in [#269](https://github.com/pourpoint-tools/hfx/issues/269).

## Execution Plan

### Phase 1: AddressSanitizer (ASan) Runtime Instrumentation

1. **Build Configuration:**
   * Configure and build both the archived fork (`src/fork`) and upstream parent (`src/parent`) under Docker `linux/amd64` using `build.sh`.
   * Set `BUILD_CFLAGS="-fsanitize=address -fsanitize-recover=address -g -O1 -fno-omit-frame-pointer"` and `BUILD_LDFLAGS="-fsanitize=address"`. Ensure debug symbols are preserved (omit `-s`).
   * Isolate builds, installs, and logs in dedicated directories (e.g., `amd64-fork-asan/`, `amd64-parent-asan/`).
2. **ASan Execution & Trace Capture:**
   * Run the test suite (`cases.py`) against both instrumented builds with `ASAN_OPTIONS=detect_leaks=0:halt_on_error=0:symbolize=1`.
   * Capture and log exact stack traces with file and line numbers demonstrating:
     - The heap buffer overflow at `seg_put` into `search_heap` during duplicate insertion in `seg/init_vars.c`.
     - The out-of-bounds heap write at `seg_put(&astar_pts, ..., 0, -1)` in `seg/do_astar.c`.
     - The exact sequence differences between upstream parent (terminating at fatal error) versus the fork (continuing past warning into accumulation and allocator abort).

### Phase 2: Root Fix Verification & Upstream Bug Report Draft

1. **Root Fix Implementation:**
   * In a clean build tree of the segmented backend, patch `raster/r.watershed/seg/init_vars.c` to mirror `ram/init_vars.c` (e.g., updating `asp_value = af.asp;` when seeding NULL-adjacent cells, or testing `af.asp == 0` for depression seeding).
2. **Validation:**
   * Execute `cases.py` on the patched build across all 12 combinations of `flat3_null`, `flat3_mask`, and `ridge5_mask` (pit/no-pit, segmented/RAM).
   * Verify that duplicate seeding is eliminated, A* count matches expected points without warnings or errors, exit codes are 0, and the drainage/accumulation oracle reports 0 diagnostics.
3. **Upstream Draft Preparation:**
   * Confirm the exact line diff against current upstream `main` (`grass-upstream` commit `aca52e2`).
   * Prepare a self-contained local draft bug report and patch diff suitable for OSGeo/GRASS submission once authorized.

### Phase 3: D8 Cycle Causation Investigation (The #269 Hypothesis)

1. **Bypassed-Execution Test Harness:**
   * Create a test fixture harness that isolates the fork's fail-open continuation by providing bounds padding or wrapping `astar_pts` and `search_heap` to prevent glibc allocator aborts.
2. **Parameterized Grid Evaluation:**
   * Construct a parameterized suite of synthetic elevation rasters:
     - Variable dimensions (3x3, 5x5, 10x10, 20x20).
     - Diverse topologies (uniform flats, directional slopes, interior ridges, multi-cell depressions).
     - Multiple sink/mask overlap placements.
3. **Cycle and Routing Analysis:**
   * Run the cycle-detection oracle from `cases.py` on all successfully generated drainage rasters.
   * Evaluate whether duplicate seed insertion or unvisited points in `astar_pts` ever induce closed D8 routing loops, or if `do_astar` routing remains strictly monotonic despite count invariant failure.

### Phase 4: Final Investigation Report

Compile all findings into a unified document (`~/Desktop/hfx-269-investigation/hfx273-final-report.md`) containing:
* Executive summary of findings for Issue #273.
* Instrumented ASan reproduction traces with line-by-line attribution.
* Root cause fix verification and upstream patch diff.
* Parameterized D8 cycle test results answering the #269 causation hypothesis.
* Prepared draft for updating GitHub Issue #273.

## Success Criteria

1. ASan execution logs cleanly identify the exact call stacks and buffer bounds for both `search_heap` and `astar_pts[-1]` corruptions.
2. The root fix in `seg/init_vars.c` achieves 100% test pass rate with 0 oracle diagnostics and parity with RAM mode across all cases.
3. The parameterized cycle evaluation provides conclusive empirical evidence whether broken-invariant continuation can generate D8 cycles.
4. Comprehensive report and upstream draft completed locally without unauthorized external publishing.
