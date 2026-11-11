# GRASS r.watershed duplicate seed: upstream report

Associated discussion: [#273](https://github.com/pourpoint-tools/hfx/issues/273). This vision does not complete #273; it remains open for the maintainer.

## Outcome

Prepare a ready-to-submit, locally verified report of a memory-safety bug in upstream GRASS `r.watershed` (segmented mode), so the maintainer can submit it to the GRASS project. The report points at the exact line causing the bug, shows it reproduces on current GRASS, and shows that a one-line change makes the reproduction pass locally.

We found the bug as a side problem while investigating GRIT D8 cycles ([#269](https://github.com/pourpoint-tools/hfx/issues/269)). The report is about GRASS only. It must not claim the bug caused GRIT's loops: a bounds-contained experiment observed no cycles (`~/Desktop/hfx-269-investigation/hfx273-native/cycle-findings.md`).

## Root cause (source reading of upstream `main` at `32785cd092089a58c9ab0fce28aeaa85066a4ee1`)

In [`raster/r.watershed/seg/init_vars.c`](https://github.com/OSGeo/grass/blob/32785cd092089a58c9ab0fce28aeaa85066a4ee1/raster/r.watershed/seg/init_vars.c#L527-L582):

- Line 527 copies the cell direction into local `asp_value` (0 for an interior cell).
- The edge branch updates both `asp_value` and `af.asp`.
- The interior NULL-neighbour branch sets `af.asp`, calls `add_pt`, but never updates `asp_value`.
- Line 578, `if (pit_flag && asp_value == 0)`, reads the stale 0, so a depression cell beside a NULL or masked cell is queued a second time.
- RAM mode ([`ram/init_vars.c:444`](https://github.com/OSGeo/grass/blob/32785cd092089a58c9ab0fce28aeaa85066a4ee1/raster/r.watershed/ram/init_vars.c#L444)) checks the stored `asp[seg_idx]` and does not duplicate.

Consequences seen on 8.4.0dev builds: the search heap is sized for valid cells plus one, so the extra insert overflows it; the extra A* extraction writes `astar_pts[-1]`; the run ends with `ERROR: bug in A* Search: doer -2 heap size 0 count 9` or memory aborts. Trigger: `-m`, a `depression=` raster, and a depression cell adjacent to a NULL or masked cell.

This is source reading plus earlier runs on `e225508f` and the `mwort/grass` fork. Current `main` and 8.5.0 have not been executed yet.

## Work

1. Reproduce on current upstream `main`, and on 8.5.0 if practical, reusing the Docker setup and `cases.py` in `~/Desktop/hfx-269-investigation/hfx273-native/`. Record `g.version -rge`, the error output, and one AddressSanitizer stack.
2. Apply the fix idea locally (check `af.asp == 0` like RAM mode, or set `asp_value` in the NULL-neighbour branch) and confirm the reproduction passes: no count error, clean AddressSanitizer run, consistent outputs. If it does not pass, report that rather than filing the suggestion.
3. Write the report as a local draft for the maintainer to review and submit.

## Report decisions

- **Channel:** private GitHub Security Advisory on `OSGeo/grass`, per their [`SECURITY.md`](https://github.com/OSGeo/grass/blob/main/SECURITY.md) (private reporting is enabled). Maintainers may convert it to a public issue.
- **Structure:** follow GRASS's [bug template](https://github.com/OSGeo/grass/blob/main/.github/ISSUE_TEMPLATE/bug_report.md) sections: describe the bug, to reproduce, expected behavior, system description, additional context.
- **Voice:** super plain human wording, following the github-writing skill and GRASS's [AI use policy](https://github.com/OSGeo/grass/blob/main/CONTRIBUTING.md). No verbose AI-style text.
- **AI disclosure (full):** we were investigating a problem in GRIT, could not find the cause, and used AI to trace upstream code bugs; this bug was found along the way and is unrelated to that problem.
- **Fix:** no patch. Give the fix idea, state that it made the reproduction pass locally, and offer to open a PR if they are interested.
- **Evidence:** only what is needed to confirm and reproduce: a minimal shell reproduction, error output, one AddressSanitizer stack, the root cause with line links, and the local fix result. Leave out ARM64 plain-char behaviour, fork details and containment experiments.
- **No duplicate found:** a search of OSGeo/grass issues found no existing report.
- Michel Wortmann is not contacted about this.

## Boundaries

- Do not submit the advisory, open an issue, or open a PR upstream. The maintainer submits after review.
- No changes to HFX production code or data.
- Keep scratch checkouts under `.worktrees/`; build outputs are not evidence.
