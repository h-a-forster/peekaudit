# peekaudit 0.1.1 (2026-10-30)

Adds the full study, a cross-check against the quick study, a like-for-like control for
spotify-confidence's group sequential test and upstream issue drafts. The 0.1.0 study
numbers are unchanged.

## Results

- Full study in `docs/results/full/`: 52 units, 20,000 A/A and 10,000 A/B replicates each,
  2,224 cells, 0 failed units, 40/40 sanity gates, 0 engine mismatches on 1,126,400
  replicate-looks. Summary in the root [REPORT.md](../REPORT.md).
- `peekaudit compare`: the quick replicates are a prefix of the full study's; 0 of 2,224
  cells differ after Holm, and all 132 quick violations are confirmed
  (`docs/results/full/quick_vs_full.md`).
- `peekaudit sensitivity`: false-positive rate counting only looks where both arms hold at
  least 150 units, which separates small-arm excesses from rare-event and skew excesses.
- `ref.gst_pow2` and `scripts/gst_reference_comparison.py`: a reference GST with spotify's
  spending function violates in the same 21 of 60 primary cells; the rest of spotify's
  largest excess comes from one-unit binary arms (`docs/results/full/gst_reference.md`).

## Tooling

- Engine-backed identity checks in `peekaudit verify` (0 mismatches on 2,000 inputs each).
- Per-unit provenance (git SHA, versions), a single completeness rule (config hash) and
  stricter `StudyConfig` validation.
- `scripts/run_study.py --units GLOB` to split long studies.

## Fixes

- In-loop verification also fails a unit on no-result or interval mismatches.
- The conformance detection limit and violation rule use a separate FWER instead of the
  method's alpha.
- `peekaudit compare` no longer treats studies as nested when replicate counts are not whole
  chunks.
- The git SHA is read only from a peekaudit checkout at its repository top level.

## Upstream drafts

Five draft issues in [docs/upstream/](upstream/README.md), none filed, each with a standalone
reproduction. All were re-checked on 2026-10-30 against gbstats 0.8.0 and
spotify-confidence 4.1.0 (still the latest PyPI releases) and still apply; the gbstats
rare-event result is identical on GrowthBook main.

Full list of changes up to 2026-10-30: [CHANGELOG.md](../CHANGELOG.md).

## Changes after the 2026-10-30 notes

These are in `polish-2026-10-30` and are not listed in the changelog.

- **Correction:** the reports said GrowthBook's product hides results below 150 units. Its
  source compares 150 with the larger arm's metric total, so `guarded` (both arms at least
  150 users) is not that rule. README, REPORT, DESIGN 4.3 and section 13, the gbstats draft
  and the generated reports now say so. No study number changed.
- `peekaudit run --units GLOB` now refuses a pattern that matches no unit, before writing
  anything. It used to exit 0 having run nothing.
- `check_conformance` and `assert_conforms` now refuse a method that supports none of the
  battery's metrics. It used to pass with 0 cells.
- `peekaudit compare` now refuses studies that hold the same replicates. It used to report 0
  differences after comparing 0 cells.
- A report no longer uses sensitivity checkpoints left by another config after `--force`.
- The report states the study's own unit total instead of a fixed 52.
- CLI usage errors (missing run folder or config, unknown method) exit 2 with a message
  instead of a traceback.
- CI actions are pinned to commit SHAs.
- DESIGN sections 2, 3 and 13 updated to match the code and the full run.
