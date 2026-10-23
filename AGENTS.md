# AGENTS.md

Guide for agents and contributors continuing this project. The authoritative spec is
[docs/DESIGN.md](docs/DESIGN.md). Read it before changing any formula, scenario or reporting
rule. If you deviate from it, update DESIGN.md in the same change and record the deviation.

## What this is

A simulation audit of the sequential-testing guarantees of open-source A/B engines
(gbstats 0.8.0, spotify-confidence 4.1.0) against reference methods (naive z, mSPRT,
asymptotic CS, Lan-DeMets OBF GST). It reports false-positive rate, coverage at stopping,
power and expected N, with exact CP intervals and Holm control. It also ships a reusable
conformance check.

## Conventions

- uv + hatchling, `src/peekaudit`, Python 3.11+ (3.11 is the canonical study environment:
  `.python-version`). Typed code, ruff (line length 100), pytest.
- Core runtime deps: numpy, scipy only. Engines are the optional `engines` extra and are
  included in the dev group. Engines are imported, never vendored; copied test vectors carry
  an attribution comment.
- `[tool.uv] override-dependencies` lifts gbstats' numpy<2 / pandas<2 caps on Python 3.12+
  only. `bokeh<3.10` is required (chartify 5.0.1 breaks on bokeh 3.10). Do not remove either
  without re-running the engine tests on 3.11 and 3.13.
- No network or LLM/model API calls at run time. CPU only. At most 4 worker processes in code;
  use `--workers 3` on the shared development machine.
- Long computations checkpoint each finished unit to disk, skip completed units on rerun and
  log progress to `runs/<study>/progress.log` (DESIGN section 7).
- Scratch files go in `tmp/` (gitignored). Never commit `runs/`, `cache/` or `data/`; only
  small derived results under `docs/results/`.
- Statistical honesty: report what was run (R, seeds, versions, wall time); intervals, not
  just point estimates; Holm within the families of DESIGN section 6.2; state the Monte
  Carlo detection limit. A suspected engine defect needs a minimal standalone repro and an
  independent check that it is not an adapter mistake (DESIGN 5.2, item 4).
- Commits: end messages with the project's co-author trailer. Never push, create remotes or
  publish without the owner's go-ahead.

## Architecture map

| Path | Role | State |
|---|---|---|
| `src/peekaudit/stats.py` | Clopper-Pearson, exact binomial test, Holm, mean CI | done |
| `src/peekaudit/gst.py` | Lan-DeMets solver (recursive integration), spending functions | done |
| `src/peekaudit/cli.py` | `run`, `status`, `sensitivity`, `report`, `compare`, `verify`, `conform`, `versions` | done |
| `src/peekaudit/scenarios.py` | metrics, true moments, effect calibration, schedules, grids | done |
| `src/peekaudit/simulate.py` | seeded, chunked generation of cumulative sufficient statistics | done |
| `src/peekaudit/methods/base.py` | `MethodOutput`, `EvalContext`, `Method`, `ScalarMethod`, stopping rule | done |
| `src/peekaudit/methods/reference.py` | naive z, mSPRT, ACS, GST-OBF | done |
| `src/peekaudit/methods/gbstats_adapter.py` | gbstats engine calls + vectorised mirror | done |
| `src/peekaudit/methods/spotify_adapter.py` | spotify-confidence engine calls + mirror | done |
| `src/peekaudit/methods/registry.py` | method names -> factories | done |
| `src/peekaudit/verify.py` | mirror-vs-engine randomised grid and in-loop agreement | done |
| `src/peekaudit/runner.py` | units, config validation, process pool, checkpoints, provenance, progress log | done |
| `src/peekaudit/sensitivity.py` | FPR counting only looks with both arms >= 150 (checkpointed) | done |
| `src/peekaudit/compare.py` | cross-study cell-by-cell check (nesting-aware Fisher + Holm) | done |
| `src/peekaudit/report.py` | tables, Holm families, gates, detection limit, findings, limitations, figures | done |
| `src/peekaudit/conformance.py` | batteries, `check_conformance`, `assert_conforms` | done |
| `scripts/reference_values.py` | recomputes the pinned oracle values and MC power table | done |
| `tests/` | unit tests; `@pytest.mark.slow` Monte Carlo tests run with `--runslow` | ongoing |
| `docs/DESIGN.md` | spec: formulas, validation, study design, milestones | done |

## How to test

```sh
uv sync
uv run pytest                # fast suite (~90 s; includes real-engine agreement tests)
uv run pytest --runslow      # adds Monte Carlo checks
uv run ruff check .
uv run ruff format --check .
```

CI (`.github/workflows/ci.yml`) runs on `pull_request` and `workflow_dispatch` only:
ubuntu py3.11, ubuntu py3.13, windows py3.13, plus a lint job.

## How to run the study

```sh
uv run peekaudit run --config quick --out runs/quick      # target < 45 min, 4 workers
uv run peekaudit status runs/quick
uv run peekaudit sensitivity runs/quick                   # ~2 min, adds the guarded FPR
uv run peekaudit report runs/quick --publish docs/results/quick
uv run python scripts/run_study.py --config full --out runs/full --workers 3   # R_AA = 20,000; 53 min on 4 workers
uv run peekaudit sensitivity runs/full --workers 3                            # ~6 min
uv run peekaudit report runs/full --publish docs/results/full
uv run peekaudit compare docs/results/quick docs/results/full --out docs/results/full/quick_vs_full
```

Before a new long run, time a one-unit pilot and project the total. If it would overrun the
estimate by more than 50%, stop and report. Keep sessions under ~45 min: split with
`--units GLOB` into the same output directory. New experiments go to a new `runs/<name>`
and `docs/results/<name>`; never write into `runs/full` or `docs/results/full`. Wrap long runs with `scripts/run_study.py`
(retry up to 3 times; checkpoints make reruns idempotent) and watch `progress.log`.
Completeness is one rule everywhere (`runner.unit_complete`: same config hash). The config
hash does not cover the code: after a code change that affects results, rerun with
`--force` (or a fresh `--out`); `run`/`status` warn about units from other code and the
report lists every git SHA found.

## Key facts already established (with evidence)

- spotify-confidence bounds equal the independent solver for power-family rho = 2 (K = 5, 10)
  to 2e-3 (`tests/test_engines_env.py`).
- The spotify sequential call on the full path equals look-by-look calls while totals stay
  <= the planned size. Past it, the engine rescales earlier fractions, so results differ;
  overrun must be verified look by look.
- spotify `ZTest.difference` costs ~0.4 s per call at 14 looks and ~8.5 s at 140 looks.
  gbstats sequential costs ~0.18 ms per call, the t-test ~0.8 ms, Bayes ~5.7 ms. Hence the
  vectorised mirrors plus in-loop verification on a replicate subsample.
- gbstats Bayes with the default improper prior and a 95% chance-to-win threshold is exactly
  a two-sided z-test at level 0.10 on the relative scale.
- gbstats 0.8.0 on PyPI (2024-06-20) is not the code GrowthBook main runs today, though both
  say 0.8.0. The audit covers PyPI only.
- GrowthBook product defaults (`packages/shared/src/constants.ts`, commit de93c5d6e1):
  sequential off, tuning 5000, confidence level 0.95, min sample size 150 (product layer).

## Status and next steps

- **Done (M0):** scaffold, `stats.py`, `gst.py`, engine environment tests, CLI `versions`,
  `scripts/reference_values.py`, DESIGN.md.
- **Done (M1):** `scenarios.py`, `simulate.py`: moment, Hill, KS-shortcut, calibration and
  chunk-determinism tests. Arm counts use round-half-up.
- **Done (M2):** `methods/base.py` (Method, MethodOutput with tie margin, EvalContext,
  ScalarMethod, stopping rule), `methods/reference.py`, `methods/registry.py`.
- **Done (M3, M4):** gbstats and spotify adapters. Mirrors agree with the engines exactly
  (bitwise-identical intervals in practice) on random grids incl. degenerate branches, on
  study data and on random GST paths with overrun. spotify's alpha is 1 - 0.95 =
  0.050000000000000044; the mirror uses it. The overrun rescaling (0.0198 full path vs
  0.0297 look by look) is a regression test. Only the scaffold's gbstats fixture is
  reproduced; the other sdist/repo fixtures are not yet copied (requires a download).
- **Done (M5):** `runner.py`, `verify.py`, CLI `run`/`status`/`verify`,
  `scripts/run_study.py`. Pilot on an idle core: `lognorm-n20000-b20-aa` 170 s (140 s of it
  real-engine verification, mostly spotify). Projection: quick ~32 min wall on 4 workers.
- **Done (M7):** `conformance.py`, `peekaudit conform`, smoke battery in CI.
- **Dropped:** `gbstats.seq.rel.guard150` (DESIGN 4.3): the compared quantity could not be
  confirmed from source.
- **Done (M6):** quick study at git `fa12471`, config hash `82e303ecf6547d58`: 52/52 units,
  0 failed, 0 unverified, 50.3 min wall / 200.6 CPU-min on 4 workers (machine loaded by
  other jobs). All 40 gates pass. Published to `docs/results/quick/`. A single-process
  rerun of one unit is byte-identical. `docs/results/verify.json`: 0 mismatches.
- **Done (M8):** README, CHANGELOG 0.1.0, version 0.1.0 (no tag). The spotify mirror's
  single-unit-control `IndexError` branch was added after the run; it cannot occur in the
  study grid (control arms start at 3 units).
- **Done (audit fixes, Unreleased in CHANGELOG):** one completeness rule (config hash) for
  run/status/report; per-unit provenance in the timing sidecar and a manifest start history
  (`runs/quick` predates this: its 52 units are attributed to the manifest's single start,
  `fa12471`, which the progress log confirms); config validation; `conformance.critical_rate`
  separates the FPR null (alpha) from the Holm FWER; wrapper log scrubbed of local paths;
  identity gates relabelled as mirror self-consistency; DESIGN 5.4 at m = 96.
  `peekaudit sensitivity` run on `runs/quick` (28 A/A units, 1.6 min, every unguarded count
  equals the run's). The generated report now has the detection limit, findings with
  Holm-adjusted p, the min-arm sensitivity, limitations and deviations; republished to
  `docs/results/quick/`. Study numbers are unchanged (no re-run needed: only reporting and
  the added sensitivity changed). README results rewritten around the small-arm mechanism.
- **Done (full study):** `runs/full` at git `e4e8dcc`, config hash `dde79da82e86a243`:
  52/52 units, 0 failed, 0 unverified, 53.2 min wall / 211.8 CPU-min on 4 workers; 40/40
  gates; 1,126,400 in-loop replicate-looks with 0 engine mismatches. Sensitivity 5.7 min on 3
  workers (every unguarded count equals the run's). Published to `docs/results/full/`
  (progress log scrubbed). Root `REPORT.md` and README lead with it.
- **Done (cross-check):** `peekaudit compare` (`compare.py`). Quick replicates are the first
  replicates of full (same seed and chunks; caches verified identical), so each quick rate is
  tested against full's remaining replicates. 0 of 2,224 cells differ after Holm (FPR/power
  and coverage); all 132 quick violations confirmed; 43 consistent -> violation, all within
  the quick 95% interval (<= 6.7% FPR). `docs/results/full/quick_vs_full.{md,json}`.
- **Not done / next:** upstream issue drafts (`docs/upstream/`); copying the gbstats sdist
  and spotify repo fixtures (needs a download); a guarded gbstats engine variant once the
  product's minimum-sample rule is confirmed from source (the sensitivity analysis is the
  interim answer); an engine-backed identity check (the two identity gates share arithmetic
  with the mirrors); variance reduction (CUPED) stretch; tagging and publishing (README
  install points at the repository without a tag until then).
