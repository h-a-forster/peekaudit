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
- No network or LLM/model API calls at run time. CPU only. At most 4 worker processes.
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
| `src/peekaudit/cli.py` | `--version`, `versions`; later `run/status/report/verify/conform` | partial |
| `src/peekaudit/scenarios.py` | metrics, true moments, effect calibration, schedules, grids | done |
| `src/peekaudit/simulate.py` | seeded, chunked generation of cumulative sufficient statistics | done |
| `src/peekaudit/methods/base.py` | `LookStats`, `MethodOutput`, `Method` protocol, `ScalarMethod` | M2 |
| `src/peekaudit/methods/reference.py` | naive z, mSPRT, ACS, GST-OBF | M2 |
| `src/peekaudit/methods/gbstats_adapter.py` | gbstats engine calls + vectorised mirror | done |
| `src/peekaudit/methods/spotify_adapter.py` | spotify-confidence engine calls + mirror | done |
| `src/peekaudit/methods/registry.py` | method names -> factories | M2 |
| `src/peekaudit/verify.py` | mirror-vs-engine randomised grid and in-loop agreement | done |
| `src/peekaudit/runner.py` | units, process pool, checkpoints, progress log | done |
| `src/peekaudit/report.py` | tables, Holm families, figures | M6 |
| `src/peekaudit/conformance.py` | batteries, `check_conformance`, `assert_conforms` | M7 |
| `scripts/reference_values.py` | recomputes the pinned oracle values and MC power table | done |
| `tests/` | unit tests; `@pytest.mark.slow` Monte Carlo tests run with `--runslow` | ongoing |
| `docs/DESIGN.md` | spec: formulas, validation, study design, milestones | done |

## How to test

```sh
uv sync
uv run pytest                # fast suite (~40 s; includes engine environment checks)
uv run pytest --runslow      # adds Monte Carlo checks
uv run ruff check .
uv run ruff format --check .
```

CI (`.github/workflows/ci.yml`) runs on `pull_request` and `workflow_dispatch` only:
ubuntu py3.11, ubuntu py3.13, windows py3.13, plus a lint job.

## How to run the study (once M5/M6 land)

```sh
uv run peekaudit run --config quick --out runs/quick      # target < 45 min, 4 workers
uv run peekaudit status runs/quick
uv run peekaudit report runs/quick --publish docs/results/quick
uv run peekaudit run --config full --out runs/full        # R_AA = 20,000; est. 3-5 h
```

Before the full run, time a one-unit pilot and project the total. If it would overrun the
estimate by more than 50%, stop and report. Wrap long runs with `scripts/run_study.py`
(retry up to 3 times; checkpoints make reruns idempotent) and watch `progress.log`.

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
  `scripts/reference_values.py`, DESIGN.md. `uv run pytest`, ruff check and format pass on
  Windows, Python 3.11.
- **Done (M1):** `scenarios.py`, `simulate.py` with moment, Hill, KS-shortcut, calibration
  and chunk-determinism tests. Arm counts use round-half-up (not banker's rounding).
- **Next (M2):** `methods/base.py`, `methods/reference.py`, `methods/registry.py`.
- **Then:** M3 gbstats adapter, M4 spotify adapter, M5 runner + pilot
  timing, M6 quick study + report, M7 conformance, M8 README findings and release notes.
- **Open questions:** the exact quantity behind GrowthBook's min-sample-size guard (DESIGN
  4.3) must be confirmed from source before `gbstats.seq.rel.guard150` is added; drop the
  variant if it cannot be confirmed.
