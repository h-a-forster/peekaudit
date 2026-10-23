# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This
project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

## 0.1.0 - 2026-10-23

First complete version: simulation engine, engine adapters, reference methods, checkpointed
study runner, report and conformance check, plus the results of the quick study.

### Added

- `peekaudit.scenarios`: seven metrics with closed-form moments (normal, Bernoulli 10% and
  1%, lognormal, Pareto tails 2.5 and 1.5, a ratio metric), lift calibration to 80% power,
  the 210-look grid and the four look schedules, the 52-unit study grid.
- `peekaudit.simulate`: seeded, chunked generation of cumulative sufficient statistics with
  exact per-interval shortcuts for normal and Bernoulli data; results do not depend on the
  worker count or chunk order.
- `peekaudit.methods`: a vectorised `Method` interface with a `ScalarMethod` wrapper for
  per-look callbacks; reference methods (naive z, mSPRT, asymptotic confidence sequence,
  Lan-DeMets O'Brien-Fleming GST); adapters for gbstats 0.8.0 (sequential relative,
  sequential relative with N* = N_max, sequential absolute, t-test, Bayesian chance to win)
  and spotify-confidence 4.1.0 (group sequential and fixed-horizon z-tests), each a
  vectorised mirror plus real engine calls.
- `peekaudit.verify` and `peekaudit verify`: mirror-versus-engine agreement on randomised
  inputs, including degenerate branches and group sequential paths that overrun.
- `peekaudit.runner` and `peekaudit run` / `peekaudit status`: per-unit checkpoints,
  resume, a process pool capped at 4 workers, a progress log, in-loop engine verification;
  `scripts/run_study.py` retries a run up to three times.
- `peekaudit.report` and `peekaudit report`: Clopper-Pearson intervals, one-sided exact
  tests with Holm families, Monte Carlo sanity gates, Markdown tables, figures and
  `summary.json`.
- `peekaudit.conformance` and `peekaudit conform`: `check_conformance` and
  `assert_conforms` for any method, with smoke, quick and standard batteries; the smoke
  battery runs in CI.
- `docs/results/quick/`: report, figures and `summary.json` of the quick study (52 units,
  4,000 / 2,000 replicates, all 40 sanity gates pass); `docs/results/verify.json`: the
  randomised mirror-versus-engine grid (0 mismatches).
- `scripts/zero_count_example.py` and `scripts/spotify_indexerror_repro.py`: standalone
  engine illustrations (z-test over-rejection with a zero-conversion arm; a
  spotify-confidence 4.1.0 `IndexError` for a single-unit control group with value 0).

### Changed

- README rewritten with install, quick start, validation steps, results and citation.
- Version 0.1.0 in `pyproject.toml`, `peekaudit.__version__` and `CITATION.cff`.

- Project scaffold: uv + hatchling, `src/` layout, CI (pull requests and manual dispatch),
  repository hygiene files.
- `docs/DESIGN.md`: architecture, exact method specifications, validation strategy, study
  design and milestones.
- `peekaudit.stats`: Clopper-Pearson intervals, one-sided exact binomial test, Holm
  adjustment, t-interval for means.
- `peekaudit.gst`: independent Lan-DeMets group sequential solver (recursive numerical
  integration) with O'Brien-Fleming-type, power-family and Pocock-type spending, validated
  against Armitage et al. (1969), Jennison and Turnbull (2000) and ldbounds/gsDesign values.
- Engine environment checks: pinned gbstats 0.8.0 and spotify-confidence 4.1.0, a gbstats
  fixture, the gbstats confidence-sequence/p-value duality, and spotify-confidence bounds
  against the reference solver.
- CLI: `peekaudit --version`, `peekaudit versions`.
