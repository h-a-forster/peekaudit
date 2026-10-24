# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This
project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Full study results in `docs/results/full/` (20,000 A/A and 10,000 A/B replicates per
  unit, 52 units, all gates pass), with the min-arm sensitivity and a scrubbed progress log.
- Root `REPORT.md`: summary of the full study with the quick study as a cross-check.
- `peekaudit compare` (`peekaudit.compare`): cell-by-cell cross-check of two published
  studies with Fisher exact tests and Holm adjustment. It detects nested studies (same seed
  and chunk size) and then compares with the larger study's remaining replicates.
- `scripts/run_study.py --units GLOB` to split a long study into shorter sessions.
- Engine-backed identity checks (`peekaudit.verify.engine_identities`, part of
  `peekaudit verify`): real gbstats engines against independent formulas and `ref.acs`;
  `docs/results/verify.json` regenerated with them.
- `docs/upstream/`: draft issue texts (not filed) for five engine behaviours, with new
  standalone reproductions `scripts/spotify_overrun_repro.py`,
  `scripts/gbstats_rare_conversion_repro.py` and `scripts/spotify_single_unit_repro.py`
  (spotify-confidence reports significance for a one-unit binary group).
- Like-for-like control for spotify's GST: `ref.gst_pow2` (power-family rho = 2 spending,
  reference variances; outside the study grid) and `peekaudit.gst_reference` /
  `scripts/gst_reference_comparison.py`, which compare it with `spotify.gst` and
  `ref.gst_obf` on a run's own A/A data (checkpointed per unit, paired McNemar tests).
  Results in `docs/results/full/gst_reference.{md,json}`.

- `peekaudit sensitivity` (`peekaudit.sensitivity`): recomputes each sequential method's
  A/A false-positive rate counting only rejections at looks where both arms hold at least
  150 units. Checkpointed per unit; the regenerated unguarded counts must equal the run's.
- The generated report now carries the Monte Carlo detection limit per method and family,
  a findings list with Holm-adjusted p-values (sequential, single-look and coverage findings
  separated, with the share of false positives at the first look and at small-arm looks),
  the min-arm sensitivity, a limitations section and the design deviations.
- Per-unit provenance (git SHA, versions) in `units/<unit>.timing.json`; the manifest keeps
  a history of starts; the report lists every distinct code version it finds.
- `StudyConfig` validation (counts, methods, schedules, unit ids, unknown or missing keys),
  run before anything is written to the output folder.

### Changed

- README results now lead with the full study; the quick study is the cross-check.

- One completeness rule for `run`, `status` and `report`: same config hash. The package
  version is recorded, not matched, so `status` no longer reports a finished run as 0 done
  after a version bump. `run` and `status` warn when finished units come from other code.
- The two identity gates are reported as mirror self-consistency checks (they compare code
  paths that share arithmetic).
- README results: states the small-arm mechanism, separates single-look from sequential
  failures, adds the stress-family results, and corrects the overrun maximum (8.2%) and the
  peeking figures (40% at N = 20,000; 46% at N = 1,000).
- REPORT and README: the O'Brien-Fleming reference GST is no longer presented as a
  like-for-like control for spotify's GST. The rho = 2 reference violates in the same cells,
  but about half of the largest spotify excess (10% conversions, 140 looks) comes from
  one-unit binary arms that the reference refuses; this is now drafted upstream. Peeking
  figures that match the exact values are restricted to N_max = 20,000, and "the engines
  report significance on 1 to 30 users" now says spotify (gbstats refuses n <= 1 arms).
- README install: `uv tool install` needs `--python 3.11` (it ignores the project's
  dependency overrides, so 3.12+ tries to build pandas 1.5.3).
- DESIGN 5.4 and `scripts/reference_values.py` use m = 96 primary cells (critical FPR 6.20%
  at R = 4,000).

### Fixed

- `critical_rate` / `ConformanceReport.detection_limit` used the method's alpha as the
  family-wise level; the Holm level is now a separate `fwer` (default 0.05), also used for
  the violation rule.
- In-loop engine verification now fails a unit on no-result or interval mismatches as well
  as decision mismatches (coverage uses the mirror's intervals); the report re-derives the
  rule from per-cell counts for older units. Published runs: 0 mismatches of every kind.
- The git SHA in manifests, timing sidecars and status was taken from whatever repository
  enclosed `<install>/../..`; it is now read only from a peekaudit checkout that is the top
  level of its repository, otherwise `unknown`.
- `peekaudit compare` treated studies as nested when the smaller replicate counts were not
  whole chunks; a partial chunk is not a prefix of a full one (arm B's stream is offset).
  The published quick/full comparison (whole chunks) is unaffected.
- `StudyConfig` validation accepted unit ids with impossible designs (N_max 0, allocation 0%
  or above 100%, arms that do not grow at every look, A/B units for infinite-variance
  metrics, non-canonical ids).
- DESIGN 4.2 said single-unit arms give `is_significant` False; for binary data they can
  give SE 0 and a significant result (the mirror already matched the engine).
- `scripts/run_study.py` no longer writes the interpreter's absolute path to the progress
  log; the published `docs/results/quick/progress.log` is scrubbed.

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
