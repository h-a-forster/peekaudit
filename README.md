# peekaudit

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](pyproject.toml)

Do open-source A/B testing engines keep their sequential-testing promises?

Experimentation platforms say you can check results as often as you like. peekaudit tests
that claim against the shipped engine code by simulation: thousands of A/A and A/B
experiments, peeked daily or near-continuously. The data include skewed and heavy-tailed
revenue, rare conversions, a ratio metric, unequal allocation, small samples, and runs past
the planned sample size. For every method, metric and look schedule it reports:

- the realised false-positive rate (nominal 5%),
- coverage of the reported interval at the stopping time,
- power and expected sample size at stopping,

with exact Clopper-Pearson intervals and Holm control across cells. It also ships a
conformance check that any sequential A/B method can run in its own test suite.

Engines audited (imported from PyPI, never vendored):

| Engine | Version | Licence | Methods |
|---|---|---|---|
| [gbstats](https://pypi.org/project/gbstats/) (GrowthBook) | 0.8.0 | MIT | asymptotic confidence sequences (relative and absolute); t-test; Bayesian chance to win |
| [spotify-confidence](https://pypi.org/project/spotify-confidence/) | 4.1.0 | Apache-2.0 | group sequential z-test; fixed-horizon z-test |

Reference implementations for comparison: naive peeked z-test, mSPRT (Johari et al.),
asymptotic confidence sequences (Waudby-Smith et al.) and a Lan-DeMets O'Brien-Fleming
group sequential test. The group sequential solver reproduces published tables (Armitage
et al. 1969; Jennison and Turnbull 2000).

Every engine result in the study comes from a vectorised mirror of the engine's arithmetic.
The mirrors are checked against real engine calls inside the run and on randomised grids
that include degenerate inputs; in v0.1 there were zero disagreements.

## Results (full study)

52 scenario units x 11 methods x 4 look schedules (2,224 cells), 20,000 A/A and 10,000 A/B
replicates per unit, 0 failed units, all 40 sanity gates passing. Summary and discussion:
[REPORT.md](REPORT.md). Every cell: [docs/results/full/REPORT.md](docs/results/full/REPORT.md)
and [docs/results/full/summary.json](docs/results/full/summary.json). Rates are FPR with
95% Clopper-Pearson intervals; "violates" means FPR > 5% after Holm adjustment. At this
size a cell can flag an FPR of about 5.5% or more.

**Where the excesses come from.** Many of the largest rates arise at the first looks of
small experiments, where one arm holds a handful of users (with N = 1,000 and 20% of traffic
to B, the first of 140 looks has 6 users in A and 1 in B). The engine libraries apply no
minimum sample size. The report therefore also gives each sequential cell's FPR counting
only rejections at looks where both arms hold at least 150 units ("guarded"). Rare
conversions are a separate problem: with 1% conversions the arms are large but the events
are few, and the guard changes nothing.

- **gbstats sequential test, default settings (relative effect):** holds 5% on normal data,
  10% conversions and the ratio metric. It violates in 14 of 96 claimed cells: lognormal
  revenue with a 20% treatment arm, 1% conversions, and one Pareto tail-2.5 cell (5.6%).
  Lognormal at N = 1,000 with 140 looks gives 16.7% [16.2, 17.2], 93% of it from looks with
  fewer than 150 users in an arm (guarded 3.7% [3.5, 4.0]); at N = 20,000 it is 11.9%
  [11.5, 12.4], guarded still 9.0% [8.6, 9.4]. 1% conversions give 6.3% to 14.1%, unchanged
  by the guard. Guarded, 10 of 96 cells still violate. Infinite-variance stress family
  (Pareto tail 1.5): 6 of 16, up to 14.3% [13.8, 14.8]. On the absolute scale: 0 of 96, and
  one stress cell, 5.7% [5.3, 6.0].
- **gbstats with N* tuned to the planned size:** worse on small experiments, 22 of 96 cells,
  up to 37.2% [36.6, 37.9] (lognormal, N = 1,000, 20% to B, 140 looks; guarded 5.7%
  [5.3, 6.0]). Guarded, 9 of 96 still violate, 7 of them with 1% conversions.
- **spotify-confidence group sequential test:** violates in 21 of 60 claimed cells, all but
  one with 20% of traffic to B. Three are single-look failures of the z-test at N = 1,000
  with 200 users in B (1% conversions 13.3% [12.9, 13.8], lognormal 7.1%, Pareto 6.1%). The
  largest sequential excesses are 20.5% [19.9, 21.1] (1% conversions, daily looks; guarded
  16.8%) and 20.4% [19.8, 20.9] (10% conversions, 140 looks; 44% of the false positives at
  the first look, guarded 5.1% [4.8, 5.4]). Guarded, 14 of 60 still violate. A reference
  GST with spotify's own spending function (power family, rho = 2) and sample variances
  violates in exactly the same 21 cells, so these excesses are how a z-based group
  sequential test behaves on small, skewed or rare-event arms. The exception is one-unit
  binary arms: spotify's variance p(1 - p) is 0 for a single user, and it reports
  significance there (SE 0, p = 0), where the reference gives no result. That accounts for
  about half of the 20.4% cell (reference: 11.5% [11.0, 11.9];
  [gst_reference.md](docs/results/full/gst_reference.md)). Past the planned sample size it reaches 8.2% [7.8, 8.5] on normal
  data (informational: the library makes no claim there).
- **Peeking at fixed-horizon tests** inflates FPR to 21-22% with 14 daily looks and about 40%
  with 140 looks at N = 20,000, matching the exact values (22.0% and 39.9%); at N = 1,000,
  23-24% (daily) and 46% (140 looks), above the exact values because the early looks are
  tiny. The gbstats Bayesian "chance to win >= 95%" rule is exactly a two-sided 10%
  z-test and reaches 37-39% when peeked daily.

**Cross-check with the quick study** (4,000 / 2,000 replicates, published in
[docs/results/quick/](docs/results/quick/REPORT.md)). The quick replicates are the first
replicates of the full study (same seeds), so `peekaudit compare` tests each quick rate
against the full study's other 16,000 replicates: 0 of 2,224 cells differ after Holm, for
FPR/power and for coverage. All 132 quick violations are confirmed. 43 cells that were
consistent in quick are flagged in full, each within the quick 95% interval and none above
6.7% FPR: small excesses the larger study can detect. See
[docs/results/full/quick_vs_full.md](docs/results/full/quick_vs_full.md).

Engine rates come from vectorised mirrors of the engines' arithmetic, checked against
the real engines on 1,126,400 replicate-looks inside the run and on randomised grids with degenerate inputs:
zero disagreements. The report lists every finding with its Holm-adjusted p-value, the
detection limit per method, the deviations from the design and the limitations, in
particular that gbstats 0.8.0 on PyPI is older than the code GrowthBook runs today and that
140 looks understate the error of truly continuous monitoring.

## Install

You need [uv](https://docs.astral.sh/uv/). peekaudit is distributed from git only.

```sh
# as a tool, from git (pin a release tag once one is published); --python 3.11 is required
uv tool install --python 3.11 "peekaudit[engines] @ git+https://github.com/h-a-forster/peekaudit"

# or for development (canonical study environment: Python 3.11)
git clone https://github.com/h-a-forster/peekaudit
cd peekaudit
uv sync --locked
uv run peekaudit versions
```

Python 3.11 satisfies gbstats 0.8.0's own pins (numpy<2, pandas<2). On Python 3.12+ the uv
project lifts those pins (gbstats' test suite passes on numpy 2 / pandas 2), but only inside
a checkout (`uv sync`, `uv run`): `uv tool install` and `pip install` ignore the project's
`[tool.uv] override-dependencies`. With the `engines` extra on 3.12+ they try to build
pandas 1.5.3 from source and fail, so pass `--python 3.11` (as above), or supply the same
overrides yourself (`uv tool install --overrides overrides.txt ...` with the four lines from
`pyproject.toml`).

## Quick start

Check a method against the conformance battery (about a minute):

```console
$ uv run peekaudit conform --method ref.gst_obf --battery smoke
ref.gst_obf on battery 'smoke': PASS (0 of 2 cells violate FPR <= 0.05; cells can flag FPR >= 6.50%)
  normal-n2000-b50-aa      daily14    FPR  4.70% [3.47, 6.20]  Holm p 1
  bern_p10-n2000-b50-aa    daily14    FPR  4.50% [3.30, 5.98]  Holm p 1
```

Exit code 0 means pass, 1 a violation, 2 a usage error. Batteries are `smoke`, `quick` and
`standard`. Your own method can be tested from Python:

```python
from peekaudit.conformance import BATTERIES, assert_conforms
from peekaudit.methods import ScalarMethod


def decide(n_a, sum_a, sumsq_a, n_b, sum_b, sumsq_b):
    """Called at every look; return (sign or 0, ci_lower, ci_upper), or None for no result."""
    ...


assert_conforms(ScalarMethod("my.method", decide), BATTERIES["smoke"])
```

Run and report the study:

```sh
uv run python scripts/run_study.py --config quick --out runs/quick --workers 3   # ~50 min
uv run peekaudit status runs/quick
uv run peekaudit sensitivity runs/quick --workers 3   # guarded FPR (both arms >= 150), ~2 min
uv run peekaudit report runs/quick --publish docs/results/quick
```

The run checkpoints one JSON file per unit, logs to `runs/quick/progress.log`, and resumes
where it stopped; the wrapper retries up to three times. `--config smoke` finishes in a few
minutes; `--config full` (20,000 / 10,000 replicates) took 53 minutes of wall time on 4
workers (212 CPU-minutes); `--units 'normal-*'` and similar globs split it into shorter
sessions writing to the same folder. `peekaudit compare SMALL LARGE --out STEM` cross-checks
two published studies cell by cell. A resume keeps finished units even if the code changed since; `run` and `status`
warn when that happens, the report lists every git SHA it finds, and `--force` recomputes.

## How to validate

```sh
uv run pytest                 # fast suite: closed forms, published constants, seeded MC
uv run pytest --runslow       # adds the slow Monte Carlo size checks
uv run peekaudit verify --n 2000 --n-gst 200 --json verify.json   # mirrors vs real engines
uv run ruff check . && uv run ruff format --check .
```

`peekaudit verify` compares every engine method's mirror with the engine on randomised inputs,
including empty and single-unit arms, zero control means and extreme variances, and exits
non-zero on any disagreement. The published result is
[docs/results/verify.json](docs/results/verify.json). The design, decision rules and every
deviation from it are documented in [docs/DESIGN.md](docs/DESIGN.md).

Standalone scripts call the engines directly, without peekaudit code:
[scripts/zero_count_example.py](scripts/zero_count_example.py) shows why z-tests over-reject
with rare conversions and a small arm;
[scripts/spotify_indexerror_repro.py](scripts/spotify_indexerror_repro.py) reproduces a
spotify-confidence exception on a single-unit control group (outside the study grid);
[scripts/spotify_overrun_repro.py](scripts/spotify_overrun_repro.py) shows spotify's
sequential bounds for earlier looks changing once a later look exceeds the planned size; and
[scripts/gbstats_rare_conversion_repro.py](scripts/gbstats_rare_conversion_repro.py)
simulates the gbstats sequential test on 1% conversions (relative vs absolute effect).
Draft upstream issue texts (not filed) are in [docs/upstream/](docs/upstream/README.md).

## Scope

- Results apply to the pinned engine versions only.
- Product-layer behaviour (for example GrowthBook's minimum sample size before results are
  shown) is not part of the engine libraries and is not audited.
- Calibration findings describe statistical behaviour under stated conditions. They are not
  claims of software defects.

## Citation

If you use peekaudit, please cite it using [CITATION.cff](CITATION.cff):

```bibtex
@software{peekaudit,
  author  = {{peekaudit contributors}},
  title   = {peekaudit: auditing sequential-testing guarantees of open-source A/B engines},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/h-a-forster/peekaudit}
}
```

## Licence

MIT. See [LICENSE](LICENSE). The audited engines keep their own licences.
