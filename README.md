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

## Results (quick study, v0.1.0)

52 scenario units x 11 methods x 4 look schedules, 4,000 A/A and 2,000 A/B replicates per
unit. Full tables and figures: [docs/results/quick/REPORT.md](docs/results/quick/REPORT.md);
raw per-cell numbers: [docs/results/quick/summary.json](docs/results/quick/summary.json).
Rates are FPR with 95% intervals; "violates" means FPR > 5% after Holm adjustment.

- **gbstats sequential test, default settings (relative effect):** holds 5% on normal data,
  10% conversions, the ratio metric and Pareto-tailed revenue. It violates in 13 of 96
  claimed cells, all lognormal revenue with a 20% treatment arm or 1% conversions, for
  example 11.2% [10.2, 12.2] (lognormal, N = 20,000, 20% to B, 140 looks) and 6.5%
  [5.7, 7.3] (1% conversions, N = 20,000, 50/50, daily looks). The same test on the
  absolute scale has no violations (0 of 96).
- **gbstats with N* tuned to the planned size:** worse on small experiments, 20 of 96 cells,
  up to 36.5% [35.0, 38.0].
- **spotify-confidence group sequential test:** holds 5% on normal data with balanced
  traffic. It violates in 17 of 60 claimed cells, all with 20% of traffic to B, for example
  20.4% [19.2, 21.7] on 10% conversions at N = 1,000 with 140 looks, against 5.5% for the
  same z-test at one look. Run past the planned sample size it reaches 8.0% on normal data
  (informational: the library makes no claim there).
- **Peeking at fixed-horizon tests** inflates FPR to about 22% with 14 daily looks and 40%
  with 140 looks, matching the exact repeated-significance value. The gbstats Bayesian
  "chance to win >= 95%" rule is exactly a two-sided 10% z-test and reaches 37% when peeked
  daily.

All 40 sanity gates pass (positive controls against exact values, power calibration,
identities, engine agreement). The study's detection limit is about 6.2% FPR per cell, so
"consistent" is not proof of calibration. See the report's limitations, in particular that
gbstats 0.8.0 on PyPI is older than the code GrowthBook runs today.

## Install

You need [uv](https://docs.astral.sh/uv/). peekaudit is distributed from git only.

```sh
# as a tool, from git
uv tool install "peekaudit[engines] @ git+https://github.com/h-a-forster/peekaudit@v0.1.0"

# or for development (canonical study environment: Python 3.11)
git clone https://github.com/h-a-forster/peekaudit
cd peekaudit
uv sync --locked
uv run peekaudit versions
```

Python 3.11 satisfies gbstats 0.8.0's own pins (numpy<2, pandas<2). On Python 3.12+ the uv
project lifts those pins (gbstats' test suite passes on numpy 2 / pandas 2); a plain pip
install with the `engines` extra on 3.12+ will not resolve without the same overrides.

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
uv run python scripts/run_study.py --config quick --out runs/quick   # ~35-50 min, 4 workers
uv run peekaudit status runs/quick
uv run peekaudit report runs/quick --publish docs/results/quick
```

The run checkpoints one JSON file per unit, logs to `runs/quick/progress.log`, and resumes
where it stopped; the wrapper retries up to three times. `--config smoke` finishes in a few
minutes; `--config full` (20,000 / 10,000 replicates) takes a few hours.

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

Two standalone scripts call the engines directly, without peekaudit code:
[scripts/zero_count_example.py](scripts/zero_count_example.py) shows why z-tests over-reject
with rare conversions and a small arm, and
[scripts/spotify_indexerror_repro.py](scripts/spotify_indexerror_repro.py) is a minimal
reproduction of a spotify-confidence exception on a single-unit control group (outside the
study grid).

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
