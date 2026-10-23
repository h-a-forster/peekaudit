# peekaudit

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](pyproject.toml)

Do open-source A/B testing engines keep their sequential-testing promises?

Experimentation platforms say you can check results as often as you like. peekaudit tests
that claim against the shipped engine code by simulation: thousands of A/A and A/B
experiments, peeked daily or near-continuously. The data include skewed and heavy-tailed
revenue, ratio metrics, unequal allocation, small samples, and runs past the planned sample
size. It reports:

- the realised false-positive rate (nominal 5%),
- coverage of the reported interval at the stopping time,
- power and expected sample size at stopping,

all with exact Clopper-Pearson intervals and Holm control across cells. It also ships a
conformance check that any sequential A/B engine can run in its own test suite.

Engines audited (imported from PyPI, never vendored):

| Engine | Version | Licence | Methods |
|---|---|---|---|
| [gbstats](https://pypi.org/project/gbstats/) (GrowthBook) | 0.8.0 | MIT | asymptotic confidence sequences; t-test; Bayesian chance to win |
| [spotify-confidence](https://pypi.org/project/spotify-confidence/) | 4.1.0 | Apache-2.0 | group sequential z-test; fixed-horizon z-test |

Reference implementations for comparison: naive peeked z-test, mSPRT (Johari et al.),
asymptotic confidence sequences (Waudby-Smith et al.) and a Lan-DeMets O'Brien-Fleming
group sequential test. The group sequential solver reproduces published tables (Armitage
et al. 1969; Jennison and Turnbull 2000; ldbounds/gsDesign).

**Status: pre-release (v0.0.1).** The design is fixed and the statistical foundations are
implemented and tested; the simulation, adapters and study are in progress. There are no
results yet. See [docs/DESIGN.md](docs/DESIGN.md) for the full specification and milestones.

## Install (development)

You need [uv](https://docs.astral.sh/uv/).

```sh
git clone https://github.com/h-a-forster/peekaudit
cd peekaudit
uv sync
uv run peekaudit versions
```

The canonical study environment is Python 3.11, which satisfies gbstats 0.8.0's own pins
(numpy<2, pandas<2). On Python 3.12+ the uv project lifts those pins (gbstats' test suite
passes on numpy 2 / pandas 2). A plain `pip install peekaudit[engines]` on 3.12+ will not
resolve without the same overrides.

## Tests

```sh
uv run pytest              # fast suite
uv run pytest --runslow    # adds Monte Carlo checks
uv run ruff check .
uv run ruff format --check .
```

## Running the study

Planned (milestone M5/M6, see [docs/DESIGN.md](docs/DESIGN.md#9-cli)):

```sh
uv run peekaudit run --config quick --out runs/quick     # ~20 min on 4 CPU workers
uv run peekaudit report runs/quick --publish docs/results/quick
uv run peekaudit run --config full --out runs/full       # several hours; resumable
```

## Scope and limitations

- gbstats is audited at its latest PyPI release (0.8.0, June 2024). GrowthBook's hosted
  product runs newer code under the same version string, so findings may not transfer.
- Product-layer behaviour (for example GrowthBook's minimum sample size before results are
  shown) is not part of the engine libraries and is treated separately.
- Calibration findings describe statistical behaviour under stated conditions. They are not
  claims of software defects. Any suspected defect is reported only with a minimal
  standalone reproduction.

## Licence

MIT. See [LICENSE](LICENSE). The audited engines keep their own licences.
