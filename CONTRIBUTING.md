# Contributing

Bug reports, engine adapters, scenarios and documentation fixes are welcome. Read
[docs/DESIGN.md](docs/DESIGN.md) first: it fixes the formulas, the validation gates and the
reporting rules.

## Setup

You need [uv](https://docs.astral.sh/uv/) and git. From a clone:

```sh
uv sync
uv run peekaudit versions
```

Python 3.11 is the canonical study environment (`.python-version`); 3.12 and 3.13 are
supported for development and CI.

## Tests and lint

```sh
uv run pytest              # fast suite, must stay under a few minutes
uv run pytest --runslow    # Monte Carlo checks
uv run ruff check .
uv run ruff format --check .
```

CI runs the fast suite on Ubuntu (Python 3.11 and 3.13) and Windows (Python 3.13), plus
lint. Tests must run offline, write only under pytest's `tmp_path`, use fixed seeds, and pass
on Windows (use `pathlib`). Mark Monte Carlo tests that take more than a few seconds with
`@pytest.mark.slow`. Statistical assertions must use an explicit exact interval (for example
a 99.9% Clopper-Pearson interval), never a hand-picked tolerance.

## Code

- Python 3.11+, typed. Core runtime dependencies are numpy and scipy only; engines are an
  optional extra. Open an issue before adding a dependency.
- Engines are imported, never vendored. Copy numerical test vectors from an engine's test
  suite only with an attribution comment (file, version or commit).
- Match the surrounding style; ruff enforces the rest (line length 100).
- Write text files as UTF-8.

## Engine adapters

An adapter provides a vectorised `evaluate` and an engine-in-the-loop `evaluate_engine`
(see docs/DESIGN.md, sections 2.1 and 5.2). Before any audit number is reported it must:

1. reproduce the engine's own published fixtures,
2. agree with the engine on the randomised input grid, including degenerate inputs, and
3. pass the in-the-loop agreement check in every study unit.

An apparent engine defect needs a minimal standalone reproduction that calls the engine
directly, plus an independent check that it is not an adapter mistake.

## Pull requests

- One change per pull request, with tests and docs.
- Add an entry under "Unreleased" in [CHANGELOG.md](CHANGELOG.md) for user-visible changes.
- Keep the "Status and next steps" section of [AGENTS.md](AGENTS.md) current.
- Do not commit study runs, caches or downloaded data; only small derived results under
  `docs/results/`.

Security issues: see [SECURITY.md](SECURITY.md). Do not open a public issue.
