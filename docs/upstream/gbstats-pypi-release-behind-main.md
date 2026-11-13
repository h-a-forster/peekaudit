# gbstats on PyPI (0.8.0) differs from the code on main, which also reports 0.8.0

**Status:** filed as [growthbook#7642](https://github.com/growthbook/growthbook/issues/7642) on 2026-11-07; text as written before filing. Target: GrowthBook repository. **Last re-checked:** 2026-10-30
against GrowthBook main at commit 830555b9c1 and gbstats 0.8.0, still the latest release on
PyPI: still applies (details under "Re-check" below). Re-check again before filing.

## Summary

The latest gbstats release on PyPI is 0.8.0 (uploaded 2024-06-20). The `packages/stats`
directory on GrowthBook's main branch still declares version 0.8.0 but contains different
code (for example, post-stratification support, which the GrowthBook product enables by
default, and a changed test API). Users who `pip install gbstats` therefore get statistics
that differ from what GrowthBook computes, under the same version number.

## Steps to check

1. `pip download gbstats==0.8.0 --no-deps --no-binary :all:` and unpack the sdist.
2. Check out GrowthBook main and compare `packages/stats/gbstats/` with the sdist
   (for example `diff -r`).
3. Compare the version strings (`pyproject.toml` / `setup.py` on main vs the sdist).

## Re-check (2026-10-30)

- `packages/stats/pyproject.toml` on main (830555b9c1) declares `version = "0.8.0"`; the
  PyPI 0.8.0 wheel metadata requires `pandas>=1.5.1,<2.0.0`, main requires `pandas ^2.1.0`.
- `gbstats/frequentist/tests.py`, `models/statistics.py`, `models/tests.py` and `utils.py`
  all differ from the PyPI wheel (`diff` of 530, 531, 1,573 and 158 lines); main adds
  `bayesian/bandits.py` and `power/midexperimentpower.py`.
- The test API changed: on main `SequentialTwoSidedTTest` takes a list of
  `(stat_a, stat_b)` pairs instead of two statistics, and results use `pValue` /
  `errorMessage` instead of `p_value` / `error_message`, so code written against the PyPI
  package fails on main with a `TypeError`.
- The library flag `post_stratification_enabled` defaults to `False`; the product default
  (`DEFAULT_POST_STRATIFICATION_ENABLED` in `packages/shared/src/constants.ts`) is `true`.
- Eight commits have touched `packages/stats` since 2026-07-01 (latest 2026-10-13,
  bddb270089), none of them with a version bump.

## Expected

A version bump on main when the code diverges from the released package, and a new PyPI
release (or a note in the package README that the PyPI package is not kept in sync).

## Impact

Independent evaluations of "the GrowthBook stats engine" that install from PyPI (as ours
did) test older code than the product runs. Results from such evaluations cannot be
attributed to the current product without this caveat.

## Environment

- gbstats 0.8.0 from PyPI
- GrowthBook main at 830555b9c1 (2026-10-30)
