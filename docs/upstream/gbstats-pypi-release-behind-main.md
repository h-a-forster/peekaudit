# gbstats on PyPI (0.8.0) differs from the code on main, which also reports 0.8.0

**Status:** draft, not filed. Target: GrowthBook repository. **Re-check before filing:** this
draft is based on a source comparison made while building peekaudit, not on a script, and
the repository may have changed since.

## Summary

The latest gbstats release on PyPI is 0.8.0 (uploaded 2024-06-20). The `packages/stats`
directory on GrowthBook's main branch still declares version 0.8.0 but contains different
code (for example, post-stratification is on by default). Users who `pip install gbstats`
therefore get statistics that differ from what GrowthBook computes, under the same version
number.

## Steps to check

1. `pip download gbstats==0.8.0 --no-deps --no-binary :all:` and unpack the sdist.
2. Check out GrowthBook main and compare `packages/stats/gbstats/` with the sdist
   (for example `diff -r`).
3. Compare the version strings (`pyproject.toml` / `setup.py` on main vs the sdist).

## Expected

A version bump on main when the code diverges from the released package, and a new PyPI
release (or a note in the package README that the PyPI package is not kept in sync).

## Impact

Independent evaluations of "the GrowthBook stats engine" that install from PyPI (as ours
did) test older code than the product runs. Results from such evaluations cannot be
attributed to the current product without this caveat.

## Environment

- gbstats 0.8.0 from PyPI
- GrowthBook main as of the comparison (record the commit when re-checking)
