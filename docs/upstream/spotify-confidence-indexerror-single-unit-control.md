# ZTest.difference raises IndexError when the control group has one unit with value 0

**Status:** draft, not filed. Target: spotify-confidence issue tracker.

## Summary

`ZTest.difference` raises `IndexError: list index out of range` from `_search_MDE_binary`
when the control group has a single unit whose value is 0 (`sum == sum_squares == 0`,
`n == 1`) and the treatment group is not binary. A clear error message or a NaN result
would be expected instead.

## Minimal reproduction

```python
import pandas as pd
import spotify_confidence as sc

df = pd.DataFrame(
    [
        {"group": "A", "sum": 0.0, "sum_squares": 0.0, "n": 1},
        {"group": "B", "sum": 120.0, "sum_squares": 900.0, "n": 50},
    ]
)
test = sc.ZTest(
    df,
    numerator_column="sum",
    numerator_sum_squares_column="sum_squares",
    denominator_column="n",
    categorical_group_columns="group",
)
test.difference(level_1="A", level_2="B")
```

## Observed

```
  File ".../spotify_confidence/analysis/frequentist/confidence_computers/z_test_computer.py", line 265, in <lambda>
    lambda row: _search_MDE_binary_local_search(
  File ".../z_test_computer.py", line 396, in _search_MDE_binary_local_search
    alt_result_x, alt_result_fun = _search_MDE_binary(
  File ".../z_test_computer.py", line 444, in _search_MDE_binary
    index = [index[int(np.ceil(len(index) / 2))]]
IndexError: list index out of range
```

## Expected

Either a result with NaN statistics (as for other zero-variance inputs) or a `ValueError`
that names the problem (for example "control group has fewer than 2 observations").

## Analysis

The control row passes the binary check used by the powered-effect search (its values are
all 0, so it looks like a proportion), while its variance is 0/0 = NaN. Every candidate
effect then evaluates to NaN, no minimum is found, and indexing the empty candidate list
fails. Related inputs behave differently:

- one control unit with value 1: no exception (an early return);
- an all-zero control with `n >= 2`: no exception.

## Environment

- spotify-confidence 4.1.0 (PyPI)
- Python 3.11.9, pandas 1.5.3, numpy 1.26.4, scipy 1.17.1

## Context

Found while auditing sequential tests by simulation; the input is degenerate (one control
unit) and did not occur in any analysed scenario. A downstream workaround is to catch the
exception and treat the comparison as "no result".
