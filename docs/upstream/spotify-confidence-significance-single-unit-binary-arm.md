# ZTest reports a significant difference when a group has one unit of binary data

**Status:** draft, not filed. Target: spotify-confidence issue tracker. Kind: robustness /
documentation request, not a computational error.

## Summary

For 0/1 data, `ZTest.difference` computes each group's variance as p(1 - p). A group with a
single unit has p in {0, 1}, so its variance is 0. The standard error of the difference then
comes from the other group alone, or is exactly 0 when the other group is constant, and the
comparison is reported as significant (p-value 0, zero-width interval in the second case).
This happens with and without `final_expected_sample_size_column`. No minimum group size is
applied or documented.

## Minimal reproduction

`scripts/spotify_single_unit_repro.py` (imports only spotify-confidence and pandas):

```python
import pandas as pd
import spotify_confidence as sc

df = pd.DataFrame(
    [
        {"group": "A", "sum": 0.0, "sum_squares": 0.0, "n": 6},
        {"group": "B", "sum": 1.0, "sum_squares": 1.0, "n": 1},
    ]
)
test = sc.ZTest(
    df,
    numerator_column="sum",
    numerator_sum_squares_column="sum_squares",
    denominator_column="n",
    categorical_group_columns="group",
)
print(test.difference(level_1="A", level_2="B"))
```

## Observed (spotify-confidence 4.1.0)

| A | B | analysis | difference | 95% CI | is_significant |
|---|---|---|---|---|---|
| 6 users, 0 conversions | 1 user, 1 conversion | fixed horizon | 1.000 | [1.000, 1.000] | True |
| 6 users, 0 conversions | 1 user, 1 conversion | sequential, planned 1,000 | 1.000 | [1.000, 1.000] | True |
| 6 users, 1 conversion | 1 user, 1 conversion | fixed horizon | 0.833 | [0.535, 1.132] | True |
| 6 users, 1 conversion | 1 user, 1 conversion | sequential, planned 1,000 | 0.833 | [0.116, 1.550] | True |

## Expected

A z-test needs a variance estimate for each group. With one observation there is none
(the sample variance is undefined), so a NaN result, a warning, or a documented minimum
group size would be expected rather than a significant result.

## Impact in a simulated A/A study

With a planned size of 1,000 users, 20% to B and 140 equally spaced looks, the first look
has 6 users in A and 1 in B. For 10% conversions, the sequential test's false-positive rate
over 20,000 A/A replicates is 20.4% [19.8, 20.9]; 9.1% of replicates (1,812) reject at that
first look. A reference group sequential test with the same spending function (power family,
rho = 2) and sample variances, which gives no result for a one-unit group, has 11.5%
[11.0, 11.9] on the same replicates: every replicate the reference rejects is also rejected
by spotify-confidence, which rejects 1,780 more (exact McNemar test, Holm-adjusted p < 1e-300).
The remaining excess over 5% is a property of the z-based group
sequential test on small arms, shared by the reference
(`docs/results/full/gst_reference.md` in the peekaudit repository).

## Environment

- spotify-confidence 4.1.0 (PyPI)
- Python 3.11.9, pandas 1.5.3, numpy 1.26.4, scipy 1.17.1
