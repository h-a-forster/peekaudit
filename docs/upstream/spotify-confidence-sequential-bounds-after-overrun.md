# Sequential ZTest: bounds of earlier looks change when later looks exceed the planned sample size

**Status:** draft, not filed. Target: spotify-confidence issue tracker.

## Summary

In a group sequential analysis (`ZTest.difference(..., final_expected_sample_size_column=...)`),
the adjusted alpha reported for a look depends on which later looks are passed in the same
call once the cumulative sample in the call exceeds the planned size. The later call
recomputes every earlier bound as if the plan had been the largest sample in the call. A
look that was significant when it was made can therefore be reported as not significant in a
later re-analysis of the same data. Within the planned size the results are consistent.

This may be intended, since the planned size was wrong, but it is not documented as far as
we can tell, and it makes past results irreproducible.

## Minimal reproduction

`scripts/spotify_overrun_repro.py` (deterministic data, no simulation): two arms, planned
size 800 units in total, a look every 100 units up to 1,200 (12 looks).

```python
import pandas as pd
import spotify_confidence as sc

PLANNED = 800


def frame(last_look):
    rows = []
    for k in range(1, last_look + 1):
        n = 50 * k  # per arm
        for group, mean in (("A", 1.0), ("B", 1.02)):
            rows.append(
                {
                    "look": k,
                    "group": group,
                    "sum": mean * n,
                    "sum_squares": (1.0 + mean**2) * n,
                    "n": n,
                    "planned": PLANNED,
                }
            )
    return pd.DataFrame(rows)


def adjusted_alpha(last_look):
    res = sc.ZTest(
        frame(last_look),
        numerator_column="sum",
        numerator_sum_squares_column="sum_squares",
        denominator_column="n",
        categorical_group_columns="group",
        ordinal_group_column="look",
    ).difference(
        level_1="A",
        level_2="B",
        groupby="look",
        final_expected_sample_size_column="planned",
        verbose=True,
    )
    return res.sort_values("look").set_index("look")["adjusted_alpha"].astype(float)


full = adjusted_alpha(12)
for k in range(1, 13):
    print(k, round(adjusted_alpha(k).loc[k], 4), round(full.loc[k], 4))
```

## Observed (spotify-confidence 4.1.0)

| look | total units | adjusted alpha, call with looks 1..k | adjusted alpha, call with looks 1..12 |
|---|---|---|---|
| 1 | 100 | 0.0008 | 0.0003 |
| 4 | 400 | 0.0083 | 0.0036 |
| 8 | 800 (planned) | 0.0311 | 0.0128 |
| 9 | 900 | 0.0304 | 0.0161 |
| 10 | 1,000 | 0.0297 | 0.0198 |
| 12 | 1,200 | 0.0287 | 0.0287 |

The right-hand column equals Lan-DeMets power-family (rho = 2) bounds at information
fractions k/12, that is, a plan of 1,200 units; we checked this with an independent
recursive-integration solver. The left-hand column equals fractions k/8 up to look 8. Past
look 8, each one-call-per-look analysis treats the current look as the final one.

## Expected

One of:

1. Bounds at looks within the planned size do not change when later rows are added (the
   alpha already spent stays spent), with overrun looks handled by an explicit rule; or
2. A documented statement that, once the sample exceeds the planned size, all bounds are
   recomputed from the largest sample in the call, so earlier results can change, ideally
   with a warning when this happens.

## Impact

Our simulation (normal data, 21 looks to 1.5 x the planned size, one call per look) gives a
false-positive rate of 7.3% to 8.2% instead of 5%, which matches the vendor's own
observation that group sequential tests inflate when oversampled. The issue here is the
retroactive change of earlier bounds, not the inflation itself.

## Environment

- spotify-confidence 4.1.0 (PyPI)
- Python 3.11.9, pandas 1.5.3, numpy 1.26.4, scipy 1.17.1
