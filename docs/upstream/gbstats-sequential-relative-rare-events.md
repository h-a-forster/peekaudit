# SequentialTwoSidedTTest: false-positive rate above alpha for relative effects with few events

**Status:** draft, not filed. Target: GrowthBook repository (gbstats).

## Summary

The documentation says sequential results can be checked as often as you like while keeping
false positives below the expected rate. In an A/A simulation with a 1% conversion rate,
20,000 users split 50/50 and 140 looks, `SequentialTwoSidedTTest` with the default relative
effect rejects in about 14% of experiments instead of at most 5%. With
`difference_type="absolute"` on the same data the rate is about 2%. The excess appears
while the control rate rests on few events, which points at the delta-method variance of
the relative effect rather than at the confidence-sequence boundary.

## Minimal reproduction

`scripts/gbstats_rare_conversion_repro.py` (gbstats, numpy and scipy only; about 35 s):

```python
import numpy as np
from gbstats.frequentist.tests import SequentialConfig, SequentialTwoSidedTTest
from gbstats.models.statistics import ProportionStatistic

rng = np.random.default_rng(1)
reps, looks, n, p = 1000, 140, 20_000, 0.01
totals = np.floor(np.arange(1, looks + 1) * n / looks + 0.5).astype(int)
n_b = totals // 2
n_a = totals - n_b
conv_a = rng.binomial(np.diff(n_a, prepend=0), p, size=(reps, looks)).cumsum(axis=1)
conv_b = rng.binomial(np.diff(n_b, prepend=0), p, size=(reps, looks)).cumsum(axis=1)

for diff_type in ("relative", "absolute"):
    cfg = SequentialConfig(difference_type=diff_type, alpha=0.05)  # tuning 5000 (default)
    fp = 0
    for i in range(reps):
        for j in range(looks):
            res = SequentialTwoSidedTTest(
                ProportionStatistic(n=int(n_a[j]), sum=float(conv_a[i, j])),
                ProportionStatistic(n=int(n_b[j]), sum=float(conv_b[i, j])),
                cfg,
            ).compute_result()
            if res.error_message is None and res.p_value < 0.05:
                fp += 1
                break
    print(diff_type, fp / reps)
```

## Observed (gbstats 0.8.0 from PyPI)

| effect | false-positive rate, 1,000 replicates (95% CI) | 20,000 replicates (95% CI) |
|---|---|---|
| relative | 15.2% [13.0, 17.6] | 14.1% [13.7, 14.6] |
| absolute | 1.5% [0.8, 2.5] | 2.1% [1.9, 2.3] |

Median control conversions at the first false rejection (relative): 8.

Other cells of the same study (20,000 replicates each, relative effect, default tuning):

- 1% conversions, 20,000 users, 50/50, 14 daily looks: 6.3% [6.0, 6.7].
- 1% conversions, 1,000 users, 50/50, 140 looks: 9.1% [8.7, 9.5].
- Lognormal revenue (sigma 1.5), 20,000 users, 20% to treatment, 140 looks: 11.9%
  [11.5, 12.4]; 9.0% [8.6, 9.4] counting only looks where both arms have at least 150 users.
- Lognormal revenue, 1,000 users, 20% to treatment, 140 looks: 16.7% [16.2, 17.2], of which
  93% start at looks with fewer than 150 users in an arm (3.7% [3.5, 4.0] counting only
  looks with at least 150 per arm). With `sequential_tuning_parameter` set to the planned
  size (1,000) the rate is 37.2% [36.6, 37.9].
- Normal data, 10% conversions and a ratio metric: no excess detected (each cell can detect
  an FPR of about 5.5% or more).

## Expected

An FPR of at most 5% (the design target of an anytime-valid test), or documentation of the
conditions under which the asymptotic guarantee applies, for example a minimum number of
events or users per arm, and a note that the relative effect needs more data than the
absolute one.

## Notes

- The confidence sequence is asymptotic (Waudby-Smith et al.), so some excess at small
  samples is expected in theory; the request is mainly to document a minimum, or to apply
  one in the library as the GrowthBook product does in its UI.
- The study's engine figures come from a vectorised re-implementation of the engine
  arithmetic that agrees with `SequentialTwoSidedTTest` exactly on the 183,040 replicate-looks
  where both were computed in the study and on randomised input grids; the reproduction
  above calls the engine only.

## Environment

- gbstats 0.8.0 (PyPI, released 2024-06-20)
- Python 3.11.9, numpy 1.26.4, scipy 1.17.1, pandas 1.5.3
