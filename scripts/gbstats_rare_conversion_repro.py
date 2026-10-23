"""Standalone Monte Carlo reproduction (gbstats 0.8.0, numpy and scipy only; no peekaudit code).

A/A test, 1% conversion rate in both arms, 20,000 users split 50/50, results checked at 140
equally spaced looks with ``SequentialTwoSidedTTest`` (default tuning 5,000, alpha 0.05). A
replicate is a false positive if any look has p < alpha. The anytime-valid design target is
an FPR of at most 5%.

Observed with gbstats 0.8.0 (1,000 replicates, seed 1, about 35 s): relative effect 15.2%
[13.0, 17.6] false positives; absolute effect 1.5% [0.8, 2.5] on the same data. The
peekaudit study (20,000 replicates, mirror checked against the engine) measured 14.1%
[13.7, 14.6] (relative) and 2.1% [1.9, 2.3] (absolute) in this cell. The script also prints
the median number of control conversions at the first false rejection (8 for the relative
effect), consistent with the relative interval's delta-method variance being unreliable while
the control rate rests on few events.

Usage: uv run python scripts/gbstats_rare_conversion_repro.py [--reps 1000] [--seed 1]
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np
from gbstats.frequentist.tests import SequentialConfig, SequentialTwoSidedTTest
from gbstats.models.statistics import ProportionStatistic
from scipy import stats as st

warnings.filterwarnings("ignore")


def clopper_pearson(k: int, n: int) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(st.beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(st.beta.ppf(0.975, k + 1, n - k))
    return lo, hi


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--looks", type=int, default=140)
    ap.add_argument("--n", type=int, default=20_000, help="total users at the last look")
    ap.add_argument("--p", type=float, default=0.01, help="conversion rate in both arms")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    totals = np.floor(np.arange(1, args.looks + 1) * args.n / args.looks + 0.5).astype(int)
    n_b = totals // 2
    n_a = totals - n_b
    # cumulative conversions per arm and look, from independent increments
    inc_a = rng.binomial(np.diff(n_a, prepend=0), args.p, size=(args.reps, args.looks))
    inc_b = rng.binomial(np.diff(n_b, prepend=0), args.p, size=(args.reps, args.looks))
    conv_a, conv_b = inc_a.cumsum(axis=1), inc_b.cumsum(axis=1)

    for diff_type in ("relative", "absolute"):
        cfg = SequentialConfig(difference_type=diff_type, alpha=0.05)
        false_pos, control_events = 0, []
        for i in range(args.reps):
            for j in range(args.looks):
                res = SequentialTwoSidedTTest(
                    ProportionStatistic(n=int(n_a[j]), sum=float(conv_a[i, j])),
                    ProportionStatistic(n=int(n_b[j]), sum=float(conv_b[i, j])),
                    cfg,
                ).compute_result()
                if res.error_message is None and res.p_value < 0.05:
                    false_pos += 1
                    control_events.append(int(conv_a[i, j]))
                    break
        lo, hi = clopper_pearson(false_pos, args.reps)
        med = f"{np.median(control_events):.0f}" if control_events else "-"
        print(
            f"{diff_type:<9} FPR {100 * false_pos / args.reps:5.1f}% "
            f"[{100 * lo:.1f}, {100 * hi:.1f}] ({false_pos}/{args.reps}); "
            f"median control conversions at first rejection: {med}"
        )


if __name__ == "__main__":
    main()
