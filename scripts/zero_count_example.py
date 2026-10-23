"""Standalone illustration (engines called directly, no peekaudit code).

Why fixed-horizon z-tests reject far too often on rare conversions with a small arm: when
the treatment arm has zero conversions, its plug-in variance p(1 - p) is 0, the standard
error comes from the control arm alone, and a difference of one control rate looks
"significant". This is a property of the Wald z-test, not a software defect.

With the study's bern_p01, N_max = 1,000, 20% to B cell (800 vs 200 units, p = 1%), B has
zero conversions with probability 0.99^200 = 13.4%, which is most of the 14.2% false-positive
rate measured at the single final look.

Usage: uv run python scripts/zero_count_example.py
"""

from __future__ import annotations

import warnings

import pandas as pd

warnings.filterwarnings("ignore")


def spotify() -> None:
    import spotify_confidence as sc

    df = pd.DataFrame(
        [
            {"group": "A", "sum": 8, "sum_squares": 8, "n": 800},
            {"group": "B", "sum": 0, "sum_squares": 0, "n": 200},
        ]
    )
    res = sc.ZTest(
        df,
        numerator_column="sum",
        numerator_sum_squares_column="sum_squares",
        denominator_column="n",
        categorical_group_columns="group",
    ).difference(level_1="A", level_2="B", verbose=True)
    row = res.iloc[0]
    print(
        f"spotify-confidence ZTest, A 8/800 vs B 0/200: difference {row['difference']:.4f}, "
        f"p = {row['p-value']:.4f}, significant = {bool(row['is_significant'])}"
    )


def gbstats() -> None:
    from gbstats.frequentist.tests import TwoSidedTTest
    from gbstats.models.statistics import ProportionStatistic

    res = TwoSidedTTest(
        ProportionStatistic(n=800, sum=8), ProportionStatistic(n=200, sum=0)
    ).compute_result()
    print(
        f"gbstats TwoSidedTTest, same data: p = {res.p_value}, "
        f"error = {res.error_message!r} (zero variance -> no result, so no rejection)"
    )


if __name__ == "__main__":
    print(f"P(B has zero conversions) = 0.99**200 = {0.99**200:.4f}")
    spotify()
    gbstats()
