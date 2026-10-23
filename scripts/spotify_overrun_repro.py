"""Minimal standalone reproduction (spotify-confidence 4.1.0 only, no peekaudit code).

In a group sequential ``ZTest.difference`` call, the bound reported for a look depends on
which *later* looks are passed in the same call, once the cumulative sample in the call
exceeds ``final_expected_sample_size``. Here a two-arm experiment is planned for 800 units
and looked at every 100 units up to 1,200 (12 looks). For each look k the script prints the
adjusted alpha from a call with looks 1..k (what an analyst sees at look k) and from a call
with all 12 looks (what a later re-analysis shows for look k).

Observed with 4.1.0: the two agree while the call stays within the planned size, but the
12-look call recomputes every earlier bound as if the plan had been 1,200 units (for example
look 8: 0.0311 at the time, 0.0128 in the re-analysis; look 10: 0.0297 vs 0.0198). A look
that was significant when it was made can therefore be reported as not significant later.
Looking past the plan one call at a time instead treats each new look as the final one.

The data are deterministic (no simulation), so the output is the same on every run.

Usage: uv run python scripts/spotify_overrun_repro.py
"""

from __future__ import annotations

import warnings

import pandas as pd
import spotify_confidence as sc

warnings.filterwarnings("ignore")

PLANNED = 800
LOOKS = list(range(1, 13))  # 100, 200, ..., 1200 units in total


def frame(last_look: int) -> pd.DataFrame:
    rows = []
    for k in LOOKS[:last_look]:
        n = 50 * k  # per arm
        for group, mean in (("A", 1.0), ("B", 1.02)):
            rows.append(
                {
                    "look": k,
                    "group": group,
                    "sum": mean * n,
                    "sum_squares": (1.0 + mean**2) * n,  # unit variance
                    "n": n,
                    "planned": PLANNED,
                }
            )
    return pd.DataFrame(rows)


def adjusted_alpha(last_look: int) -> pd.Series:
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


def main() -> None:
    full = adjusted_alpha(LOOKS[-1])
    print(f"planned size {PLANNED}; adjusted alpha per look")
    print("look  total  call with looks 1..k  call with all 12 looks")
    for k in LOOKS:
        now = adjusted_alpha(k).loc[k]
        print(f"{k:>4}  {100 * k:>5}  {now:>20.4f}  {full.loc[k]:>22.4f}")


if __name__ == "__main__":
    main()
