"""Minimal standalone reproduction (spotify-confidence 4.1.0 only, no peekaudit code).

``ZTest.difference`` reports a significant difference when one group holds a single unit of
0/1 data. For binary data the variance is computed as p(1 - p), which is 0 for one unit, so
the standard error is 0 (when the other group is also constant) or rests on the other group
alone, and the test declares significance. Sample variances (n - 1 denominator) are
undefined for one unit, so a test built on them gives no result.

In a group sequential analysis with 1,000 planned users, 20% to B and frequent looks, the
first look has 6 users in A and 1 in B; this behaviour then accounts for about half of the
false positives of an A/A test with 10% conversions (see docs/results/full/gst_reference.md).

Usage: uv run python scripts/spotify_single_unit_repro.py
"""

from __future__ import annotations

import warnings

import pandas as pd
import spotify_confidence as sc

warnings.filterwarnings("ignore")

CASES = (
    ("A: 6 users, 0 conversions; B: 1 user, 1 conversion", (6, 0), (1, 1)),
    ("A: 6 users, 1 conversion;  B: 1 user, 1 conversion", (6, 1), (1, 1)),
)


def frame(a: tuple[int, int], b: tuple[int, int], planned: int | None) -> pd.DataFrame:
    rows = []
    for group, (n, conversions) in (("A", a), ("B", b)):
        row = {"group": group, "sum": float(conversions), "sum_squares": float(conversions)}
        row.update({"n": n, "look": 1})
        if planned is not None:
            row["planned"] = planned
        rows.append(row)
    return pd.DataFrame(rows)


def difference(a: tuple[int, int], b: tuple[int, int], planned: int | None) -> pd.Series:
    test = sc.ZTest(
        frame(a, b, planned),
        numerator_column="sum",
        numerator_sum_squares_column="sum_squares",
        denominator_column="n",
        categorical_group_columns="group",
        ordinal_group_column="look",
    )
    kw = {"final_expected_sample_size_column": "planned"} if planned else {}
    res = test.difference(level_1="A", level_2="B", groupby="look", **kw)
    return res.iloc[0]


def main() -> None:
    for label, a, b in CASES:
        for planned in (None, 1000):
            row = difference(a, b, planned)
            kind = "sequential, planned 1000" if planned else "fixed horizon"
            pre = "adjusted " if planned else ""
            ci = (row[f"{pre}ci_lower"], row[f"{pre}ci_upper"])
            print(
                f"{label} | {kind}: difference {row['difference']:.4f}, "
                f"CI [{ci[0]:.4f}, {ci[1]:.4f}], is_significant {bool(row['is_significant'])}"
            )


if __name__ == "__main__":
    main()
