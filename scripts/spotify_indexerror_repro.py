"""Minimal standalone reproduction (spotify-confidence 4.1.0 only, no peekaudit code).

``ZTest.difference`` raises ``IndexError: list index out of range`` when the control group
has a single unit whose value is 0 (sum == sum_squares == 0) and the treatment group is not
binary. The control row then passes the "binary" check used by the powered-effect search
(``_search_MDE_binary``), while its variance is 0/0 = NaN, so every candidate effect gives
NaN, no minimum is found and indexing the empty list fails. A single control unit with
value 1 does not crash (an early return), and nor does an all-zero control with n >= 2.

This is a robustness issue on a degenerate input (one control unit), outside every cell of
the peekaudit study (control arms start at 3 units). It does not affect the reported
statistics. peekaudit's adapter records the exception as "no result" and its mirror
reproduces the branch.

Usage: uv run python scripts/spotify_indexerror_repro.py
"""

from __future__ import annotations

import warnings

import pandas as pd
import spotify_confidence as sc

warnings.filterwarnings("ignore")

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
try:
    test.difference(level_1="A", level_2="B")
    print("no exception (behaviour differs from spotify-confidence 4.1.0)")
except IndexError as exc:
    print(f"IndexError: {exc}")
