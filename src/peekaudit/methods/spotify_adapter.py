"""spotify-confidence 4.1.0 adapter: ``ZTest`` group sequential and fixed-horizon z-tests.

Mirrors follow ``analysis/frequentist/confidence_computers/z_test_computer.py`` and
``multiple_comparison.py``: variance Q/n - (S/n)^2 when every row of the call is binary
(sum_squares == sum), else (Q - S^2/n)/(n - 1); SE = sqrt(v_A/n_A + v_B/n_B); two-sided
p = 2 Phi(-|D|/SE). The engine's alpha is ``1 - interval_size`` = 1 - 0.95, which in floating
point is 0.050000000000000044; the mirror uses the same value.

Group sequential semantics (DESIGN 4.2): at look k the analyst passes looks 1..k with the
planned total N_max. Information fractions are t_j = N_j / max(N_max, max_{i<=k} N_i), so
past the planned size every earlier fraction is rescaled. Bounds come from the engine's own
``sequential_bounds`` (power-family spending, rho = 2), cached per fraction vector: they do
not depend on the data because arm counts are deterministic.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
from scipy.stats import norm

from peekaudit.methods.base import EvalContext, Method, MethodOutput
from peekaudit.scenarios import Metric
from peekaudit.simulate import LookStats

#: The engine's default interval size; its alpha is computed as 1 - INTERVAL_SIZE.
INTERVAL_SIZE = 0.95
SPOTIFY_ALPHA = 1 - INTERVAL_SIZE


@lru_cache(maxsize=256)
def engine_bounds(t: tuple[float, ...], alpha: float = SPOTIFY_ALPHA) -> np.ndarray:
    """The engine's two-sided sequential z bounds for information fractions ``t`` (cached)."""
    from spotify_confidence.analysis.frequentist.confidence_computers.z_test_computer import (
        sequential_bounds,
    )

    zb = sequential_bounds(np.asarray(t, dtype=np.float64), alpha=alpha, sides=2).df["zb"]
    return zb.to_numpy(dtype=np.float64).copy()


def look_bounds(totals: np.ndarray, n_max: int, alpha: float = SPOTIFY_ALPHA) -> np.ndarray:
    """Bound used at each look under look-by-look semantics (see module docstring)."""
    totals = np.asarray(totals, dtype=np.int64)
    k = totals.size
    zb = np.empty(k)
    running_max = np.maximum.accumulate(totals)
    m = int(np.sum(running_max <= n_max))  # running max is non-decreasing: a prefix
    if m:
        t = tuple((totals[:m] / max(int(running_max[m - 1]), n_max)).tolist())
        zb[:m] = engine_bounds(t, alpha)
    for j in range(m, k):
        denom = max(int(running_max[j]), n_max)
        t = tuple((totals[: j + 1] / denom).tolist())
        zb[j] = engine_bounds(t, alpha)[-1]
    return zb


def _moments(stats: LookStats, cumulative: bool) -> dict[str, np.ndarray]:
    """Engine variance branch, D, SE, z and p per (replicate, look)."""
    na = np.broadcast_to(stats.n_a.astype(np.float64), stats.sum_a.shape)
    nb = np.broadcast_to(stats.n_b.astype(np.float64), stats.sum_b.shape)
    binary = (stats.sumsq_a == stats.sum_a) & (stats.sumsq_b == stats.sum_b)
    if cumulative:  # look-by-look calls pass looks 1..k: the binary check spans them all
        binary = np.logical_and.accumulate(binary, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = {}
        bad = (na == 0) | (nb == 0)
        for arm, n in (("a", na), ("b", nb)):
            s = getattr(stats, f"sum_{arm}")
            q = getattr(stats, f"sumsq_{arm}")
            pe = s / n
            v = np.where(binary, q / n - pe**2, (q - np.power(s, 2) / n) / (n - 1))
            bad |= v < 0
            out[f"pe_{arm}"], out[f"v_{arm}"] = pe, v
        if cumulative:
            bad = np.logical_or.accumulate(bad, axis=1)  # a refused row fails later calls
        se = np.sqrt(out["v_a"] / na + out["v_b"] / nb)
        d = out["pe_b"] - out["pe_a"]
        z = (out["pe_b"] - out["pe_a"] - 0) / se
        p = norm.sf(np.abs(z)) * 2
    return {"d": d, "se": se, "p": p, "bad": bad}


def _finish(m: dict[str, np.ndarray], adj_alpha: np.ndarray | float) -> MethodOutput:
    bad = m["bad"]
    with np.errstate(invalid="ignore"):
        sig = (m["p"] < adj_alpha) & ~bad
        zcrit = norm.ppf(1 - np.asarray(adj_alpha) / 2.0)
        lo = m["d"] - zcrit * m["se"]
        hi = m["d"] + zcrit * m["se"]
        margin = (adj_alpha - m["p"]) / adj_alpha
    sign = np.where(sig, np.sign(m["d"]), 0).astype(np.int8)
    return MethodOutput(
        sign,
        np.where(bad, np.nan, lo),
        np.where(bad, np.nan, hi),
        bad.copy(),
        np.where(bad, np.nan, margin),
    )


def _frame(stats: LookStats, i: int, looks: range, n_max: int | None) -> Any:
    import pandas as pd

    rows = []
    for j in looks:
        for g in ("a", "b"):
            row = {
                "look": j + 1,
                "group": g.upper(),
                "sum": float(getattr(stats, f"sum_{g}")[i, j]),
                "sum_squares": float(getattr(stats, f"sumsq_{g}")[i, j]),
                "n": int(getattr(stats, f"n_{g}")[j]),
            }
            if n_max is not None:
                row["planned"] = int(n_max)
            rows.append(row)
    return pd.DataFrame(rows)


def _difference(frame: Any, sequential: bool) -> Any:
    import spotify_confidence as sc

    test = sc.ZTest(
        frame,
        numerator_column="sum",
        numerator_sum_squares_column="sum_squares",
        denominator_column="n",
        categorical_group_columns="group",
        ordinal_group_column="look",
        interval_size=INTERVAL_SIZE,
    )
    kw = {"final_expected_sample_size_column": "planned"} if sequential else {}
    return (
        test.difference(level_1="A", level_2="B", groupby="look", verbose=True, **kw)
        .sort_values("look")
        .reset_index(drop=True)
    )


class _SpotifyMethod(Method):
    def supports(self, metric: Metric) -> bool:
        return not metric.ratio


class SpotifyGST(_SpotifyMethod):
    """``ZTest`` with ``final_expected_sample_size`` = N_max (group sequential test)."""

    def __init__(self) -> None:
        super().__init__(
            "spotify.gst", "absolute", "planned", "spotify GST (power spending, rho 2)", "spotify"
        )

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        totals = stats.n_a + stats.n_b
        zb = look_bounds(totals, ctx.n_max)
        adj = 2 * (1 - norm.cdf(zb))
        return _finish(_moments(stats, cumulative=True), adj[None, :])

    def evaluate_engine(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        """Real engine calls with look-by-look semantics.

        Looks whose running maximum stays within N_max share one full-path call (the engine
        gives identical rows); every later look gets its own call on looks 1..k.
        """
        r, k = stats.n_reps, stats.n_looks
        out = MethodOutput.empty((r, k))
        running = np.maximum.accumulate(stats.n_a + stats.n_b)
        m = int(np.sum(running <= ctx.n_max))
        for i in range(r):
            calls: list[tuple[range, list[int]]] = []
            if m:
                calls.append((range(m), list(range(m))))
            calls += [(range(j + 1), [j]) for j in range(m, k)]
            for looks, targets in calls:
                try:
                    res = _difference(_frame(stats, i, looks, ctx.n_max), sequential=True)
                except ValueError:
                    continue  # engine refuses the input: recorded as "no result"
                for j in targets:
                    row = res.iloc[j]
                    out.no_result[i, j] = False
                    out.lower[i, j] = float(row["adjusted ci_lower"])
                    out.upper[i, j] = float(row["adjusted ci_upper"])
                    sig = bool(row["is_significant"])
                    out.reject[i, j] = int(np.sign(row["difference"])) if sig else 0
        return out


class SpotifyZ(_SpotifyMethod):
    """``ZTest`` without sequential correction, re-run at every look (peeked)."""

    def __init__(self) -> None:
        super().__init__("spotify.z", "absolute", "fixed", "spotify z-test (peeked)", "spotify")

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:  # noqa: ARG002
        return _finish(_moments(stats, cumulative=False), SPOTIFY_ALPHA)

    def evaluate_engine(
        self, stats: LookStats, ctx: EvalContext, per_look: bool = False
    ) -> MethodOutput:
        """Real engine calls.

        ``per_look=True`` makes one single-look call per look (what a peeking analyst runs).
        The default makes one call per replicate on all looks and reads the unadjusted
        p-value and 95% CI per row; with one comparison per call these equal the single-look
        decision ``p < alpha`` and interval (tested), at a fraction of the cost.
        """
        r, k = stats.n_reps, stats.n_looks
        out = MethodOutput.empty((r, k))
        for i in range(r):
            groups = [range(j, j + 1) for j in range(k)] if per_look else [range(k)]
            for looks in groups:
                try:
                    res = _difference(_frame(stats, i, looks, None), sequential=False)
                except ValueError:
                    if per_look or k == 1:
                        continue
                    # One bad row fails the whole call; fall back to single-look calls.
                    sub = self.evaluate_engine(stats.rows(slice(i, i + 1)), ctx, per_look=True)
                    for f in ("reject", "lower", "upper", "no_result"):
                        getattr(out, f)[i] = getattr(sub, f)[0]
                    break
                for pos, j in enumerate(looks):
                    row = res.iloc[pos]
                    p = row["p-value"]
                    out.no_result[i, j] = False
                    out.lower[i, j] = float(row["ci_lower"])
                    out.upper[i, j] = float(row["ci_upper"])
                    sig = p is not None and float(p) < SPOTIFY_ALPHA
                    out.reject[i, j] = int(np.sign(row["difference"])) if sig else 0
        return out
