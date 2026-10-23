"""Reference methods on the absolute scale (DESIGN 4.4): naive z, mSPRT, ACS, LD-OBF GST.

All use the plug-in variance V = s_A^2/n_A + s_B^2/n_B with n - 1 variances (delta method
for ratio metrics, see :func:`peekaudit.methods.base.sample_moments`). A look with n <= 1 in
an arm, or V <= 0, gives "no result" (no decision, no interval).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.stats import norm

from peekaudit.gst import lan_demets_bounds, spend_obf
from peekaudit.methods.base import EvalContext, Method, MethodOutput, sample_moments
from peekaudit.scenarios import fine_grid
from peekaudit.simulate import LookStats


def acs_rho(alpha: float, n_star: float) -> float:
    """Waudby-Smith et al. tuning rho (arXiv:2103.06476v7, eq. 161), as in gbstats."""
    return np.sqrt((-2 * np.log(alpha) + np.log(-2 * np.log(alpha) + 1)) / n_star)


def acs_pvalue_halfwidth(
    diff: np.ndarray, var: np.ndarray, n_total: np.ndarray, n_star: float, alpha: float
) -> tuple[np.ndarray, np.ndarray]:
    """Asymptotic confidence sequence p-value (eq. 155) and half-width (eq. 9).

    ``var`` is the variance of the estimated difference (not per unit); ``n_total`` is
    n_A + n_B. The operation order follows gbstats 0.8.0 so the gbstats mirror can reuse it.
    """
    rho = acs_rho(alpha, n_star)
    n = np.asarray(n_total, dtype=np.float64)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        s2 = var * n
        width = np.sqrt(s2) * np.sqrt(
            (2 * (n * np.power(rho, 2) + 1))
            * np.log(np.sqrt(n * np.power(rho, 2) + 1) / alpha)
            / (np.power(n * rho, 2))
        )
        st2 = np.power(diff - 0, 2) * n / var
        tr2p1 = n * np.power(rho, 2) + 1
        evalue = np.exp(np.power(rho, 2) * st2 / (2 * tr2p1)) / np.sqrt(tr2p1)
        p = np.minimum(1 / evalue, 1)
    return p, width


def _abs_diff_var(stats: LookStats) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ma, va, mb, vb = sample_moments(stats)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = va / stats.n_a + vb / stats.n_b
    d = mb - ma
    bad = ~np.isfinite(v) | (v <= 0) | ~np.isfinite(d)
    return d, v, bad


def _finish(
    d: np.ndarray, reject: np.ndarray, half: np.ndarray, margin: np.ndarray, bad: np.ndarray
) -> MethodOutput:
    sign = np.where(reject & ~bad, np.sign(d), 0).astype(np.int8)
    lo = np.where(bad, np.nan, d - half)
    hi = np.where(bad, np.nan, d + half)
    return MethodOutput(sign, lo, hi, bad.copy(), np.where(bad, np.nan, margin))


class NaiveZ(Method):
    """Fixed-horizon two-sided z-test, re-run at every look (the classic peeking error)."""

    def __init__(self) -> None:
        super().__init__("ref.naive_z", "absolute", "fixed", "naive z-test (peeked)")

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        d, v, bad = _abs_diff_var(stats)
        zc = float(norm.isf(ctx.alpha / 2))
        with np.errstate(divide="ignore", invalid="ignore"):
            se = np.sqrt(v)
            z = np.abs(d) / se
        return _finish(d, z > zc, zc * se, (z - zc) / zc, bad)


class MSPRT(Method):
    """Normal-mixture mSPRT (Johari et al. 2022) with plug-in variance and mixing sd tau."""

    def __init__(self) -> None:
        super().__init__("ref.msprt", "absolute", "anytime", "mSPRT (normal mixture)")

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        d, v, bad = _abs_diff_var(stats)
        t2 = ctx.tau**2
        thr = np.log(1.0 / ctx.alpha)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            log_lam = 0.5 * np.log(v / (v + t2)) + t2 * d * d / (2.0 * v * (v + t2))
            half = np.sqrt(v * (v + t2) / t2 * (2.0 * thr + np.log((v + t2) / v)))
        return _finish(d, log_lam >= thr, half, (log_lam - thr) / thr, bad)


class ACS(Method):
    """Asymptotic confidence sequence (Waudby-Smith et al.) with tuning N* = N_max."""

    def __init__(self) -> None:
        super().__init__("ref.acs", "absolute", "anytime", "asymptotic CS (N* = N_max)")

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        d, v, bad = _abs_diff_var(stats)
        n = (stats.n_a + stats.n_b).astype(np.float64)
        p, half = acs_pvalue_halfwidth(d, v, n, float(ctx.n_max), ctx.alpha)
        return _finish(d, p < ctx.alpha, half, (ctx.alpha - p) / ctx.alpha, bad)


@lru_cache(maxsize=64)
def obf_bounds(n_max: int, cols: tuple[int, ...]) -> np.ndarray:
    """LD O'Brien-Fleming-type bounds (0.025 per side) at t_k = N_k / N_max, cached.

    Spending is capped at t = 1: looks past N_max spend nothing, so their bound is +inf (no
    rejection and an unbounded interval), which is the correct protocol for an overrun.
    """
    t = fine_grid(n_max)[list(cols)] / n_max
    b = lan_demets_bounds(t, spend_obf(0.025))
    b[t > 1.0] = np.inf
    return b


class GstOBF(Method):
    """Lan-DeMets O'Brien-Fleming-type group sequential z-test planned for N_max."""

    def __init__(self) -> None:
        super().__init__("ref.gst_obf", "absolute", "planned", "GST, LD O'Brien-Fleming")

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        if ctx.alpha != 0.05:
            raise ValueError("ref.gst_obf is implemented for alpha = 0.05")
        b = obf_bounds(ctx.n_max, tuple(int(c) for c in ctx.cols))
        d, v, bad = _abs_diff_var(stats)
        with np.errstate(divide="ignore", invalid="ignore"):
            se = np.sqrt(v)
            z = np.abs(d) / se
        return _finish(d, z >= b, b * se, (z - b) / b, bad)
