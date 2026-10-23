"""gbstats 0.8.0 adapter: vectorised mirrors of the engine formulas plus real engine calls.

Mirrors follow the engine source (``gbstats/frequentist/tests.py``,
``gbstats/bayesian/tests.py``, ``gbstats/models/statistics.py``) operation by operation,
including the "default output" branches (control mean 0, a variance <= 0). Engine calls
build the engine's own statistic objects per (replicate, look) and call ``compute_result``.
The engine is imported lazily so the core package does not depend on it. DESIGN 4.1.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.stats import norm
from scipy.stats import t as student_t

from peekaudit.methods.base import EvalContext, Method, MethodOutput
from peekaudit.methods.reference import acs_pvalue_halfwidth
from peekaudit.scenarios import Metric
from peekaudit.simulate import LookStats

#: GrowthBook's default sequential tuning parameter (expected total sample size N*).
DEFAULT_TUNING = 5000.0
#: GrowthBook's default "chance to win" decision threshold.
CTW_THRESHOLD = 0.95


def gbstats_arm(stats: LookStats, arm: str) -> tuple[np.ndarray, np.ndarray]:
    """Engine statistic mean and variance for one arm, (R, K) each.

    Bernoulli -> ``ProportionStatistic`` (variance p(1-p)); continuous ->
    ``SampleMeanStatistic`` (n - 1 variance, 0 when n <= 1); ratio -> ``RatioStatistic``
    (delta method, 0 when the denominator mean is 0 or n <= 1). A mean of an empty arm is 0.
    """
    n_int = getattr(stats, f"n_{arm}")
    n = np.broadcast_to(n_int.astype(np.float64), stats.sum_a.shape)
    s = getattr(stats, f"sum_{arm}")
    q = getattr(stats, f"sumsq_{arm}")
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_s = np.where(n == 0, 0.0, s / n)

        def sm_var(sum_: np.ndarray, sumsq: np.ndarray) -> np.ndarray:
            return np.where(n <= 1, 0.0, (sumsq - np.power(sum_, 2) / n) / (n - 1))

        if stats.ratio:
            d = getattr(stats, f"dsum_{arm}")
            dq = getattr(stats, f"dsumsq_{arm}")
            md = getattr(stats, f"mdsum_{arm}")
            mean = np.where(d == 0, 0.0, s / d)
            d_mean = np.where(n == 0, 0.0, d / n)
            m_var, d_var = sm_var(s, q), sm_var(d, dq)
            cov = np.where(n <= 1, 0.0, (md - s * d / n) / (n - 1))
            var = (
                m_var / np.power(d_mean, 2)
                - 2 * cov * mean_s / np.power(d_mean, 3)
                + np.power(mean_s, 2) * d_var / np.power(d_mean, 4)
            )
            var = np.where((d_mean == 0) | (n <= 1), 0.0, var)
            return mean, var
        if stats.binary:
            return mean_s, mean_s * (1 - mean_s)
        return mean_s, sm_var(s, q)


def _diff_var(
    stats: LookStats, relative: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
    """Point estimate, its variance, the default-output mask and the arm moments."""
    ma, va = gbstats_arm(stats, "a")
    mb, vb = gbstats_arm(stats, "b")
    na = stats.n_a.astype(np.float64)
    nb = stats.n_b.astype(np.float64)
    default = (ma == 0) | (va <= 0) | (vb <= 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        if relative:
            d = (mb - ma) / ma
            v = vb / (np.power(ma, 2) * nb) + va * np.power(mb, 2) / (np.power(ma, 4) * na)
        else:
            d = mb - ma
            v = vb / nb + va / na
    return d, v, default, (ma, va, mb, vb, na, nb)


def _output(
    d: np.ndarray,
    reject: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    margin: np.ndarray,
    default: np.ndarray,
) -> MethodOutput:
    # The engine's default output is p = 1 (chance to win 0.5) and CI [0, 0].
    sign = np.where(reject & ~default, np.sign(d), 0).astype(np.int8)
    return MethodOutput(
        sign,
        np.where(default, 0.0, lo),
        np.where(default, 0.0, hi),
        default.copy(),
        np.where(default, np.nan, margin),
    )


def _engine_stat(stats: LookStats, arm: str, i: int, j: int) -> Any:
    from gbstats.models.statistics import (
        ProportionStatistic,
        RatioStatistic,
        SampleMeanStatistic,
    )

    n = int(getattr(stats, f"n_{arm}")[j])
    s = float(getattr(stats, f"sum_{arm}")[i, j])
    q = float(getattr(stats, f"sumsq_{arm}")[i, j])
    if stats.ratio:
        return RatioStatistic(
            n=n,
            m_statistic=SampleMeanStatistic(n=n, sum=s, sum_squares=q),
            d_statistic=SampleMeanStatistic(
                n=n,
                sum=float(getattr(stats, f"dsum_{arm}")[i, j]),
                sum_squares=float(getattr(stats, f"dsumsq_{arm}")[i, j]),
            ),
            m_d_sum_of_products=float(getattr(stats, f"mdsum_{arm}")[i, j]),
        )
    if stats.binary:
        return ProportionStatistic(n=n, sum=s)
    return SampleMeanStatistic(n=n, sum=s, sum_squares=q)


class _GbstatsMethod(Method):
    def supports(self, metric: Metric) -> bool:  # noqa: ARG002
        return True

    def _engine_result(self, stat_a: Any, stat_b: Any, ctx: EvalContext) -> Any:
        raise NotImplementedError

    def _engine_decision(self, res: Any, ctx: EvalContext) -> int:
        raise NotImplementedError

    def evaluate_engine(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        r, k = stats.n_reps, stats.n_looks
        out = MethodOutput.empty((r, k))
        for i in range(r):
            for j in range(k):
                res = self._engine_result(
                    _engine_stat(stats, "a", i, j), _engine_stat(stats, "b", i, j), ctx
                )
                out.no_result[i, j] = res.error_message is not None
                out.lower[i, j] = float(res.ci[0])
                out.upper[i, j] = float(res.ci[1])
                out.reject[i, j] = self._engine_decision(res, ctx)
        return out


class GbstatsSequential(_GbstatsMethod):
    """``SequentialTwoSidedTTest``: asymptotic confidence sequence; reject iff p < alpha."""

    def __init__(self, relative: bool = True, tuned: bool = False) -> None:
        name = "gbstats.seq." + ("rel" if relative else "abs") + (".tuned" if tuned else "")
        label = f"gbstats sequential, {'relative' if relative else 'absolute'}, " + (
            "N* = N_max" if tuned else "N* = 5000"
        )
        super().__init__(name, "relative" if relative else "absolute", "anytime", label, "gbstats")
        self.relative = relative
        self.tuned = tuned

    def tuning(self, ctx: EvalContext) -> float:
        return float(ctx.n_max) if self.tuned else DEFAULT_TUNING

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        d, v, default, _ = _diff_var(stats, self.relative)
        n_total = stats.n_a + stats.n_b
        p, width = acs_pvalue_halfwidth(d, v, n_total, self.tuning(ctx), ctx.alpha)
        return _output(d, p < ctx.alpha, d - width, d + width, (ctx.alpha - p) / ctx.alpha, default)

    def _engine_result(self, stat_a: Any, stat_b: Any, ctx: EvalContext) -> Any:
        from gbstats.frequentist.tests import SequentialConfig, SequentialTwoSidedTTest

        cfg = SequentialConfig(
            difference_type="relative" if self.relative else "absolute",
            alpha=ctx.alpha,
            sequential_tuning_parameter=self.tuning(ctx),
        )
        return SequentialTwoSidedTTest(stat_a, stat_b, cfg).compute_result()

    def _engine_decision(self, res: Any, ctx: EvalContext) -> int:
        if res.error_message is None and res.p_value < ctx.alpha:
            return int(np.sign(res.expected))
        return 0


class GbstatsTTest(_GbstatsMethod):
    """``TwoSidedTTest`` (Welch df, relative), a fixed-horizon test re-run at every look."""

    def __init__(self) -> None:
        super().__init__(
            "gbstats.ttest.rel", "relative", "fixed", "gbstats t-test (peeked)", "gbstats"
        )

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        d, v, default, (_, va, _, vb, na, nb) = _diff_var(stats, relative=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            dof = np.power(vb / nb + va / na, 2) / (
                np.power(vb, 2) / (np.power(nb, 2) * (nb - 1))
                + np.power(va, 2) / (np.power(na, 2) * (na - 1))
            )
            crit = d / np.sqrt(v)
            dof_ok = np.where(default, 10.0, dof)
            p = 2 * (1 - student_t.cdf(np.abs(crit), dof_ok))
            width = student_t.ppf(1 - ctx.alpha / 2, dof_ok) * np.sqrt(v)
        return _output(d, p < ctx.alpha, d - width, d + width, (ctx.alpha - p) / ctx.alpha, default)

    def _engine_result(self, stat_a: Any, stat_b: Any, ctx: EvalContext) -> Any:
        from gbstats.frequentist.tests import FrequentistConfig, TwoSidedTTest

        cfg = FrequentistConfig(difference_type="relative", alpha=ctx.alpha)
        return TwoSidedTTest(stat_a, stat_b, cfg).compute_result()

    def _engine_decision(self, res: Any, ctx: EvalContext) -> int:
        if res.error_message is None and res.p_value < ctx.alpha:
            return int(np.sign(res.expected))
        return 0


class GbstatsBayes(_GbstatsMethod):
    """``EffectBayesianABTest`` with the default improper prior; ship when chance to win
    >= 95% (or <= 5% for a loss). Equivalent to a two-sided z-test at level 0.10."""

    def __init__(self) -> None:
        super().__init__(
            "gbstats.bayes.ctw95",
            "relative",
            "none",
            "gbstats Bayesian, chance to win 95%",
            "gbstats",
        )

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        d, v, default, _ = _diff_var(stats, relative=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            post_prec = 1 / v + 0 / 1
            std = np.sqrt(1 / post_prec)
            std_ok = np.where(default, 1.0, std)
            ctw = norm.sf(0, d, std_ok)
            lo = norm.ppf(ctx.alpha / 2, d, std_ok)
            hi = norm.ppf(1 - ctx.alpha / 2, d, std_ok)
        win, lose = ctw >= CTW_THRESHOLD, ctw <= 1 - CTW_THRESHOLD
        half = CTW_THRESHOLD - 0.5
        margin = (np.abs(ctw - 0.5) - half) / half
        out = _output(d, win | lose, lo, hi, margin, default)
        return out

    def _engine_result(self, stat_a: Any, stat_b: Any, ctx: EvalContext) -> Any:
        from gbstats.bayesian.tests import EffectBayesianABTest, EffectBayesianConfig

        cfg = EffectBayesianConfig(difference_type="relative", alpha=ctx.alpha)
        return EffectBayesianABTest(stat_a, stat_b, cfg).compute_result()

    def _engine_decision(self, res: Any, ctx: EvalContext) -> int:  # noqa: ARG002
        if res.error_message is not None:
            return 0
        if res.chance_to_win >= CTW_THRESHOLD:
            return 1
        if res.chance_to_win <= 1 - CTW_THRESHOLD:
            return -1
        return 0
