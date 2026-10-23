"""Method interface, outputs, evaluation context and the stopping rule.

A method maps cumulative sufficient statistics at K looks for R replicates to a decision
and an interval per (replicate, look). Methods are vectorised (``evaluate``); engine
adapters additionally call the real engine for a few replicates (``evaluate_engine``) so
the vectorised "mirror" can be checked against it (DESIGN 2.1 and 5.2).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np

from peekaudit.scenarios import ALPHA, SCHEDULES, Metric, fine_grid, schedule_within_nmax
from peekaudit.simulate import LookStats

Scale = Literal["absolute", "relative"]
#: What a method promises (DESIGN 6.3): "anytime" = every schedule; "planned" = schedules
#: within N_max; "fixed" = a single look at N_max; "none" = no frequentist claim.
Claim = Literal["anytime", "planned", "fixed", "none"]


@dataclass(frozen=True)
class MethodOutput:
    """Per-(replicate, look) output, arrays of shape (R, K).

    ``reject``: +1 reject with B > A, -1 reject with B < A, 0 no rejection.
    ``lower``/``upper``: the interval the method reports at that look (NaN when none).
    ``no_result``: the method returned its "cannot analyse" output.
    ``margin``: signed relative distance of the test statistic from its decision threshold
    (positive = reject side; NaN when not applicable). Used only to classify exact ties
    when a mirror is compared with its engine.
    """

    reject: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    no_result: np.ndarray
    margin: np.ndarray

    def rows(self, idx: np.ndarray | slice) -> MethodOutput:
        return MethodOutput(
            self.reject[idx], self.lower[idx], self.upper[idx], self.no_result[idx], self.margin[idx]
        )

    @staticmethod
    def empty(shape: tuple[int, int]) -> MethodOutput:
        nan = np.full(shape, np.nan)
        return MethodOutput(
            np.zeros(shape, dtype=np.int8), nan, nan.copy(), np.ones(shape, bool), nan.copy()
        )


@dataclass(frozen=True)
class EvalContext:
    """What a method may know besides the data: the design, never the truth."""

    metric: Metric
    n_max: int
    schedule: str
    cols: np.ndarray
    tau: float
    alpha: float = ALPHA

    @property
    def totals(self) -> np.ndarray:
        """Total sample size N_k at each look of the schedule."""
        return fine_grid(self.n_max)[self.cols]

    @staticmethod
    def for_schedule(metric: Metric, n_max: int, schedule: str, tau: float) -> EvalContext:
        return EvalContext(metric, n_max, schedule, SCHEDULES[schedule], tau)


class Method(ABC):
    """Base class of every audited or reference method."""

    def __init__(
        self,
        name: str,
        scale: Scale,
        claim: Claim,
        label: str,
        engine: str | None = None,
    ) -> None:
        self.name = name
        self.scale: Scale = scale
        self.claim: Claim = claim
        self.label = label
        self.engine = engine

    def supports(self, metric: Metric) -> bool:  # noqa: ARG002 - overridden
        """Whether the method can analyse this metric at all."""
        return True

    def claims(self, schedule: str) -> bool:
        """Whether the method promises FPR <= alpha (and coverage) under this schedule."""
        if self.claim == "anytime":
            return True
        if self.claim == "planned":
            return schedule_within_nmax(schedule)
        if self.claim == "fixed":
            return schedule == "fixed1"
        return False

    @abstractmethod
    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        """Vectorised decisions and intervals for every replicate and look."""

    def evaluate_engine(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        """Call the real engine for every replicate of ``stats`` (slow; adapters only)."""
        raise NotImplementedError(f"{self.name} has no engine")

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.name}>"


Decision = tuple[int, float, float] | None
DecideFn = Callable[[int, float, float, int, float, float], Decision]


class ScalarMethod(Method):
    """Wrap a per-look callback so any engine can be audited without vectorising it.

    ``decide(n_a, sum_a, sumsq_a, n_b, sum_b, sumsq_b)`` returns ``(reject, lower, upper)``
    with reject in {-1, 0, +1}, or ``None`` for "no result". It is called for every
    replicate and look, so it should be fast (it is used by the conformance check).
    Ratio metrics are not supported.
    """

    def __init__(
        self,
        name: str,
        decide: DecideFn,
        scale: Scale = "absolute",
        claim: Claim = "anytime",
        label: str = "",
    ) -> None:
        super().__init__(name, scale, claim, label or name)
        self._decide = decide

    def supports(self, metric: Metric) -> bool:
        return not metric.ratio

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:  # noqa: ARG002
        r, k = stats.n_reps, stats.n_looks
        out = MethodOutput.empty((r, k))
        for i in range(r):
            for j in range(k):
                d = self._decide(
                    int(stats.n_a[j]),
                    float(stats.sum_a[i, j]),
                    float(stats.sumsq_a[i, j]),
                    int(stats.n_b[j]),
                    float(stats.sum_b[i, j]),
                    float(stats.sumsq_b[i, j]),
                )
                if d is None:
                    continue
                rej, lo, hi = d
                if rej not in (-1, 0, 1):
                    raise ValueError(f"{self.name}: reject must be -1, 0 or +1, got {rej}")
                out.reject[i, j] = rej
                out.lower[i, j] = lo
                out.upper[i, j] = hi
                out.no_result[i, j] = False
        return out


def sample_moments(stats: LookStats) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-arm means and per-unit variances on the absolute scale, (R, K) each.

    Means are S/n. Variances use the n - 1 denominator (NaN when n <= 1). For ratio metrics
    the mean is sum(m)/sum(d) and the variance is the delta-method per-unit variance with
    n - 1 (co)variances, the same formula gbstats uses.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        out = []
        for arm in ("a", "b"):
            n = getattr(stats, f"n_{arm}").astype(np.float64)
            s = getattr(stats, f"sum_{arm}")
            q = getattr(stats, f"sumsq_{arm}")
            nm1 = np.where(n > 1, n - 1.0, np.nan)
            var_m = (q - s * s / n) / nm1
            if not stats.ratio:
                out += [s / n, var_m]
                continue
            d = getattr(stats, f"dsum_{arm}")
            dq = getattr(stats, f"dsumsq_{arm}")
            md = getattr(stats, f"mdsum_{arm}")
            var_d = (dq - d * d / n) / nm1
            cov = (md - s * d / n) / nm1
            mm, dm = s / n, d / n
            var = var_m / dm**2 - 2.0 * cov * mm / dm**3 + mm**2 * var_d / dm**4
            out += [s / d, var]
    return out[0], out[1], out[2], out[3]


def stopping_outcome(out: MethodOutput) -> dict[str, np.ndarray]:
    """Apply the stopping rule: stop at the first rejecting look, else at the last look.

    Returns per-replicate arrays: ``stop`` (look index), ``sign`` (0 if never rejected),
    ``lower``/``upper`` and ``no_result`` at the stopping look.
    """
    rej = out.reject != 0
    r, k = rej.shape
    any_rej = rej.any(axis=1)
    stop = np.where(any_rej, rej.argmax(axis=1), k - 1)
    idx = np.arange(r)
    return {
        "stop": stop.astype(np.int16),
        "sign": np.where(any_rej, out.reject[idx, stop], 0).astype(np.int8),
        "lower": out.lower[idx, stop],
        "upper": out.upper[idx, stop],
        "no_result": out.no_result[idx, stop],
    }
