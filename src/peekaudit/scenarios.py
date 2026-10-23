"""Study scenarios: metrics and their true moments, effect calibration, look schedules, units.

Everything here is deterministic and cheap. ``peekaudit.simulate`` turns a :class:`Unit` into
data. See docs/DESIGN.md, sections 3.1-3.4.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import optimize
from scipy.stats import norm

Kind = Literal["normal", "bernoulli", "lognormal", "pareto", "ratio"]
Hypothesis = Literal["aa", "ab"]

#: Two-sided level of every audited test.
ALPHA = 0.05
#: Power at N_max of a fixed-horizon two-sided z-test used to calibrate the A/B lift.
DESIGN_POWER = 0.80
#: Fine look grid: N_k = round(k * N_max / GRID_PER_NMAX), k = 1..GRID_LOOKS (to 1.5 N_max).
GRID_PER_NMAX = 140
GRID_LOOKS = 210
#: Master seed of the study (DESIGN 3.4).
MASTER_SEED = 20261023


@dataclass(frozen=True)
class Metric:
    """A per-unit outcome distribution with its true mean and variance (control arm)."""

    id: str
    kind: Kind
    label: str
    mean: float
    variance: float  # math.inf when infinite
    binary: bool = False
    ratio: bool = False
    p: float = 0.0  # Bernoulli rate
    mu: float = 0.0  # normal mean / lognormal log-mean
    sigma: float = 0.0  # normal sd / lognormal log-sd
    tail: float = 0.0  # Pareto tail index a (x_m = 1)
    sessions_lambda: float = 0.0  # ratio: sessions = 1 + Poisson(lambda)
    click_q: float = 0.0  # ratio: clicks | sessions ~ Binomial(sessions, q)

    @property
    def finite_variance(self) -> bool:
        return math.isfinite(self.variance)

    @property
    def ab_allowed(self) -> bool:
        """A/B units need a finite variance to calibrate the lift."""
        return self.finite_variance

    def variance_b(self, lift: float) -> float:
        """True per-unit variance of the treatment arm under multiplicative lift ``lift``.

        For ratio metrics this is the linearised (delta-method) per-unit variance
        q_B (1 - q_B) / E[sessions].
        """
        if self.kind == "bernoulli":
            pb = self.p * (1.0 + lift)
            return pb * (1.0 - pb)
        if self.kind == "ratio":
            qb = self.click_q * (1.0 + lift)
            return qb * (1.0 - qb) / (1.0 + self.sessions_lambda)
        return (1.0 + lift) ** 2 * self.variance

    def max_lift(self) -> float:
        """Largest admissible lift (rates must stay below 1)."""
        if self.kind == "bernoulli":
            return 1.0 / self.p - 1.0
        if self.kind == "ratio":
            return 1.0 / self.click_q - 1.0
        return math.inf


def _lognormal(mu: float, sigma: float) -> tuple[float, float]:
    m = math.exp(mu + sigma**2 / 2.0)
    v = (math.exp(sigma**2) - 1.0) * math.exp(2.0 * mu + sigma**2)
    return m, v


def _pareto(a: float) -> tuple[float, float]:
    m = a / (a - 1.0)
    v = a / ((a - 1.0) ** 2 * (a - 2.0)) if a > 2.0 else math.inf
    return m, v


_LN_M, _LN_V = _lognormal(0.0, 1.5)
_P25_M, _P25_V = _pareto(2.5)
_P15_M, _P15_V = _pareto(1.5)
_RATIO_Q, _RATIO_LAMBDA = 0.2, 1.5

METRICS: dict[str, Metric] = {
    m.id: m
    for m in (
        Metric("normal", "normal", "Normal(1, 1)", 1.0, 1.0, mu=1.0, sigma=1.0),
        Metric("bern_p10", "bernoulli", "conversion 10%", 0.10, 0.09, binary=True, p=0.10),
        Metric("bern_p01", "bernoulli", "conversion 1%", 0.01, 0.0099, binary=True, p=0.01),
        Metric("lognorm", "lognormal", "LogNormal(0, 1.5)", _LN_M, _LN_V, mu=0.0, sigma=1.5),
        Metric("pareto_a25", "pareto", "Pareto tail 2.5", _P25_M, _P25_V, tail=2.5),
        Metric("pareto_a15", "pareto", "Pareto tail 1.5 (stress)", _P15_M, _P15_V, tail=1.5),
        Metric(
            "ratio",
            "ratio",
            "ratio clicks/sessions",
            _RATIO_Q,
            _RATIO_Q * (1.0 - _RATIO_Q) / (1.0 + _RATIO_LAMBDA),
            ratio=True,
            sessions_lambda=_RATIO_LAMBDA,
            click_q=_RATIO_Q,
        ),
    )
}

N_MAX_VALUES = (1_000, 20_000)
ALLOCATIONS = (0.5, 0.2)


def arm_counts(n_total: np.ndarray | int, pi: float) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic arm sizes: n_B = round_half_up(pi * N), n_A = N - n_B."""
    n = np.asarray(n_total, dtype=np.int64)
    n_b = np.floor(pi * n + 0.5).astype(np.int64)
    return n - n_b, n_b


def fine_grid(n_max: int) -> np.ndarray:
    """Total sample size at each fine-grid look k = 1..210: round_half_up(k N_max / 140)."""
    k = np.arange(1, GRID_LOOKS + 1, dtype=np.int64)
    return np.floor(k * n_max / GRID_PER_NMAX + 0.5).astype(np.int64)


#: Schedules as 0-based columns of the fine grid (DESIGN 3.3).
SCHEDULES: dict[str, np.ndarray] = {
    "fixed1": np.array([139]),
    "daily14": np.arange(9, 140, 10),
    "stream140": np.arange(0, 140),
    "overrun21": np.arange(9, 210, 10),
}


def schedule_within_nmax(schedule: str) -> bool:
    """True when the schedule never looks past the planned maximum N_max."""
    return int(SCHEDULES[schedule][-1]) <= GRID_PER_NMAX - 1


def calibrate_lift(
    metric: Metric, n_max: int, pi: float, alpha: float = ALPHA, power: float = DESIGN_POWER
) -> float:
    """Multiplicative lift r giving a fixed-horizon two-sided z-test power ``power`` at N_max.

    Solves |r * mean| = (z_{1-alpha/2} + z_{power}) sqrt(var_A / n_A + var_B(r) / n_B) with the
    true variances (DESIGN 3.2). Raises ValueError for infinite-variance metrics.
    """
    if not metric.finite_variance:
        raise ValueError(f"metric {metric.id} has infinite variance; no lift can be calibrated")
    n_a_arr, n_b_arr = arm_counts(n_max, pi)
    n_a, n_b = float(n_a_arr), float(n_b_arr)
    zsum = float(norm.isf(alpha / 2.0) + norm.ppf(power))

    def f(r: float) -> float:
        se = math.sqrt(metric.variance / n_a + metric.variance_b(r) / n_b)
        return r * metric.mean - zsum * se

    hi = min(metric.max_lift() * (1.0 - 1e-9), 1e6)
    if f(hi) < 0:
        raise ValueError(f"no admissible lift reaches power {power} for {metric.id}, n={n_max}")
    return float(optimize.brentq(f, 1e-12, hi, xtol=1e-14, rtol=1e-12))


@dataclass(frozen=True)
class Unit:
    """One work unit: (metric, N_max, allocation to B, hypothesis)."""

    metric_id: str
    n_max: int
    pi: float
    hypothesis: Hypothesis

    @property
    def metric(self) -> Metric:
        return METRICS[self.metric_id]

    @property
    def id(self) -> str:
        return f"{self.metric_id}-n{self.n_max}-b{round(self.pi * 100):02d}-{self.hypothesis}"

    @property
    def cell(self) -> str:
        """The (metric, N_max, allocation) cell shared by the A/A and A/B units."""
        return f"{self.metric_id}-n{self.n_max}-b{round(self.pi * 100):02d}"

    @property
    def design_lift(self) -> float:
        """Calibrated lift of the cell (used as the A/B effect and the mSPRT mixing scale).

        Infinite-variance metrics borrow the relative lift of ``pareto_a25`` for the same
        N_max and allocation; they have A/A units only, so this only sets the mSPRT scale.
        """
        m = self.metric
        if m.finite_variance:
            return calibrate_lift(m, self.n_max, self.pi)
        return calibrate_lift(METRICS["pareto_a25"], self.n_max, self.pi)

    @property
    def lift(self) -> float:
        """True multiplicative lift of B over A (0 for A/A)."""
        return 0.0 if self.hypothesis == "aa" else self.design_lift

    @property
    def tau(self) -> float:
        """mSPRT mixing standard deviation: the absolute design effect r * mean."""
        return self.design_lift * self.metric.mean

    def true_effect(self, scale: Literal["absolute", "relative"]) -> float:
        """True effect of B vs A on the given scale (multiplicative lift model)."""
        if scale == "relative":
            return self.lift
        return self.lift * self.metric.mean

    def seed_key(self) -> int:
        """First 8 bytes of SHA-256 of the unit id, as an unsigned integer."""
        return int.from_bytes(hashlib.sha256(self.id.encode()).digest()[:8], "big")


def study_units() -> list[Unit]:
    """The 52 units of the study grid (no A/B unit for infinite-variance metrics)."""
    units = []
    for metric_id, metric in METRICS.items():
        for n_max in N_MAX_VALUES:
            for pi in ALLOCATIONS:
                units.append(Unit(metric_id, n_max, pi, "aa"))
                if metric.ab_allowed:
                    units.append(Unit(metric_id, n_max, pi, "ab"))
    return units


def unit_from_id(unit_id: str) -> Unit:
    """Parse a unit id produced by :attr:`Unit.id`."""
    metric_id, n_part, b_part, hyp = unit_id.rsplit("-", 3)
    if metric_id not in METRICS or hyp not in ("aa", "ab"):
        raise ValueError(f"not a unit id: {unit_id}")
    return Unit(metric_id, int(n_part[1:]), int(b_part[1:]) / 100.0, hyp)  # type: ignore[arg-type]
