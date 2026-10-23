"""Inference on Monte Carlo output: exact binomial intervals and tests, Holm, mean intervals.

Every rate the audit reports (false-positive rate, coverage, power) is a binomial proportion
over independent replicates, so it gets an exact Clopper-Pearson interval. Calibration
claims are tested one-sided with an exact binomial test and adjusted with Holm's step-down
procedure across the cells of a family. See docs/DESIGN.md, section "Inference on the
simulation output".
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import stats as st


def clopper_pearson(k: int, n: int, conf: float = 0.95) -> tuple[float, float]:
    """Exact two-sided Clopper-Pearson interval for a binomial proportion k/n.

    lower = Beta(a/2; k, n-k+1) and upper = Beta(1-a/2; k+1, n-k), with a = 1 - conf, and
    the conventions lower = 0 when k = 0 and upper = 1 when k = n.
    """
    _check_counts(k, n)
    if not 0.0 < conf < 1.0:
        raise ValueError(f"conf must be in (0, 1), got {conf}")
    a = 1.0 - conf
    lo = 0.0 if k == 0 else float(st.beta.ppf(a / 2.0, k, n - k + 1))
    hi = 1.0 if k == n else float(st.beta.ppf(1.0 - a / 2.0, k + 1, n - k))
    return lo, hi


def binom_sf_test(k: int, n: int, p0: float) -> float:
    """One-sided exact binomial p-value for H0: p <= p0 against H1: p > p0.

    Returns P(X >= k) for X ~ Binomial(n, p0). Used for "the false-positive rate exceeds
    alpha" and, on miss counts, for "coverage falls below 1 - alpha".
    """
    _check_counts(k, n)
    if not 0.0 <= p0 <= 1.0:
        raise ValueError(f"p0 must be in [0, 1], got {p0}")
    if k == 0:
        return 1.0
    return float(st.binom.sf(k - 1, n, p0))


def holm(pvalues: Sequence[float] | np.ndarray) -> np.ndarray:
    """Holm step-down adjusted p-values (strong family-wise error control).

    With p-values sorted ascending p_(1) <= ... <= p_(m), the adjusted value is
    max_{j <= i} min(1, (m - j + 1) p_(j)), returned in the original order.
    """
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1:
        raise ValueError("pvalues must be one-dimensional")
    if np.any((p < 0) | (p > 1)) or np.any(np.isnan(p)):
        raise ValueError("pvalues must lie in [0, 1]")
    m = p.size
    if m == 0:
        return p.copy()
    order = np.argsort(p, kind="stable")
    scaled = np.minimum(1.0, (m - np.arange(m)) * p[order])
    adjusted_sorted = np.maximum.accumulate(scaled)
    out = np.empty(m)
    out[order] = adjusted_sorted
    return out


def mean_ci(x: Sequence[float] | np.ndarray, conf: float = 0.95) -> tuple[float, float, float]:
    """Sample mean with a two-sided Student-t interval: (mean, lower, upper).

    Used for expected sample size at stopping, which is bounded, so the t-interval is
    adequate at the replicate counts used here (thousands).
    """
    a = np.asarray(x, dtype=float)
    if a.ndim != 1 or a.size < 2:
        raise ValueError("need a one-dimensional sample with at least two values")
    m = float(a.mean())
    se = float(a.std(ddof=1)) / np.sqrt(a.size)
    q = float(st.t.ppf(0.5 + conf / 2.0, a.size - 1))
    return m, m - q * se, m + q * se


def _check_counts(k: int, n: int) -> None:
    if n <= 0:
        raise ValueError(f"n must be positive, got {n}")
    if not 0 <= k <= n:
        raise ValueError(f"k must be in [0, n], got k={k}, n={n}")
