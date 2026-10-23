"""Group sequential boundaries by recursive numerical integration (reference implementation).

Under H0 the standardised statistics Z_k = S(t_k) / sqrt(t_k) at information fractions
0 < t_1 < ... < t_K are those of a standard Brownian motion S sampled at t_k (canonical joint
distribution; Jennison and Turnbull 2000, ch. 3). The probability of a first crossing of the
symmetric two-sided boundary |Z_k| >= b_k at look k is computed with the recursion of
Armitage, McPherson and Rowe (1969): carry the sub-density of S on the continuation region
|S| < b_k sqrt(t_k) from look to look with one Gaussian transition per look (trapezoidal
rule on a uniform grid).

Lan-DeMets boundaries (Lan and DeMets 1983) spend a cumulative error alpha(t_k) by look k;
b_k solves P(first crossing at look k) = alpha(t_k) - alpha(t_{k-1}). Spending functions are
given per side (one-sided alpha); the two-sided total is twice the per-side spend.

This module is independent of every audited engine. It reproduces the Armitage et al.
repeated-significance error rates, the Pocock and O'Brien-Fleming constants of Jennison and
Turnbull (2000, Tables 2.1 and 2.3), the Lan-DeMets O'Brien-Fleming-type bounds published
with the R packages ldbounds and gsDesign, and the power-family (rho = 2) bounds computed by
spotify-confidence 4.1.0 (see tests/test_gst.py).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
from scipy import optimize
from scipy.stats import norm

SpendingFunction = Callable[[float], float]

#: Bounds are capped here when a look is allotted (numerically) no error to spend.
Z_CAP = 12.0


def spend_obf(alpha_side: float) -> SpendingFunction:
    """Lan-DeMets O'Brien-Fleming-type spending per side: 2 - 2 Phi(z_{1-a/2} / sqrt(t)).

    This is the gsDesign ``sfLDOF`` / ldbounds ``iuse=1`` convention with a = alpha_side.
    """
    z = float(norm.isf(alpha_side / 2.0))
    return lambda t: float(2.0 - 2.0 * norm.cdf(z / np.sqrt(t))) if t > 0 else 0.0


def spend_power(alpha_side: float, rho: float) -> SpendingFunction:
    """Kim-DeMets power-family spending per side: a * t**rho (rho = 2 in spotify-confidence)."""
    return lambda t: float(alpha_side * t**rho)


def spend_pocock(alpha_side: float) -> SpendingFunction:
    """Lan-DeMets Pocock-type spending per side: a * log(1 + (e - 1) t)."""
    return lambda t: float(alpha_side * np.log1p((np.e - 1.0) * t))


def _validate_fractions(t: Sequence[float]) -> np.ndarray:
    ta = np.asarray(t, dtype=float)
    if ta.ndim != 1 or ta.size == 0:
        raise ValueError("t must be a non-empty one-dimensional sequence")
    if ta[0] <= 0 or np.any(np.diff(ta) <= 0):
        raise ValueError("information fractions must be positive and strictly increasing")
    return ta


class _Recursion:
    """Sub-density of S on the continuation region, advanced one look at a time."""

    def __init__(self, ngrid: int) -> None:
        if ngrid < 101 or ngrid % 2 == 0:
            raise ValueError("ngrid must be an odd integer >= 101")
        self.ngrid = ngrid
        self.t_prev = 0.0
        self.grid: np.ndarray | None = None  # S values; None means S(0) = 0 exactly
        self.weights: np.ndarray | None = None  # density * trapezoid weights

    def crossing(self, t: float, b: float) -> float:
        """P(continue to the previous look, then |Z| >= b at information t)."""
        dt = t - self.t_prev
        sd = np.sqrt(dt)
        s_hi = b * np.sqrt(t)
        if self.grid is None:
            return float(2.0 * norm.sf(s_hi / sd))
        up = norm.sf((s_hi - self.grid) / sd)
        lo = norm.cdf((-s_hi - self.grid) / sd)
        return float(np.sum(self.weights * (up + lo)))

    def advance(self, t: float, b: float) -> None:
        """Condition on not crossing |Z| >= b at information t and move to that look."""
        dt = t - self.t_prev
        sd = np.sqrt(dt)
        s_hi = min(b, Z_CAP) * np.sqrt(t)
        new = np.linspace(-s_hi, s_hi, self.ngrid)
        if self.grid is None:
            dens = norm.pdf(new, 0.0, sd)
        else:
            kernel = norm.pdf((new[:, None] - self.grid[None, :]) / sd) / sd
            dens = kernel @ self.weights
        h = new[1] - new[0]
        w = dens * h
        w[0] *= 0.5
        w[-1] *= 0.5
        self.grid, self.weights, self.t_prev = new, w, t


def crossing_probability(
    t: Sequence[float], bounds: Sequence[float], ngrid: int = 1001
) -> np.ndarray:
    """Per-look first-crossing probabilities of |Z_k| >= b_k under H0.

    The sum is the overall type I error of the two-sided design. With b_k = z_{1-a/2} at
    every look it gives the error rate of naively repeated significance tests.
    """
    ta = _validate_fractions(t)
    ba = np.asarray(bounds, dtype=float)
    if ba.shape != ta.shape:
        raise ValueError("t and bounds must have the same length")
    rec = _Recursion(ngrid)
    out = np.empty(ta.size)
    for k, (tk, bk) in enumerate(zip(ta, ba, strict=True)):
        out[k] = rec.crossing(tk, bk)
        rec.advance(tk, bk)
    return out


def lan_demets_bounds(
    t: Sequence[float], spending: SpendingFunction, ngrid: int = 1001
) -> np.ndarray:
    """Two-sided symmetric Lan-DeMets bounds b_k for per-side spending function ``spending``.

    ``t`` are information fractions (the final planned look has t = 1; looks past t = 1 are
    allowed and spend nothing beyond alpha(1) if the spending function is capped by the
    caller). Each b_k solves 2 * [alpha(t_k) - alpha(t_{k-1})] = P(first crossing at look k).
    """
    ta = _validate_fractions(t)
    rec = _Recursion(ngrid)
    bounds = np.empty(ta.size)
    spent_prev = 0.0
    for k, tk in enumerate(ta):
        spent = spending(min(tk, 1.0))
        target = 2.0 * (spent - spent_prev)
        spent_prev = spent
        if target <= 1e-15 or rec.crossing(tk, Z_CAP) >= target:
            bk = Z_CAP
        else:
            bk = optimize.brentq(
                lambda b, tk=tk, target=target: rec.crossing(tk, b) - target, 0.0, Z_CAP, xtol=1e-10
            )
        bounds[k] = bk
        rec.advance(tk, bk)
    return bounds
