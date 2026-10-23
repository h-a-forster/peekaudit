"""Recompute the pinned reference values and the Monte Carlo power table in docs/DESIGN.md.

Usage: uv run python scripts/reference_values.py

Prints (1) error rates of naively repeated 5% z-tests, (2) Pocock and O'Brien-Fleming
constants, (3) Lan-DeMets bounds for OBF-type and power-family spending, all from the
independent solver in peekaudit.gst, and (4) the critical false-positive rate and power of
the Holm-adjusted one-sided exact binomial test used to flag violations.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize
from scipy.stats import binom, norm

from peekaudit.gst import crossing_probability, lan_demets_bounds, spend_obf, spend_power

ALPHA = 0.05


def equal(k: int) -> np.ndarray:
    return np.arange(1, k + 1) / k


def naive_rates() -> None:
    z = float(norm.isf(ALPHA / 2))
    print("Naive repeated two-sided 5% z-tests, K equally spaced looks")
    for k in (1, 2, 3, 4, 5, 10, 14, 20, 50):
        print(f"  K={k:<3d} P(any rejection) = {crossing_probability(equal(k), [z] * k).sum():.4f}")


def constants() -> None:
    print("Pocock / O'Brien-Fleming constants, two-sided alpha = 0.05")
    for k in (2, 3, 4, 5):
        t = equal(k)
        c_p = optimize.brentq(
            lambda c, t=t, k=k: crossing_probability(t, [c] * k).sum() - ALPHA, 1.5, 4
        )
        shape = np.sqrt(k / np.arange(1, k + 1))
        c_b = optimize.brentq(
            lambda c, t=t, s=shape: crossing_probability(t, c * s).sum() - ALPHA, 1.5, 4
        )
        print(f"  K={k}: Pocock {c_p:.3f}  OBF {c_b:.3f}")


def lan_demets() -> None:
    print("Lan-DeMets bounds, 0.025 per side")
    print("  OBF-type K=5:     ", np.round(lan_demets_bounds(equal(5), spend_obf(0.025)), 4))
    print("  power rho=2 K=5:  ", np.round(lan_demets_bounds(equal(5), spend_power(0.025, 2)), 4))
    print("  power rho=2 K=10: ", np.round(lan_demets_bounds(equal(10), spend_power(0.025, 2)), 4))


def mc_power() -> None:
    print("Violation test: one-sided exact binomial, H0 FPR <= 0.05, Holm over m cells")
    for reps in (2000, 4000, 20000):
        for m in (1, 96):  # 96 = primary A/A cells of an anytime method
            level = ALPHA / m
            k = int(binom.isf(level, reps, ALPHA))
            while binom.sf(k - 1, reps, ALPHA) > level:
                k += 1
            powers = "  ".join(
                f"{p:.1%}: {binom.sf(k - 1, reps, p):.3f}" for p in (0.06, 0.065, 0.07, 0.08)
            )
            print(f"  R={reps:<6d} m={m:<3d} critical FPR {k / reps:.2%}  power {powers}")


if __name__ == "__main__":
    naive_rates()
    constants()
    lan_demets()
    mc_power()
