"""Reference group sequential solver against published values.

Sources of the reference numbers:
- Armitage, McPherson and Rowe (1969): repeated significance tests at nominal two-sided 5%
  with K equally spaced looks (reproduced in Jennison and Turnbull 2000, Table 1.1).
- Pocock and O'Brien-Fleming constants, two-sided alpha = 0.05 (Jennison and Turnbull 2000,
  Tables 2.1 and 2.3).
- Lan-DeMets O'Brien-Fleming-type bounds, K = 5 equally spaced, one-sided 0.025 per side
  (R packages ldbounds and gsDesign).
- Power-family rho = 2 bounds as computed by spotify-confidence 4.1.0 (also cross-checked
  directly against the engine in tests/test_engines_env.py).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from scipy import optimize
from scipy.stats import norm

from peekaudit.gst import (
    crossing_probability,
    lan_demets_bounds,
    spend_obf,
    spend_pocock,
    spend_power,
)

Z975 = float(norm.isf(0.025))


def _equal(k: int) -> np.ndarray:
    return np.arange(1, k + 1) / k


@pytest.mark.parametrize(
    ("k", "expected"),
    [(1, 0.05), (2, 0.083), (3, 0.107), (4, 0.126), (5, 0.142), (10, 0.193), (20, 0.248)],
)
def test_naive_repeated_testing_error(k: int, expected: float) -> None:
    total = crossing_probability(_equal(k), [Z975] * k).sum()
    assert total == pytest.approx(expected, abs=1e-3)


@pytest.mark.parametrize(
    ("k", "pocock", "obf"),
    [(2, 2.178, 1.977), (3, 2.289, 2.004), (4, 2.361, 2.024), (5, 2.413, 2.040)],
)
def test_pocock_and_obrien_fleming_constants(k: int, pocock: float, obf: float) -> None:
    t = _equal(k)
    c_p = optimize.brentq(lambda c: crossing_probability(t, [c] * k).sum() - 0.05, 1.5, 4)
    shape = np.sqrt(k / np.arange(1, k + 1))
    c_b = optimize.brentq(lambda c: crossing_probability(t, c * shape).sum() - 0.05, 1.5, 4)
    assert c_p == pytest.approx(pocock, abs=1.5e-3)
    assert c_b == pytest.approx(obf, abs=1.5e-3)


def test_lan_demets_obf_matches_ldbounds() -> None:
    b = lan_demets_bounds(_equal(5), spend_obf(0.025))
    np.testing.assert_allclose(b, [4.8769, 3.3570, 2.6803, 2.2898, 2.0310], atol=2e-3)


def test_power_family_matches_spotify_values() -> None:
    b = lan_demets_bounds(_equal(5), spend_power(0.025, 2.0))
    np.testing.assert_allclose(b, [3.0902, 2.7141, 2.4728, 2.2799, 2.1140], atol=2e-3)


@pytest.mark.parametrize(
    "spending", [spend_obf(0.025), spend_power(0.025, 2.0), spend_pocock(0.025)]
)
def test_bounds_spend_exactly_alpha(spending) -> None:
    t = np.array([0.15, 0.3, 0.55, 0.8, 1.0])  # unequal spacing
    b = lan_demets_bounds(t, spending)
    probs = crossing_probability(t, b)
    np.testing.assert_allclose(np.cumsum(probs), [2 * spending(x) for x in t], atol=2e-5)


def test_many_looks_is_fast_and_spends_alpha() -> None:
    t = _equal(140)
    start = time.perf_counter()
    b = lan_demets_bounds(t, spend_power(0.025, 2.0), ngrid=501)
    assert time.perf_counter() - start < 30
    assert crossing_probability(t, b, ngrid=501).sum() == pytest.approx(0.05, abs=5e-4)


def test_input_validation() -> None:
    with pytest.raises(ValueError):
        crossing_probability([0.5, 0.4], [2.0, 2.0])
    with pytest.raises(ValueError):
        crossing_probability([0.5, 1.0], [2.0])
    with pytest.raises(ValueError):
        lan_demets_bounds([0.0, 1.0], spend_obf(0.025))


@pytest.mark.slow
def test_monte_carlo_type_i_error_of_obf_design() -> None:
    """Simulated Brownian motion crosses the LD-OBF K=5 bounds at rate 0.05 (exact CI)."""
    from peekaudit.stats import clopper_pearson

    rng = np.random.default_rng(20261023)
    t = _equal(5)
    b = lan_demets_bounds(t, spend_obf(0.025))
    reps = 400_000
    steps = rng.standard_normal((reps, 5)) * np.sqrt(np.diff(t, prepend=0.0))
    z = np.cumsum(steps, axis=1) / np.sqrt(t)
    hits = int(np.any(np.abs(z) >= b, axis=1).sum())
    lo, hi = clopper_pearson(hits, reps, conf=0.999)
    assert lo <= 0.05 <= hi
