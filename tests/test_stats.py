"""Exact-inference helpers against closed forms and statsmodels."""

from __future__ import annotations

import numpy as np
import pytest

from peekaudit.stats import binom_sf_test, clopper_pearson, holm, mean_ci


def test_clopper_pearson_closed_forms() -> None:
    # k = 0: upper = 1 - (a/2)^(1/n); k = n: lower = (a/2)^(1/n).
    lo, hi = clopper_pearson(0, 10)
    assert lo == 0.0
    assert hi == pytest.approx(1 - 0.025 ** (1 / 10), abs=1e-12)
    lo, hi = clopper_pearson(10, 10)
    assert lo == pytest.approx(0.025 ** (1 / 10), abs=1e-12)
    assert hi == 1.0


def test_clopper_pearson_tabulated_value() -> None:
    # 5/10 -> (0.18709, 0.81291), the standard textbook value.
    lo, hi = clopper_pearson(5, 10)
    assert lo == pytest.approx(0.187086, abs=1e-6)
    assert hi == pytest.approx(0.812914, abs=1e-6)


def test_clopper_pearson_matches_statsmodels() -> None:
    from statsmodels.stats.proportion import proportion_confint

    for k, n in [(1, 50), (200, 4000), (37, 1000), (999, 1000)]:
        ref = proportion_confint(k, n, alpha=0.05, method="beta")
        assert clopper_pearson(k, n) == pytest.approx(ref, abs=1e-10)


@pytest.mark.parametrize("bad", [(-1, 10), (11, 10), (0, 0)])
def test_counts_validated(bad: tuple[int, int]) -> None:
    with pytest.raises(ValueError):
        clopper_pearson(*bad)


def test_binom_sf_test() -> None:
    assert binom_sf_test(0, 100, 0.05) == 1.0
    # P(X >= 1) = 1 - 0.95^100
    assert binom_sf_test(1, 100, 0.05) == pytest.approx(1 - 0.95**100, rel=1e-12)
    # Far above the null: tiny p-value; at the null mean: about one half.
    assert binom_sf_test(300, 4000, 0.05) < 1e-9
    assert 0.4 < binom_sf_test(200, 4000, 0.05) < 0.6


def test_holm_matches_statsmodels() -> None:
    from statsmodels.stats.multitest import multipletests

    rng = np.random.default_rng(0)
    p = np.concatenate([rng.uniform(0, 0.01, 5), rng.uniform(0, 1, 20)])
    ref = multipletests(p, method="holm")[1]
    np.testing.assert_allclose(holm(p), ref, atol=1e-12)


def test_holm_small_example() -> None:
    # Sorted p = 0.01, 0.02, 0.03 with m = 3 -> 0.03, 0.04, 0.04 (monotone step-down).
    np.testing.assert_allclose(holm([0.03, 0.01, 0.02]), [0.04, 0.03, 0.04])
    assert holm([]).size == 0
    with pytest.raises(ValueError):
        holm([0.5, 1.5])


def test_mean_ci() -> None:
    m, lo, hi = mean_ci([1.0, 2.0, 3.0, 4.0])
    assert m == 2.5
    # t_{0.975, 3} = 3.182446, s = 1.290994, se = 0.645497
    assert hi - m == pytest.approx(3.182446 * 0.645497, rel=1e-5)
    assert m - lo == pytest.approx(hi - m)
