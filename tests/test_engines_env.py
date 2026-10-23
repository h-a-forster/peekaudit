"""The audited engines are installed at the pinned versions and behave as documented.

These are environment and adapter-precondition checks, not audit results: they reproduce a
fixture from an engine's own test suite, check a duality the engine's formulas imply, and
cross-check the spotify-confidence bound solver against the independent reference solver in
peekaudit.gst.
"""

from __future__ import annotations

from importlib import metadata

import numpy as np
import pytest

from peekaudit.gst import lan_demets_bounds, spend_power

GBSTATS_VERSION = "0.8.0"
SPOTIFY_VERSION = "4.1.0"


def test_pinned_engine_versions() -> None:
    assert metadata.version("gbstats") == GBSTATS_VERSION
    assert metadata.version("spotify-confidence") == SPOTIFY_VERSION


def test_gbstats_sequential_fixture() -> None:
    """Fixture from gbstats 0.8.0 tests/frequentist/test_tests.py
    (TestSequentialTTest.test_sequential_test_runs_prop), rounded to 5 decimals there."""
    from gbstats.frequentist.tests import SequentialTwoSidedTTest
    from gbstats.models.statistics import ProportionStatistic

    res = SequentialTwoSidedTTest(
        ProportionStatistic(sum=1396, n=3000), ProportionStatistic(sum=2422, n=3461)
    ).compute_result()
    assert round(res.expected, 5) == 0.50386
    assert [round(x, 5) for x in res.ci] == [0.40098, 0.60675]
    assert round(res.uplift.stddev, 5) == 0.03386


def test_gbstats_sequential_ci_and_pvalue_are_dual() -> None:
    """Reject (p < alpha) exactly when the confidence sequence excludes zero."""
    from gbstats.frequentist.tests import SequentialConfig, SequentialTwoSidedTTest
    from gbstats.models.statistics import SampleMeanStatistic

    rng = np.random.default_rng(1)
    checked = 0
    for _ in range(200):
        n_a, n_b = (int(x) for x in rng.integers(20, 5000, size=2))
        xa = rng.lognormal(0.0, 1.0, n_a)
        xb = rng.lognormal(rng.normal(0, 0.2), 1.0, n_b)
        res = SequentialTwoSidedTTest(
            SampleMeanStatistic(n=n_a, sum=xa.sum(), sum_squares=(xa**2).sum()),
            SampleMeanStatistic(n=n_b, sum=xb.sum(), sum_squares=(xb**2).sum()),
            SequentialConfig(difference_type="absolute"),
        ).compute_result()
        excludes_zero = res.ci[0] > 0 or res.ci[1] < 0
        if abs(res.p_value - 0.05) > 1e-9:
            assert (res.p_value < 0.05) == excludes_zero
            checked += 1
    assert checked > 150


@pytest.mark.parametrize("k", [5, 10])
def test_spotify_bounds_match_reference_solver(k: int) -> None:
    from spotify_confidence.analysis.frequentist.confidence_computers.z_test_computer import (
        sequential_bounds,
    )

    t = np.arange(1, k + 1) / k
    engine = sequential_bounds(t, alpha=0.05, sides=2).df["zb"].to_numpy()
    ours = lan_demets_bounds(t, spend_power(0.025, 2.0))
    np.testing.assert_allclose(engine, ours, atol=2e-3)
