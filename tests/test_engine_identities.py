"""Engine-backed identity checks: real gbstats engines vs independent computations."""

from __future__ import annotations

import numpy as np
import pytest

from peekaudit.simulate import LookStats
from peekaudit.verify import _independent_relative_z, engine_identities

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_independent_relative_z_hand_value() -> None:
    # A: 1000 users, 100 conversions; B: 1000 users, 130 conversions.
    st = LookStats(
        np.array([1000]), np.array([[100.0]]), np.array([[100.0]]),
        np.array([1000]), np.array([[130.0]]), np.array([[130.0]]), True,
    )  # fmt: skip
    z, valid = _independent_relative_z(st)
    var = 0.13 * 0.87 / (0.1**2 * 1000) + 0.1 * 0.9 * 0.13**2 / (0.1**4 * 1000)
    assert valid[0]
    assert z[0] == pytest.approx(0.3 / np.sqrt(var), rel=1e-12)


def test_independent_relative_z_refusals() -> None:
    st = LookStats(
        np.array([10, 10, 1]), np.array([[0.0, 10.0, 3.0]]), np.array([[0.0, 10.0, 9.0]]),
        np.array([10, 10, 10]), np.array([[2.0, 2.0, 2.0]]), np.array([[4.0, 4.0, 4.0]]),
        False,
    )  # fmt: skip
    _, valid = _independent_relative_z(st)
    # zero control mean; zero control variance; single control unit (variance 0)
    assert not valid.any()


def test_engine_identities_agree() -> None:
    res = engine_identities(n=300, seed=3)
    assert set(res) == {"bayes_vs_z10", "acs_vs_gbstats_abs"}
    for agr in res.values():
        assert agr.ok, agr
        assert agr.engine_rejections > 0
        assert agr.cells - agr.engine_no_results > 50
