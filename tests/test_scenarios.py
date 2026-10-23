"""Scenario definitions: true moments, calibration, grids, schedules and units."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import norm

from peekaudit.scenarios import (
    METRICS,
    SCHEDULES,
    Unit,
    arm_counts,
    calibrate_lift,
    fine_grid,
    schedule_within_nmax,
    study_units,
    unit_from_id,
)


def test_closed_form_moments() -> None:
    ln = METRICS["lognorm"]
    assert ln.mean == pytest.approx(3.080216848918031, rel=1e-12)
    assert ln.variance == pytest.approx((math.e**2.25 - 1) * math.e**2.25, rel=1e-12)  # 80.53
    p25 = METRICS["pareto_a25"]
    assert p25.mean == pytest.approx(5 / 3)
    assert p25.variance == pytest.approx(20 / 9)
    assert not METRICS["pareto_a15"].finite_variance
    assert METRICS["pareto_a15"].mean == pytest.approx(3.0)
    assert METRICS["ratio"].variance == pytest.approx(0.064)
    assert METRICS["bern_p01"].variance == pytest.approx(0.0099)


def test_arm_counts_round_half_up() -> None:
    n_a, n_b = arm_counts(np.array([7, 5, 1000, 3]), 0.5)
    assert n_b.tolist() == [4, 3, 500, 2]  # 3.5 -> 4, 2.5 -> 3, 1.5 -> 2 (not banker's)
    assert (n_a + n_b).tolist() == [7, 5, 1000, 3]
    n_a, n_b = arm_counts(np.array([7]), 0.2)
    assert (n_a[0], n_b[0]) == (6, 1)


def test_fine_grid_and_schedules() -> None:
    g = fine_grid(1000)
    assert g.size == 210 and g[0] == 7 and g[139] == 1000 and g[209] == 1500
    assert np.all(np.diff(g) > 0)
    assert np.all(np.diff(fine_grid(20000)) > 0)
    for n_max in (1000, 20000):
        for pi in (0.5, 0.2):
            n_a, n_b = arm_counts(fine_grid(n_max), pi)
            assert np.all(np.diff(n_a) > 0) and np.all(np.diff(n_b) > 0) and n_b[0] >= 1
    assert [len(v) for v in SCHEDULES.values()] == [1, 14, 140, 21]
    assert SCHEDULES["daily14"][-1] == 139 and SCHEDULES["overrun21"][-1] == 209
    assert schedule_within_nmax("stream140") and not schedule_within_nmax("overrun21")


@pytest.mark.parametrize("metric_id", [m for m in METRICS if METRICS[m].finite_variance])
@pytest.mark.parametrize(("n_max", "pi"), [(1000, 0.5), (1000, 0.2), (20000, 0.5), (20000, 0.2)])
def test_calibrated_lift_gives_80_percent_power(metric_id: str, n_max: int, pi: float) -> None:
    m = METRICS[metric_id]
    r = calibrate_lift(m, n_max, pi)
    assert r > 0
    n_a, n_b = (float(x) for x in arm_counts(n_max, pi))
    se = math.sqrt(m.variance / n_a + m.variance_b(r) / n_b)
    shift = r * m.mean / se
    z = norm.isf(0.025)
    power = norm.sf(z - shift) + norm.cdf(-z - shift)
    assert shift == pytest.approx(z + norm.ppf(0.8), rel=1e-10)
    assert power == pytest.approx(0.80, abs=2e-4)  # the far tail adds < 2e-4


def test_calibration_rejects_infinite_variance() -> None:
    with pytest.raises(ValueError, match="infinite variance"):
        calibrate_lift(METRICS["pareto_a15"], 1000, 0.5)


def test_study_units() -> None:
    units = study_units()
    assert len(units) == 52
    assert len({u.id for u in units}) == 52
    assert not any(u.metric_id == "pareto_a15" and u.hypothesis == "ab" for u in units)
    for u in units:
        assert unit_from_id(u.id) == u


def test_unit_effects_and_seed_keys() -> None:
    aa = Unit("lognorm", 20000, 0.2, "aa")
    ab = Unit("lognorm", 20000, 0.2, "ab")
    assert aa.lift == 0.0 and ab.lift == pytest.approx(aa.design_lift)
    assert aa.tau == pytest.approx(ab.design_lift * METRICS["lognorm"].mean)
    assert ab.true_effect("relative") == pytest.approx(ab.lift)
    assert ab.true_effect("absolute") == pytest.approx(ab.lift * METRICS["lognorm"].mean)
    assert aa.seed_key() != ab.seed_key()
    assert aa.id == "lognorm-n20000-b20-aa" and aa.cell == ab.cell
    stress = Unit("pareto_a15", 1000, 0.5, "aa")
    assert stress.tau > 0  # borrows the pareto_a25 relative lift


def test_unit_from_id_rejects_garbage() -> None:
    with pytest.raises(ValueError):
        unit_from_id("nope-n10-b50-aa")
