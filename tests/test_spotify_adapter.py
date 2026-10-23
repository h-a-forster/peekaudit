"""spotify-confidence 4.1.0 adapter: bounds, call semantics, overrun rescaling, single-look
equivalence and mirror == engine."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from peekaudit.gst import lan_demets_bounds, spend_power
from peekaudit.methods import EvalContext, get_method
from peekaudit.methods.spotify_adapter import (
    SPOTIFY_ALPHA,
    _difference,
    _frame,
    engine_bounds,
    look_bounds,
)
from peekaudit.scenarios import METRICS, SCHEDULES, Unit
from peekaudit.simulate import LookStats, generate_chunk
from peekaudit.verify import compare, random_grid, random_paths

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _ctx(n_max: int, k: int) -> EvalContext:
    return EvalContext(METRICS["normal"], n_max, "custom", np.arange(k), 0.0)


def _path(totals: list[int], seed: int = 0, lift: float = 0.0) -> LookStats:
    rng = np.random.default_rng(seed)
    t = np.asarray(totals)
    n_b = t // 2
    n_a = t - n_b
    xa, xb = rng.normal(1, 1, n_a[-1]), rng.normal(1 + lift, 1, n_b[-1])
    cs = lambda x, n: np.cumsum(x)[n - 1][None, :]  # noqa: E731
    return LookStats(n_a, cs(xa, n_a), cs(xa * xa, n_a), n_b, cs(xb, n_b), cs(xb * xb, n_b))


def test_alpha_constant() -> None:
    assert SPOTIFY_ALPHA == 1 - 0.95 and SPOTIFY_ALPHA != 0.05


def test_engine_bounds_match_reference_solver_and_cache() -> None:
    t = tuple((np.arange(1, 6) / 5).tolist())
    np.testing.assert_allclose(
        engine_bounds(t), [3.0902, 2.7141, 2.4728, 2.2799, 2.1140], atol=1e-4
    )
    np.testing.assert_allclose(
        engine_bounds(t), lan_demets_bounds(np.array(t), spend_power(0.025, 2.0)), atol=2e-3
    )
    assert engine_bounds(t) is engine_bounds(t)


def test_engine_adjusted_alpha_and_decision_rule() -> None:
    st = _path([200, 400, 600, 800], lift=0.25)
    res = _difference(_frame(st, 0, range(4), 800), sequential=True)
    zb = engine_bounds(tuple((np.array([200, 400, 600, 800]) / 800).tolist()))
    np.testing.assert_allclose(res["adjusted_alpha"].to_numpy(float), 2 * (1 - norm.cdf(zb)))
    mirror = get_method("spotify.gst").evaluate(st, _ctx(800, 4))
    assert np.array_equal(mirror.reject[0], np.where(res["is_significant"], 1, 0))


def test_full_path_equals_look_by_look_within_planned_size() -> None:
    st = _path([150, 300, 450, 600, 750], seed=3, lift=0.15)
    full = _difference(_frame(st, 0, range(5), 1000), sequential=True)
    for k in range(5):
        one = _difference(_frame(st, 0, range(k + 1), 1000), sequential=True).iloc[k]
        assert one["adjusted ci_lower"] == full.iloc[k]["adjusted ci_lower"]
        assert one["adjusted_alpha"] == full.iloc[k]["adjusted_alpha"]


def test_overrun_rescales_earlier_fractions() -> None:
    """Regression: a 12-look path planned for 8 looks. At look 10 the full-path call uses
    t = N/N_12 while look-by-look uses t = N/N_10, so the adjusted alpha differs."""
    totals = [100 * k for k in range(1, 13)]
    st = _path(totals, seed=4)
    full = _difference(_frame(st, 0, range(12), 800), sequential=True)
    lbl = _difference(_frame(st, 0, range(10), 800), sequential=True)
    a_full = float(full.iloc[9]["adjusted_alpha"])
    a_lbl = float(lbl.iloc[9]["adjusted_alpha"])
    assert a_full == pytest.approx(0.0198, abs=5e-4)
    assert a_lbl == pytest.approx(0.0297, abs=5e-4)
    zb = look_bounds(np.array(totals), 800)
    assert 2 * (1 - norm.cdf(zb[9])) == pytest.approx(a_lbl, rel=1e-12)


def test_single_look_call_equals_full_path_unadjusted_columns() -> None:
    st = _path([100, 200, 300], seed=5, lift=0.3)
    full = _difference(_frame(st, 0, range(3), None), sequential=False)
    for k in range(3):
        one = _difference(_frame(st, 0, range(k, k + 1), None), sequential=False).iloc[0]
        assert bool(one["is_significant"]) == (float(full.iloc[k]["p-value"]) < SPOTIFY_ALPHA)
        assert one["ci_lower"] == full.iloc[k]["ci_lower"]
        assert one["adjusted ci_lower"] == full.iloc[k]["ci_lower"]


@pytest.mark.parametrize("kind", ["binary", "continuous"])
def test_z_mirror_equals_single_look_engine_on_random_grid(kind: str) -> None:
    m = get_method("spotify.z")
    st = random_grid(kind, 40, seed=21)
    ctx = _ctx(10_000, 40)
    agr = compare(m.evaluate(st, ctx), m.evaluate_engine(st, ctx, per_look=True))  # type: ignore[call-arg]
    assert agr.ok, agr
    assert agr.engine_no_results > 0  # n = 0 inputs are refused by the engine


def test_gst_mirror_equals_engine_on_random_paths() -> None:
    m = get_method("spotify.gst")
    total_rej = 0
    for st, n_max in random_paths(15, seed=8, max_looks=4):
        ctx = _ctx(n_max, st.n_looks)
        agr = compare(m.evaluate(st, ctx), m.evaluate_engine(st, ctx))
        assert agr.ok, agr
        total_rej += agr.engine_rejections
    assert total_rej > 0


def test_mirror_equals_engine_on_study_data_including_overrun() -> None:
    unit = Unit("bern_p10", 1000, 0.2, "ab")
    st = generate_chunk(unit, 0, reps=1)
    for sched, name in (("overrun21", "spotify.gst"), ("daily14", "spotify.z")):
        s = st.take(SCHEDULES[sched])
        ctx = EvalContext.for_schedule(unit.metric, unit.n_max, sched, unit.tau)
        m = get_method(name)
        agr = compare(m.evaluate(s, ctx), m.evaluate_engine(s, ctx))
        assert agr.ok, (name, agr)


def test_edge_cases() -> None:
    m = get_method("spotify.z")
    # n = 1 continuous arm: variance 0/0 = NaN, so no rejection (engine and mirror).
    st = LookStats(
        np.array([5]), np.array([[5.0]]), np.array([[7.0]]),
        np.array([1]), np.array([[3.0]]), np.array([[9.0]]),
    )  # fmt: skip
    agr = compare(m.evaluate(st, _ctx(100, 1)), m.evaluate_engine(st, _ctx(100, 1)))
    assert agr.ok and agr.engine_rejections == 0
    # Zero conversions in both arms: SE = 0 and D = 0, so z is NaN and nothing is significant.
    st = LookStats(np.array([50]), np.zeros((1, 1)), np.zeros((1, 1)),
                   np.array([50]), np.zeros((1, 1)), np.zeros((1, 1)), True)  # fmt: skip
    out = m.evaluate(st, _ctx(100, 1))
    assert out.reject[0, 0] == 0 and not out.no_result[0, 0]
    assert not m.supports(METRICS["ratio"])


def test_engine_indexerror_branch_is_mirrored_as_no_result() -> None:
    """spotify-confidence 4.1.0 raises IndexError for one control unit with value 0 and a
    non-binary treatment (scripts/spotify_indexerror_repro.py); both sides give no result."""
    m = get_method("spotify.z")
    st = LookStats(
        np.array([1, 1, 50]), np.array([[0.0, 1.0, 0.0]]), np.array([[0.0, 1.0, 0.0]]),
        np.array([50, 50, 50]), np.array([[120.0] * 3]), np.array([[900.0] * 3]),
    )  # fmt: skip
    ctx = _ctx(1000, 3)
    mir = m.evaluate(st, ctx)
    eng = m.evaluate_engine(st, ctx, per_look=True)  # type: ignore[call-arg]
    assert mir.no_result[0].tolist() == [True, False, False]
    assert compare(mir, eng).ok
