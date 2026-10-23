"""Reference methods: closed forms, identities, branches, stopping rule and Monte Carlo size."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.stats import norm

from peekaudit.gst import crossing_probability, lan_demets_bounds, spend_obf
from peekaudit.methods import EvalContext, ScalarMethod, get_method, stopping_outcome
from peekaudit.methods.base import MethodOutput, sample_moments
from peekaudit.methods.reference import ACS, MSPRT, GstOBF, NaiveZ, acs_rho, obf_bounds
from peekaudit.scenarios import METRICS, SCHEDULES, Unit, fine_grid
from peekaudit.simulate import LookStats, generate
from peekaudit.stats import clopper_pearson

Z975 = float(norm.isf(0.025))


def _stats(n_a, s_a, q_a, n_b, s_b, q_b) -> LookStats:  # type: ignore[no-untyped-def]
    """One replicate, K looks, from per-look lists."""
    row = lambda x: np.asarray(x, dtype=float)[None, :]  # noqa: E731
    return LookStats(
        np.asarray(n_a), row(s_a), row(q_a), np.asarray(n_b), row(s_b), row(q_b), binary=False
    )


def _ctx(schedule: str = "fixed1", n_max: int = 1000, tau: float = 0.1) -> EvalContext:
    return EvalContext.for_schedule(METRICS["normal"], n_max, schedule, tau)


def _one_look(n_a=100, mean_a=1.0, var_a=1.0, n_b=100, mean_b=1.3, var_b=2.0) -> LookStats:  # type: ignore[no-untyped-def]
    s_a, s_b = n_a * mean_a, n_b * mean_b
    q_a = var_a * (n_a - 1) + s_a**2 / n_a
    q_b = var_b * (n_b - 1) + s_b**2 / n_b
    return _stats([n_a], [s_a], [q_a], [n_b], [s_b], [q_b])


def test_sample_moments_closed_form() -> None:
    ma, va, mb, vb = sample_moments(_one_look())
    assert ma[0, 0] == pytest.approx(1.0) and mb[0, 0] == pytest.approx(1.3)
    assert va[0, 0] == pytest.approx(1.0) and vb[0, 0] == pytest.approx(2.0)


def test_naive_z_closed_form() -> None:
    out = NaiveZ().evaluate(_one_look(), _ctx())
    se = math.sqrt(1 / 100 + 2 / 100)
    z = 0.3 / se
    assert out.reject[0, 0] == (1 if z > Z975 else 0)
    assert out.lower[0, 0] == pytest.approx(0.3 - Z975 * se)
    assert out.upper[0, 0] == pytest.approx(0.3 + Z975 * se)
    assert out.margin[0, 0] == pytest.approx((z - Z975) / Z975)
    neg = NaiveZ().evaluate(_one_look(mean_b=0.5), _ctx())
    assert neg.reject[0, 0] == -1


def test_msprt_closed_form() -> None:
    tau = 0.2
    v = 1 / 100 + 2 / 100
    d = 0.3
    lam = math.sqrt(v / (v + tau**2)) * math.exp(tau**2 * d**2 / (2 * v * (v + tau**2)))
    half = math.sqrt(v * (v + tau**2) / tau**2 * (2 * math.log(20) + math.log((v + tau**2) / v)))
    out = MSPRT().evaluate(_one_look(), _ctx(tau=tau))
    assert out.reject[0, 0] == (1 if lam >= 20 else 0)
    assert out.upper[0, 0] - out.lower[0, 0] == pytest.approx(2 * half)
    assert out.margin[0, 0] == pytest.approx((math.log(lam) - math.log(20)) / math.log(20))


def test_acs_equals_msprt_with_data_dependent_tau() -> None:
    """ref.acs is ref.msprt with tau^2 = rho^2 N V (DESIGN 4.4)."""
    rng = np.random.default_rng(0)
    for _ in range(50):
        n_a, n_b = (int(x) for x in rng.integers(20, 3000, 2))
        st = _one_look(n_a, 1.0, rng.uniform(0.5, 2), n_b, rng.normal(1, 0.1), rng.uniform(0.5, 2))
        ctx = _ctx(n_max=int(rng.integers(500, 30000)))
        acs = ACS().evaluate(st, ctx)
        ma, va, mb, vb = sample_moments(st)
        v = va[0, 0] / n_a + vb[0, 0] / n_b
        tau = math.sqrt(acs_rho(0.05, ctx.n_max) ** 2 * (n_a + n_b) * v)
        ms = MSPRT().evaluate(st, _ctx(tau=tau))
        np.testing.assert_allclose(acs.lower, ms.lower, rtol=1e-10)
        np.testing.assert_allclose(acs.upper, ms.upper, rtol=1e-10)
        if abs(acs.margin[0, 0]) > 1e-9:
            assert acs.reject[0, 0] == ms.reject[0, 0]


def test_acs_rho_matches_paper_formula() -> None:
    a = 0.05
    assert acs_rho(a, 5000) == pytest.approx(
        math.sqrt((-2 * math.log(a) + math.log(-2 * math.log(a) + 1)) / 5000)
    )


def test_gst_obf_bounds_and_overrun() -> None:
    cols = tuple(int(c) for c in SCHEDULES["daily14"])
    t = fine_grid(20000)[list(cols)] / 20000
    np.testing.assert_allclose(obf_bounds(20000, cols), lan_demets_bounds(t, spend_obf(0.025)))
    assert obf_bounds(1000, (139,))[0] == pytest.approx(Z975, abs=1e-6)
    over = obf_bounds(1000, tuple(int(c) for c in SCHEDULES["overrun21"]))
    assert np.all(np.isfinite(over[:14])) and np.all(np.isinf(over[14:]))


def test_gst_obf_decision_uses_bound() -> None:
    st = _one_look(10000, 1.0, 1.0, 10000, 1.03, 1.0)  # z ~ 2.12
    out = GstOBF().evaluate(st, _ctx("fixed1", n_max=20000))
    assert out.reject[0, 0] == 1
    with pytest.raises(ValueError):
        ctx = EvalContext(METRICS["normal"], 1000, "fixed1", SCHEDULES["fixed1"], 0.1, 0.1)
        GstOBF().evaluate(st, ctx)


@pytest.mark.parametrize("cls", [NaiveZ, MSPRT, ACS, GstOBF])
def test_no_result_branches(cls) -> None:  # type: ignore[no-untyped-def]
    # n = 1 in B, and zero variance in both arms: no decision and no interval.
    st = _stats([5, 5], [5.0, 5.0], [5.0, 5.0], [1, 5], [2.0, 5.0], [4.0, 5.0])
    ctx = EvalContext(METRICS["normal"], 1000, "x", np.array([0, 1]), 0.1)
    if cls is GstOBF:
        ctx = EvalContext(METRICS["normal"], 1000, "x", np.array([69, 139]), 0.1)
    out = cls().evaluate(st, ctx)
    assert out.no_result.all() and (out.reject == 0).all() and np.isnan(out.lower).all()


def test_ratio_moments_match_delta_formula() -> None:
    unit = Unit("ratio", 1000, 0.5, "aa")
    st = generate(unit, 3).take(np.array([50]))
    ma, va, _, _ = sample_moments(st)
    n = st.n_a[0]
    s, q, d, dq, md = (
        getattr(st, f)[0, 0] for f in ("sum_a", "sumsq_a", "dsum_a", "dsumsq_a", "mdsum_a")
    )
    vm, vd = (q - s * s / n) / (n - 1), (dq - d * d / n) / (n - 1)
    cov = (md - s * d / n) / (n - 1)
    mm, dm = s / n, d / n
    assert ma[0, 0] == pytest.approx(s / d)
    assert va[0, 0] == pytest.approx(vm / dm**2 - 2 * cov * mm / dm**3 + mm**2 * vd / dm**4)


def test_stopping_outcome() -> None:
    rej = np.array([[0, 1, -1], [0, 0, 0], [-1, 0, 0]], dtype=np.int8)
    lo = np.arange(9, dtype=float).reshape(3, 3)
    out = MethodOutput(rej, lo, lo + 1, np.zeros((3, 3), bool), np.zeros((3, 3)))
    s = stopping_outcome(out)
    assert s["stop"].tolist() == [1, 2, 0]
    assert s["sign"].tolist() == [1, 0, -1]
    assert s["lower"].tolist() == [1.0, 5.0, 6.0]


def test_scalar_method_matches_vectorised_naive_z() -> None:
    def decide(n_a, s_a, q_a, n_b, s_b, q_b):  # type: ignore[no-untyped-def]
        if n_a < 2 or n_b < 2:
            return None
        va = (q_a - s_a**2 / n_a) / (n_a - 1)
        vb = (q_b - s_b**2 / n_b) / (n_b - 1)
        se = math.sqrt(va / n_a + vb / n_b)
        d = s_b / n_b - s_a / n_a
        rej = int(np.sign(d)) if abs(d) / se > Z975 else 0
        return rej, d - Z975 * se, d + Z975 * se

    st = generate(Unit("normal", 1000, 0.2, "ab"), 20).take(SCHEDULES["daily14"])
    ctx = _ctx("daily14")
    a = ScalarMethod("scalar.z", decide, claim="fixed").evaluate(st, ctx)
    b = NaiveZ().evaluate(st, ctx)
    assert np.array_equal(a.reject, b.reject)
    np.testing.assert_allclose(a.lower, b.lower, rtol=1e-12)
    assert not ScalarMethod("x", decide).supports(METRICS["ratio"])


def test_claims() -> None:
    assert NaiveZ().claims("fixed1") and not NaiveZ().claims("daily14")
    assert GstOBF().claims("stream140") and not GstOBF().claims("overrun21")
    assert ACS().claims("overrun21") and MSPRT().claims("stream140")
    assert get_method("ref.acs").name == "ref.acs"
    with pytest.raises(KeyError):
        get_method("nope")


def _fpr(method_name: str, schedule: str, reps: int, n_max: int = 1000) -> tuple[int, int]:
    unit = Unit("normal", n_max, 0.5, "aa")
    st = generate(unit, reps).take(SCHEDULES[schedule])
    ctx = EvalContext.for_schedule(unit.metric, n_max, schedule, unit.tau)
    out = get_method(method_name).evaluate(st, ctx)
    return int((stopping_outcome(out)["sign"] != 0).sum()), reps


def test_mc_naive_z_daily14_matches_exact_crossing_probability() -> None:
    k, n = _fpr("ref.naive_z", "daily14", 4000)
    t = fine_grid(1000)[SCHEDULES["daily14"]] / 1000
    exact = crossing_probability(t, [Z975] * 14).sum()
    lo, hi = clopper_pearson(k, n, 0.999)
    assert lo <= exact <= hi, (k / n, exact)
    assert exact > 0.10


@pytest.mark.parametrize("name", ["ref.acs", "ref.msprt", "ref.gst_obf"])
def test_mc_sequential_methods_hold_size_on_normal(name: str) -> None:
    k, n = _fpr(name, "daily14", 4000)
    lo, _ = clopper_pearson(k, n, 0.999)
    assert lo <= 0.05, (name, k / n)


@pytest.mark.slow
@pytest.mark.parametrize("name", ["ref.acs", "ref.msprt", "ref.gst_obf", "ref.naive_z"])
def test_mc_slow_stream140_n20000(name: str) -> None:
    k, n = _fpr(name, "stream140", 20000, n_max=20000)
    lo, hi = clopper_pearson(k, n, 0.999)
    if name == "ref.naive_z":
        t = fine_grid(20000)[SCHEDULES["stream140"]] / 20000
        exact = crossing_probability(t, [Z975] * 140).sum()
        assert lo <= exact <= hi
    elif name == "ref.gst_obf":
        assert lo <= 0.05 <= hi
    else:
        assert lo <= 0.05
