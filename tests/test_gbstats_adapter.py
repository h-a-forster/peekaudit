"""gbstats 0.8.0 adapter: engine fixture, mirror == engine (random grid, branches, study
data), the Bayes identity and the ACS identity."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from peekaudit.methods import EvalContext, get_method
from peekaudit.methods.gbstats_adapter import GbstatsSequential, _diff_var
from peekaudit.methods.reference import ACS
from peekaudit.scenarios import METRICS, SCHEDULES, Unit
from peekaudit.simulate import LookStats, generate_chunk
from peekaudit.verify import compare, random_grid

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

GB = [
    "gbstats.seq.rel",
    "gbstats.seq.rel.tuned",
    "gbstats.seq.abs",
    "gbstats.ttest.rel",
    "gbstats.bayes.ctw95",
]


def _ctx(n_max: int = 10_000, k: int = 1) -> EvalContext:
    return EvalContext(METRICS["normal"], n_max, "custom", np.arange(k), 0.1)


def _single(n_a, s_a, q_a, n_b, s_b, q_b, binary=False, **ratio) -> LookStats:  # type: ignore[no-untyped-def]
    r = lambda x: np.array([[float(x)]])  # noqa: E731
    extra = {k: r(v) for k, v in ratio.items()}
    return LookStats(
        np.array([n_a]), r(s_a), r(q_a), np.array([n_b]), r(s_b), r(q_b), binary, **extra
    )


def test_engine_fixture_through_the_mirror() -> None:
    """gbstats 0.8.0 tests/frequentist/test_tests.py::test_sequential_test_runs_prop
    (values rounded to 5 decimals there)."""
    st = _single(3000, 1396, 1396, 3461, 2422, 2422, binary=True)
    out = get_method("gbstats.seq.rel").evaluate(st, _ctx())
    d, _, _, _ = _diff_var(st, relative=True)
    assert round(float(d[0, 0]), 5) == 0.50386
    assert [round(float(out.lower[0, 0]), 5), round(float(out.upper[0, 0]), 5)] == [
        0.40098,
        0.60675,
    ]
    assert out.reject[0, 0] == 1


@pytest.mark.parametrize("name", GB)
@pytest.mark.parametrize("kind", ["binary", "continuous", "ratio"])
def test_mirror_equals_engine_on_random_grid(name: str, kind: str) -> None:
    m = get_method(name)
    st = random_grid(kind, 300, seed=11)
    ctx = _ctx(k=300)
    agr = compare(m.evaluate(st, ctx), m.evaluate_engine(st, ctx))
    assert agr.ok, agr
    assert agr.engine_rejections > 0 and agr.engine_no_results > 0


@pytest.mark.parametrize("name", GB)
def test_default_output_branches(name: str) -> None:
    m = get_method(name)
    cases = [
        _single(100, 0, 0, 100, 5, 5, binary=True),  # zero control mean
        _single(100, 100, 100, 100, 50, 50, binary=True),  # all-ones control: variance 0
        _single(1, 2.0, 4.0, 100, 50, 60),  # n = 1 control: SampleMean variance 0
        _single(100, 50.0, 60.0, 0, 0.0, 0.0),  # empty treatment arm
        _single(50, 10.0, 30.0, 50, 12, 40, dsum_a=0, dsumsq_a=0, mdsum_a=0, dsum_b=60,
                dsumsq_b=100, mdsum_b=20),  # ratio with zero denominator
    ]  # fmt: skip
    for st in cases:
        mir, eng = m.evaluate(st, _ctx()), m.evaluate_engine(st, _ctx())
        assert mir.no_result[0, 0] and eng.no_result[0, 0]
        assert mir.reject[0, 0] == 0 == eng.reject[0, 0]
        assert mir.lower[0, 0] == 0.0 == eng.lower[0, 0]


@pytest.mark.parametrize(
    "unit_id", ["lognorm-n1000-b20-ab", "ratio-n1000-b50-aa", "bern_p01-n1000-b20-ab"]
)
def test_mirror_equals_engine_on_study_data(unit_id: str) -> None:
    from peekaudit.scenarios import unit_from_id

    unit = unit_from_id(unit_id)
    st = generate_chunk(unit, 0, reps=2).take(SCHEDULES["stream140"])
    ctx = EvalContext.for_schedule(unit.metric, unit.n_max, "stream140", unit.tau)
    for name in GB:
        m = get_method(name)
        agr = compare(m.evaluate(st, ctx), m.evaluate_engine(st, ctx))
        assert agr.ok, (name, agr)


def test_bayes_rule_is_a_level_010_relative_z_test() -> None:
    st = generate_chunk(Unit("lognorm", 1000, 0.5, "aa"), 0, reps=250).take(SCHEDULES["daily14"])
    ctx = EvalContext.for_schedule(METRICS["lognorm"], 1000, "daily14", 1.0)
    out = get_method("gbstats.bayes.ctw95").evaluate(st, ctx)
    d, v, default, _ = _diff_var(st, relative=True)
    z = d / np.sqrt(v)
    zt = norm.ppf(0.95)
    expected = np.where(np.abs(z) >= zt, np.sign(z), 0)
    tie = np.abs(np.abs(z) - zt) < 1e-9
    assert np.array_equal(out.reject[~tie & ~default], expected[~tie & ~default])
    assert (out.reject != 0).mean() > 0.1  # peeking a 10% test 14 times rejects a lot


@pytest.mark.parametrize("kind", ["continuous", "ratio"])
def test_ref_acs_equals_engine_absolute_sequential_with_tuning_nmax(kind: str) -> None:
    """DESIGN 4.4 identity: ref.acs == gbstats.seq.abs with N* = N_max on non-binary data
    (away from gbstats' zero-control-mean branch)."""
    st = random_grid(kind, 300, seed=5, allow_empty=False)
    ctx = _ctx(n_max=7_000, k=300)
    ref = ACS().evaluate(st, ctx)
    eng = GbstatsSequential(relative=False, tuned=True).evaluate_engine(st, ctx)
    both = ~ref.no_result & ~eng.no_result
    assert both.sum() > 150
    assert np.array_equal(ref.reject[both], eng.reject[both])
    np.testing.assert_allclose(ref.lower[both], eng.lower[both], rtol=1e-12)


def test_sequential_tuning_changes_width() -> None:
    st = _single(3000, 1396, 1396, 3461, 2422, 2422, binary=True)
    a = GbstatsSequential(relative=True).evaluate(st, _ctx(n_max=5000))
    b = GbstatsSequential(relative=True, tuned=True).evaluate(st, _ctx(n_max=5000))
    c = GbstatsSequential(relative=True, tuned=True).evaluate(st, _ctx(n_max=50_000))
    assert a.lower[0, 0] == b.lower[0, 0] and c.lower[0, 0] != a.lower[0, 0]
    assert get_method("gbstats.seq.rel.tuned").name == "gbstats.seq.rel.tuned"
