"""Data generation: generator moments, exact shortcuts, determinism and chunk independence."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats as st

from peekaudit.scenarios import METRICS, Unit, arm_counts, fine_grid
from peekaudit.simulate import (
    LookStats,
    _arm_per_unit,
    chunk_rng,
    chunk_sizes,
    generate,
    generate_chunk,
)

N_DRAWS = 1_000_000


def _draw(metric_id: str, seed: int = 0, lift: float = 0.0) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    return _arm_per_unit(rng, METRICS[metric_id], lift, 1, np.array([N_DRAWS]))


@pytest.mark.parametrize("metric_id", ["normal", "bern_p10", "bern_p01", "lognorm", "pareto_a25"])
def test_generator_mean_within_clt_tolerance(metric_id: str) -> None:
    m = METRICS[metric_id]
    d = _draw(metric_id)
    mean = d["sum"][0, 0] / N_DRAWS
    assert abs(mean - m.mean) < 5 * math.sqrt(m.variance / N_DRAWS)


@pytest.mark.parametrize("metric_id", ["normal", "bern_p10", "bern_p01"])
def test_generator_variance(metric_id: str) -> None:
    m = METRICS[metric_id]
    d = _draw(metric_id, seed=1)
    s, q = d["sum"][0, 0], d["sumsq"][0, 0]
    var = (q - s * s / N_DRAWS) / (N_DRAWS - 1)
    assert var == pytest.approx(m.variance, rel=0.01)


def test_lognormal_log_scale() -> None:
    rng = np.random.default_rng(2)
    x = rng.lognormal(0.0, 1.5, N_DRAWS)
    logx = np.log(x)
    assert abs(logx.mean()) < 5 * 1.5 / math.sqrt(N_DRAWS)
    assert logx.std() == pytest.approx(1.5, rel=0.005)


@pytest.mark.parametrize("tail", [2.5, 1.5])
def test_pareto_tail_index_hill(tail: float) -> None:
    """Classical Pareto x_m (1 + Lomax): minimum 1 and Hill tail estimate close to a."""
    rng = np.random.default_rng(3)
    x = 1.0 + rng.pareto(tail, N_DRAWS)
    assert x.min() >= 1.0
    k = 20_000
    top = np.sort(x)[-(k + 1) :]
    hill = 1.0 / np.mean(np.log(top[1:] / top[0]))
    # Exact Pareto: the Hill estimator is unbiased-ish with sd a / sqrt(k).
    assert abs(hill - tail) < 5 * tail / math.sqrt(k)


def test_ratio_generator_moments() -> None:
    d = _draw("ratio", seed=4)
    m = METRICS["ratio"]
    q_hat = d["sum"][0, 0] / d["dsum"][0, 0]
    sessions_mean = d["dsum"][0, 0] / N_DRAWS
    assert sessions_mean == pytest.approx(2.5, abs=5 * math.sqrt(1.5 / N_DRAWS))
    assert abs(q_hat - m.mean) < 5 * math.sqrt(m.variance / N_DRAWS)
    # Linearised per-unit variance Var(C - q S) / E[S]^2 = q(1-q)/E[S] = 0.064
    n = N_DRAWS
    s_m, s_d = d["sum"][0, 0], d["dsum"][0, 0]
    var_m = (d["sumsq"][0, 0] - s_m**2 / n) / (n - 1)
    var_d = (d["dsumsq"][0, 0] - s_d**2 / n) / (n - 1)
    cov = (d["mdsum"][0, 0] - s_m * s_d / n) / (n - 1)
    mm, md = s_m / n, s_d / n
    delta = var_m / md**2 - 2 * cov * mm / md**3 + mm**2 * var_d / md**4
    assert delta == pytest.approx(m.variance, rel=0.01)


def test_lift_scales_mean() -> None:
    d = _draw("lognorm", seed=5, lift=0.5)
    assert d["sum"][0, 0] / N_DRAWS == pytest.approx(1.5 * METRICS["lognorm"].mean, rel=0.02)


@pytest.mark.parametrize("metric_id", ["normal", "bern_p10"])
def test_exact_shortcut_matches_per_unit_in_distribution(metric_id: str) -> None:
    unit = Unit(metric_id, 1000, 0.5, "aa")
    exact = generate_chunk(unit, 0, reps=10_000, exact=True, looks=3)
    per_unit = generate_chunk(unit, 1, reps=10_000, exact=False, looks=3)
    for name in ("sum_a", "sumsq_b"):
        for k in (0, 2):
            a = getattr(exact, name)[:, k]
            b = getattr(per_unit, name)[:, k]
            assert st.ks_2samp(a, b, method="asymp").pvalue > 1e-3, (name, k)


def test_bernoulli_sumsq_equals_sum_and_counts() -> None:
    unit = Unit("bern_p01", 1000, 0.2, "ab")
    s = generate_chunk(unit, 0, reps=50)
    assert s.binary and np.array_equal(s.sum_a, s.sumsq_a)
    n_a, n_b = arm_counts(fine_grid(1000), 0.2)
    assert np.array_equal(s.n_a, n_a) and np.array_equal(s.n_b, n_b)
    assert np.all(np.diff(s.sum_a, axis=1) >= 0) and np.all(s.sum_b <= s.n_b)


def test_single_unit_arm_has_zero_variance_exactly() -> None:
    """With one observation, sumsq == sum**2 exactly (engines' n=1 branches depend on it)."""
    unit = Unit("normal", 1000, 0.2, "aa")
    s = generate_chunk(unit, 0, reps=100, looks=1)
    assert s.n_b[0] == 1
    assert np.array_equal(s.sumsq_b[:, 0], s.sum_b[:, 0] ** 2)
    unit = Unit("lognorm", 1000, 0.2, "aa")
    s = generate_chunk(unit, 0, reps=100, looks=1)
    assert np.array_equal(s.sumsq_b[:, 0], s.sum_b[:, 0] ** 2)


def test_same_seed_identical_for_any_chunk_order() -> None:
    unit = Unit("lognorm", 1000, 0.5, "ab")
    a = generate(unit, 600, chunk=250)
    b = generate(unit, 600, chunk=250, order=[2, 0, 1])
    for name in ("sum_a", "sumsq_a", "sum_b", "sumsq_b"):
        assert getattr(a, name).tobytes() == getattr(b, name).tobytes()
    assert a.n_reps == 600 and a.n_looks == 210


def test_chunk_streams_independent_of_other_chunks_and_units() -> None:
    u1 = Unit("ratio", 1000, 0.5, "aa")
    u2 = Unit("ratio", 1000, 0.5, "ab")
    c0 = generate_chunk(u1, 0, reps=20)
    c0_again = generate_chunk(u1, 0, reps=20)
    assert c0.ratio and np.array_equal(c0.mdsum_a, c0_again.mdsum_a)
    assert not np.array_equal(generate_chunk(u1, 1, reps=20).sum_a, c0.sum_a)
    assert not np.array_equal(generate_chunk(u2, 0, reps=20).sum_a, c0.sum_a)
    x = chunk_rng(u1, 3).random(3)
    assert np.array_equal(x, chunk_rng(u1, 3).random(3))


def test_lookstats_take_rows_concat() -> None:
    unit = Unit("ratio", 1000, 0.2, "aa")
    s = generate_chunk(unit, 0, reps=10)
    t = s.take(np.array([0, 5])).rows(slice(2, 4))
    assert t.sum_a.shape == (2, 2) and t.n_a.shape == (2,) and t.dsum_b.shape == (2, 2)
    both = LookStats.concat([s.rows(slice(0, 4)), s.rows(slice(4, 10))])
    assert np.array_equal(both.mdsum_b, s.mdsum_b)


def test_chunk_sizes() -> None:
    assert chunk_sizes(600, 250) == [250, 250, 100]
    assert chunk_sizes(500, 250) == [250, 250]
    with pytest.raises(ValueError):
        generate(Unit("normal", 1000, 0.5, "aa"), 500, order=[0, 0])
