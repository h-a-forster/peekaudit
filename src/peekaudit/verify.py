"""Mirror-versus-engine agreement: randomised input grids and in-loop replicate checks.

DESIGN 5.2: decisions must agree exactly (a decision is a tie when the mirror's statistic
is within 1e-9 relative of its threshold; ties are counted, not failed), "no result" flags
must agree, and intervals must agree to rtol 1e-9 / atol 1e-12.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from peekaudit.methods.base import EvalContext, Method, MethodOutput
from peekaudit.scenarios import METRICS
from peekaudit.simulate import LookStats

TIE_TOL = 1e-9
RTOL = 1e-9
ATOL = 1e-12


@dataclass
class Agreement:
    """Counts from comparing a mirror with its engine on the same inputs."""

    cells: int = 0
    decision_mismatches: int = 0
    ties: int = 0
    tie_mismatches: int = 0
    no_result_mismatches: int = 0
    interval_mismatches: int = 0
    max_interval_rel_diff: float = 0.0
    engine_rejections: int = 0
    engine_no_results: int = 0

    @property
    def ok(self) -> bool:
        return (
            self.decision_mismatches == 0
            and self.no_result_mismatches == 0
            and self.interval_mismatches == 0
        )

    def add(self, other: Agreement) -> Agreement:
        out = Agreement()
        for f in asdict(self):
            a, b = getattr(self, f), getattr(other, f)
            setattr(out, f, max(a, b) if f == "max_interval_rel_diff" else a + b)
        return out

    def to_dict(self) -> dict[str, float | int | bool]:
        d: dict[str, float | int | bool] = dict(asdict(self))
        d["ok"] = self.ok
        return d


def compare(mirror: MethodOutput, engine: MethodOutput) -> Agreement:
    """Compare two outputs of identical shape."""
    tie = np.abs(mirror.margin) < TIE_TOL
    mism = mirror.reject != engine.reject
    nr_mism = mirror.no_result != engine.no_result
    both = ~mirror.no_result & ~engine.no_result
    rel = 0.0
    bad_int = 0
    for m_arr, e_arr in ((mirror.lower, engine.lower), (mirror.upper, engine.upper)):
        m, e = m_arr[both], e_arr[both]
        same_nan = np.isnan(m) & np.isnan(e)
        same_inf = np.isinf(m) & np.isinf(e) & (np.sign(m) == np.sign(e))
        fin = ~(same_nan | same_inf)
        with np.errstate(invalid="ignore"):
            close = np.isclose(m[fin], e[fin], rtol=RTOL, atol=ATOL)
            diff = np.abs(m[fin] - e[fin]) / np.maximum(np.abs(e[fin]), ATOL)
        bad_int += int((~close).sum())
        if diff.size:
            rel = max(rel, float(np.nanmax(diff)) if np.isfinite(diff).any() else np.inf)
    return Agreement(
        cells=int(mirror.reject.size),
        decision_mismatches=int((mism & ~tie).sum()),
        ties=int(tie.sum()),
        tie_mismatches=int((mism & tie).sum()),
        no_result_mismatches=int(nr_mism.sum()),
        interval_mismatches=bad_int,
        max_interval_rel_diff=rel,
        engine_rejections=int((engine.reject != 0).sum()),
        engine_no_results=int(engine.no_result.sum()),
    )


def verify_rows(
    method: Method, stats: LookStats, ctx: EvalContext, mirror: MethodOutput
) -> Agreement:
    """Run the real engine on every replicate of ``stats`` and compare with ``mirror``."""
    engine = method.evaluate_engine(stats, ctx)
    return compare(mirror, engine)


# --------------------------------------------------------------------------- random grids


def _continuous_sums(
    rng: np.random.Generator, n: np.ndarray, mean: np.ndarray, var: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    s = n * mean
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.where(n > 0, var * np.maximum(n - 1, 0) + s * s / np.maximum(n, 1), 0.0)
    return s, q


def random_grid(kind: str, size: int, seed: int = 0, allow_empty: bool = True) -> LookStats:
    """``size`` independent single-look inputs as one replicate with ``size`` looks.

    Includes degenerate inputs: n in {0, 1, 2}, zero sums, all-ones Bernoulli, zero control
    mean, zero and huge variances, very unequal n. ``kind`` is "binary", "continuous" or
    "ratio".
    """
    rng = np.random.default_rng(seed)
    lo = 0 if allow_empty else 1

    def counts() -> np.ndarray:
        n = rng.integers(2, 20_000, size)
        u = rng.random(size)
        n = np.where(u < 0.15, rng.integers(lo, 3, size), n)  # n in {0,1,2}
        n = np.where((u >= 0.15) & (u < 0.25), rng.integers(3, 30, size), n)
        return n.astype(np.int64)

    n_a, n_b = counts(), counts()
    if kind == "binary":
        p_a, p_b = rng.uniform(0, 0.3, size), rng.uniform(0, 0.3, size)
        s_a = rng.binomial(n_a, p_a).astype(float)
        s_b = rng.binomial(n_b, p_b).astype(float)
        u = rng.random(size)
        s_a = np.where(u < 0.08, 0.0, s_a)  # zero control mean
        s_b = np.where((u >= 0.08) & (u < 0.14), 0.0, s_b)
        s_a = np.where((u >= 0.14) & (u < 0.18), n_a.astype(float), s_a)  # all ones
        s_b = np.where((u >= 0.18) & (u < 0.22), n_b.astype(float), s_b)
        row = lambda x: x[None, :]  # noqa: E731
        return LookStats(n_a, row(s_a), row(s_a.copy()), n_b, row(s_b), row(s_b.copy()), True)

    mean_a = rng.lognormal(0, 1, size)
    mean_b = mean_a * rng.lognormal(0, 0.1, size)
    var_a = rng.lognormal(0, 2, size)
    var_b = var_a * rng.lognormal(0, 0.5, size)
    u = rng.random(size)
    mean_a = np.where(u < 0.06, 0.0, mean_a)
    var_a = np.where((u >= 0.06) & (u < 0.11), 0.0, var_a)
    var_b = np.where((u >= 0.11) & (u < 0.16), 0.0, var_b)
    var_a = np.where((u >= 0.16) & (u < 0.2), var_a * 1e8, var_a)  # huge variance
    s_a, q_a = _continuous_sums(rng, n_a, mean_a, var_a)
    s_b, q_b = _continuous_sums(rng, n_b, mean_b, var_b)
    row = lambda x: x[None, :]  # noqa: E731
    if kind == "continuous":
        return LookStats(n_a, row(s_a), row(q_a), n_b, row(s_b), row(q_b), False)
    if kind != "ratio":
        raise ValueError(kind)
    extra = {}
    for arm, n, s_m, q_m in (("a", n_a, s_a, q_a), ("b", n_b, s_b, q_b)):
        d_mean = rng.uniform(1, 4, size)
        d_mean = np.where(rng.random(size) < 0.05, 0.0, d_mean)
        d_var = rng.uniform(0, 3, size)
        s_d, q_d = _continuous_sums(rng, n, d_mean, d_var)
        var_m = np.where(n > 1, (q_m - s_m * s_m / np.maximum(n, 1)) / np.maximum(n - 1, 1), 0)
        corr = rng.uniform(-0.9, 0.9, size)
        cov = corr * np.sqrt(np.maximum(var_m, 0) * d_var)
        md = np.where(n > 0, cov * np.maximum(n - 1, 0) + s_m * s_d / np.maximum(n, 1), 0.0)
        extra[f"dsum_{arm}"], extra[f"dsumsq_{arm}"], extra[f"mdsum_{arm}"] = (
            row(s_d),
            row(q_d),
            row(md),
        )
    return LookStats(n_a, row(s_a), row(q_a), n_b, row(s_b), row(q_b), False, **extra)


def random_paths(size: int, seed: int = 0, max_looks: int = 6) -> list[tuple[LookStats, int]]:
    """Random multi-look paths (with planned N_max, sometimes overrun) for group sequential
    adapters. Each item is (stats with one replicate, planned N_max)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(size):
        k = int(rng.integers(1, max_looks + 1))
        binary = bool(rng.random() < 0.4)
        n_max = int(rng.integers(20, 5000))
        frac = np.sort(rng.uniform(0.02, 1.4 if rng.random() < 0.3 else 1.0, k))
        totals = np.maximum.accumulate(np.maximum((frac * n_max).astype(np.int64), 4))
        totals = totals + np.arange(k)  # strictly increasing
        pi = float(rng.choice([0.5, 0.2]))
        n_b = np.maximum(np.floor(pi * totals + 0.5).astype(np.int64), 1)
        n_a = totals - n_b
        inc_a, inc_b = np.diff(n_a, prepend=0), np.diff(n_b, prepend=0)
        lift = float(rng.uniform(0.0, 0.6)) if rng.random() < 0.5 else 0.0
        if binary:
            p = rng.uniform(0.005, 0.3)
            s_a = np.cumsum(rng.binomial(inc_a, p)).astype(float)
            s_b = np.cumsum(rng.binomial(inc_b, p * (1 + lift))).astype(float)
            q_a, q_b = s_a.copy(), s_b.copy()
        else:
            xa = [rng.lognormal(0, 1, int(m)) for m in inc_a]
            xb = [rng.lognormal(0, 1, int(m)) * (1 + lift) for m in inc_b]
            s_a = np.cumsum([x.sum() for x in xa])
            s_b = np.cumsum([x.sum() for x in xb])
            q_a = np.cumsum([(x * x).sum() for x in xa])
            q_b = np.cumsum([(x * x).sum() for x in xb])
        row = lambda x: np.asarray(x, dtype=float)[None, :]  # noqa: E731
        st = LookStats(n_a, row(s_a), row(q_a), n_b, row(s_b), row(q_b), binary)
        out.append((st, n_max))
    return out


def grid_agreement(method: Method, n: int = 2000, seed: int = 0) -> dict[str, Agreement]:
    """Mirror vs engine on random single-look inputs (or random paths for spotify.gst)."""
    results: dict[str, Agreement] = {}
    if method.name == "spotify.gst":
        total = Agreement()
        for st, n_max in random_paths(n, seed):
            ctx = EvalContext(METRICS["normal"], n_max, "custom", np.arange(st.n_looks), 0.0)
            total = total.add(compare(method.evaluate(st, ctx), method.evaluate_engine(st, ctx)))
        results["paths"] = total
        return results
    kinds = ["binary", "continuous"] + (["ratio"] if method.supports(METRICS["ratio"]) else [])
    per = max(1, n // len(kinds))
    for i, kind in enumerate(kinds):
        st = random_grid(kind, per, seed + i)
        ctx = EvalContext(METRICS["normal"], 10_000, "custom", np.arange(per), 0.1)
        mirror = method.evaluate(st, ctx)
        if method.name == "spotify.z":
            engine = method.evaluate_engine(st, ctx, per_look=True)  # type: ignore[call-arg]
        else:
            engine = method.evaluate_engine(st, ctx)
        results[kind] = compare(mirror, engine)
    return results
