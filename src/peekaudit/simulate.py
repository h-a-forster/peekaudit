"""Seeded, chunked generation of cumulative per-arm sufficient statistics.

For a :class:`~peekaudit.scenarios.Unit`, replicates are generated in fixed chunks
(DESIGN 3.4). Chunk ``c`` uses its own stream
``SeedSequence(master_seed, spawn_key=(unit.seed_key(), c))``, so a chunk's data does not
depend on which other chunks exist, the worker count or the order of evaluation.

Within a chunk the control arm A is drawn first, then the treatment arm B. Normal and
Bernoulli metrics use exact per-interval shortcuts (the interval sums have exactly the
right joint distribution); other metrics are drawn per unit and cumulated.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from peekaudit.scenarios import GRID_LOOKS, MASTER_SEED, Metric, Unit, arm_counts, fine_grid

DEFAULT_CHUNK = 250


@dataclass(frozen=True)
class LookStats:
    """Cumulative per-arm sufficient statistics at each look.

    ``n_a``/``n_b`` have shape (K,) (deterministic arm counts); every other array has shape
    (R, K). For ratio metrics, ``sum``/``sumsq`` describe the numerator (e.g. clicks),
    ``dsum``/``dsumsq`` the denominator (e.g. sessions) and ``mdsum`` their cross-products.
    ``n_total_max`` is the planned maximum N_max (used by group sequential methods).
    """

    n_a: np.ndarray
    sum_a: np.ndarray
    sumsq_a: np.ndarray
    n_b: np.ndarray
    sum_b: np.ndarray
    sumsq_b: np.ndarray
    binary: bool = False
    dsum_a: np.ndarray | None = None
    dsumsq_a: np.ndarray | None = None
    mdsum_a: np.ndarray | None = None
    dsum_b: np.ndarray | None = None
    dsumsq_b: np.ndarray | None = None
    mdsum_b: np.ndarray | None = None

    @property
    def ratio(self) -> bool:
        return self.dsum_a is not None

    @property
    def n_reps(self) -> int:
        return int(self.sum_a.shape[0])

    @property
    def n_looks(self) -> int:
        return int(self.sum_a.shape[1])

    def _map(self, f_rk, f_k) -> LookStats:  # type: ignore[no-untyped-def]
        kw = {}
        for name in ("sum_a", "sumsq_a", "sum_b", "sumsq_b"):
            kw[name] = f_rk(getattr(self, name))
        for name in ("dsum_a", "dsumsq_a", "mdsum_a", "dsum_b", "dsumsq_b", "mdsum_b"):
            v = getattr(self, name)
            kw[name] = None if v is None else f_rk(v)
        kw["n_a"] = f_k(self.n_a)
        kw["n_b"] = f_k(self.n_b)
        return replace(self, **kw)

    def take(self, cols: np.ndarray | slice) -> LookStats:
        """Restrict to a subset of looks (columns)."""
        return self._map(lambda a: a[:, cols], lambda a: a[cols])

    def rows(self, idx: np.ndarray | slice) -> LookStats:
        """Restrict to a subset of replicates (rows)."""
        return self._map(lambda a: a[idx], lambda a: a)

    @staticmethod
    def concat(parts: list[LookStats]) -> LookStats:
        """Stack replicate chunks (all parts must share the look structure)."""
        first = parts[0]
        kw = {}
        for name in (
            "sum_a",
            "sumsq_a",
            "sum_b",
            "sumsq_b",
            "dsum_a",
            "dsumsq_a",
            "mdsum_a",
            "dsum_b",
            "dsumsq_b",
            "mdsum_b",
        ):
            if getattr(first, name) is None:
                kw[name] = None
            else:
                kw[name] = np.concatenate([getattr(p, name) for p in parts], axis=0)
        return replace(first, **kw)


def chunk_rng(unit: Unit, chunk: int, master_seed: int = MASTER_SEED) -> np.random.Generator:
    """Generator for one chunk of one unit (independent of every other chunk and unit)."""
    ss = np.random.SeedSequence(master_seed, spawn_key=(unit.seed_key(), chunk))
    return np.random.Generator(np.random.PCG64(ss))


def _cum_at(values: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Cumulative row sums of ``values`` (R, n) evaluated after counts[k] units -> (R, K)."""
    cs = np.cumsum(values, axis=1)
    out = cs[:, counts - 1]
    return out


def _arm_per_unit(
    rng: np.random.Generator, metric: Metric, lift: float, reps: int, counts: np.ndarray
) -> dict[str, np.ndarray]:
    """Draw every unit of one arm and cumulate at the look counts."""
    n = int(counts[-1])
    scale = 1.0 + lift
    if metric.kind == "lognormal":
        x = rng.lognormal(metric.mu, metric.sigma, size=(reps, n)) * scale
    elif metric.kind == "pareto":
        x = (1.0 + rng.pareto(metric.tail, size=(reps, n))) * scale
    elif metric.kind == "normal":
        x = rng.normal(metric.mu * scale, metric.sigma * scale, size=(reps, n))
    elif metric.kind == "bernoulli":
        x = (rng.random(size=(reps, n)) < metric.p * scale).astype(np.float64)
    elif metric.kind == "ratio":
        sessions = 1.0 + rng.poisson(metric.sessions_lambda, size=(reps, n))
        clicks = rng.binomial(sessions.astype(np.int64), metric.click_q * scale).astype(np.float64)
        return {
            "sum": _cum_at(clicks, counts),
            "sumsq": _cum_at(clicks * clicks, counts),
            "dsum": _cum_at(sessions, counts),
            "dsumsq": _cum_at(sessions * sessions, counts),
            "mdsum": _cum_at(clicks * sessions, counts),
        }
    else:  # pragma: no cover - exhaustive
        raise ValueError(metric.kind)
    return {"sum": _cum_at(x, counts), "sumsq": _cum_at(x * x, counts)}


def _arm_exact(
    rng: np.random.Generator, metric: Metric, lift: float, reps: int, counts: np.ndarray
) -> dict[str, np.ndarray]:
    """Exact per-interval shortcut for normal and Bernoulli metrics."""
    inc = np.diff(counts, prepend=0)
    if np.any(inc <= 0):
        raise ValueError("look counts must be strictly increasing")
    scale = 1.0 + lift
    if metric.kind == "bernoulli":
        s = rng.binomial(inc, metric.p * scale, size=(reps, inc.size)).astype(np.float64)
        cs = np.cumsum(s, axis=1)
        return {"sum": cs, "sumsq": cs.copy()}
    if metric.kind == "normal":
        mu, sd = metric.mu * scale, metric.sigma * scale
        m = inc.astype(np.float64)
        s = rng.normal(m * mu, np.sqrt(m) * sd, size=(reps, inc.size))
        # Interval sum of squares = S^2/m + sd^2 * chi2_{m-1}, independent of S.
        chi = 2.0 * rng.standard_gamma((m - 1.0) / 2.0, size=(reps, inc.size))
        q = s * s / m + sd * sd * chi
        return {"sum": np.cumsum(s, axis=1), "sumsq": np.cumsum(q, axis=1)}
    raise ValueError(f"no exact shortcut for {metric.kind}")


def generate_chunk(
    unit: Unit,
    chunk: int,
    reps: int = DEFAULT_CHUNK,
    master_seed: int = MASTER_SEED,
    exact: bool = True,
    looks: int = GRID_LOOKS,
) -> LookStats:
    """Cumulative statistics for ``reps`` replicates of one chunk on the fine look grid.

    ``exact=False`` forces per-unit draws for normal and Bernoulli metrics (used to test the
    shortcuts). ``looks`` truncates the fine grid (tests only).
    """
    metric = unit.metric
    totals = fine_grid(unit.n_max)[:looks]
    n_a, n_b = arm_counts(totals, unit.pi)
    rng = chunk_rng(unit, chunk, master_seed)
    use_exact = exact and metric.kind in ("normal", "bernoulli")
    draw = _arm_exact if use_exact else _arm_per_unit
    a = draw(rng, metric, 0.0, reps, n_a)
    b = draw(rng, metric, unit.lift, reps, n_b)
    return LookStats(
        n_a=n_a,
        sum_a=a["sum"],
        sumsq_a=a["sumsq"],
        n_b=n_b,
        sum_b=b["sum"],
        sumsq_b=b["sumsq"],
        binary=metric.binary,
        dsum_a=a.get("dsum"),
        dsumsq_a=a.get("dsumsq"),
        mdsum_a=a.get("mdsum"),
        dsum_b=b.get("dsum"),
        dsumsq_b=b.get("dsumsq"),
        mdsum_b=b.get("mdsum"),
    )


def chunk_sizes(n_reps: int, chunk: int = DEFAULT_CHUNK) -> list[int]:
    """Replicates per chunk: full chunks, then a remainder."""
    full, rem = divmod(n_reps, chunk)
    return [chunk] * full + ([rem] if rem else [])


def generate(
    unit: Unit,
    n_reps: int,
    chunk: int = DEFAULT_CHUNK,
    master_seed: int = MASTER_SEED,
    order: list[int] | None = None,
) -> LookStats:
    """All replicates of a unit, assembled in chunk order (whatever order they are drawn in)."""
    sizes = chunk_sizes(n_reps, chunk)
    idx = list(range(len(sizes))) if order is None else list(order)
    if sorted(idx) != list(range(len(sizes))):
        raise ValueError("order must be a permutation of the chunk indices")
    parts = {c: generate_chunk(unit, c, sizes[c], master_seed) for c in idx}
    return LookStats.concat([parts[c] for c in range(len(sizes))])
