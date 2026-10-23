"""Conformance check: does a sequential A/B method hold its false-positive rate under peeking?

Any project can call :func:`assert_conforms` from its own test suite with a
:class:`~peekaudit.methods.base.Method` (for example a :class:`ScalarMethod` wrapping a
per-look callback). The check simulates A/A experiments, applies the method under the
schedules it claims, and tests each cell's false-positive rate against alpha with a one-sided
exact binomial test, Holm-adjusted across cells (DESIGN 6.2 and 8).

A pass means "no excess detected at this battery's size", never "proven calibrated": see
:attr:`ConformanceReport.detection_limit`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from scipy.stats import binom

from peekaudit.methods.base import EvalContext, Method, stopping_outcome
from peekaudit.scenarios import ALPHA, METRICS, SCHEDULES, Unit
from peekaudit.simulate import generate
from peekaudit.stats import binom_sf_test, clopper_pearson, holm

#: Family-wise error rate of the Holm-adjusted violation test (DESIGN 6.2). It is a property
#: of the audit, not of the method: ``alpha`` is the FPR a method promises (the null rate of
#: the binomial test), ``FWER`` the level at which an excess over it is declared.
FWER = 0.05


@dataclass(frozen=True)
class Battery:
    name: str
    metrics: tuple[str, ...]
    n_max: tuple[int, ...]
    allocations: tuple[float, ...]
    schedules: tuple[str, ...]
    reps: int


BATTERIES: dict[str, Battery] = {
    "smoke": Battery("smoke", ("normal", "bern_p10"), (2000,), (0.5,), ("daily14",), 1000),
    "quick": Battery(
        "quick",
        ("normal", "bern_p10", "lognorm", "pareto_a25"),
        (1000, 20000),
        (0.5,),
        ("daily14", "stream140"),
        2000,
    ),
    "standard": Battery(
        "standard",
        tuple(m for m, v in METRICS.items() if v.finite_variance),
        (1000, 20000),
        (0.5, 0.2),
        tuple(SCHEDULES),
        4000,
    ),
}


@dataclass
class ConformanceCell:
    unit: str
    schedule: str
    n: int
    false_positives: int
    rate: float
    ci95: tuple[float, float]
    p_value: float
    holm_p: float = 1.0
    violation: bool = False


@dataclass
class ConformanceReport:
    method: str
    battery: str
    alpha: float
    cells: list[ConformanceCell] = field(default_factory=list)
    fwer: float = FWER

    @property
    def violations(self) -> list[ConformanceCell]:
        return [c for c in self.cells if c.violation]

    @property
    def passed(self) -> bool:
        return not self.violations

    @property
    def detection_limit(self) -> float:
        """Smallest FPR that would be flagged in a cell at the strictest Holm level."""
        if not self.cells:
            return float("nan")
        return critical_rate(self.cells[0].n, len(self.cells), self.alpha, self.fwer)

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "battery": self.battery,
            "alpha": self.alpha,
            "fwer": self.fwer,
            "passed": self.passed,
            "detection_limit": self.detection_limit,
            "cells": [asdict(c) for c in self.cells],
        }

    def summary(self) -> str:
        head = (
            f"{self.method} on battery '{self.battery}': "
            f"{'PASS' if self.passed else 'FAIL'} ({len(self.violations)} of {len(self.cells)}"
            f" cells violate FPR <= {self.alpha}; cells can flag FPR >= "
            f"{100 * self.detection_limit:.2f}%)"
        )
        lines = [head]
        for c in self.cells:
            flag = "  VIOLATION" if c.violation else ""
            lines.append(
                f"  {c.unit:<24} {c.schedule:<10} FPR {100 * c.rate:5.2f}% "
                f"[{100 * c.ci95[0]:.2f}, {100 * c.ci95[1]:.2f}]  Holm p {c.holm_p:.3g}{flag}"
            )
        return "\n".join(lines)


def critical_rate(reps: int, m: int, alpha: float = ALPHA, fwer: float = FWER) -> float:
    """Smallest k/reps that the exact one-sided test of H0: FPR <= alpha rejects at level
    fwer/m (the strictest Holm step over m cells)."""
    level = fwer / max(m, 1)
    k = int(binom.isf(level, reps, alpha))
    while binom.sf(k - 1, reps, alpha) > level:
        k += 1
    while k > 0 and binom.sf(k - 2, reps, alpha) <= level:
        k -= 1
    return k / reps


def check_conformance(
    method: Method, battery: str | Battery = "quick", alpha: float = ALPHA, fwer: float = FWER
) -> ConformanceReport:
    """Simulate the battery's A/A units and test FPR <= alpha in every claimed cell.

    ``alpha`` is the false-positive rate the method promises; ``fwer`` is the family-wise
    level of the Holm-adjusted violation test across cells.
    """
    bat = BATTERIES[battery] if isinstance(battery, str) else battery
    schedules = [s for s in bat.schedules if method.claims(s)]
    if not schedules:
        raise ValueError(
            f"{method.name} claims none of the battery's schedules {bat.schedules}"
            f" (claim: {method.claim!r})"
        )
    rep = ConformanceReport(method.name, bat.name, alpha, fwer=fwer)
    for metric_id in bat.metrics:
        metric = METRICS[metric_id]
        if not method.supports(metric):
            continue
        for n_max in bat.n_max:
            for pi in bat.allocations:
                unit = Unit(metric_id, n_max, pi, "aa")
                data = generate(unit, bat.reps)
                for s in schedules:
                    ctx = EvalContext(metric, n_max, s, SCHEDULES[s], unit.tau, alpha)
                    out = method.evaluate(data.take(ctx.cols), ctx)
                    k = int((stopping_outcome(out)["sign"] != 0).sum())
                    rep.cells.append(
                        ConformanceCell(
                            unit.id,
                            s,
                            bat.reps,
                            k,
                            k / bat.reps,
                            clopper_pearson(k, bat.reps),
                            binom_sf_test(k, bat.reps, alpha),
                        )
                    )
    if rep.cells:
        adj = holm([c.p_value for c in rep.cells])
        for c, a in zip(rep.cells, adj, strict=True):
            c.holm_p = float(a)
            c.violation = bool(a < fwer)
    return rep


def assert_conforms(
    method: Method, battery: str | Battery = "smoke", alpha: float = ALPHA, fwer: float = FWER
) -> None:
    """Raise AssertionError listing the violating cells if the method fails the battery."""
    rep = check_conformance(method, battery, alpha, fwer)
    if not rep.passed:
        raise AssertionError(rep.summary())
