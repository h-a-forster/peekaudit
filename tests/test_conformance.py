"""Conformance check: positive and negative controls, scalar methods, CLI exit codes."""

from __future__ import annotations

import math

import pytest

from peekaudit.cli import main
from peekaudit.conformance import (
    BATTERIES,
    Battery,
    assert_conforms,
    check_conformance,
    critical_rate,
)
from peekaudit.methods import ScalarMethod, get_method
from peekaudit.methods.reference import NaiveZ

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

SMALL = Battery("small", ("normal", "bern_p10"), (2000,), (0.5,), ("daily14",), 600)


def test_reference_sequential_methods_pass() -> None:
    for name in ("ref.gst_obf", "ref.acs"):
        rep = check_conformance(get_method(name), SMALL)
        assert rep.passed, rep.summary()
        assert len(rep.cells) == 2


def test_peeked_naive_z_claimed_anytime_fails() -> None:
    m = NaiveZ()
    m.claim = "anytime"
    with pytest.raises(AssertionError, match="VIOLATION"):
        assert_conforms(m, SMALL)


def test_method_claiming_no_battery_schedule_is_an_error() -> None:
    with pytest.raises(ValueError, match="claims none"):
        check_conformance(NaiveZ(), SMALL)  # fixed-horizon only; battery uses daily14


def test_scalar_method_bonferroni_conforms() -> None:
    """A user-supplied per-look callback: z-test at alpha/14 per look (Bonferroni)."""
    from scipy.stats import norm

    zc = float(norm.isf(0.05 / 14 / 2))

    def decide(n_a, s_a, q_a, n_b, s_b, q_b):  # type: ignore[no-untyped-def]
        va = (q_a - s_a**2 / n_a) / (n_a - 1)
        vb = (q_b - s_b**2 / n_b) / (n_b - 1)
        se = math.sqrt(va / n_a + vb / n_b)
        d = s_b / n_b - s_a / n_a
        if se <= 0:
            return None
        return (int(math.copysign(1, d)) if abs(d) / se > zc else 0, d - zc * se, d + zc * se)

    bat = Battery("tiny", ("normal",), (1000,), (0.5,), ("daily14",), 300)
    assert_conforms(ScalarMethod("bonferroni.z", decide), bat)


def test_critical_rate_matches_design_table() -> None:
    assert critical_rate(4000, 84) == pytest.approx(0.0617, abs=5e-4)
    assert critical_rate(20000, 84) == pytest.approx(0.0551, abs=5e-4)
    assert critical_rate(4000, 96) == pytest.approx(0.0620, abs=5e-4)  # anytime methods
    assert critical_rate(1000, 2) == pytest.approx(0.065, abs=1e-3)


def test_batteries_defined() -> None:
    assert BATTERIES["smoke"].reps == 1000 and BATTERIES["standard"].reps == 4000


def test_cli_conform_exit_codes(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["conform", "--method", "nope.method"]) == 2
    assert main(["conform", "--method", "ref.naive_z", "--claim", "anytime"]) == 1
    assert "FAIL" in capsys.readouterr().out


def test_critical_rate_separates_null_rate_from_fwer() -> None:
    # alpha is the FPR under test (binomial null); the Holm level is fwer / m, fixed at 0.05.
    from scipy.stats import binom

    reps, m, alpha = 1000, 2, 0.10
    rate = critical_rate(reps, m, alpha)
    k = round(rate * reps)
    assert binom.sf(k - 1, reps, alpha) <= 0.05 / m < binom.sf(k - 2, reps, alpha)
    assert rate == pytest.approx(0.120, abs=1e-3)
    assert critical_rate(reps, m, alpha, fwer=0.10) < rate


def test_detection_limit_uses_fwer_not_alpha() -> None:
    bat = Battery("tiny", ("normal",), (1000,), (0.5,), ("daily14",), 300)
    rep = check_conformance(get_method("ref.acs"), bat, alpha=0.10)
    assert rep.fwer == 0.05
    assert rep.detection_limit == critical_rate(300, len(rep.cells), 0.10, 0.05)
