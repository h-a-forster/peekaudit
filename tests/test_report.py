"""Report: families, Holm verdicts, gates and output files on a tiny run."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from peekaudit.report import build_rows, classify, report
from peekaudit.runner import StudyConfig, run_study

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

TINY = StudyConfig(
    "tiny-report",
    r_aa=400,
    r_ab=200,
    k_verify=1,
    k_verify_slow=0,
    chunk=200,
    methods=("ref.naive_z", "ref.acs", "ref.gst_obf", "gbstats.seq.rel"),
    schedules=("fixed1", "daily14", "overrun21"),
    units=("normal-n1000-b50-aa", "normal-n1000-b50-ab", "pareto_a15-n1000-b50-aa"),
)


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    d = tmp_path_factory.mktemp("run")
    run_study(TINY, d, workers=1, echo=lambda _: None)
    return d


def test_families_and_verdicts(run_dir: Path) -> None:
    from peekaudit.report import load_units

    _, units = load_units(run_dir)
    rows = build_rows(units)
    classify(rows)

    def row(method: str, unit: str, sched: str) -> dict:  # type: ignore[type-arg]
        return next(
            r for r in rows if (r["method"], r["unit"], r["schedule"]) == (method, unit, sched)
        )

    peeked = row("ref.naive_z", "normal-n1000-b50-aa", "daily14")
    assert peeked["family"] == "informational" and peeked["verdict"] == "informational"
    assert peeked["rate"] > 0.15
    fixed = row("ref.naive_z", "normal-n1000-b50-aa", "fixed1")
    assert fixed["family"] == "fpr_primary" and fixed["holm_m"] == 1
    assert row("ref.acs", "pareto_a15-n1000-b50-aa", "daily14")["family"] == "fpr_stress"
    assert row("ref.gst_obf", "normal-n1000-b50-aa", "overrun21")["family"] == "informational"
    assert row("ref.acs", "normal-n1000-b50-aa", "daily14")["holm_m"] == 3
    ab = row("ref.acs", "normal-n1000-b50-ab", "daily14")
    assert ab["family"] == "coverage" and 0 <= ab["power"] <= 1 and "type_s_ci95" in ab
    assert row("ref.acs", "normal-n1000-b50-aa", "daily14")["verdict"] in (
        "consistent",
        "conservative",
    )


def test_report_outputs(run_dir: Path, tmp_path: Path) -> None:
    summary = report(run_dir, publish=tmp_path / "pub")
    assert (tmp_path / "pub" / "REPORT.md").exists()
    assert (tmp_path / "pub" / "summary.json").exists()
    data = json.loads((tmp_path / "pub" / "summary.json").read_text(encoding="utf-8"))
    assert len(data["cells"]) == len(summary["cells"]) > 0
    gates = {g["gate"]: g["pass"] for g in summary["gates"]}
    assert gates["engine in-loop agreement"] is True
    assert gates["identity: bayes_vs_z10"] is True  # zero cells: trivially no mismatch
    assert gates["positive control: naive z daily14"] is None  # N_max 20,000 cell not run
    text = (tmp_path / "pub" / "REPORT.md").read_text(encoding="utf-8")
    assert "## Sanity gates" in text and "ref.acs" in text
