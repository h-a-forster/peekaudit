"""Minimum-arm sensitivity: guarded counts, regeneration matches the run, report section."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from peekaudit.methods.base import MethodOutput
from peekaudit.runner import StudyConfig, run_study
from peekaudit.sensitivity import (
    guarded_counts,
    guarded_looks,
    load_sensitivity,
    run_sensitivity,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

TINY = StudyConfig(
    "tiny-sens",
    r_aa=300,
    r_ab=150,
    k_verify=0,
    k_verify_slow=0,
    chunk=100,
    methods=("ref.naive_z", "ref.msprt", "gbstats.seq.rel"),
    schedules=("fixed1", "stream140"),
    units=("bern_p10-n1000-b20-aa", "normal-n1000-b50-ab"),
)


def test_guarded_looks_uses_the_smaller_arm() -> None:
    # 20% to B: n_B = round(0.2 N); both arms >= 150 needs N >= 748.
    mask = guarded_looks(np.array([7, 700, 748, 1000]), 0.2, 150)
    assert mask.tolist() == [False, False, True, True]


def test_guarded_counts_ignore_early_rejections_only() -> None:
    reject = np.array([[1, 0, 0], [0, 0, -1], [1, 0, 1], [0, 0, 0]], dtype=np.int8)
    nan = np.full(reject.shape, np.nan)
    out = MethodOutput(reject, nan, nan, np.zeros(reject.shape, bool), nan)
    got = guarded_counts(out, np.array([False, True, True]))
    assert got == {"rejections": 3, "rejections_guarded": 2}


def test_sensitivity_regenerates_the_run_and_feeds_the_report(tmp_path: Path) -> None:
    from peekaudit.report import report

    run_study(TINY, tmp_path, workers=1, echo=lambda _: None)
    done = run_sensitivity(tmp_path, workers=1, echo=lambda _: None)
    assert done == ["bern_p10-n1000-b20-aa"]  # A/A units only
    (res,) = load_sensitivity(tmp_path)
    assert res["matches_run"]
    cell = res["cells"]["ref.msprt"]["stream140"]
    assert cell["rejections_guarded"] <= cell["rejections"]
    assert cell["guarded_looks"] == 36 and cell["looks"] == 140
    assert "ref.naive_z" not in res["cells"]  # not a sequential method
    # Resume skips finished units.
    assert run_sensitivity(tmp_path, workers=1, echo=lambda _: None) == []
    summary = report(tmp_path)
    assert summary["sensitivity_min_arm"]
    text = (tmp_path / "report" / "REPORT.md").read_text(encoding="utf-8")
    assert "## Sensitivity: rejections only where both arms hold >= 150 units" in text
    assert "Not computed" not in text
