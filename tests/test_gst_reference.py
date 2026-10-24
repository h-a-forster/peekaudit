"""Like-for-like reference GST (rho = 2) and the spotify.gst comparison on a tiny run."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from peekaudit.gst_reference import (
    cell_counts,
    mcnemar_p,
    render_markdown,
    run_gst_reference,
    summarise,
)
from peekaudit.methods import EvalContext
from peekaudit.methods.base import MethodOutput
from peekaudit.methods.reference import GstPower2, power2_bounds
from peekaudit.runner import StudyConfig, run_study
from peekaudit.scenarios import METRICS, fine_grid

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_power2_bounds_match_spotify_engine_bounds() -> None:
    from peekaudit.methods.spotify_adapter import engine_bounds

    cols = tuple(range(9, 140, 10))
    t = tuple((fine_grid(1000)[list(cols)] / 1000).tolist())
    ours = power2_bounds(1000, cols)
    np.testing.assert_allclose(ours, engine_bounds(t), atol=3e-3)
    over = power2_bounds(1000, (9, 139, 149))
    assert np.isinf(over[-1]) and np.isfinite(over[:-1]).all()


def test_power2_refuses_single_unit_arms() -> None:
    from peekaudit.simulate import LookStats

    ctx = EvalContext(METRICS["bern_p10"], 1000, "custom", np.arange(1), 0.0)
    st = LookStats(np.array([6]), np.zeros((1, 1)), np.zeros((1, 1)),
                   np.array([1]), np.ones((1, 1)), np.ones((1, 1)), True)  # fmt: skip
    out = GstPower2().evaluate(st, ctx)
    assert out.no_result[0, 0] and out.reject[0, 0] == 0


def _out(reject: list[list[int]], no_result: list[list[bool]], zero: list[bool]) -> MethodOutput:
    rej = np.asarray(reject, dtype=np.int8)
    lo = np.where(np.asarray(zero)[:, None], 0.0, -1.0) * np.ones(rej.shape)
    hi = np.ones(rej.shape) * np.where(np.asarray(zero)[:, None], 0.0, 1.0)
    return MethodOutput(rej, lo, hi, np.asarray(no_result), np.zeros(rej.shape))


def test_cell_counts_pairs_and_attribution() -> None:
    nr = [[False, False]] * 3
    sp = _out([[1, 0], [0, 1], [0, 0]], nr, [True, False, False])
    pw = _out([[0, 0], [0, 1], [1, 0]], [[True, False], [False, False], [False, False]],
              [False] * 3)  # fmt: skip
    obf = _out([[0, 0], [0, 0], [0, 0]], nr, [False] * 3)
    c = cell_counts({"spotify.gst": sp, "ref.gst_pow2": pw, "ref.gst_obf": obf})
    assert c["n"] == 3
    assert c["spotify.gst"] == {"rejections": 2, "rejections_look1": 1}
    assert c["paired"] == {"both": 1, "spotify_only": 1, "pow2_only": 1}
    assert c["spotify_reject_where_ref_no_result"] == 1
    assert c["spotify_reject_zero_se"] == 1
    assert mcnemar_p(0, 0) == 1.0 and mcnemar_p(10, 0) < 0.01


def test_comparison_reproduces_run_counts(tmp_path: Path) -> None:
    cfg = StudyConfig(
        "tiny-gstref", r_aa=300, r_ab=100, k_verify=0, k_verify_slow=0, chunk=100,
        methods=("spotify.gst", "ref.gst_obf"), schedules=("daily14", "stream140", "overrun21"),
        units=("bern_p10-n1000-b20-aa", "ratio-n1000-b50-aa", "normal-n1000-b50-ab"),
    )  # fmt: skip
    run_study(cfg, tmp_path, workers=1, echo=lambda _: None)
    assert run_gst_reference(tmp_path, workers=1, echo=lambda _: None) == ["bern_p10-n1000-b20-aa"]
    assert run_gst_reference(tmp_path, workers=1, echo=lambda _: None) == []  # resumed
    s = summarise(tmp_path)
    assert s["all_match_run"] and {r["schedule"] for r in s["cells"]} == {"daily14", "stream140"}
    stream = next(r for r in s["cells"] if r["schedule"] == "stream140")
    assert stream["look1_arms"] == [6, 1]
    assert stream["spotify.gst"]["rejections"] >= stream["ref.gst_pow2"]["rejections"]
    assert "spotify.gst against like-for-like" in render_markdown(s, min_rate=0.0)
    ck = json.loads((tmp_path / "gst_reference" / "bern_p10-n1000-b20-aa.json").read_text())
    assert ck["matches_run"]
