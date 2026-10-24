"""Tests for the cross-study comparison."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import pytest

from peekaudit.cli import main
from peekaudit.compare import compare_summaries, is_nested, render_markdown, write_comparison


def _summary(name: str, r_aa: int, cells: list[dict[str, Any]], seed: int = 1) -> dict[str, Any]:
    return {
        "study": name,
        "manifest": {
            "config_hash": name + "hash",
            "config": {"master_seed": seed, "chunk": 250, "r_aa": r_aa, "r_ab": r_aa // 2},
        },
        "cells": cells,
    }


def _cell(k: int, n: int, verdict: str = "consistent", family: str = "fpr_primary") -> dict:
    return {
        "method": "m",
        "unit": "u-aa",
        "schedule": "daily14",
        "family": family,
        "verdict": verdict,
        "rejections": k,
        "n": n,
        "rate": k / n,
        "rate_ci95": [max(0.0, k / n - 0.01), k / n + 0.01],
        "cover_hits": n - k,
        "cover_n": n,
        "coverage": 1 - k / n,
        "coverage_ci95": [1 - k / n - 0.01, min(1.0, 1 - k / n + 0.01)],
    }


def test_nesting_detection() -> None:
    small = _summary("q", 4000, [])
    assert is_nested(small, _summary("f", 20000, []))
    assert not is_nested(small, _summary("f", 20000, [], seed=2))
    assert not is_nested(_summary("f", 20000, []), small)
    # A partial last chunk is not a prefix of the larger study's full chunk.
    assert not is_nested(_summary("q", 4100, []), _summary("f", 20000, []))
    assert is_nested(_summary("q", 4100, []), _summary("f", 4100, []))


def test_partial_chunk_is_not_a_prefix() -> None:
    """Why is_nested needs whole chunks: arm B of a shorter chunk is drawn from another offset."""
    import numpy as np

    from peekaudit.scenarios import unit_from_id
    from peekaudit.simulate import generate_chunk

    u = unit_from_id("lognorm-n1000-b50-aa")
    a, b = generate_chunk(u, 0, 100), generate_chunk(u, 0, 250)
    assert np.array_equal(a.sum_a, b.sum_a[:100])
    assert not np.array_equal(a.sum_b, b.sum_b[:100])


def test_nested_counts_are_differenced() -> None:
    small = _summary("q", 4000, [_cell(200, 4000)])
    large = _summary("f", 20000, [_cell(1000, 20000)])
    res = compare_summaries(small, large)
    row = next(r for r in res["comparisons"] if r["kind"] == "rate")
    assert (row["other_k"], row["other_n"]) == (800, 16000)
    assert row["p"] == pytest.approx(1.0)
    assert res["by_kind"]["rate"]["holm_p_lt_0.05"] == 0


def test_detects_a_real_difference_and_verdict_change() -> None:
    small = _summary("q", 4000, [_cell(200, 4000)])
    large = _summary("f", 20000, [_cell(2200, 20000, verdict="violation")])
    res = compare_summaries(small, large)
    assert res["by_kind"]["rate"]["holm_p_lt_0.05"] == 1
    assert len(res["violation_changes"]) == 1
    assert "violation" in render_markdown(res)


def test_inconsistent_nesting_raises() -> None:
    small = _summary("q", 4000, [_cell(300, 4000)])
    large = _summary("f", 20000, [_cell(100, 20000)])
    with pytest.raises(ValueError, match="not nested"):
        compare_summaries(small, large)
    # treated as independent samples it is a valid comparison
    assert compare_summaries(small, large, nested=False)["nested"] is False


def test_no_shared_cells_raises() -> None:
    other = copy.deepcopy(_cell(1, 10))
    other["unit"] = "v-aa"
    with pytest.raises(ValueError, match="share no cells"):
        compare_summaries(_summary("q", 10, [_cell(1, 10)]), _summary("f", 20, [other]))


def test_write_and_cli(tmp_path: Path) -> None:
    import json

    for name, cell, r in (("q", _cell(200, 4000), 4000), ("f", _cell(1000, 20000), 20000)):
        d = tmp_path / name
        d.mkdir()
        (d / "summary.json").write_text(json.dumps(_summary(name, r, [cell])), encoding="utf-8")
    res = write_comparison(tmp_path / "q", tmp_path / "f", tmp_path / "out" / "cmp")
    assert res["nested"]
    assert (tmp_path / "out" / "cmp.md").exists()
    assert main(["compare", str(tmp_path / "q"), str(tmp_path / "f"), "--out",
                 str(tmp_path / "cmp2")]) == 0  # fmt: skip
