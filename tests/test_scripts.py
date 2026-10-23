"""Standalone scripts: wrapper logging and the upstream reproductions."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_run_study_logs_no_absolute_interpreter_path() -> None:
    mod = _load("run_study")
    cmd = [str(Path("/home/someone/project/.venv/bin/python")), "-m", "peekaudit", "run"]
    assert mod.loggable_command(cmd) == "python -m peekaudit run"


def test_published_progress_logs_hold_no_local_paths() -> None:
    for log in (ROOT / "docs" / "results").rglob("progress.log"):
        text = log.read_text(encoding="utf-8")
        assert ":\\" not in text and "/home/" not in text and "/Users/" not in text, log


def test_spotify_overrun_repro_values() -> None:
    mod = _load("spotify_overrun_repro")
    full = mod.adjusted_alpha(12)
    within = mod.adjusted_alpha(8)
    at10 = mod.adjusted_alpha(10).loc[10]
    assert abs(at10 - 0.0297) < 5e-4
    assert abs(full.loc[10] - 0.0198) < 5e-4
    assert abs(within.loc[8] - 0.0311) < 5e-4
    assert abs(full.loc[8] - 0.0128) < 5e-4


def test_gbstats_rare_conversion_repro_runs(capsys) -> None:  # type: ignore[no-untyped-def]
    mod = _load("gbstats_rare_conversion_repro")
    mod.main(["--reps", "3", "--looks", "5", "--n", "2000"])
    out = capsys.readouterr().out
    assert "relative" in out and "absolute" in out
