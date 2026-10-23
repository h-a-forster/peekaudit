"""Standalone scripts: the resume wrapper never logs the interpreter's absolute path."""

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
