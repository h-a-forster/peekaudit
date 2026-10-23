"""CLI smoke tests."""

from __future__ import annotations

import pytest

from peekaudit import __version__
from peekaudit.cli import environment_report, main


def test_version_flag(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_versions_command(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["versions"]) == 0
    out = capsys.readouterr().out
    assert "numpy" in out
    assert "gbstats" in out


def test_no_command_prints_help() -> None:
    assert main([]) == 2


def test_environment_report_keys() -> None:
    report = environment_report()
    assert report["peekaudit"] == __version__
    assert {"python", "numpy", "scipy", "gbstats", "spotify-confidence"} <= set(report)
