"""Living audit (scripts/living.py, peekaudit.living): battery, comparison and the issue step.

The weekly job is replayed end to end without network: target results already on disk
(as a finished run leaves them) with one verdict changed against a baseline, then the
issue step against a fake ``gh`` that records its calls.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from peekaudit import living

ROOT = Path(__file__).resolve().parents[1]


def _script():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("living_script", ROOT / "scripts/living.py")
    mod = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


FAKE_GH = """
import json, sys
from pathlib import Path
state = Path(sys.argv[1]); args = sys.argv[2:]
calls = json.loads(state.read_text()) if state.exists() else []
calls.append(args); state.write_text(json.dumps(calls))
if args[:2] == ["issue", "list"]:
    opened = [c for c in calls if c[:2] == ["issue", "create"]]
    print(json.dumps([{"number": 7, "title": c[c.index("--title") + 1]} for c in opened]))
"""


def test_committed_baseline_covers_every_target() -> None:
    base = json.loads((ROOT / "docs/results/living/baseline.json").read_text("utf-8"))
    assert set(base["targets"]) == set(living.TARGETS)
    for name, t in base["targets"].items():
        assert t["verdicts"], name
        # the battery never errors on the audited release
        assert not any(k.endswith("/error") for k in t["verdicts"]), name
    commit = base["targets"]["gbstats-main"]["environment"]["gbstats commit"]
    assert commit == living.GROWTHBOOK_COMMIT


def test_engine_battery_verdicts_and_resume(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An engine section and a section that raises, through the checkpointed battery."""
    from peekaudit.conformance import Battery

    calls = []
    tiny = Battery("tiny", ("normal",), (1000,), (0.5,), ("daily14",), 200)

    class Engine(living.EngineMethod):
        def evaluate(self, stats, ctx):  # type: ignore[no-untyped-def]
            calls.append(1)
            return self.inner.evaluate(stats, ctx)

    def section():  # type: ignore[no-untyped-def]
        from peekaudit.conformance import check_conformance
        from peekaudit.methods.registry import get_method

        rep = check_conformance(Engine(get_method("ref.msprt")), tiny)
        return {
            f"fpr/{c.unit}": {"verdict": "violates" if c.violation else "ok", "detail": ""}
            for c in rep.cells
        }

    def boom():  # type: ignore[no-untyped-def]
        raise ImportError("engine moved a module")

    monkeypatch.setattr(living, "sections", lambda t: {"fpr": section, "mirror": boom})
    res = living.battery("gbstats", tmp_path)
    assert res["verdicts"] == {"fpr/normal-n1000-b50-aa": "ok", "mirror/error": "ERROR"}
    assert "ImportError: engine moved a module" in res["details"]["mirror/error"]
    n_calls = len(calls)
    assert living.battery("gbstats", tmp_path)["verdicts"] == res["verdicts"]
    assert len(calls) == n_calls  # finished sections are read back, not rerun


def test_scrub_removes_machine_paths() -> None:
    home = str(Path.home())
    out = living.scrub({"a": [f"{sys.prefix}/x.py", f"{home}\\y"]})
    assert home not in json.dumps(out) and sys.prefix not in json.dumps(out)


def test_living_job_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    script = _script()
    baseline = {
        "targets": {
            "tea-tasting": {
                "environment": {"tea-tasting": "2.0.0"},
                "verdicts": {"coverage/prop_pooled": "violates", "fpr/welch_t": "ok"},
                "details": {},
            }
        }
    }
    base_path = tmp_path / "baseline.json"
    base_path.write_text(json.dumps(baseline))
    monkeypatch.setattr(script, "BASELINE", base_path)
    out = tmp_path / "run"

    def finished_run(verdict: str) -> None:
        (out / "results").mkdir(parents=True, exist_ok=True)
        (out / "run.json").write_text(json.dumps({"targets": ["tea-tasting"], "pins": []}))
        result = {
            "target": "tea-tasting",
            "environment": {"tea-tasting": "2.1.0"},
            "verdicts": {"coverage/prop_pooled": verdict, "fpr/welch_t": "ok"},
            "details": {"coverage/prop_pooled": "covers 0.951"},
        }
        (out / "results" / "tea-tasting.json").write_text(json.dumps(result))

    def no_install(*a, **k):  # type: ignore[no-untyped-def]
        raise AssertionError("a finished target must not be reinstalled")

    monkeypatch.setattr(script, "run_target", no_install)
    gh_state = tmp_path / "gh.json"
    fake = tmp_path / "fake_gh.py"
    fake.write_text(FAKE_GH)
    monkeypatch.setenv("LIVING_GH", json.dumps([sys.executable, str(fake), str(gh_state)]))

    # Same verdicts as the baseline: exit 0 and no issue.
    finished_run("violates")
    assert script.main(["--out", str(out), "--targets", "tea-tasting"]) == 0
    assert script.main(["notify", str(out), "--repo", "o/r"]) == 0
    assert not gh_state.exists()

    # The release fixed the pooled interval: exit 1, an issue, then a comment on it.
    finished_run("ok")
    assert script.main(["--out", str(out), "--targets", "tea-tasting"]) == 1
    summary = (out / "summary.md").read_text("utf-8")
    assert "`coverage/prop_pooled`: violates -> ok" in summary
    assert "| tea-tasting | 2.0.0 | 2.1.0 | 2 | 0 | 1 |" in summary
    assert script.main(["notify", str(out), "--repo", "o/r", "--run-url", "U"]) == 0
    assert "Run and artifact: U" in (out / "issue.md").read_text("utf-8")
    assert script.main(["notify", str(out), "--repo", "o/r"]) == 0
    calls = json.loads(gh_state.read_text())
    verbs = [c[:2] for c in calls]
    assert verbs == [
        ["issue", "list"],
        ["issue", "create"],
        ["issue", "list"],
        ["issue", "comment"],
    ]
    assert calls[1][calls[1].index("--title") + 1] == script.ISSUE_TITLE
    assert calls[3][2] == "7"


def test_bad_pin_is_refused_before_anything_runs(tmp_path: Path) -> None:
    script = _script()
    out = tmp_path / "run"
    assert (
        script.main(["--out", str(out), "--targets", "tea-tasting", "--pin", "gbstats==0.8.0"]) == 2
    )
    assert not out.exists()
