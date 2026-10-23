"""Runner: checkpoints, resume, determinism across worker counts, worker cap, status."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from peekaudit.runner import (
    MAX_WORKERS,
    StudyConfig,
    effective_workers,
    is_done,
    load_config,
    run_study,
    run_unit,
    study_status,
    unit_path,
)
from peekaudit.scenarios import unit_from_id

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

TINY = StudyConfig(
    "tiny",
    r_aa=300,
    r_ab=150,
    k_verify=1,
    k_verify_slow=0,
    chunk=100,
    methods=("ref.naive_z", "ref.acs", "gbstats.seq.rel", "gbstats.bayes.ctw95"),
    schedules=("fixed1", "daily14"),
    units=("normal-n1000-b50-aa", "bern_p10-n1000-b20-ab", "ratio-n1000-b50-aa"),
)


def _quiet(_: str) -> None:
    pass


def test_run_unit_summary_shape() -> None:
    s, arrays = run_unit(unit_from_id("bern_p10-n1000-b20-ab"), TINY)
    assert s["reps"] == 150 and s["verified_ok"]
    cell = s["cells"]["ref.naive_z"]["fixed1"]
    assert cell["n"] == 150 and cell["claimed"]
    assert cell["rejections"] == cell["reject_pos"] + cell["reject_neg"]
    assert not s["cells"]["ref.naive_z"]["daily14"]["claimed"]
    assert s["cells"]["gbstats.seq.rel"]["daily14"]["verification"]["cells"] == 14
    assert s["identities"]["bayes_vs_z10"]["mismatches"] == 0
    assert arrays["ref.acs|daily14|stop"].shape == (150,)
    # A/B power of the calibrated z-test at N_max is in the right region.
    assert 0.6 < cell["reject_pos"] / cell["n"] < 0.95


def test_checkpoint_resume_and_determinism(tmp_path: Path) -> None:
    a, b = tmp_path / "a", tmp_path / "b"
    first = run_study(TINY, a, workers=1, units_glob="normal-*", echo=_quiet)
    assert first.completed == ["normal-n1000-b50-aa"]
    mtime = unit_path(a, "normal-n1000-b50-aa").stat().st_mtime_ns
    second = run_study(TINY, a, workers=1, echo=_quiet)  # resume: one skipped, two new
    assert second.skipped == 1 and len(second.completed) == 2
    assert unit_path(a, "normal-n1000-b50-aa").stat().st_mtime_ns == mtime
    # Simulate a crash that lost one unit: only that unit is recomputed.
    unit_path(a, "ratio-n1000-b50-aa").unlink()
    third = run_study(TINY, a, workers=1, echo=_quiet)
    assert third.completed == ["ratio-n1000-b50-aa"] and third.skipped == 2
    # A fresh run with a process pool gives byte-identical unit JSON.
    run_study(TINY, b, workers=2, echo=_quiet)
    for uid in TINY.units:
        assert unit_path(a, uid).read_bytes() == unit_path(b, uid).read_bytes()
    st = study_status(a)
    assert st["done"] == 3 and st["total"] == 3 and not st["failed"]
    log = (a / "progress.log").read_text(encoding="utf-8")
    assert "[3/3]" in log and "ETA" in log
    manifest = json.loads((a / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["config_hash"] == TINY.hash() and "gbstats" in manifest["versions"]


def test_config_change_is_refused_and_invalidates(tmp_path: Path) -> None:
    run_study(TINY, tmp_path, workers=1, units_glob="normal-*", echo=_quiet)
    other = StudyConfig.from_dict({**TINY.to_dict(), "r_aa": 200})
    with pytest.raises(RuntimeError, match="different config"):
        run_study(other, tmp_path, workers=1, echo=_quiet)
    assert not is_done(tmp_path, unit_from_id("normal-n1000-b50-aa"), other)


def test_worker_cap() -> None:
    assert MAX_WORKERS == 4
    assert effective_workers(16) <= 4
    assert effective_workers(None) <= 4
    assert effective_workers(0) == 1


def test_named_configs() -> None:
    q, f = load_config("quick"), load_config("full")
    assert (q.r_aa, q.r_ab, q.k_verify) == (4000, 2000, 10)
    assert (f.r_aa, f.r_ab) == (20000, 10000)
    assert q.hash() != f.hash()
    assert len(q.unit_list()) == 52
    assert StudyConfig.from_dict(q.to_dict()) == q


def test_failed_unit_is_logged_and_run_continues(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import peekaudit.runner as runner

    real = runner.run_unit

    def boom(unit, cfg):  # type: ignore[no-untyped-def]
        if unit.id.startswith("normal"):
            raise RuntimeError("simulated failure")
        return real(unit, cfg)

    monkeypatch.setattr(runner, "run_unit", boom)
    res = run_study(TINY, tmp_path, workers=1, units_glob="*-aa", echo=_quiet)
    assert res.failed == ["normal-n1000-b50-aa"] and res.completed == ["ratio-n1000-b50-aa"]
    assert "FAILED normal-n1000-b50-aa" in (tmp_path / "progress.log").read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"r_aa": 0}, "r_aa"),
        ({"chunk": 0}, "chunk"),
        ({"k_verify": -1}, "k_verify"),
        ({"methods": ["no.such.method"]}, "unknown methods"),
        ({"methods": []}, "methods must not be empty"),
        ({"schedules": ["nope"]}, "unknown schedules"),
        ({"units": ["bogus"]}, "not a unit id"),
        ({"methods": "ref.acs"}, "must be a list"),
        ({"extra_key": 1}, "unknown keys"),
    ],
)
def test_invalid_config_is_rejected_before_anything_is_written(
    tmp_path: Path, change: dict[str, object], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        cfg = StudyConfig.from_dict({**TINY.to_dict(), **change})
        run_study(cfg, tmp_path / "out", workers=1, echo=_quiet)
    assert not (tmp_path / "out").exists()


def test_missing_config_key_is_rejected(tmp_path: Path) -> None:
    d = TINY.to_dict()
    del d["r_ab"]
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(d), encoding="utf-8")
    with pytest.raises(ValueError, match="missing keys"):
        load_config(str(path))


def test_one_completeness_rule_and_provenance(tmp_path: Path) -> None:
    from peekaudit.report import load_units, provenance_groups

    run_study(TINY, tmp_path, workers=1, units_glob="normal-*", echo=_quiet)
    uid = "normal-n1000-b50-aa"
    tpath = tmp_path / "units" / f"{uid}.timing.json"
    timing = json.loads(tpath.read_text(encoding="utf-8"))
    assert timing["provenance"]["git_sha"] and "gbstats" in timing["provenance"]["versions"]
    # A unit written by another package version is still complete for run, status and report.
    p = unit_path(tmp_path, uid)
    d = json.loads(p.read_text(encoding="utf-8"))
    d["peekaudit_version"] = "0.0.0-other"
    p.write_text(json.dumps(d), encoding="utf-8")
    assert is_done(tmp_path, unit_from_id(uid), TINY)
    assert study_status(tmp_path)["done"] == 1
    assert [u["unit"] for u in load_units(tmp_path)[1]] == [uid]
    res = run_study(TINY, tmp_path, workers=1, units_glob="normal-*", echo=_quiet)
    assert res.skipped == 1 and not res.completed
    # Each start is kept in the manifest history; the first start stays at the top level.
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert len(manifest["starts"]) == 2
    assert manifest["started"] == manifest["starts"][0]["started"]
    # A unit computed by other code is flagged on resume and reported as its own group.
    timing["provenance"]["git_sha"] = "0123abc"
    tpath.write_text(json.dumps(timing), encoding="utf-8")
    lines: list[str] = []
    run_study(TINY, tmp_path, workers=1, units_glob="normal-*", echo=lines.append)
    assert any("computed by other code" in ln for ln in lines)
    groups = provenance_groups(tmp_path, *load_units(tmp_path))
    assert [g["git_sha"] for g in groups] == ["0123abc"]
    assert study_status(tmp_path)["units_by_git_sha"] == {"0123abc": 1}


def test_legacy_units_take_manifest_provenance_only_for_a_single_start(tmp_path: Path) -> None:
    from peekaudit.report import load_units, provenance_groups

    run_study(TINY, tmp_path, workers=1, units_glob="normal-*", echo=_quiet)
    uid = "normal-n1000-b50-aa"
    (tmp_path / "units" / f"{uid}.timing.json").write_text('{"seconds": 1.0}', encoding="utf-8")
    manifest, units = load_units(tmp_path)
    assert provenance_groups(tmp_path, manifest, units)[0]["source"].startswith("run manifest")
    manifest["starts"].append(manifest["starts"][0])
    assert provenance_groups(tmp_path, manifest, units)[0]["git_sha"] == "unknown"
