"""Study runner: work units, a capped process pool, per-unit checkpoints and a progress log.

Layout of ``runs/<study>/`` (DESIGN 7):

- ``manifest.json``: config, config hash, versions, git SHA, start time.
- ``units/<unit_id>.json``: aggregates for every (method, schedule) cell, engine
  verification counts and identity checks. Written atomically (temp file + ``os.replace``).
  The content depends only on the config and the code, never on timing or worker count.
- ``units/<unit_id>.timing.json``: wall time of the unit (kept out of the unit JSON).
- ``cache/<unit_id>.npz``: per-replicate stopping look, sign and interval at the stop.
- ``progress.log``: one line per finished unit with done/total, elapsed, rate and ETA.

A unit is skipped when its JSON exists with the same config hash and peekaudit version.
"""

from __future__ import annotations

import fnmatch
import hashlib
import json
import os
import subprocess
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import norm

from peekaudit import __version__
from peekaudit.methods.base import EvalContext, Method, MethodOutput, stopping_outcome
from peekaudit.methods.registry import METHOD_NAMES, all_methods
from peekaudit.scenarios import (
    MASTER_SEED,
    SCHEDULES,
    Unit,
    fine_grid,
    study_units,
    unit_from_id,
)
from peekaudit.simulate import DEFAULT_CHUNK, chunk_sizes, generate_chunk
from peekaudit.verify import TIE_TOL, Agreement, compare

#: Hard cap on worker processes (shared machine; DESIGN 7).
MAX_WORKERS = 4
_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


@dataclass(frozen=True)
class StudyConfig:
    """Replicate counts, verification depth and the method/schedule grid of a study."""

    name: str
    r_aa: int
    r_ab: int
    k_verify: int
    k_verify_slow: int  # spotify.gst on stream140 (about 8 s per engine call)
    chunk: int = DEFAULT_CHUNK
    master_seed: int = MASTER_SEED
    methods: tuple[str, ...] = METHOD_NAMES
    schedules: tuple[str, ...] = tuple(SCHEDULES)
    units: tuple[str, ...] = ()  # empty = the full study grid

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["methods"] = list(self.methods)
        d["schedules"] = list(self.schedules)
        d["units"] = list(self.units)
        return d

    def hash(self) -> str:
        """Hash of everything that changes unit results (the name does not)."""
        d = self.to_dict()
        d.pop("name")
        d.pop("units")
        return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:16]

    @staticmethod
    def from_dict(d: dict[str, Any]) -> StudyConfig:
        d = dict(d)
        for key in ("methods", "schedules", "units"):
            if key in d:
                d[key] = tuple(d[key])
        return StudyConfig(**d)

    def unit_list(self) -> list[Unit]:
        if not self.units:
            return study_units()
        return [unit_from_id(u) for u in self.units]


CONFIGS: dict[str, StudyConfig] = {
    "quick": StudyConfig("quick", r_aa=4000, r_ab=2000, k_verify=10, k_verify_slow=2),
    "full": StudyConfig("full", r_aa=20000, r_ab=10000, k_verify=20, k_verify_slow=4),
    "smoke": StudyConfig("smoke", r_aa=500, r_ab=250, k_verify=2, k_verify_slow=1),
}


def load_config(spec: str) -> StudyConfig:
    """A named config (quick, full, smoke) or a path to a JSON config."""
    if spec in CONFIGS:
        return CONFIGS[spec]
    return StudyConfig.from_dict(json.loads(Path(spec).read_text(encoding="utf-8")))


# ------------------------------------------------------------------------------- one unit


def _k_for(method: Method, schedule: str, cfg: StudyConfig) -> int:
    if method.engine is None:
        return 0
    if method.name == "spotify.gst" and schedule == "stream140":
        return cfg.k_verify_slow
    return cfg.k_verify


def _aggregate(
    outcome: dict[str, np.ndarray], totals: np.ndarray, truth: float, n_looks: int
) -> dict[str, Any]:
    sign = outcome["sign"]
    stop = outcome["stop"].astype(np.int64)
    lo, hi, nr = outcome["lower"], outcome["upper"], outcome["no_result"]
    has_int = ~nr & ~np.isnan(lo) & ~np.isnan(hi)
    hits = has_int & (lo <= truth) & (truth <= hi)
    stop_n = totals[stop].astype(np.float64)
    first = np.bincount(stop[sign != 0], minlength=n_looks)
    return {
        "n": int(sign.size),
        "rejections": int((sign != 0).sum()),
        "reject_pos": int((sign > 0).sum()),
        "reject_neg": int((sign < 0).sum()),
        "no_result_at_stop": int(nr.sum()),
        "cover_n": int(has_int.sum()),
        "cover_hits": int(hits.sum()),
        "truth": truth,
        "stop_n_mean": float(stop_n.mean()),
        "stop_n_sd": float(stop_n.std(ddof=1)) if stop_n.size > 1 else 0.0,
        "first_rejection_hist": first.tolist(),
    }


def _identity_counts(a: np.ndarray, b: np.ndarray, tie: np.ndarray) -> dict[str, int]:
    mism = a != b
    return {
        "cells": int(a.size),
        "mismatches": int((mism & ~tie).sum()),
        "ties": int(tie.sum()),
    }


def run_unit(unit: Unit, cfg: StudyConfig) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Simulate one unit: every supported method under every schedule, plus engine checks.

    Returns the unit summary (JSON-serialisable, deterministic) and per-replicate arrays.
    """
    from peekaudit.methods.gbstats_adapter import GbstatsSequential, _diff_var
    from peekaudit.methods.reference import ACS

    metric = unit.metric
    methods = [m for m in all_methods(cfg.methods) if m.supports(metric)]
    reps = cfg.r_aa if unit.hypothesis == "aa" else cfg.r_ab
    sizes = chunk_sizes(reps, cfg.chunk)
    keys = [(m.name, s) for m in methods for s in cfg.schedules]
    parts: dict[tuple[str, str], list[dict[str, np.ndarray]]] = {k: [] for k in keys}
    verification: dict[tuple[str, str], Agreement] = {}
    ident_counts: dict[str, dict[str, int]] = {
        "bayes_vs_z10": {"cells": 0, "mismatches": 0, "ties": 0},
        "acs_vs_gbstats_abs_nmax": {"cells": 0, "mismatches": 0, "ties": 0},
    }
    contexts = {s: EvalContext.for_schedule(metric, unit.n_max, s, unit.tau) for s in cfg.schedules}
    acs_ref, gb_abs = ACS(), GbstatsSequential(relative=False, tuned=True)
    bayes = next((m for m in methods if m.name == "gbstats.bayes.ctw95"), None)
    z10 = float(norm.ppf(0.95))  # chance to win 95% <=> |z| >= z_0.95 (level 0.10, two-sided)

    for c, size in enumerate(sizes):
        data = generate_chunk(unit, c, size, cfg.master_seed)
        # Identity checks on the whole fine grid (DESIGN 5.3).
        grid_ctx = EvalContext(metric, unit.n_max, "grid", np.arange(data.n_looks), unit.tau)
        if bayes is not None:
            out = bayes.evaluate(data, grid_ctx)
            d, v, default, _ = _diff_var(data, relative=True)
            with np.errstate(divide="ignore", invalid="ignore"):
                z = d / np.sqrt(v)
                expected = np.where(~default & (np.abs(z) >= z10), np.sign(z), 0)
                tie = np.abs(np.abs(z) - z10) < TIE_TOL * z10
            cnt = _identity_counts(out.reject, expected, tie)
            for key in cnt:
                ident_counts["bayes_vs_z10"][key] += cnt[key]
        if not metric.binary:
            a = acs_ref.evaluate(data, grid_ctx)
            b = gb_abs.evaluate(data, grid_ctx)
            ok = ~a.no_result & ~b.no_result
            tie = np.abs(a.margin) < TIE_TOL
            cnt = _identity_counts(a.reject[ok], b.reject[ok], tie[ok])
            for key in cnt:
                ident_counts["acs_vs_gbstats_abs_nmax"][key] += cnt[key]

        for s in cfg.schedules:
            ctx = contexts[s]
            sub = data.take(ctx.cols)
            for m in methods:
                out: MethodOutput = m.evaluate(sub, ctx)
                parts[(m.name, s)].append(stopping_outcome(out))
                k = min(_k_for(m, s, cfg), size) if c == 0 else 0
                if k > 0:
                    rows = slice(0, k)
                    engine = m.evaluate_engine(sub.rows(rows), ctx)
                    verification[(m.name, s)] = compare(out.rows(rows), engine)

    cells: dict[str, dict[str, Any]] = {}
    arrays: dict[str, np.ndarray] = {}
    for m in methods:
        truth = unit.true_effect(m.scale)
        cells[m.name] = {}
        for s in cfg.schedules:
            outcome = {
                f: np.concatenate([p[f] for p in parts[(m.name, s)]])
                for f in ("stop", "sign", "lower", "upper", "no_result")
            }
            totals = contexts[s].totals
            agg = _aggregate(outcome, totals, truth, len(totals))
            agg["claimed"] = m.claims(s)
            if (m.name, s) in verification:
                agg["verification"] = verification[(m.name, s)].to_dict()
            cells[m.name][s] = agg
            for f, arr in outcome.items():
                arrays[f"{m.name}|{s}|{f}"] = arr
    verified_ok = all(a.decision_mismatches == 0 for a in verification.values())
    summary = {
        "unit": unit.id,
        "cell": unit.cell,
        "metric": unit.metric_id,
        "n_max": unit.n_max,
        "pi": unit.pi,
        "hypothesis": unit.hypothesis,
        "lift": unit.lift,
        "design_lift": unit.design_lift,
        "tau": unit.tau,
        "reps": reps,
        "config_hash": cfg.hash(),
        "peekaudit_version": __version__,
        "methods": {
            m.name: {"scale": m.scale, "claim": m.claim, "label": m.label} for m in methods
        },
        "totals": {s: contexts[s].totals.tolist() for s in cfg.schedules},
        "cells": cells,
        "identities": ident_counts,
        "verified_ok": verified_ok,
    }
    return summary, arrays


# --------------------------------------------------------------------------- checkpoints


def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def unit_path(out: Path, unit_id: str) -> Path:
    return out / "units" / f"{unit_id}.json"


def is_done(out: Path, unit: Unit, cfg: StudyConfig) -> bool:
    p = unit_path(out, unit.id)
    if not p.exists():
        return False
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return d.get("config_hash") == cfg.hash() and d.get("peekaudit_version") == __version__


def _run_and_save(unit_id: str, cfg_dict: dict[str, Any], out_str: str) -> tuple[str, float, bool]:
    t0 = time.perf_counter()
    cfg = StudyConfig.from_dict(cfg_dict)
    out = Path(out_str)
    summary, arrays = run_unit(unit_from_id(unit_id), cfg)
    cache = out / "cache" / f"{unit_id}.npz"
    tmp_cache = cache.with_name(cache.stem + ".tmp.npz")
    np.savez_compressed(tmp_cache, **arrays)
    os.replace(tmp_cache, cache)
    _atomic_write(unit_path(out, unit_id), json.dumps(summary, indent=1, sort_keys=True))
    elapsed = time.perf_counter() - t0
    _atomic_write(
        out / "units" / f"{unit_id}.timing.json", json.dumps({"seconds": round(elapsed, 2)})
    )
    return unit_id, elapsed, bool(summary["verified_ok"])


def _git_sha() -> str:
    try:
        root = Path(__file__).resolve().parents[2]
        r = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, timeout=10
        )
        sha = r.stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=root,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        return (sha + ("-dirty" if dirty else "")) if sha else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def effective_workers(requested: int | None) -> int:
    """Worker count after the hard cap (never more than MAX_WORKERS or the CPU count)."""
    cpu = os.cpu_count() or 1
    req = MAX_WORKERS if requested is None else requested
    return max(1, min(MAX_WORKERS, cpu, req))


def _cost(unit: Unit) -> float:
    """Rough relative cost, used only to start expensive units first."""
    spotify = 0.0 if unit.metric.ratio else 1.0
    gen = {"ratio": 3.0, "lognormal": 1.5, "pareto": 1.5}.get(unit.metric.kind, 0.2)
    size = 1.0 if unit.hypothesis == "aa" else 0.5
    return 100.0 * spotify + gen * size * unit.n_max / 1000.0


@dataclass
class RunResult:
    total: int
    skipped: int
    completed: list[str] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)
    unverified: list[str] = field(default_factory=list)
    seconds: float = 0.0


def run_study(
    cfg: StudyConfig,
    out_dir: str | Path,
    workers: int | None = None,
    units_glob: str | None = None,
    force: bool = False,
    echo: Callable[[str], None] = print,
) -> RunResult:
    """Run (or resume) a study. Completed units are skipped unless ``force``."""
    out = Path(out_dir)
    (out / "units").mkdir(parents=True, exist_ok=True)
    (out / "cache").mkdir(exist_ok=True)
    from peekaudit.cli import environment_report

    manifest_path = out / "manifest.json"
    manifest = {
        "study": cfg.name,
        "config": cfg.to_dict(),
        "config_hash": cfg.hash(),
        "versions": environment_report(),
        "git_sha": _git_sha(),
        "started": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old.get("config_hash") != cfg.hash() and not force:
            raise RuntimeError(
                f"{out} holds a study with a different config (hash {old.get('config_hash')});"
                " use --force or another --out"
            )
        manifest["first_started"] = old.get("first_started", old.get("started"))
    else:
        manifest["first_started"] = manifest["started"]
    _atomic_write(manifest_path, json.dumps(manifest, indent=1, sort_keys=True))

    units = cfg.unit_list()
    if units_glob:
        units = [u for u in units if fnmatch.fnmatch(u.id, units_glob)]
    todo = [u for u in units if force or not is_done(out, u, cfg)]
    todo.sort(key=_cost, reverse=True)
    res = RunResult(total=len(units), skipped=len(units) - len(todo))
    n_workers = effective_workers(workers)
    if workers is not None and workers > n_workers:
        echo(f"workers capped at {n_workers} (requested {workers})")
    log_path = out / "progress.log"

    def log(line: str) -> None:
        stamp = datetime.now(UTC).isoformat(timespec="seconds")
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"{stamp} {line}\n")
        echo(line)

    log(
        f"start study={cfg.name} hash={cfg.hash()} units={len(units)} todo={len(todo)} "
        f"skipped={res.skipped} workers={n_workers}"
    )
    t0 = time.perf_counter()
    done = res.skipped

    def finished(uid: str, secs: float, ok: bool) -> None:
        nonlocal done
        done += 1
        res.completed.append(uid)
        if not ok:
            res.unverified.append(uid)
        el = time.perf_counter() - t0
        n_new = len(res.completed) + len(res.failed)
        rate = n_new / el * 60 if el > 0 else 0.0
        remaining = len(units) - done - len(res.failed)
        eta = remaining / rate if rate > 0 else float("nan")
        flag = "" if ok else " ENGINE-MISMATCH"
        log(
            f"[{done}/{len(units)}] done {uid} in {secs:.1f}s{flag} | elapsed {el / 60:.1f} min"
            f" | rate {rate:.2f} units/min | ETA {eta:.1f} min"
        )

    cfg_dict = cfg.to_dict()
    if n_workers == 1 or len(todo) <= 1:
        for u in todo:
            try:
                finished(*_run_and_save(u.id, cfg_dict, str(out)))
            except Exception as exc:  # noqa: BLE001 - logged and reported, run continues
                res.failed.append(u.id)
                log(f"FAILED {u.id}: {type(exc).__name__}: {exc}")
    else:
        for var in _THREAD_VARS:  # inherited by spawned workers before they import numpy
            os.environ[var] = "1"
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = {ex.submit(_run_and_save, u.id, cfg_dict, str(out)): u for u in todo}
            for fut in as_completed(futs):
                u = futs[fut]
                try:
                    finished(*fut.result())
                except Exception as exc:  # noqa: BLE001
                    res.failed.append(u.id)
                    log(f"FAILED {u.id}: {type(exc).__name__}: {exc}")
    res.seconds = time.perf_counter() - t0
    log(
        f"end completed={len(res.completed)} failed={len(res.failed)} "
        f"unverified={len(res.unverified)} wall={res.seconds / 60:.1f} min"
    )
    return res


def study_status(out_dir: str | Path) -> dict[str, Any]:
    """Done/total, failed and unverified units, and the last progress line."""
    out = Path(out_dir)
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    cfg = StudyConfig.from_dict(manifest["config"])
    units = cfg.unit_list()
    done, unverified = [], []
    for u in units:
        if is_done(out, u, cfg):
            done.append(u.id)
            d = json.loads(unit_path(out, u.id).read_text(encoding="utf-8"))
            if not d.get("verified_ok", False):
                unverified.append(u.id)
    log = out / "progress.log"
    lines = log.read_text(encoding="utf-8").splitlines() if log.exists() else []
    failed = sorted({ln.split("FAILED ", 1)[1].split(":")[0] for ln in lines if "FAILED " in ln})
    failed = [f for f in failed if f not in done]
    secs = []
    for u in done:
        t = out / "units" / f"{u}.timing.json"
        if t.exists():
            secs.append(json.loads(t.read_text(encoding="utf-8"))["seconds"])
    return {
        "study": manifest["study"],
        "config_hash": manifest["config_hash"],
        "done": len(done),
        "total": len(units),
        "failed": failed,
        "unverified": unverified,
        "cpu_minutes_done": round(sum(secs) / 60, 1),
        "last_log": lines[-1] if lines else "",
    }


def grid_totals(n_max: int, schedule: str) -> np.ndarray:
    return fine_grid(n_max)[SCHEDULES[schedule]]
