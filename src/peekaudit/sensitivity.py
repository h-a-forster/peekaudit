"""Minimum-arm-size sensitivity: how much of each false-positive rate comes from tiny arms?

Many FPR excesses in the study arise at the first looks of small experiments, where one arm
holds a handful of users: a zero-variance arm gives a standard error of 0, or a z-statistic
rests on a few observations, and the engine reports significance. This module recomputes,
for the A/A units of a finished run, the FPR when rejections are only counted at looks where
both arms hold at least ``min_arm`` units (the analysis is unchanged; earlier results are
ignored, as a dashboard with a minimum sample size would hide them). GrowthBook's product
layer uses 150 (DESIGN 4.3); the engine libraries audited here apply no minimum.

The data are regenerated from the run's config and seeds, so the unguarded rejection count
of every cell must equal the run's own count; cells that do not match are flagged and the
report does not use them. One JSON checkpoint per unit under
``<run>/sensitivity/min_arm<k>/``, a progress log next to them, at most 4 worker processes.
"""

from __future__ import annotations

import json
import os
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from peekaudit.methods.base import EvalContext, MethodOutput, stopping_outcome
from peekaudit.methods.registry import all_methods
from peekaudit.runner import StudyConfig, effective_workers, read_unit, unit_complete
from peekaudit.scenarios import Unit, arm_counts, unit_from_id
from peekaudit.simulate import chunk_sizes, generate_chunk

#: GrowthBook's product-layer minimum sample size (packages/shared/src/constants.ts).
MIN_ARM = 150
#: Methods with a sequential claim (DESIGN 6.3).
SEQUENTIAL_METHODS: tuple[str, ...] = (
    "gbstats.seq.rel",
    "gbstats.seq.rel.tuned",
    "gbstats.seq.abs",
    "spotify.gst",
    "ref.msprt",
    "ref.acs",
    "ref.gst_obf",
)


def guarded_looks(totals: np.ndarray, pi: float, min_arm: int) -> np.ndarray:
    """Boolean mask of looks at which both arms hold at least ``min_arm`` units."""
    n_a, n_b = arm_counts(np.asarray(totals), pi)
    return np.minimum(n_a, n_b) >= min_arm


def guarded_counts(out: MethodOutput, mask: np.ndarray) -> dict[str, int]:
    """Replicates rejecting at any look, and at any look inside ``mask``."""
    rej = out.reject != 0
    return {
        "rejections": int((stopping_outcome(out)["sign"] != 0).sum()),
        "rejections_guarded": int(rej[:, mask].any(axis=1).sum()),
    }


def sensitivity_dir(run_dir: str | Path, min_arm: int = MIN_ARM) -> Path:
    return Path(run_dir) / "sensitivity" / f"min_arm{min_arm}"


def unit_sensitivity(
    unit: Unit, cfg: StudyConfig, methods: tuple[str, ...], min_arm: int
) -> dict[str, Any]:
    """Guarded and unguarded rejection counts for every claimed (method, schedule) cell."""
    metric = unit.metric
    ms = [m for m in all_methods(methods) if m.supports(metric)]
    contexts = {s: EvalContext.for_schedule(metric, unit.n_max, s, unit.tau) for s in cfg.schedules}
    masks = {s: guarded_looks(ctx.totals, unit.pi, min_arm) for s, ctx in contexts.items()}
    counts: dict[str, dict[str, dict[str, int]]] = {
        m.name: {
            s: {"rejections": 0, "rejections_guarded": 0} for s in cfg.schedules if m.claims(s)
        }
        for m in ms
    }
    reps = cfg.r_aa if unit.hypothesis == "aa" else cfg.r_ab
    for c, size in enumerate(chunk_sizes(reps, cfg.chunk)):
        data = generate_chunk(unit, c, size, cfg.master_seed)
        for s, ctx in contexts.items():
            sub = data.take(ctx.cols)
            for m in ms:
                if s not in counts[m.name]:
                    continue
                got = guarded_counts(m.evaluate(sub, ctx), masks[s])
                for key, v in got.items():
                    counts[m.name][s][key] += v
    cells: dict[str, dict[str, dict[str, Any]]] = {}
    for name, by_s in counts.items():
        cells[name] = {}
        for s, cnt in by_s.items():
            cells[name][s] = {
                "n": reps,
                **cnt,
                "guarded_looks": int(masks[s].sum()),
                "looks": int(masks[s].size),
            }
    return {
        "unit": unit.id,
        "config_hash": cfg.hash(),
        "min_arm": min_arm,
        "methods": list(methods),
        "cells": cells,
    }


def _check_against_run(result: dict[str, Any], run_unit: dict[str, Any] | None) -> bool:
    """True when every unguarded count equals the run's rejection count for the cell."""
    if run_unit is None:
        return False
    ok = True
    for name, by_s in result["cells"].items():
        for s, c in by_s.items():
            ref = run_unit["cells"].get(name, {}).get(s)
            match = ref is not None and ref["rejections"] == c["rejections"]
            c["matches_run"] = match
            ok &= match
    return ok


def _work(unit_id: str, cfg_dict: dict[str, Any], run_dir: str, methods: list[str], k: int) -> str:
    cfg = StudyConfig.from_dict(cfg_dict)
    res = unit_sensitivity(unit_from_id(unit_id), cfg, tuple(methods), k)
    res["matches_run"] = _check_against_run(res, read_unit(Path(run_dir), unit_id))
    dest = sensitivity_dir(run_dir, k) / f"{unit_id}.json"
    tmp = dest.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(res, indent=1, sort_keys=True), encoding="utf-8")
    os.replace(tmp, dest)
    return unit_id


def _done(path: Path, cfg: StudyConfig, methods: tuple[str, ...], min_arm: int) -> bool:
    if not path.exists():
        return False
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        d.get("config_hash") == cfg.hash()
        and d.get("min_arm") == min_arm
        and d.get("methods") == list(methods)
    )


def run_sensitivity(
    run_dir: str | Path,
    min_arm: int = MIN_ARM,
    methods: tuple[str, ...] = SEQUENTIAL_METHODS,
    workers: int | None = None,
    echo: Callable[[str], None] = print,
) -> list[str]:
    """Compute (or resume) the sensitivity checkpoints for the completed A/A units of a run."""
    if min_arm < 1:
        raise ValueError("min_arm must be >= 1")
    run = Path(run_dir)
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    cfg = StudyConfig.from_dict(manifest["config"])
    methods = tuple(m for m in methods if m in cfg.methods)
    out = sensitivity_dir(run, min_arm)
    out.mkdir(parents=True, exist_ok=True)
    units = [
        u
        for u in cfg.unit_list()
        if u.hypothesis == "aa" and unit_complete(read_unit(run, u.id), cfg)
    ]
    todo = [u for u in units if not _done(out / f"{u.id}.json", cfg, methods, min_arm)]
    log_path = out / "progress.log"

    def log(line: str) -> None:
        stamp = datetime.now(UTC).isoformat(timespec="seconds")
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"{stamp} {line}\n")
        echo(line)

    n_workers = effective_workers(workers)
    log(f"start min_arm={min_arm} units={len(units)} todo={len(todo)} workers={n_workers}")
    t0 = time.perf_counter()
    done: list[str] = []
    args = (cfg.to_dict(), str(run), list(methods), min_arm)

    def finished(uid: str) -> None:
        done.append(uid)
        el = time.perf_counter() - t0
        log(f"[{len(done)}/{len(todo)}] done {uid} | elapsed {el / 60:.1f} min")

    if n_workers == 1 or len(todo) <= 1:
        for u in todo:
            finished(_work(u.id, *args))
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = [ex.submit(_work, u.id, *args) for u in todo]
            for fut in as_completed(futs):
                finished(fut.result())
    log(f"end completed={len(done)} wall={(time.perf_counter() - t0) / 60:.1f} min")
    return done


def load_sensitivity(run_dir: str | Path, min_arm: int = MIN_ARM) -> list[dict[str, Any]]:
    """Checkpoints of a finished sensitivity run (empty if none)."""
    d = sensitivity_dir(run_dir, min_arm)
    if not d.exists():
        return []
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted(d.glob("*.json"))]
