"""Like-for-like control for spotify's group sequential test (spotify.gst).

The study's reference GST (``ref.gst_obf``) differs from spotify.gst in two ways: it uses
O'Brien-Fleming-type spending where spotify uses power-family spending with rho = 2, and it
uses n - 1 variances, so an arm with n <= 1 (or a zero plug-in variance) gives no result.
spotify's z-test uses the binary variance p(1 - p)/n for 0/1 data, which is 0 for a
single-unit arm, so a look with 1 user in B can still give a finite (or zero) standard error
and a significant result.

This module separates the two differences on the A/A units of a finished run. For every
cell where spotify.gst claims FPR <= alpha it evaluates, on the run's own data (regenerated
from its config and seeds):

- ``spotify.gst`` (the verified mirror),
- ``ref.gst_pow2``: the same rho = 2 spending, peekaudit's own bounds, reference variance,
- ``ref.gst_obf``: the study's reference,

and counts rejections, rejections at the first look, paired discordance between spotify.gst
and ref.gst_pow2 (exact McNemar test, Holm across cells), and spotify rejections at looks
where ref.gst_pow2 has no result (an arm with n <= 1 or zero variance) or where spotify's
own standard error is 0. The unguarded spotify.gst and ref.gst_obf counts must equal the
run's; cells that do not are flagged.

One JSON checkpoint per unit under ``<run>/gst_reference/``, a progress log next to them,
completed units are skipped on a rerun, at most ``MAX_WORKERS`` worker processes.
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
from scipy import stats as st

from peekaudit.methods.base import EvalContext, MethodOutput, stopping_outcome
from peekaudit.methods.reference import GstOBF, GstPower2
from peekaudit.runner import StudyConfig, effective_workers, read_unit, unit_complete
from peekaudit.scenarios import Unit, arm_counts, unit_from_id
from peekaudit.simulate import chunk_sizes, generate_chunk
from peekaudit.stats import binom_sf_test, clopper_pearson, holm

#: Version of the checkpoint content; bump when the computed fields change.
SCHEMA = 1
METHODS = ("spotify.gst", "ref.gst_pow2", "ref.gst_obf")
ALPHA = 0.05
FWER = 0.05


def _methods() -> dict[str, Any]:
    from peekaudit.methods.spotify_adapter import SpotifyGST

    return {"spotify.gst": SpotifyGST(), "ref.gst_pow2": GstPower2(), "ref.gst_obf": GstOBF()}


def cell_counts(outs: dict[str, MethodOutput]) -> dict[str, Any]:
    """Counts for one chunk of one cell, from the three methods' outputs on the same data."""
    stops = {name: stopping_outcome(o) for name, o in outs.items()}
    res: dict[str, Any] = {"n": int(next(iter(stops.values()))["sign"].size)}
    for name, so in stops.items():
        rej = so["sign"] != 0
        res[name] = {
            "rejections": int(rej.sum()),
            "rejections_look1": int((rej & (so["stop"] == 0)).sum()),
        }
    sp, pw = stops["spotify.gst"], stops["ref.gst_pow2"]
    sp_rej, pw_rej = sp["sign"] != 0, pw["sign"] != 0
    idx = np.arange(sp_rej.size)
    stop = sp["stop"].astype(np.int64)
    ref_refuses = outs["ref.gst_pow2"].no_result[idx, stop]
    zero_se = sp["lower"] == sp["upper"]
    res["paired"] = {
        "both": int((sp_rej & pw_rej).sum()),
        "spotify_only": int((sp_rej & ~pw_rej).sum()),
        "pow2_only": int((~sp_rej & pw_rej).sum()),
    }
    res["spotify_reject_where_ref_no_result"] = int((sp_rej & ref_refuses).sum())
    res["spotify_reject_zero_se"] = int((sp_rej & zero_se).sum())
    return res


def _add(acc: dict[str, Any], new: dict[str, Any]) -> dict[str, Any]:
    if not acc:
        return json.loads(json.dumps(new))
    for k, v in new.items():
        acc[k] = _add(acc[k], v) if isinstance(v, dict) else acc[k] + v
    return acc


def unit_gst_reference(unit: Unit, cfg: StudyConfig) -> dict[str, Any]:
    """Counts for every schedule where spotify.gst claims FPR <= alpha, for one A/A unit."""
    ms = _methods()
    if not ms["spotify.gst"].supports(unit.metric):
        raise ValueError(f"spotify.gst does not support {unit.metric_id}")
    schedules = [s for s in cfg.schedules if ms["spotify.gst"].claims(s)]
    contexts = {
        s: EvalContext.for_schedule(unit.metric, unit.n_max, s, unit.tau) for s in schedules
    }
    cells: dict[str, dict[str, Any]] = {s: {} for s in schedules}
    for c, size in enumerate(chunk_sizes(cfg.r_aa, cfg.chunk)):
        data = generate_chunk(unit, c, size, cfg.master_seed)
        for s, ctx in contexts.items():
            sub = data.take(ctx.cols)
            outs = {name: m.evaluate(sub, ctx) for name, m in ms.items()}
            cells[s] = _add(cells[s], cell_counts(outs))
    for s, ctx in contexts.items():
        n_a, n_b = arm_counts(ctx.totals[:1], unit.pi)
        cells[s]["look1_arms"] = [int(n_a[0]), int(n_b[0])]
    return {
        "unit": unit.id,
        "config_hash": cfg.hash(),
        "schema": SCHEMA,
        "cells": cells,
    }


def _check_against_run(result: dict[str, Any], run_unit: dict[str, Any] | None) -> bool:
    """True when the spotify.gst and ref.gst_obf counts equal the run's for every cell."""
    if run_unit is None:
        return False
    ok = True
    for s, c in result["cells"].items():
        match = True
        for name in ("spotify.gst", "ref.gst_obf"):
            ref = run_unit["cells"].get(name, {}).get(s)
            match &= ref is not None and ref["rejections"] == c[name]["rejections"]
        c["matches_run"] = match
        ok &= match
    return ok


def gst_reference_dir(run_dir: str | Path) -> Path:
    return Path(run_dir) / "gst_reference"


def _work(unit_id: str, cfg_dict: dict[str, Any], run_dir: str) -> str:
    cfg = StudyConfig.from_dict(cfg_dict)
    res = unit_gst_reference(unit_from_id(unit_id), cfg)
    res["matches_run"] = _check_against_run(res, read_unit(Path(run_dir), unit_id))
    dest = gst_reference_dir(run_dir) / f"{unit_id}.json"
    tmp = dest.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(res, indent=1, sort_keys=True), encoding="utf-8")
    os.replace(tmp, dest)
    return unit_id


def _done(path: Path, cfg: StudyConfig) -> bool:
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return d.get("config_hash") == cfg.hash() and d.get("schema") == SCHEMA


def eligible_units(cfg: StudyConfig, run: Path) -> list[Unit]:
    """Completed A/A units of the run whose metric spotify.gst supports."""
    from peekaudit.methods.spotify_adapter import SpotifyGST

    sp = SpotifyGST()
    return [
        u
        for u in cfg.unit_list()
        if u.hypothesis == "aa"
        and sp.supports(u.metric)
        and "spotify.gst" in cfg.methods
        and "ref.gst_obf" in cfg.methods
        and unit_complete(read_unit(run, u.id), cfg)
    ]


def run_gst_reference(
    run_dir: str | Path,
    workers: int | None = None,
    units: str | None = None,
    echo: Callable[[str], None] = print,
) -> list[str]:
    """Compute (or resume) the checkpoints for the eligible units of a finished run."""
    import fnmatch

    run = Path(run_dir)
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    cfg = StudyConfig.from_dict(manifest["config"])
    out = gst_reference_dir(run)
    out.mkdir(parents=True, exist_ok=True)
    todo_units = eligible_units(cfg, run)
    if units:
        todo_units = [u for u in todo_units if fnmatch.fnmatch(u.id, units)]
    todo = [u for u in todo_units if not _done(out / f"{u.id}.json", cfg)]
    log_path = out / "progress.log"

    def log(line: str) -> None:
        stamp = datetime.now(UTC).isoformat(timespec="seconds")
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"{stamp} {line}\n")
        echo(line)

    n_workers = effective_workers(workers)
    log(f"start units={len(todo_units)} todo={len(todo)} workers={n_workers}")
    t0 = time.perf_counter()
    done: list[str] = []
    args = (cfg.to_dict(), str(run))

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


def mcnemar_p(b: int, c: int) -> float:
    """Exact two-sided McNemar p-value from the discordant counts."""
    if b + c == 0:
        return 1.0
    return float(st.binomtest(b, b + c, 0.5).pvalue)


def summarise(run_dir: str | Path) -> dict[str, Any]:
    """Rows per cell with rates, 95% Clopper-Pearson intervals and Holm-adjusted McNemar p."""
    run = Path(run_dir)
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    cfg = StudyConfig.from_dict(manifest["config"])
    rows: list[dict[str, Any]] = []
    for path in sorted(gst_reference_dir(run).glob("*.json")):
        d = json.loads(path.read_text(encoding="utf-8"))
        if d.get("config_hash") != cfg.hash() or d.get("schema") != SCHEMA:
            continue
        for s, c in sorted(d["cells"].items()):
            n = c["n"]
            row: dict[str, Any] = {
                "unit": d["unit"],
                "schedule": s,
                "n": n,
                "look1_arms": c["look1_arms"],
                "matches_run": c.get("matches_run", False),
            }
            for name in METHODS:
                k = c[name]["rejections"]
                row[name] = {
                    "rejections": k,
                    "rate": k / n,
                    "rate_ci95": list(clopper_pearson(k, n)),
                    "rejections_look1": c[name]["rejections_look1"],
                }
            row["paired"] = c["paired"]
            row["mcnemar_p"] = mcnemar_p(c["paired"]["spotify_only"], c["paired"]["pow2_only"])
            row["spotify_reject_where_ref_no_result"] = c["spotify_reject_where_ref_no_result"]
            row["spotify_reject_zero_se"] = c["spotify_reject_zero_se"]
            rows.append(row)
    adj = holm([r["mcnemar_p"] for r in rows]) if rows else []
    for r, p in zip(rows, adj, strict=True):
        r["mcnemar_p_holm"] = float(p)
    # Verdicts as in the study (DESIGN 6.2): one-sided binomial test of FPR > alpha, Holm per
    # (method, family); infinite-variance metrics form the stress family.
    violations: dict[str, dict[str, int]] = {}
    for r in rows:
        metric = unit_from_id(r["unit"]).metric
        r["family"] = "fpr_primary" if metric.finite_variance else "fpr_stress"
    for name in METHODS:
        violations[name] = {}
        for fam in ("fpr_primary", "fpr_stress"):
            rs = [r for r in rows if r["family"] == fam]
            ps = holm([binom_sf_test(r[name]["rejections"], r["n"], ALPHA) for r in rs])
            for r, p in zip(rs, ps, strict=True):
                r[name]["holm_p"] = float(p)
                r[name]["violation"] = bool(p < FWER)
            violations[name][fam] = sum(r[name]["violation"] for r in rs)
            violations[name][f"{fam}_cells"] = len(rs)
    return {
        "study": manifest.get("study"),
        "config_hash": cfg.hash(),
        "r_aa": cfg.r_aa,
        "master_seed": cfg.master_seed,
        "cells": rows,
        "violations": violations,
        "all_match_run": bool(rows) and all(r["matches_run"] for r in rows),
    }


def render_markdown(summary: dict[str, Any], min_rate: float = 0.06) -> str:
    """A table of the cells where any of the three methods exceeds ``min_rate``."""

    def pct(x: float) -> str:
        return f"{100 * x:.1f}"

    def _p(p: float) -> str:
        return f"{p:.2g}" if p > 0 else "< 1e-300"

    def rate(m: dict[str, Any]) -> str:
        lo, hi = m["rate_ci95"]
        return f"{pct(m['rate'])} [{pct(lo)}, {pct(hi)}]"

    rows = [r for r in summary["cells"] if max(r[name]["rate"] for name in METHODS) >= min_rate]
    lines = [
        "# spotify.gst against like-for-like reference GSTs (A/A units)",
        "",
        f"Study `{summary['study']}`, config hash `{summary['config_hash']}`, "
        f"{summary['r_aa']:,} replicates per unit, master seed {summary['master_seed']}. "
        "Same data for all three methods. `ref.gst_pow2` uses spotify's spending function "
        "(power family, rho = 2) with reference variances (n - 1; no result for an arm "
        "with n <= 1 or zero variance). FPR in % with 95% Clopper-Pearson intervals. "
        "McNemar: exact test of spotify.gst vs ref.gst_pow2 on the same replicates, "
        "Holm-adjusted over all "
        f"{len(summary['cells'])} cells. 'ref n/a' = spotify rejections at looks where "
        "ref.gst_pow2 gives no result; 'SE 0' = spotify rejections with a zero standard "
        f"error. Rows: any method >= {pct(min_rate)}%. Run counts reproduced: "
        f"{'yes, every cell' if summary['all_match_run'] else 'NO (see JSON)'}.",
        "",
        "Violations (FPR > 5%, one-sided binomial test, Holm per method and family; "
        "primary / stress cells): "
        + "; ".join(
            f"{name} {v['fpr_primary']} of {v['fpr_primary_cells']} / "
            f"{v['fpr_stress']} of {v['fpr_stress_cells']}"
            for name, v in summary["violations"].items()
        )
        + ".",
        "",
        "| unit | schedule | look-1 arms (A, B) | spotify.gst | ref.gst_pow2 | ref.gst_obf "
        "| spotify at look 1 | spotify only / pow2 only | McNemar Holm p | ref n/a | SE 0 |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        sp = r["spotify.gst"]
        k = sp["rejections"]
        l1 = f"{100 * sp['rejections_look1'] / k:.0f}%" if k else "-"
        lines.append(
            f"| {r['unit']} | {r['schedule']} | {r['look1_arms'][0]}, {r['look1_arms'][1]} "
            f"| {rate(sp)} | {rate(r['ref.gst_pow2'])} | {rate(r['ref.gst_obf'])} | {l1} "
            f"| {r['paired']['spotify_only']} / {r['paired']['pow2_only']} "
            f"| {_p(r['mcnemar_p_holm'])} | {r['spotify_reject_where_ref_no_result']} "
            f"| {r['spotify_reject_zero_se']} |"
        )
    return "\n".join(lines) + "\n"


def write_summary(run_dir: str | Path, out_dir: str | Path) -> dict[str, Any]:
    summary = summarise(run_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "gst_reference.json").write_text(
        json.dumps(summary, indent=1, sort_keys=True), encoding="utf-8"
    )
    (out / "gst_reference.md").write_text(render_markdown(summary), encoding="utf-8")
    return summary
