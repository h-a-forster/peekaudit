"""Turn a study run into inference: exact intervals, Holm families, gates, tables, figures.

DESIGN 6: every rate is a binomial proportion over replicates with a Clopper-Pearson
interval. FPR violations are one-sided exact binomial tests of H0: FPR <= alpha, Holm-
adjusted per method within the primary family (claimed schedules, finite-variance metrics);
the infinite-variance stress cells form their own family; coverage (A/B cells within claim)
is a third family. Cells outside a method's claim are informational.
"""

from __future__ import annotations

import json
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats as st

from peekaudit.conformance import FWER, critical_rate
from peekaudit.gst import crossing_probability
from peekaudit.methods.registry import METHOD_NAMES
from peekaudit.runner import StudyConfig, read_unit, unit_complete, unit_provenance
from peekaudit.scenarios import ALPHA, METRICS, arm_counts
from peekaudit.sensitivity import MIN_ARM, load_sensitivity
from peekaudit.stats import binom_sf_test, clopper_pearson, holm

#: Nominal fixed-look level of each method (Bayes ctw95 is a level-0.10 test, DESIGN 4.1).
FIXED_LEVEL = {"gbstats.bayes.ctw95": 0.10}
FIXED_HORIZON = {
    "gbstats.ttest.rel",
    "spotify.z",
    "ref.naive_z",
    "spotify.gst",
    "ref.gst_obf",
    "gbstats.bayes.ctw95",
}
ANYTIME = {"gbstats.seq.rel", "gbstats.seq.rel.tuned", "gbstats.seq.abs", "ref.acs", "ref.msprt"}
SCHEDULE_ORDER = ("fixed1", "daily14", "stream140", "overrun21")
GATE_CONF = 0.999


def load_units(run_dir: str | Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    out = Path(run_dir)
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    cfg = StudyConfig.from_dict(manifest["config"])
    units = []
    for u in cfg.unit_list():
        d = read_unit(out, u.id)
        # Same completeness rule as `run` and `status` (runner.unit_complete).
        if d is not None and unit_complete(d, cfg):
            units.append(d)
    return manifest, units


def provenance_groups(
    run_dir: str | Path, manifest: dict[str, Any], units: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Distinct (git SHA, versions) sets of the loaded units, with their unit counts.

    Units computed before per-unit provenance was recorded are attributed to the manifest
    only when the run has exactly one recorded start (manifest history or progress log);
    otherwise their provenance is reported as unknown.
    """
    out = Path(run_dir)
    starts = manifest.get("starts")
    if starts is None:
        log = out / "progress.log"
        text = log.read_text(encoding="utf-8") if log.exists() else ""
        n_starts = text.count(" start study=")
    else:
        n_starts = len(starts)
    groups: dict[str, dict[str, Any]] = {}
    for u in units:
        prov = unit_provenance(out, u["unit"])
        source = "unit record"
        if prov is None:
            if n_starts == 1:
                prov = {"git_sha": manifest.get("git_sha"), "versions": manifest.get("versions")}
                source = "run manifest (single start)"
            else:
                prov = {"git_sha": "unknown", "versions": {}}
                source = "not recorded"
        key = json.dumps([prov.get("git_sha"), prov.get("versions"), source], sort_keys=True)
        g = groups.setdefault(key, {**prov, "source": source, "units": 0})
        g["units"] += 1
    return list(groups.values())


def _cp(k: int, n: int, conf: float = 0.95) -> list[float]:
    return list(clopper_pearson(k, n, conf)) if n > 0 else [math.nan, math.nan]


def build_rows(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One row per (method, unit, schedule) with rates, intervals and family."""
    rows = []
    for u in units:
        metric = METRICS[u["metric"]]
        for method, by_sched in u["cells"].items():
            for sched, c in by_sched.items():
                n, k = c["n"], c["rejections"]
                row: dict[str, Any] = {
                    "method": method,
                    "unit": u["unit"],
                    "metric": u["metric"],
                    "n_max": u["n_max"],
                    "pi": u["pi"],
                    "hypothesis": u["hypothesis"],
                    "schedule": sched,
                    "claimed": c["claimed"],
                    "verified_ok": u["verified_ok"],
                    "n": n,
                    "rejections": k,
                    "rate": k / n,
                    "rate_ci95": _cp(k, n),
                    "reject_pos": c["reject_pos"],
                    "reject_neg": c["reject_neg"],
                    "no_result_at_stop": c["no_result_at_stop"],
                    "cover_n": c["cover_n"],
                    "cover_hits": c["cover_hits"],
                    "coverage": c["cover_hits"] / c["cover_n"] if c["cover_n"] else math.nan,
                    "coverage_ci95": _cp(c["cover_hits"], c["cover_n"]),
                    "stop_frac_mean": c["stop_n_mean"] / u["n_max"],
                    "stop_frac_ci95": _mean_ci(c["stop_n_mean"], c["stop_n_sd"], n, u["n_max"]),
                    "verification": c.get("verification"),
                }
                if u["hypothesis"] == "aa":
                    row["p_fpr_gt_alpha"] = binom_sf_test(k, n, ALPHA)
                    hist = c["first_rejection_hist"]
                    n_a, n_b = arm_counts(np.asarray(u["totals"][sched]), u["pi"])
                    small = np.minimum(n_a, n_b) < MIN_ARM
                    row["rejections_first_look"] = int(hist[0])
                    row["rejections_small_arm"] = int(np.asarray(hist)[small].sum())
                    if not c["claimed"]:
                        row["family"] = "informational"
                    elif metric.finite_variance:
                        row["family"] = "fpr_primary"
                    else:
                        row["family"] = "fpr_stress"
                else:
                    row["power"] = c["reject_pos"] / n
                    row["power_ci95"] = _cp(c["reject_pos"], n)
                    row["type_s"] = c["reject_neg"] / n
                    row["type_s_ci95"] = _cp(c["reject_neg"], n)
                    row["family"] = "coverage" if c["claimed"] else "informational"
                if c["cover_n"]:
                    misses = c["cover_n"] - c["cover_hits"]
                    row["p_miss_gt_alpha"] = binom_sf_test(misses, c["cover_n"], ALPHA)
                rows.append(row)
    return rows


def _mean_ci(mean: float, sd: float, n: int, scale: float) -> list[float]:
    if n < 2:
        return [mean / scale, mean / scale]
    q = float(st.t.ppf(0.975, n - 1)) * sd / math.sqrt(n)
    return [(mean - q) / scale, (mean + q) / scale]


def classify(rows: list[dict[str, Any]]) -> None:
    """Holm per (method, family) and the verdict per row (DESIGN 6.2), in place."""
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if not r["verified_ok"]:
            r["verdict"] = "excluded (engine mismatch)"
            continue
        if r["family"] in ("fpr_primary", "fpr_stress"):
            groups[(r["method"], r["family"])].append(r)
        elif r["family"] == "coverage":
            groups[(r["method"], "coverage")].append(r)
    for (_, fam), rs in groups.items():
        key = "p_miss_gt_alpha" if fam == "coverage" else "p_fpr_gt_alpha"
        rs_ok = [r for r in rs if key in r]
        adj = holm([r[key] for r in rs_ok]) if rs_ok else []
        for r, a in zip(rs_ok, adj, strict=True):
            r["holm_p"] = float(a)
            r["holm_m"] = len(rs_ok)
            if a < 0.05:
                r["verdict"] = "violation"
            elif fam != "coverage" and r["rate_ci95"][1] < 0.025:
                r["verdict"] = "conservative"
            else:
                r["verdict"] = "consistent"
    for r in rows:
        r.setdefault("verdict", "informational")


def detection_limits(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Smallest FPR a cell can flag at the strictest Holm step, per (method, FPR family)."""
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if r["family"] in ("fpr_primary", "fpr_stress") and "holm_m" in r:
            groups[(r["method"], r["family"])].append(r)
    out = []
    for (method, fam), rs in sorted(groups.items(), key=lambda kv: (_order(kv[0][0]), kv[0][1])):
        n, m = min(r["n"] for r in rs), rs[0]["holm_m"]
        out.append({"method": method, "family": fam, "m": m, "n": n,
                    "critical_rate": critical_rate(n, m, ALPHA, FWER)})  # fmt: skip
    return out


def findings(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Every violation, most severe first within method, with its kind (DESIGN 6.4 item 5).

    Kinds: "sequential" (FPR excess under repeated looks), "single look" (FPR excess of the
    one analysis at N_max, so not a sequential failure), "coverage" (A/B interval misses).
    """
    viol = [r for r in rows if r["verdict"] == "violation"]
    out = []
    for r in viol:
        if r["family"] == "coverage":
            kind = "coverage"
        elif r["schedule"] == "fixed1":
            kind = "single look"
        else:
            kind = "sequential"
        out.append({**r, "kind": kind})
    kinds = ("sequential", "single look", "coverage")
    out.sort(key=lambda r: (_order(r["method"]), kinds.index(r["kind"]), -r.get("rate", 0.0)))
    return out


def sensitivity_rows(
    rows: list[dict[str, Any]], sens: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Guarded FPR per claimed A/A cell, Holm per (method, family) as in the main analysis."""
    by_key = {(r["method"], r["unit"], r["schedule"]): r for r in rows if r["hypothesis"] == "aa"}
    out = []
    for unit in sens:
        if not unit.get("matches_run"):
            continue
        for method, by_s in unit["cells"].items():
            for sched, c in by_s.items():
                r = by_key.get((method, unit["unit"], sched))
                if r is None or r["family"] not in ("fpr_primary", "fpr_stress"):
                    continue
                k, n = c["rejections_guarded"], c["n"]
                out.append({
                    "method": method, "unit": unit["unit"], "metric": r["metric"],
                    "n_max": r["n_max"], "pi": r["pi"], "schedule": sched, "family": r["family"],
                    "min_arm": unit["min_arm"], "guarded_looks": c["guarded_looks"],
                    "looks": c["looks"], "n": n, "rejections": c["rejections"], "rate": r["rate"],
                    "rate_ci95": r["rate_ci95"], "verdict": r["verdict"],
                    "rejections_guarded": k, "rate_guarded": k / n,
                    "rate_guarded_ci95": _cp(k, n), "p_guarded": binom_sf_test(k, n, ALPHA),
                })  # fmt: skip
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in out:
        groups[(r["method"], r["family"])].append(r)
    for rs in groups.values():
        for r, a in zip(rs, holm([r["p_guarded"] for r in rs]), strict=True):
            r["holm_p_guarded"] = float(a)
            r["verdict_guarded"] = "violation" if a < FWER else "not flagged"
    return out


def design_deviations() -> list[str]:
    """Bullet lines of DESIGN.md section 13 when the source tree is available."""
    path = Path(__file__).resolve().parents[2] / "docs" / "DESIGN.md"
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8")
    start = text.find("## 13.")
    if start < 0:
        return []
    end = text.find("\n## ", start + 1)
    body = text[start : end if end > 0 else None].splitlines()[1:]
    return [ln for ln in body if ln.strip()]


# ---------------------------------------------------------------------------------- gates


def _find(rows: list[dict[str, Any]], **kw: Any) -> dict[str, Any] | None:
    hits = [r for r in rows if all(r[k] == v for k, v in kw.items())]
    if len(hits) > 1:
        raise KeyError(f"expected one row for {kw}, found {len(hits)}")
    return hits[0] if hits else None


def gates(rows: list[dict[str, Any]], units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Monte Carlo sanity gates of DESIGN 5.3 (99.9% CP intervals)."""
    out: list[dict[str, Any]] = []
    z = float(st.norm.isf(0.025))
    totals = {u["unit"]: u["totals"] for u in units}

    def add(name: str, ok: bool | None, detail: str) -> None:
        out.append({"gate": name, "pass": None if ok is None else bool(ok), "detail": detail})

    def missing(name: str) -> None:
        add(name, None, "cell not in this run")

    base = {"metric": "normal", "n_max": 20000, "pi": 0.5, "hypothesis": "aa"}
    for sched in ("daily14", "stream140"):
        r = _find(rows, method="ref.naive_z", schedule=sched, **base)
        if r is None:
            missing(f"positive control: naive z {sched}")
            continue
        t = np.asarray(totals[r["unit"]][sched], dtype=float)
        exact = float(crossing_probability(t / t[-1], [z] * t.size).sum())
        lo, hi = clopper_pearson(r["rejections"], r["n"], GATE_CONF)
        add(
            f"positive control: naive z {sched}",
            lo <= exact <= hi,
            f"FPR {r['rate']:.4f}, 99.9% CI [{lo:.4f}, {hi:.4f}], exact {exact:.4f}",
        )
    for sched in ("fixed1", "daily14", "stream140"):
        r = _find(rows, method="ref.gst_obf", schedule=sched, **base)
        if r is None:
            missing(f"negative control: ref.gst_obf {sched}")
            continue
        lo, hi = clopper_pearson(r["rejections"], r["n"], GATE_CONF)
        add(f"negative control: ref.gst_obf {sched}", lo <= 0.05 <= hi,
            f"FPR {r['rate']:.4f}, 99.9% CI [{lo:.4f}, {hi:.4f}] vs 0.05")  # fmt: skip
        for m in ("ref.acs", "ref.msprt"):
            r = _find(rows, method=m, schedule=sched, **base)
            if r is None:
                missing(f"negative control: {m} {sched}")
                continue
            lo, hi = clopper_pearson(r["rejections"], r["n"], GATE_CONF)
            add(f"negative control: {m} {sched}", lo <= 0.05,
                f"FPR {r['rate']:.4f}, 99.9% CI [{lo:.4f}, {hi:.4f}], lower <= 0.05")  # fmt: skip
    for pi in (0.5, 0.2):
        for r in rows:
            if not (
                r["metric"] == "normal" and r["n_max"] == 20000 and r["pi"] == pi
                and r["hypothesis"] == "aa" and r["schedule"] == "fixed1"
            ):  # fmt: skip
                continue
            lo, hi = clopper_pearson(r["rejections"], r["n"], GATE_CONF)
            if r["method"] in FIXED_HORIZON:
                level = FIXED_LEVEL.get(r["method"], 0.05)
                ok, rule = lo <= level <= hi, f"contains {level}"
            else:
                ok, rule = lo <= 0.05, "lower <= 0.05 (anytime-valid: conservative at one look)"
            add(f"fixed-horizon calibration: {r['method']} (pi={pi})", ok,
                f"FPR {r['rate']:.4f}, 99.9% CI [{lo:.4f}, {hi:.4f}], {rule}")  # fmt: skip
    for metric in ("normal", "bern_p10"):
        for pi in (0.5, 0.2):
            r = _find(rows, method="ref.naive_z", schedule="fixed1", metric=metric,
                      n_max=20000, pi=pi, hypothesis="ab")  # fmt: skip
            if r is None:
                missing(f"power calibration: naive z fixed1 {metric} pi={pi}")
                continue
            lo, hi = clopper_pearson(r["reject_pos"], r["n"], GATE_CONF)
            add(f"power calibration: naive z fixed1 {metric} pi={pi}", lo <= 0.80 <= hi,
                f"power {r['power']:.4f}, 99.9% CI [{lo:.4f}, {hi:.4f}] vs 0.80")  # fmt: skip
    # These compare two code paths that share arithmetic (the Bayes mirror and a z built from
    # the same variance helper; ref.acs and the gbstats mirror share acs_pvalue_halfwidth),
    # so they check self-consistency, not engine correctness. Engine correctness rests on
    # the in-loop and randomised-grid comparisons with the real engines.
    for key, what in (
        ("bayes_vs_z10", "Bayes mirror vs z at level 0.10, same variance helper"),
        ("acs_vs_gbstats_abs_nmax", "ref.acs vs gbstats absolute mirror, shared CS formula"),
    ):
        cells = sum(u["identities"][key]["cells"] for u in units)
        mism = sum(u["identities"][key]["mismatches"] for u in units)
        ties = sum(u["identities"][key]["ties"] for u in units)
        add(f"mirror self-consistency: {key}", mism == 0,
            f"{what}: {mism} mismatches, {ties} ties in {cells:,} replicate-looks")  # fmt: skip
    bad = [u["unit"] for u in units if not u["verified_ok"]]
    add("engine in-loop agreement", not bad, "all units agree" if not bad else ", ".join(bad))
    return out


def verification_totals(rows: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    tot: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for r in rows:
        v = r.get("verification")
        if not v:
            continue
        t = tot[r["method"]]
        t["replicate_looks"] += int(v["cells"])
        t["decision_mismatches"] += int(v["decision_mismatches"])
        t["ties"] += int(v["ties"])
        t["no_result_mismatches"] += int(v["no_result_mismatches"])
        t["interval_mismatches"] += int(v["interval_mismatches"])
        t["engine_rejections"] += int(v["engine_rejections"])
    return {k: dict(v) for k, v in tot.items()}


# --------------------------------------------------------------------------------- output


def _order(method: str) -> int:
    return METHOD_NAMES.index(method) if method in METHOD_NAMES else len(METHOD_NAMES)


def _pct(x: float) -> str:
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{100 * x:.1f}"


def _ci(ci: list[float]) -> str:
    return f"[{_pct(ci[0])}, {_pct(ci[1])}]"


def fpr_table(rows: list[dict[str, Any]], method: str) -> list[str]:
    lines = [
        "| metric | N_max | B share | " + " | ".join(SCHEDULE_ORDER) + " |",
        "|---|---|---|" + "---|" * len(SCHEDULE_ORDER),
    ]
    sel = [r for r in rows if r["method"] == method and r["hypothesis"] == "aa"]
    keys = sorted({(r["metric"], r["n_max"], r["pi"]) for r in sel},
                  key=lambda k: (list(METRICS).index(k[0]), k[1], -k[2]))  # fmt: skip
    for metric, n_max, pi in keys:
        cells = []
        for s in SCHEDULE_ORDER:
            hit = [r for r in sel if (r["metric"], r["n_max"], r["pi"], r["schedule"])
                   == (metric, n_max, pi, s)]  # fmt: skip
            if not hit:
                cells.append("")
                continue
            r = hit[0]
            mark = {"violation": " **V**", "conservative": " c", "informational": " i"}.get(
                r["verdict"], ""
            )
            cells.append(f"{_pct(r['rate'])} {_ci(r['rate_ci95'])}{mark}")
        lines.append(f"| {metric} | {n_max:,} | {int(pi * 100)}% | " + " | ".join(cells) + " |")
    return lines


def ab_table(rows: list[dict[str, Any]], method: str) -> list[str]:
    lines = [
        "| metric | N_max | B share | schedule | power % [95% CI] | type S % | coverage % "
        "[95% CI] | mean stop N / N_max |",
        "|---|---|---|---|---|---|---|---|",
    ]
    sel = [r for r in rows if r["method"] == method and r["hypothesis"] == "ab"]
    sel.sort(key=lambda r: (list(METRICS).index(r["metric"]), r["n_max"], -r["pi"],
                            SCHEDULE_ORDER.index(r["schedule"])))  # fmt: skip
    for r in sel:
        mark = " **V**" if r["verdict"] == "violation" else ""
        lines.append(
            f"| {r['metric']} | {r['n_max']:,} | {int(r['pi'] * 100)}% | {r['schedule']} | "
            f"{_pct(r['power'])} {_ci(r['power_ci95'])} | {_pct(r['type_s'])} | "
            f"{_pct(r['coverage'])} {_ci(r['coverage_ci95'])}{mark} | "
            f"{r['stop_frac_mean']:.3f} |"
        )
    return lines


def write_report(
    manifest: dict[str, Any],
    rows: list[dict[str, Any]],
    gate_list: list[dict[str, Any]],
    units: list[dict[str, Any]],
    dest: Path,
    figures: list[str],
    provenance: list[dict[str, Any]] | None = None,
    sens_rows: list[dict[str, Any]] | None = None,
    cpu_minutes: float | None = None,
) -> None:
    labels: dict[str, str] = {}
    for u in units:
        for m, info in u["methods"].items():
            labels[m] = info["label"]
    methods = [m for m in METHOD_NAMES if m in labels] + sorted(set(labels) - set(METHOD_NAMES))
    cfg = manifest["config"]
    vt = verification_totals(rows)
    lines = [
        f"# peekaudit results: {manifest['study']} study",
        "",
        "Generated by `peekaudit report` from the unit checkpoints of this run. Every number",
        "below is recomputable from `summary.json`. Rates are percentages with exact 95%",
        "Clopper-Pearson intervals. Marks: **V** = violation (Holm-adjusted one-sided exact",
        "test, p < 0.05), c = conservative (upper bound < 2.5%), i = informational (outside",
        "the method's claim; no pass/fail). Unmarked claimed cells are consistent.",
        "",
        "## What was run",
        "",
        f"- Replicates: {cfg['r_aa']:,} per A/A unit, {cfg['r_ab']:,} per A/B unit; chunks of "
        f"{cfg['chunk']}; master seed {cfg['master_seed']}; config hash "
        f"`{manifest['config_hash']}`.",
        f"- Units: {len(units)} of 52 completed; methods: {len(methods)}; schedules: "
        + ", ".join(cfg["schedules"])
        + ".",
        f"- Engine verification: first {cfg['k_verify']} replicates per unit and schedule "
        f"(spotify.gst on stream140: {cfg['k_verify_slow']}), real engine called look by look.",
        f"- Started {manifest.get('first_started')}; "
        f"{len(manifest.get('starts') or [None])} recorded start(s)"
        + (f"; {cpu_minutes:.1f} CPU-minutes over finished units." if cpu_minutes else "."),
    ]
    for g in provenance or [
        {"git_sha": manifest.get("git_sha"), "versions": manifest.get("versions") or {},
         "source": "run manifest", "units": len(units)}
    ]:  # fmt: skip
        vers = ", ".join(f"{k} {v}" for k, v in (g["versions"] or {}).items())
        lines.append(
            f"- Code: git `{g['git_sha']}` for {g['units']} units (source: {g['source']}); "
            f"versions: {vers or 'not recorded'}."
        )
    if provenance and len(provenance) > 1:
        lines.append("- **Units come from more than one code version or environment.**")
    lines += [
        "- Deviations from the design: see the last section of this report.",
        "",
        "## Sanity gates (DESIGN 5.3, 99.9% intervals)",
        "",
        "| gate | result | detail |",
        "|---|---|---|",
    ]
    verdict = {True: "pass", False: "**FAIL**", None: "not run"}
    lines += [f"| {g['gate']} | {verdict[g['pass']]} | {g['detail']} |"
              for g in gate_list]  # fmt: skip
    lines += [
        "",
        "## Engine mirror verification (in-loop, real engine vs vectorised mirror)",
        "",
        "| method | replicate-looks checked | decision mismatches | ties | no-result "
        "mismatches | interval mismatches | engine rejections |",
        "|---|---|---|---|---|---|---|",
    ]
    for m, t in sorted(vt.items(), key=lambda kv: _order(kv[0])):
        lines.append(
            f"| {m} | {t['replicate_looks']:,} | {t['decision_mismatches']} | {t['ties']} | "
            f"{t['no_result_mismatches']} | {t['interval_mismatches']} | {t['engine_rejections']} |"
        )
    lines += ["", "## Verdict counts (claimed cells)", "",
              "| method | family | cells | violations | consistent | conservative |",
              "|---|---|---|---|---|---|"]  # fmt: skip
    for m in methods:
        for fam in ("fpr_primary", "fpr_stress", "coverage"):
            rs = [r for r in rows if r["method"] == m and r["family"] == fam]
            if not rs:
                continue
            c = defaultdict(int)
            for r in rs:
                c[r["verdict"]] += 1
            lines.append(f"| {m} | {fam} | {len(rs)} | {c['violation']} | {c['consistent']} | "
                         f"{c['conservative']} |")  # fmt: skip
    lines += _detection_section(rows)
    lines += _findings_section(rows)
    lines += _small_arm_section(sens_rows or [])
    lines += _limitations_section(rows)
    lines += ["", "## False-positive rate (A/A), % [95% CI]", ""]
    for m in methods:
        lines += [f"### {m}: {labels.get(m, '')}", ""] + fpr_table(rows, m) + [""]
    lines += ["## A/B cells: power, type S, coverage at stopping, stopping N", ""]
    for m in methods:
        lines += [f"### {m}", ""] + ab_table(rows, m) + [""]
    if figures:
        lines += ["## Figures", ""] + [f"![{f}]({f})" for f in figures] + [""]
    dev = design_deviations()
    lines += ["## Deviations from the design (DESIGN.md section 13)", ""]
    lines += dev or ["See docs/DESIGN.md, section 13, in the source repository."]
    lines += [""]
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _cell_label(r: dict[str, Any]) -> str:
    return (f"{r['metric']}, N_max = {r['n_max']:,}, {int(r['pi'] * 100)}% to B, "
            f"{r['schedule']}")  # fmt: skip


def _detection_section(rows: list[dict[str, Any]]) -> list[str]:
    lim = detection_limits(rows)
    lines = [
        "",
        "## Monte Carlo detection limit",
        "",
        "A cell is flagged when the one-sided exact binomial test of FPR <= 5% rejects after",
        f"Holm adjustment at family-wise level {FWER} within the method's family. The table",
        "gives the smallest observed FPR that can be flagged at the strictest Holm step (m",
        "cells). Smaller excesses can go undetected, so 'consistent' means 'no excess",
        "detected at this size', not 'calibrated'.",
        "",
        "| method | family | cells m | replicates | smallest flaggable FPR % |",
        "|---|---|---|---|---|",
    ]
    lines += [f"| {d['method']} | {d['family']} | {d['m']} | {d['n']:,} | "
              f"{100 * d['critical_rate']:.2f} |" for d in lim]  # fmt: skip
    return lines


def _findings_section(rows: list[dict[str, Any]]) -> list[str]:
    fs = findings(rows)
    lines = [
        "",
        "## Findings (every violation)",
        "",
        "Kinds: 'sequential' = FPR excess under repeated looks; 'single look' = excess of the",
        "one analysis at N_max (the single z-test failing with few events or heavy skew in",
        "a small arm, not a sequential failure); 'coverage' = the interval at stopping",
        "misses the true effect more than 5% of the time. For A/A cells, 'first look' and",
        f"'arm < {MIN_ARM}' give the share of false",
        "positives whose first rejection came at the first look, or at a look where the",
        f"smaller arm held fewer than {MIN_ARM} units.",
        "",
    ]
    if not fs:
        return lines + ["No violations."]
    for r in fs:
        if r["kind"] == "coverage":
            lines.append(
                f"- {r['method']} (coverage): {_pct(r['coverage'])}% "
                f"{_ci(r['coverage_ci95'])} on {_cell_label(r)} "
                f"(Holm-adjusted p = {r['holm_p']:.2g}, m = {r['holm_m']})."
            )
            continue
        k = max(r["rejections"], 1)
        stress = ", stress" if r["family"] == "fpr_stress" else ""
        lines.append(
            f"- {r['method']} ({r['kind']}{stress}): FPR {_pct(r['rate'])}% "
            f"{_ci(r['rate_ci95'])} on {_cell_label(r)} "
            f"(Holm-adjusted p = {r['holm_p']:.2g}, m = {r['holm_m']}); first look "
            f"{100 * r['rejections_first_look'] / k:.0f}%, arm < {MIN_ARM} "
            f"{100 * r['rejections_small_arm'] / k:.0f}%."
        )
    return lines


def _small_arm_section(sens_rows: list[dict[str, Any]]) -> list[str]:
    lines = [
        "",
        f"## Sensitivity: rejections only where both arms hold >= {MIN_ARM} units",
        "",
        "Many excesses arise at the first looks of small experiments. With 20% of traffic to",
        "B and N_max = 1,000, the first stream140 look has 6 users in A and 1 in B. A zero-",
        "variance arm gives a standard error of 0, or a z-statistic rests on a few",
        "observations, and the engine reports significance on 1 to 30 users in B. The engines",
        "apply no minimum; GrowthBook's product layer hides results below a minimum sample",
        "size of 150 (not part of the audited library). The rows below recompute FPR counting",
        f"only rejections at looks where both arms hold at least {MIN_ARM} units (same data,",
        "same analysis, earlier results ignored), Holm-adjusted per method and family as in",
        "the main analysis. This is a sensitivity analysis, not a claim about any product.",
        "",
    ]
    if not sens_rows:
        return lines + ["Not computed for this run (`peekaudit sensitivity RUN_DIR` adds it)."]
    lines += [
        "| method | claimed A/A cells (primary + stress) | flagged (all looks) | "
        f"flagged (arms >= {MIN_ARM}) |",
        "|---|---|---|---|",
    ]
    by_m: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in sens_rows:
        by_m[r["method"]].append(r)
    for m in sorted(by_m, key=_order):
        rs = by_m[m]
        lines.append(
            f"| {m} | {len(rs)} | {sum(r['verdict'] == 'violation' for r in rs)} | "
            f"{sum(r['verdict_guarded'] == 'violation' for r in rs)} |"
        )
    lines += [
        "",
        "Cells flagged in the main analysis:",
        "",
        "| method | cell | FPR % all looks | FPR % arms >= 150 | looks kept | flagged |",
        "|---|---|---|---|---|---|",
    ]
    flagged = [r for r in sens_rows if r["verdict"] == "violation"]
    flagged.sort(key=lambda r: (_order(r["method"]), -r["rate"]))
    for r in flagged:
        lines.append(
            f"| {r['method']} | {_cell_label(r)} | {_pct(r['rate'])} {_ci(r['rate_ci95'])} | "
            f"{_pct(r['rate_guarded'])} {_ci(r['rate_guarded_ci95'])} | "
            f"{r['guarded_looks']}/{r['looks']} | "
            f"{'**V**' if r['verdict_guarded'] == 'violation' else 'no'} |"
        )
    return lines


def _limitations_section(rows: list[dict[str, Any]]) -> list[str]:
    lim = [d["critical_rate"] for d in detection_limits(rows) if d["family"] == "fpr_primary"]
    lim_txt = f"{100 * min(lim):.2f}-{100 * max(lim):.2f}%" if lim else "n/a"
    return [
        "",
        "## Limitations",
        "",
        "- **PyPI release, not the hosted code.** gbstats 0.8.0 on PyPI (2024-06-20) is older",
        "  than the code GrowthBook runs today, which still reports version 0.8.0 (for example,",
        "  post-stratification is now on by default). This audit covers the PyPI release only.",
        f"- **Detection limit.** Per-cell excesses below about {lim_txt} FPR (primary",
        "  families, see above) can go undetected; 'consistent' is not proof of calibration.",
        "- **stream140 understates continuous monitoring.** stream140 looks 140 times up to",
        "  N_max. Looking after every observation can only add rejections, so the error under",
        "  truly continuous monitoring is at least as high as reported here for every method.",
        "- **Small arms drive many excesses.** Rejections at looks with a handful of users in",
        "  one arm account for most of several headline rates (see the sensitivity section).",
        "  They are real outputs of the engines on these inputs, not adapter errors (the",
        "  mirrors match the engines on them), but a product with a minimum sample size would",
        "  hide them.",
        "- **Single-look failures are not sequential failures.** A spotify.gst violation at",
        "  `fixed1` is the single z-test at N_max with few events or heavy skew in the small arm",
        "  (for example 1% conversions with 200 users in B); the findings list labels them",
        "  separately.",
        "- **Mirrors, verified on a subsample.** Engine rates come from vectorised mirrors,",
        "  compared with real engine calls on the first replicates of every unit and schedule",
        "  and on randomised grids (zero disagreements). The two identity gates are",
        "  self-consistency checks between code paths that share arithmetic, not independent",
        "  checks of engine correctness.",
        "- **Scope.** Product-layer behaviour (minimum sample sizes, data-quality checks) is not",
        "  audited. Results apply to the pinned engine versions and the simulated conditions.",
        "  Calibration findings describe statistical behaviour, not software defects.",
    ]


def make_figures(rows: list[dict[str, Any]], dest: Path) -> list[str]:
    """FPR by schedule per method (faceted by metric) and power vs stopping N. PNG files."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    colors = {"fixed1": "#2a78d6", "daily14": "#eb6834", "stream140": "#1baf7a",
              "overrun21": "#6250d6"}  # fmt: skip
    markers = {"fixed1": "o", "daily14": "s", "stream140": "^", "overrun21": "D"}
    methods = sorted(dict.fromkeys(r["method"] for r in rows), key=_order)
    files = []
    for n_max in (1000, 20000):
        for pi in (0.5, 0.2):
            metrics = [m for m in METRICS if any(
                r["metric"] == m and r["n_max"] == n_max and r["pi"] == pi
                and r["hypothesis"] == "aa" for r in rows)]  # fmt: skip
            if not metrics:
                continue
            fig, axes = plt.subplots(len(metrics), 1, figsize=(10, 2.1 * len(metrics)),
                                     sharex=True, squeeze=False)  # fmt: skip
            for ax, metric in zip(axes[:, 0], metrics, strict=True):
                for si, s in enumerate(SCHEDULE_ORDER):
                    xs, ys, lo, hi = [], [], [], []
                    for mi, m in enumerate(methods):
                        hit = [r for r in rows if r["method"] == m and r["metric"] == metric
                               and r["n_max"] == n_max and r["pi"] == pi
                               and r["hypothesis"] == "aa" and r["schedule"] == s]  # fmt: skip
                        if hit:
                            r = hit[0]
                            xs.append(mi + (si - 1.5) * 0.18)
                            ys.append(100 * r["rate"])
                            lo.append(100 * (r["rate"] - r["rate_ci95"][0]))
                            hi.append(100 * (r["rate_ci95"][1] - r["rate"]))
                    ax.errorbar(xs, ys, yerr=[lo, hi], fmt=markers[s], color=colors[s], ms=4,
                                lw=1, capsize=0, label=s)  # fmt: skip
                ax.axhline(5, color="#52514e", lw=1, ls="--")
                ax.set_yscale("symlog", linthresh=10)
                ax.set_ylim(0, 100)
                ax.set_yticks([0, 5, 10, 20, 50, 100])
                ax.set_yticklabels(["0", "5", "10", "20", "50", "100"])
                ax.set_ylabel(f"{metric}\nFPR %", fontsize=8)
                ax.tick_params(labelsize=7)
                ax.grid(axis="y", color="#e5e4e0", lw=0.5)
                for side in ("top", "right"):
                    ax.spines[side].set_visible(False)
            axes[-1, 0].set_xticks(range(len(methods)))
            axes[-1, 0].set_xticklabels(methods, rotation=30, ha="right", fontsize=7)
            axes[0, 0].legend(ncol=4, fontsize=7, frameon=False, loc="upper left")
            fig.suptitle(f"A/A false-positive rate, N_max = {n_max:,}, B share {int(pi * 100)}%"
                         " (dashed: nominal 5%; bars: 95% CI)", fontsize=9)  # fmt: skip
            fig.tight_layout()
            name = f"fpr_n{n_max}_b{int(pi * 100)}.png"
            fig.savefig(dest / name, dpi=110)
            plt.close(fig)
            files.append(name)
    # Power vs mean stopping N for claimed sequential methods (A/B, daily14 and stream140).
    seq = ["gbstats.seq.rel", "gbstats.seq.rel.tuned", "gbstats.seq.abs", "spotify.gst",
           "ref.msprt", "ref.acs", "ref.gst_obf"]  # fmt: skip
    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
    shapes = ["o", "s", "^", "D", "v", "P", "X"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, s in zip(axes, ("daily14", "stream140"), strict=True):
        for m, col, mk in zip(seq, palette, shapes, strict=True):
            hit = [r for r in rows if r["method"] == m and r["hypothesis"] == "ab"
                   and r["schedule"] == s and r["n_max"] == 20000]  # fmt: skip
            if hit:
                ax.scatter([r["stop_frac_mean"] for r in hit], [100 * r["power"] for r in hit],
                           color=col, marker=mk, s=24, label=m)  # fmt: skip
        ax.set_title(f"{s}, N_max = 20,000 (all metrics and allocations)", fontsize=9)
        ax.set_xlabel("mean stopping N / N_max", fontsize=8)
        ax.grid(color="#e5e4e0", lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("power % (fixed-horizon z-test at N_max: 80%)", fontsize=8)
    axes[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(dest / "power_vs_stopping_n.png", dpi=110)
    plt.close(fig)
    files.append("power_vs_stopping_n.png")
    return files


def report(run_dir: str | Path, publish: str | Path | None = None) -> dict[str, Any]:
    """Write summary.json, REPORT.md and figures into ``<run_dir>/report`` (and ``publish``)."""
    manifest, units = load_units(run_dir)
    rows = build_rows(units)
    classify(rows)
    gate_list = gates(rows, units) if units else []
    dest = Path(run_dir) / "report"
    dest.mkdir(parents=True, exist_ok=True)
    figures = make_figures(rows, dest)
    summary = {
        "study": manifest["study"],
        "manifest": manifest,
        "gates": gate_list,
        "verification_totals": verification_totals(rows),
        "units": [
            {
                k: u[k]
                for k in ("unit", "lift", "design_lift", "tau", "reps", "identities", "verified_ok")
            }  # fmt: skip
            for u in units
        ],
        "cells": rows,
    }
    provenance = provenance_groups(run_dir, manifest, units)
    sens = sensitivity_rows(rows, load_sensitivity(run_dir))
    summary["provenance"] = provenance
    summary["sensitivity_min_arm"] = sens
    (dest / "summary.json").write_text(json.dumps(summary, separators=(",", ":"), default=float),
                                       encoding="utf-8")  # fmt: skip
    cpu = 0.0
    for u in units:
        t = Path(run_dir) / "units" / f"{u['unit']}.timing.json"
        if t.exists():
            cpu += json.loads(t.read_text(encoding="utf-8")).get("seconds", 0.0) / 60
    write_report(manifest, rows, gate_list, units, dest / "REPORT.md", figures,
                 provenance, sens, cpu or None)  # fmt: skip
    if publish:
        pub = Path(publish)
        pub.mkdir(parents=True, exist_ok=True)
        for f in ["summary.json", "REPORT.md", *figures]:
            shutil.copy2(dest / f, pub / f)
    return summary
