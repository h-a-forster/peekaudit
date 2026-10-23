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

from peekaudit.gst import crossing_probability
from peekaudit.runner import StudyConfig, is_done, unit_path
from peekaudit.scenarios import ALPHA, METRICS
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
        if is_done(out, u, cfg):
            units.append(json.loads(unit_path(out, u.id).read_text(encoding="utf-8")))
    return manifest, units


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
    for key in ("bayes_vs_z10", "acs_vs_gbstats_abs_nmax"):
        cells = sum(u["identities"][key]["cells"] for u in units)
        mism = sum(u["identities"][key]["mismatches"] for u in units)
        ties = sum(u["identities"][key]["ties"] for u in units)
        add(f"identity: {key}", mism == 0,
            f"{mism} mismatches, {ties} ties in {cells:,} replicate-looks")  # fmt: skip
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
) -> None:
    methods = list(units[0]["methods"]) if units else []
    labels: dict[str, str] = {}
    for u in units:
        for m, info in u["methods"].items():
            labels[m] = info["label"]
            if m not in methods:
                methods.append(m)
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
        f"- Code: git `{manifest['git_sha']}`; started {manifest.get('first_started')}.",
        "- Versions: " + ", ".join(f"{k} {v}" for k, v in manifest["versions"].items()) + ".",
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
    for m, t in vt.items():
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
    lines += ["", "## False-positive rate (A/A), % [95% CI]", ""]
    for m in methods:
        lines += [f"### {m}: {labels.get(m, '')}", ""] + fpr_table(rows, m) + [""]
    lines += ["## A/B cells: power, type S, coverage at stopping, stopping N", ""]
    for m in methods:
        lines += [f"### {m}", ""] + ab_table(rows, m) + [""]
    if figures:
        lines += ["## Figures", ""] + [f"![{f}]({f})" for f in figures] + [""]
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    methods = list(dict.fromkeys(r["method"] for r in rows))
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
    (dest / "summary.json").write_text(json.dumps(summary, indent=1, default=float),
                                       encoding="utf-8")  # fmt: skip
    write_report(manifest, rows, gate_list, units, dest / "REPORT.md", figures)
    if publish:
        pub = Path(publish)
        pub.mkdir(parents=True, exist_ok=True)
        for f in ["summary.json", "REPORT.md", *figures]:
            shutil.copy2(dest / f, pub / f)
    return summary
