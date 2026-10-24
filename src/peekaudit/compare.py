"""Cross-check two published studies cell by cell (for example quick against full).

Both studies draw replicates from ``SeedSequence(master_seed, spawn_key=(unit, chunk))``.
With the same master seed and chunk size, and replicate counts that are whole chunks, the
smaller study's replicates are exactly the first chunks of the larger one, so the two are
*nested*, not independent. Comparing the small study with the large one would then compare
a sample with a superset of itself.
:func:`compare_summaries` detects nesting from the manifests and, when it holds, compares
the small study with the *remaining* replicates of the large one (counts differenced),
which are independent of it.

For every cell present in both studies and every rate (FPR or power in ``rate``, coverage
at stopping in ``coverage``), a two-sided Fisher exact test asks whether the two samples
have the same underlying rate. P-values are Holm-adjusted within each rate kind across all
compared cells. Under correct, deterministic code the null holds everywhere, so a rejection
indicates a code or seeding difference between the studies, not a statistical finding.
Verdicts (DESIGN 6.2) are also cross-tabulated; a larger study flags smaller excesses, so
"consistent -> violation" is expected where the excess lies between the two detection limits.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from scipy import stats as st

from peekaudit.stats import holm

RATE_KINDS = (("rate", "rejections", "n"), ("coverage", "cover_hits", "cover_n"))
_CFG_KEYS_FOR_NESTING = ("master_seed", "chunk")


def load_summary(path: str | Path) -> dict[str, Any]:
    """Read a ``summary.json`` (or a results directory that contains one)."""
    p = Path(path)
    if p.is_dir():
        p = p / "summary.json"
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)  # type: ignore[no-any-return]


def is_nested(small: dict[str, Any], large: dict[str, Any]) -> bool:
    """True when the small study's replicates are the leading replicates of the large one.

    Needs the same master seed and chunk size, and for each replicate count either the same
    count or a multiple of the chunk size no larger than the large study's. A partial last
    chunk is not a prefix of the full chunk: arm B is drawn after a shorter arm-A draw, so
    its random stream is offset and its data differ.
    """
    cs, cl = small["manifest"]["config"], large["manifest"]["config"]
    if any(cs.get(k) != cl.get(k) for k in _CFG_KEYS_FOR_NESTING):
        return False
    chunk = int(cs["chunk"])
    for key in ("r_aa", "r_ab"):
        r_s, r_l = int(cs[key]), int(cl[key])
        if not (r_s == r_l or (r_s < r_l and r_s % chunk == 0)):
            return False
    return True


def _key(cell: dict[str, Any]) -> tuple[str, str, str]:
    return (cell["method"], cell["unit"], cell["schedule"])


def compare_summaries(
    small: dict[str, Any], large: dict[str, Any], *, nested: bool | None = None
) -> dict[str, Any]:
    """Cell-by-cell comparison; see the module docstring for the method."""
    if nested is None:
        nested = is_nested(small, large)
    big = {_key(c): c for c in large["cells"]}
    pairs = [(c, big[_key(c)]) for c in small["cells"] if _key(c) in big]
    if not pairs:
        raise ValueError("the two studies share no cells")

    comparisons: list[dict[str, Any]] = []
    for kind, k_field, n_field in RATE_KINDS:
        rows = []
        for a, b in pairs:
            if a.get(n_field) is None or b.get(n_field) is None:
                continue
            k1, n1 = int(a[k_field]), int(a[n_field])
            k2, n2 = int(b[k_field]), int(b[n_field])
            if nested:
                k2, n2 = k2 - k1, n2 - n1
                if k2 < 0 or n2 < 0 or k2 > n2:
                    raise ValueError(
                        f"counts are not nested for {_key(a)} ({kind}): the studies differ"
                    )
            if n1 == 0 or n2 == 0:
                continue
            p = float(st.fisher_exact([[k1, n1 - k1], [k2, n2 - k2]])[1])
            rows.append(
                {
                    "method": a["method"],
                    "unit": a["unit"],
                    "schedule": a["schedule"],
                    "kind": kind,
                    "small_k": k1,
                    "small_n": n1,
                    "other_k": k2,
                    "other_n": n2,
                    "small_rate": k1 / n1,
                    "other_rate": k2 / n2,
                    "p": min(p, 1.0),
                }
            )
        adj = holm([r["p"] for r in rows]) if rows else []
        for r, q in zip(rows, adj, strict=True):
            r["holm_p"] = float(q)
        comparisons.extend(rows)

    verdicts: Counter[tuple[str, str]] = Counter()
    claimed_changes: list[dict[str, Any]] = []
    for a, b in pairs:
        if a.get("family") in (None, "informational"):
            continue
        va, vb = a.get("verdict", "?"), b.get("verdict", "?")
        verdicts[(va, vb)] += 1
        if (va == "violation") != (vb == "violation"):
            claimed_changes.append(
                {
                    "method": a["method"],
                    "unit": a["unit"],
                    "schedule": a["schedule"],
                    "family": a["family"],
                    "small_verdict": va,
                    "large_verdict": vb,
                    "small_rate": a.get("coverage" if a["family"] == "coverage" else "rate"),
                    "large_rate": b.get("coverage" if b["family"] == "coverage" else "rate"),
                    "small_ci95": a.get(
                        "coverage_ci95" if a["family"] == "coverage" else "rate_ci95"
                    ),
                    "large_ci95": b.get(
                        "coverage_ci95" if b["family"] == "coverage" else "rate_ci95"
                    ),
                }
            )

    def _overlap(a: dict[str, Any], b: dict[str, Any], f: str) -> bool:
        lo1, hi1 = a[f]
        lo2, hi2 = b[f]
        return bool(lo1 <= hi2 and lo2 <= hi1)

    overlap = {
        "rate": sum(_overlap(a, b, "rate_ci95") for a, b in pairs),
        "coverage": sum(
            _overlap(a, b, "coverage_ci95")
            for a, b in pairs
            if a.get("coverage_ci95") and b.get("coverage_ci95")
        ),
    }
    by_kind = {
        kind: {
            "compared": sum(r["kind"] == kind for r in comparisons),
            "raw_p_lt_0.05": sum(r["kind"] == kind and r["p"] < 0.05 for r in comparisons),
            "holm_p_lt_0.05": sum(r["kind"] == kind and r["holm_p"] < 0.05 for r in comparisons),
            "min_holm_p": min(
                (r["holm_p"] for r in comparisons if r["kind"] == kind), default=float("nan")
            ),
        }
        for kind, _, _ in RATE_KINDS
    }
    return {
        "small_study": small["study"],
        "large_study": large["study"],
        "small_config_hash": small["manifest"]["config_hash"],
        "large_config_hash": large["manifest"]["config_hash"],
        "nested": nested,
        "cells_compared": len(pairs),
        "by_kind": by_kind,
        "ci95_overlap": overlap,
        "verdict_crosstab": [
            {"small": va, "large": vb, "cells": n} for (va, vb), n in sorted(verdicts.items())
        ],
        "violation_changes": claimed_changes,
        "comparisons": comparisons,
    }


def _pct(x: float) -> str:
    return f"{100 * x:.1f}"


def render_markdown(res: dict[str, Any]) -> str:
    """A short Markdown section summarising :func:`compare_summaries`."""
    s, lg = res["small_study"], res["large_study"]
    other = f"replicates of {lg} not in {s}" if res["nested"] else lg
    lines = [
        f"# Cross-check: {s} vs {lg}",
        "",
        f"Config hashes `{res['small_config_hash']}` ({s}) and `{res['large_config_hash']}` "
        f"({lg}); {res['cells_compared']} cells in both.",
        "",
    ]
    if res["nested"]:
        lines += [
            f"The two studies share the master seed and chunk size, so the {s} replicates are "
            f"the first replicates of {lg} (nested, not independent). Each {s} rate is "
            f"therefore compared with the {other} (counts differenced), which is independent "
            "of it.",
            "",
        ]
    lines += [
        "Two-sided Fisher exact test per cell, Holm-adjusted within each rate kind. Under "
        "unchanged, deterministic code every null is true, so rejections would point to a "
        "code or seeding difference.",
        "",
        "| rate | cells compared | raw p < 0.05 | Holm p < 0.05 | smallest Holm p |",
        "|---|---|---|---|---|",
    ]
    for kind, d in res["by_kind"].items():
        label = "FPR / power" if kind == "rate" else "coverage at stopping"
        lines.append(
            f"| {label} | {d['compared']:,} | {d['raw_p_lt_0.05']:,} | {d['holm_p_lt_0.05']:,} "
            f"| {d['min_holm_p']:.3g} |"
        )
    ov = res["ci95_overlap"]
    lines += [
        "",
        f"95% CP intervals of the two studies overlap in {ov['rate']:,} of "
        f"{res['cells_compared']:,} cells for FPR / power and in {ov['coverage']:,} for "
        "coverage (the intervals are not independent when nested; shown for orientation).",
        "",
        "Verdicts on claimed cells (DESIGN 6.2):",
        "",
        f"| {s} | {lg} | cells |",
        "|---|---|---|",
    ]
    for row in res["verdict_crosstab"]:
        lines.append(f"| {row['small']} | {row['large']} | {row['cells']} |")
    ch = res["violation_changes"]
    lines += ["", f"Cells whose violation status differs ({len(ch)}):", ""]
    if ch:
        lines += [
            f"| method | cell | family | {s} | {lg} |",
            "|---|---|---|---|---|",
        ]
        for c in sorted(ch, key=lambda c: (c["method"], c["unit"], c["schedule"])):
            lines.append(
                f"| {c['method']} | {c['unit']} {c['schedule']} | {c['family']} | "
                f"{_pct(c['small_rate'])} [{_pct(c['small_ci95'][0])}, "
                f"{_pct(c['small_ci95'][1])}] {c['small_verdict']} | "
                f"{_pct(c['large_rate'])} [{_pct(c['large_ci95'][0])}, "
                f"{_pct(c['large_ci95'][1])}] {c['large_verdict']} |"
            )
    return "\n".join(lines) + "\n"


def write_comparison(
    small_path: str | Path, large_path: str | Path, out: str | Path
) -> dict[str, Any]:
    """Compare two summaries and write ``<out>.md`` and ``<out>.json``."""
    res = compare_summaries(load_summary(small_path), load_summary(large_path))
    out_p = Path(out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.with_suffix(".md").write_text(render_markdown(res), encoding="utf-8", newline="\n")
    with open(out_p.with_suffix(".json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh)
    return res
