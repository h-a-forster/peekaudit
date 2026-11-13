"""Living audit: a quick battery against whatever engine release is installed.

``scripts/living.py`` builds one fresh environment per target with the newest release (or,
for ``gbstats-main``, GrowthBook's repository at a pinned commit), runs
``python -m peekaudit.living battery <target>`` inside it, and compares the verdicts with
``docs/results/living/baseline.json`` (:func:`compare`). The weekly workflow
(``.github/workflows/living.yml``) opens an issue when a verdict changes.

Every number is a real engine call; no mirror is trusted:

* ``mirror/<method>``: the study's vectorised mirror against the engine on randomised
  inputs (:func:`peekaudit.verify.grid_agreement`). ``differs`` means the engine's formulas
  changed and the study's numbers no longer describe this release.
* ``fpr/<method>/<unit>/<schedule>``: :func:`peekaudit.conformance.check_conformance`
  with the method's engine calls; ``violates`` after Holm within the method.
* tea-tasting and expectation rows: the 2026-11-09 study's checks
  (``docs/results/newlib-2026-11-09/``) at fixed replicate counts, Holm within the target.

Seeds are fixed, so the same environment always gives the same verdicts. The engine is not
the only thing that moves: numpy, scipy and pandas also resolve to their newest releases,
and the platform can differ from the baseline's, so the summary lists every environment
entry that changed next to the verdicts. A battery change regenerates the baseline in the
same commit.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import sys
import time
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from peekaudit.conformance import Battery, check_conformance
from peekaudit.methods.base import EvalContext, Method, MethodOutput, ScalarMethod
from peekaudit.simulate import LookStats
from peekaudit.stats import binom_sf_test, clopper_pearson, holm

#: GrowthBook main at the commit the ``gbstats-main`` target installs. Bump it on purpose
#: (and regenerate the baseline), never implicitly: main is not a release.
GROWTHBOOK_COMMIT = "c4cd05e0d725d075c0f0b33e744e7894b3b3b640"
GROWTHBOOK_GBSTATS = (
    f"gbstats @ git+https://github.com/growthbook/growthbook@{GROWTHBOOK_COMMIT}"
    "#subdirectory=packages/stats"
)
FWER = 0.05


@dataclass(frozen=True)
class Target:
    """One environment of the living audit."""

    name: str
    dist: str  # distribution whose version is reported
    requirements: tuple[str, ...]  # newest release unless pinned
    python: str
    build_from_source: bool = False  # only for a pinned repository commit


TARGETS: dict[str, Target] = {
    t.name: t
    for t in (
        # 3.11: gbstats 0.8.0 declares numpy<2 and pandas<2 (the study's environment).
        Target("gbstats", "gbstats", ("gbstats",), "3.11"),
        Target("gbstats-main", "gbstats", (GROWTHBOOK_GBSTATS,), "3.11", True),
        # bokeh<3.10: chartify 5.0.1 fails to import with bokeh 3.10 (pyproject.toml).
        Target(
            "spotify-confidence", "spotify-confidence", ("spotify-confidence", "bokeh<3.10"), "3.11"
        ),
        Target("tea-tasting", "tea-tasting", ("tea-tasting",), "3.12"),
        Target("expectation", "expectation", ("expectation",), "3.12"),
    )
}


def environment(target: Target) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for dist in (target.dist, "numpy", "scipy", "pandas", "peekaudit"):
        try:
            out[dist] = importlib.metadata.version(dist)
        except importlib.metadata.PackageNotFoundError:
            out[dist] = None
    try:
        direct = importlib.metadata.distribution(target.dist).read_text("direct_url.json")
        commit = json.loads(direct or "{}").get("vcs_info", {}).get("commit_id")
        if commit:
            out[f"{target.dist} commit"] = commit
    except importlib.metadata.PackageNotFoundError:
        pass
    out["python"] = platform.python_version()
    out["platform"] = sys.platform
    return out


def scrub(obj: Any) -> Any:
    """``obj`` with this environment's and the home directory's paths replaced, so results
    (and the committed baseline) never carry a machine's paths."""
    text = json.dumps(obj)
    for path, label in ((sys.prefix, "<env>"), (str(Path.home()), "<home>")):
        for form in {path, path.replace("\\", "/")}:
            text = text.replace(json.dumps(form)[1:-1], label)
    return json.loads(text)


def _text(exc: BaseException, limit: int) -> str:
    """``exc``'s message, scrubbed before it is cut to ``limit`` characters: a cut never
    splits a path the scrub would then miss, and the kept text does not depend on how long
    this machine's paths are."""
    return str(scrub(str(exc)))[:limit]


# --- engine-backed FPR battery ------------------------------------------------------------


class EngineMethod(Method):
    """A study method whose ``evaluate`` is its real engine call (no mirror)."""

    def __init__(self, inner: Method) -> None:
        super().__init__(inner.name, inner.scale, inner.claim, inner.label, inner.engine)
        self.inner = inner

    def supports(self, metric: Any) -> bool:
        return self.inner.supports(metric)

    def evaluate(self, stats: LookStats, ctx: EvalContext) -> MethodOutput:
        return self.inner.evaluate_engine(stats, ctx)


#: Engine FPR batteries: the study's headline cells (rare conversions, small skewed arms)
#: and clean controls, at a size each engine's call cost allows.
GBSTATS_BATTERY = Battery(
    "living-gbstats",
    ("normal", "bern_p01", "lognorm"),
    (1000, 20000),
    (0.5, 0.2),
    ("daily14", "stream140"),
    500,
)
SPOTIFY_BATTERY = Battery(
    "living-spotify", ("normal", "bern_p10"), (1000,), (0.5, 0.2), ("daily14",), 150
)


def _mirror_rows(names: Sequence[str], n: int) -> Callable[[], dict[str, Any]]:
    def run() -> dict[str, Any]:
        from peekaudit.methods.registry import get_method
        from peekaudit.verify import grid_agreement

        out: dict[str, Any] = {}
        for name in names:
            res = grid_agreement(get_method(name), n=n, seed=0)
            ok = all(a.ok for a in res.values())
            out[f"mirror/{name}"] = {
                "verdict": "agrees" if ok else "differs",
                "detail": {k: a.to_dict() for k, a in res.items()},
            }
        return out

    return run


def _fpr_rows(name: str, battery: Battery) -> Callable[[], dict[str, Any]]:
    def run() -> dict[str, Any]:
        from peekaudit.methods.registry import get_method

        rep = check_conformance(EngineMethod(get_method(name)), battery)
        return {
            f"fpr/{name}/{c.unit}/{c.schedule}": {
                "verdict": "violates" if c.violation else "ok",
                "detail": f"{c.false_positives}/{c.n} = {100 * c.rate:.2f}% "
                f"[{100 * c.ci95[0]:.2f}, {100 * c.ci95[1]:.2f}], Holm p {c.holm_p:.3g}",
            }
            for c in rep.cells
        }

    return run


# --- tea-tasting ----------------------------------------------------------------------------

TEA_BATTERY = Battery(
    "living-tea",
    ("normal", "bern_p01", "lognorm"),
    (1000, 20000),
    (0.5, 0.2),
    ("fixed1",),
    2000,
)
#: The Proportion interval defect (tea-tasting#327): (control rate, treatment rate,
#: control n, treatment n), from ``repro_tea.py``.
TEA_COVERAGE = ((0.10, 0.30, 8000, 2000), (0.30, 0.10, 8000, 2000), (0.10, 0.12, 5000, 5000))


def _tea_methods() -> dict[str, tuple[Any, bool]]:
    from tea_tasting.metrics.mean import Mean
    from tea_tasting.metrics.proportion import Proportion

    return {
        "welch_t": (Mean("x"), False),
        "z": (Mean("x", use_t=False), False),
        "student_t": (Mean("x", equal_var=True), False),
        "prop_pooled": (Proportion("x", method="norm"), True),
        "prop_unpooled": (Proportion("x", method="norm", equal_var=False), True),
    }


def _tea_scalar(name: str, metric: Any, binary: bool, errors: list[str]) -> ScalarMethod:
    from tea_tasting.aggr import Aggregates

    def aggs(n: int, s: float, q: float) -> Any:
        mean = s / n
        if binary:
            var = mean * (1 - mean) * n / (n - 1)
        else:
            var = max((q - n * mean * mean) / (n - 1), 0.0) if n > 1 else 0.0
        return Aggregates(count_=n, mean_={"x": mean}, var_={"x": var})

    def decide(na: int, sa: float, qa: float, nb: int, sb: float, qb: float) -> Any:
        try:
            r = metric.analyze_aggregates(aggs(na, sa, qa), aggs(nb, sb, qb))
        except Exception as exc:  # noqa: BLE001 - counted and reported, never hidden
            errors.append(f"{type(exc).__name__}: {_text(exc, 80)}")
            return None
        if not np.isfinite(r.pvalue):
            return None
        sign = int(np.sign(r.effect_size)) if r.pvalue < 0.05 else 0
        return sign, float(r.effect_size_ci_lower), float(r.effect_size_ci_upper)

    supports_binary = binary

    class _Tea(ScalarMethod):
        def supports(self, m: Any) -> bool:
            return not m.ratio and (m.binary if supports_binary else True)

    return _Tea(f"tea.{name}", decide, "absolute", "fixed")


def _tea_rows() -> dict[str, Any]:
    from tea_tasting.aggr import Aggregates

    rows: dict[str, Any] = {}
    stat: list[tuple[str, int, int, float]] = []  # row, k, n, null rate of k
    for name, (metric, binary) in _tea_methods().items():
        errors: list[str] = []
        rep = check_conformance(_tea_scalar(name, metric, binary, errors), TEA_BATTERY)
        for c in rep.cells:
            key = f"fpr/{name}/{c.unit}/{c.schedule}"
            stat.append((key, c.false_positives, c.n, 0.05))
            rows[key] = {"detail": f"{c.false_positives}/{c.n} = {100 * c.rate:.2f}%"}
        if errors:
            rows[f"errors/{name}"] = {"verdict": "raises", "detail": sorted(set(errors))[:5]}
    rng = np.random.default_rng(20261113)
    reps = 2000
    for pc, pt, nc, nt in TEA_COVERAGE:
        sc, st = rng.binomial(nc, pc, reps), rng.binomial(nt, pt, reps)
        for name in ("prop_pooled", "prop_unpooled"):
            metric = _tea_methods()[name][0]
            miss_abs = miss_rel = 0
            for a, b in zip(sc, st, strict=True):
                r = metric.analyze_aggregates(
                    Aggregates(count_=nc, mean_={"x": a / nc}),
                    Aggregates(count_=nt, mean_={"x": b / nt}),
                )
                miss_abs += not r.effect_size_ci_lower <= pt - pc <= r.effect_size_ci_upper
                miss_rel += not (
                    r.rel_effect_size_ci_lower <= pt / pc - 1 <= r.rel_effect_size_ci_upper
                )
            for scale, miss in (("abs", miss_abs), ("rel", miss_rel)):
                key = f"coverage/{name}/{scale}/{pc}->{pt}/n{nc}+{nt}"
                stat.append((key, miss, reps, 0.05))
                rows[key] = {"detail": f"covers {1 - miss / reps:.3f} (nominal 0.950)"}
    _holm_verdicts(rows, stat)
    return rows


# --- expectation ----------------------------------------------------------------------------

EXP_REPS = 200
EXP_HORIZON = 500
#: (row, data, SequentialTesting keyword arguments, batch size), from ``pk_exp.py``.
EXP_CELLS: tuple[tuple[str, str, dict[str, Any], int], ...] = (
    ("mean/default", "normal", {"test_type": "mean", "null_value": 0.0}, 1),
    ("mean/default", "uniform", {"test_type": "mean", "null_value": 0.0}, 1),
    ("mean/default", "normal", {"test_type": "mean", "null_value": 0.0}, 10),
    (
        "mean/known_var",
        "normal",
        {"test_type": "mean", "null_value": 0.0, "known_variance": 1.0},
        1,
    ),
    ("proportion/default", "bern10", {"test_type": "proportion", "null_value": 0.1}, 1),
    ("proportion/default", "bern10", {"test_type": "proportion", "null_value": 0.1}, 10),
)


def _exp_data(name: str, rng: np.random.Generator, n: int) -> np.ndarray:
    if name == "normal":
        return rng.standard_normal(n)
    if name == "uniform":
        return rng.uniform(-1, 1, n)
    if name == "bern10":
        return (rng.random(n) < 0.1).astype(float)
    raise ValueError(name)


def _exp_rows() -> dict[str, Any]:
    rows: dict[str, Any] = {"readme_import": {"verdict": "works", "detail": ""}}
    try:
        from expectation.seqtest import SequentialTesting  # noqa: F401 - the README's import
    except ImportError as exc:
        rows["readme_import"] = {"verdict": "fails", "detail": _text(exc, 200)}
    from expectation.seqtest.sequential_e_testing import SequentialTesting as ST

    stat: list[tuple[str, int, int, float]] = []
    for ci, (row, data, kw, batch) in enumerate(EXP_CELLS):
        rej = ran = 0
        errors: dict[str, int] = {}
        for rep in range(EXP_REPS):
            x = _exp_data(data, np.random.default_rng([20261113, ci, rep]), EXP_HORIZON)
            test = ST(**kw)
            try:
                hit = False
                for i in range(0, EXP_HORIZON, batch):
                    if test.update(x[i : i + batch]).reject_null:
                        hit = True
                        break
            except Exception as exc:  # noqa: BLE001 - counted and reported
                msg = f"{type(exc).__name__}: {_text(exc, 60)}"
                errors[msg] = errors.get(msg, 0) + 1
                continue
            ran += 1
            rej += hit
        key = f"fpr/{row}/{data}/batch{batch}"
        rows[key] = {"detail": f"{rej}/{ran} rejected, {EXP_REPS - ran} raised {errors}"}
        if ran:
            stat.append((key, rej, ran, 0.05))
        if errors:
            rows[f"errors/{row}/{data}/batch{batch}"] = {
                "verdict": "all raise" if not ran else "some raise",
                "detail": sorted(errors),
            }
    _holm_verdicts(rows, stat)
    for rec in rows.values():
        rec.setdefault("verdict", "not run (every replicate raised)")
    return rows


def _holm_verdicts(rows: dict[str, Any], stat: list[tuple[str, int, int, float]]) -> None:
    """One-sided exact binomial test of each excess rate, Holm across the target."""
    if not stat:
        return
    adj = holm([binom_sf_test(k, n, p0) for _, k, n, p0 in stat])
    for (key, k, n, _), a in zip(stat, adj, strict=True):
        lo, hi = clopper_pearson(k, n)
        rows[key]["verdict"] = "violates" if a < FWER else "ok"
        rows[key]["detail"] += f"; rate {k / n:.4f} [{lo:.4f}, {hi:.4f}], Holm p {a:.3g}"


# --- targets --------------------------------------------------------------------------------

GBSTATS_METHODS = (
    "gbstats.seq.rel",
    "gbstats.seq.rel.tuned",
    "gbstats.seq.abs",
    "gbstats.ttest.rel",
    "gbstats.bayes.ctw95",
)


def sections(target: str) -> dict[str, Callable[[], dict[str, Any]]]:
    """The battery's checkpoint units for a target, in order."""
    if target in ("gbstats", "gbstats-main"):
        return {
            "mirror": _mirror_rows(GBSTATS_METHODS, 600),
            "fpr-seq-rel": _fpr_rows("gbstats.seq.rel", GBSTATS_BATTERY),
            "fpr-seq-abs": _fpr_rows("gbstats.seq.abs", GBSTATS_BATTERY),
        }
    if target == "spotify-confidence":
        return {
            "mirror-z": _mirror_rows(("spotify.z",), 300),
            "mirror-gst": _mirror_rows(("spotify.gst",), 40),
            "fpr-gst": _fpr_rows("spotify.gst", SPOTIFY_BATTERY),
        }
    if target == "tea-tasting":
        return {"tea": _tea_rows}
    if target == "expectation":
        return {"expectation": _exp_rows}
    raise KeyError(f"unknown target {target!r}; known: {', '.join(TARGETS)}")


def battery(target: str, work: str | Path) -> dict[str, Any]:
    """Run a target's battery in this environment; each finished section is checkpointed
    to ``work/sections.jsonl`` and skipped on rerun."""
    secs = sections(target)
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    done_path = work / "sections.jsonl"
    done: dict[str, Any] = {}
    if done_path.exists():
        for line in done_path.read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            done[rec["section"]] = rec
    t0 = time.time()
    warnings.filterwarnings("ignore")
    for name, run in secs.items():
        if name in done:
            continue
        t1 = time.time()
        try:
            rows = run()
        except Exception as exc:  # noqa: BLE001 - an engine that cannot run is a verdict
            rows = {
                f"{name}/error": {
                    "verdict": "ERROR",
                    "detail": f"{type(exc).__name__}: {_text(exc, 300)}",
                }
            }
        rec = {"section": name, "rows": rows, "secs": round(time.time() - t1, 1)}
        with done_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        done[name] = rec
    rows = {k: v for name in secs for k, v in done[name]["rows"].items()}
    return scrub(
        {
            "target": target,
            "environment": environment(TARGETS[target]),
            "verdicts": {k: v["verdict"] for k, v in rows.items()},
            "details": {k: v["detail"] for k, v in rows.items()},
            "section_secs": {name: done[name]["secs"] for name in secs},
            "elapsed_s": round(time.time() - t0, 1),
        }
    )


# --- comparison -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Change:
    target: str
    row: str
    old: str | None  # None: the row is new
    new: str | None  # None: the row is gone


def compare(results: dict[str, dict[str, Any]], baseline: dict[str, Any]) -> list[Change]:
    """Every row whose verdict differs from the baseline, over the targets in ``results``."""
    changes: list[Change] = []
    base_targets = baseline.get("targets", {})
    for target in sorted(results):
        new = results[target]["verdicts"]
        old = base_targets.get(target, {}).get("verdicts", {})
        for row in sorted(set(old) | set(new)):
            if old.get(row) != new.get(row):
                changes.append(Change(target, row, old.get(row), new.get(row)))
    return changes


def make_baseline(results: dict[str, dict[str, Any]], note: str = "") -> dict[str, Any]:
    return {
        "note": note,
        "targets": {
            t: {k: r[k] for k in ("environment", "verdicts", "details")}
            for t, r in sorted(results.items())
        },
    }


def _label(env: dict[str, Any], dist: str) -> str:
    v = env.get(dist) or "-"
    commit = env.get(f"{dist} commit")
    return f"{v} @ {commit[:10]}" if commit else v


def summary_markdown(
    results: dict[str, dict[str, Any]],
    baseline: dict[str, Any],
    changes: list[Change],
    failed: dict[str, str] | None = None,
) -> str:
    """The run summary; it is also the issue body when a verdict changed."""
    failed = failed or {}
    base_targets = baseline.get("targets", {})
    lines = [
        "# peekaudit living audit",
        "",
        f"{len(changes)} verdict change(s) against `docs/results/living/baseline.json`"
        + (f"; {len(failed)} target(s) could not run" if failed else "")
        + ".",
        "",
        "| Target | Baseline | Now | Rows | Not ok | Changed |",
        "|---|---|---|---|---|---|",
    ]
    for target in sorted(set(results) | set(failed)):
        dist = TARGETS[target].dist if target in TARGETS else target
        old = _label(base_targets.get(target, {}).get("environment", {}), dist)
        if target in failed:
            lines.append(f"| {target} | {old} | not run | - | - | - |")
            continue
        r = results[target]
        bad = sum(v not in ("ok", "agrees", "works") for v in r["verdicts"].values())
        n_changed = sum(c.target == target for c in changes)
        lines.append(
            f"| {target} | {old} | {_label(r['environment'], dist)} | {len(r['verdicts'])} | "
            f"{bad} | {n_changed} |"
        )
    env_lines = []
    for target in sorted(results):
        old_env = base_targets.get(target, {}).get("environment", {})
        new_env = results[target].get("environment", {})
        diff = [
            f"{k} {old_env.get(k, '-')} -> {new_env.get(k, '-')}"
            for k in sorted(set(old_env) | set(new_env))
            if old_env.get(k) != new_env.get(k)
        ]
        if diff:
            env_lines.append(f"- **{target}**: {', '.join(diff)}")
    if env_lines:
        lines += ["", "## Environment changes", "", "Any of these can move a verdict.", ""]
        lines += env_lines
    if changes:
        lines += ["", "## Changed verdicts", ""]
        for c in changes:
            lines.append(f"- **{c.target}** `{c.row}`: {c.old or 'new row'} -> {c.new or 'gone'}")
            detail = results.get(c.target, {}).get("details", {}).get(c.row)
            if detail:
                lines.append(f"  - now: {json.dumps(detail)[:400]}")
    if failed:
        lines += ["", "## Not run", ""]
        lines += [f"- **{t}**: {msg}" for t, msg in sorted(failed.items())]
    lines += [
        "",
        "A change means the newest release behaves differently from the audited one. Check",
        "it against the study before updating REPORT.md, then regenerate the baseline",
        "(`uv run python scripts/living.py --update-baseline`).",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python -m peekaudit.living")
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("battery", help="run one target's battery in this environment")
    b.add_argument("target", choices=sorted(TARGETS))
    b.add_argument("--work", required=True, help="checkpoint directory")
    b.add_argument("--result", required=True, help="JSON file to write")
    args = p.parse_args(argv)
    res = battery(args.target, args.work)
    Path(args.result).parent.mkdir(parents=True, exist_ok=True)
    Path(args.result).write_text(json.dumps(res, indent=1, sort_keys=True) + "\n", "utf-8")
    bad = sum(v not in ("ok", "agrees", "works") for v in res["verdicts"].values())
    print(f"{args.target}: {len(res['verdicts'])} rows, {bad} not ok, {res['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
