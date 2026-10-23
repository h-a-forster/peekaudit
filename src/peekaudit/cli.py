"""Command-line interface (docs/DESIGN.md, section "CLI").

``run`` / ``status`` run and monitor a checkpointed study, ``report`` turns a run into
tables and figures, ``verify`` compares engine mirrors with the real engines on random
inputs, ``conform`` runs the conformance battery on any method, ``versions`` prints the
environment.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from collections.abc import Sequence
from importlib import metadata

from peekaudit import __version__

#: Distributions recorded in every study manifest and printed by ``peekaudit versions``.
TRACKED_DISTRIBUTIONS = (
    "numpy",
    "scipy",
    "pandas",
    "gbstats",
    "spotify-confidence",
    "statsmodels",
    "pydantic",
)


def environment_report() -> dict[str, str]:
    """Versions of Python and the tracked distributions ("not installed" when absent)."""
    report = {"python": platform.python_version(), "peekaudit": __version__}
    for name in TRACKED_DISTRIBUTIONS:
        try:
            report[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            report[name] = "not installed"
    return report


def _cmd_versions(_: argparse.Namespace) -> int:
    for name, version in environment_report().items():
        print(f"{name:<20} {version}")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    from peekaudit.runner import load_config, run_study

    cfg = load_config(args.config)
    out = args.out or f"runs/{cfg.name}"
    res = run_study(cfg, out, workers=args.workers, units_glob=args.units, force=args.force)
    return 1 if res.failed else 0


def _cmd_status(args: argparse.Namespace) -> int:
    from peekaudit.runner import study_status

    st = study_status(args.run_dir)
    print(
        f"study {st['study']} (config {st['config_hash']}): {st['done']}/{st['total']} units done"
    )
    print(f"CPU time of finished units: {st['cpu_minutes_done']} min")
    shas = st["units_by_git_sha"]
    if shas and (len(shas) > 1 or st["current_git_sha"] not in shas):
        print(
            "finished units by code version (a resume keeps them; --force recomputes): "
            + ", ".join(f"{sha} x{n}" for sha, n in shas.items())
            + f"; current code {st['current_git_sha']}"
        )
    if st["failed"]:
        print("failed: " + ", ".join(st["failed"]))
    if st["unverified"]:
        print("engine mismatch (not reported): " + ", ".join(st["unverified"]))
    if st["last_log"]:
        print("last: " + st["last_log"])
    return 0


def _cmd_verify(args: argparse.Namespace) -> int:
    from peekaudit.methods.registry import get_method
    from peekaudit.verify import Agreement, engine_identities, grid_agreement

    if args.method == "identities":
        names: list[str] = []
    elif args.method:
        names = [args.method]
    else:
        names = [
            "gbstats.seq.rel", "gbstats.seq.rel.tuned", "gbstats.seq.abs", "gbstats.ttest.rel",
            "gbstats.bayes.ctw95", "spotify.z", "spotify.gst",
        ]  # fmt: skip
    report: dict[str, dict[str, object]] = {}
    ok = True

    def show(name: str, res: dict[str, Agreement]) -> None:
        nonlocal ok
        report[name] = {kind: agr.to_dict() for kind, agr in res.items()}
        for kind, agr in res.items():
            ok &= agr.ok
            print(
                f"{name:<22} {kind:<10} cells={agr.cells:<5} decision_mismatches="
                f"{agr.decision_mismatches} ties={agr.ties} no_result_mismatches="
                f"{agr.no_result_mismatches} interval_mismatches={agr.interval_mismatches} "
                f"max_rel_diff={agr.max_interval_rel_diff:.1e} {'OK' if agr.ok else 'FAIL'}"
            )

    for name in names:
        n = args.n_gst if name == "spotify.gst" else args.n
        show(name, grid_agreement(get_method(name), n=n, seed=args.seed))
    if args.method in (None, "identities"):
        show("engine_identities", engine_identities(n=args.n, seed=args.seed))
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump({"seed": args.seed, "n": args.n, "n_gst": args.n_gst, "results": report},
                      fh, indent=1)  # fmt: skip
    return 0 if ok else 1


def _cmd_report(args: argparse.Namespace) -> int:
    from peekaudit.report import report

    summary = report(args.run_dir, publish=args.publish)
    failed = [g["gate"] for g in summary["gates"] if g["pass"] is False]
    print(f"{len(summary['cells'])} cells, {len(summary['gates'])} gates, {len(failed)} failed")
    for g in failed:
        print(f"  gate FAILED: {g}")
    return 1 if failed else 0


def _cmd_sensitivity(args: argparse.Namespace) -> int:
    from peekaudit.sensitivity import load_sensitivity, run_sensitivity

    run_sensitivity(args.run_dir, min_arm=args.min_arm, workers=args.workers)
    bad = [r["unit"] for r in load_sensitivity(args.run_dir, args.min_arm) if not r["matches_run"]]
    if bad:
        print("unguarded counts differ from the run (not used by the report): " + ", ".join(bad))
    return 1 if bad else 0


def _cmd_compare(args: argparse.Namespace) -> int:
    from peekaudit.compare import write_comparison

    res = write_comparison(args.small, args.large, args.out)
    for kind, d in res["by_kind"].items():
        print(f"{kind}: {d['compared']} cells, {d['holm_p_lt_0.05']} differ (Holm p < 0.05)")
    print(f"nested={res['nested']} violation changes={len(res['violation_changes'])}")
    return 1 if any(d["holm_p_lt_0.05"] for d in res["by_kind"].values()) else 0


def _cmd_conform(args: argparse.Namespace) -> int:
    from peekaudit.conformance import check_conformance
    from peekaudit.methods.registry import get_method

    try:
        method = get_method(args.method)
        if args.claim:
            method.claim = args.claim
        rep = check_conformance(method, args.battery)
    except Exception as exc:  # noqa: BLE001 - CLI reports errors with exit code 2
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    print(rep.summary())
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(rep.to_dict(), fh, indent=1)
    return 0 if rep.passed else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="peekaudit",
        description="Audit the sequential-testing guarantees of A/B testing engines.",
    )
    parser.add_argument("--version", action="version", version=f"peekaudit {__version__}")
    sub = parser.add_subparsers(dest="command")
    p_versions = sub.add_parser("versions", help="print Python, library and engine versions")
    p_versions.set_defaults(func=_cmd_versions)

    p_run = sub.add_parser("run", help="run or resume a checkpointed study")
    p_run.add_argument("--config", default="quick", help="quick | full | smoke | PATH.json")
    p_run.add_argument("--out", help="output directory (default runs/<config name>)")
    p_run.add_argument("--workers", type=int, help="worker processes (capped at 4)")
    p_run.add_argument("--units", help="glob on unit ids, e.g. 'normal-*'")
    p_run.add_argument("--force", action="store_true", help="recompute finished units")
    p_run.set_defaults(func=_cmd_run)

    p_status = sub.add_parser("status", help="progress of a study run")
    p_status.add_argument("run_dir")
    p_status.set_defaults(func=_cmd_status)

    p_verify = sub.add_parser("verify", help="mirror vs engine on random inputs")
    p_verify.add_argument(
        "--method", help="one engine method, or 'identities' (default: all, plus identities)"
    )
    p_verify.add_argument("--n", type=int, default=2000, help="single-look inputs per method")
    p_verify.add_argument("--n-gst", type=int, default=200, help="random paths for spotify.gst")
    p_verify.add_argument("--seed", type=int, default=0)
    p_verify.add_argument("--json", help="write the agreement counts to this file")
    p_verify.set_defaults(func=_cmd_verify)

    p_report = sub.add_parser("report", help="tables, gates and figures from a study run")
    p_report.add_argument("run_dir")
    p_report.add_argument("--publish", help="also copy the results here (e.g. docs/results/quick)")
    p_report.set_defaults(func=_cmd_report)

    p_sens = sub.add_parser(
        "sensitivity", help="FPR counting only looks where both arms reach a minimum size"
    )
    p_sens.add_argument("run_dir")
    p_sens.add_argument("--min-arm", type=int, default=150, help="minimum units per arm")
    p_sens.add_argument("--workers", type=int, help="worker processes (capped at 4)")
    p_sens.set_defaults(func=_cmd_sensitivity)

    p_cmp = sub.add_parser("compare", help="cell-by-cell cross-check of two published studies")
    p_cmp.add_argument("small", help="summary.json (or its directory) of the smaller study")
    p_cmp.add_argument("large", help="summary.json (or its directory) of the larger study")
    p_cmp.add_argument("--out", required=True, help="output stem; writes .md and .json")
    p_cmp.set_defaults(func=_cmd_compare)

    p_conform = sub.add_parser("conform", help="conformance battery for one method")
    p_conform.add_argument("--method", required=True, help="registry name or module:factory")
    p_conform.add_argument("--battery", default="smoke", choices=["smoke", "quick", "standard"])
    p_conform.add_argument(
        "--claim", choices=["anytime", "planned", "fixed"], help="override the method's claim"
    )
    p_conform.add_argument("--json", help="write the report to this file")
    p_conform.set_defaults(func=_cmd_conform)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        return 2
    return int(args.func(args))
