"""Living audit: the quick battery against the newest release of every audited engine.

One command (about 40 min; needs uv and network for the installs):

    uv run python scripts/living.py                  # every target, newest releases
    uv run python scripts/living.py --targets tea-tasting  # a subset
    uv run python scripts/living.py --pin gbstats==0.8.0 --targets gbstats

Each target gets a fresh wheels-only venv (``<out>/venvs/<target>``, the target's Python)
holding a wheel of this checkout plus the engine (``gbstats-main``: GrowthBook's repository
at the commit pinned in ``peekaudit.living``, built from source without its dependencies,
which then come as wheels); ``python -m peekaudit.living battery`` runs
inside it and checkpoints to ``<out>/work/<target>``. A finished target's
``<out>/results/<target>.json`` is kept, so a rerun skips it (``--fresh`` starts over);
``<out>`` defaults to ``tmp/living/<today>``.
The verdicts are then compared with ``docs/results/living/baseline.json``; ``<out>`` gets
``summary.md`` and ``changes.json``.

Exit code: 0 no verdict changed, 1 a verdict changed, 2 a target could not be installed or
run (nothing is compared then for that target).

``--update-baseline`` writes the baseline from this run instead (every target must run).
``notify <out> --repo OWNER/NAME`` opens an issue, or comments on the open one, when
``changes.json`` lists a change; it calls ``gh`` (or ``LIVING_GH``: a command, or a JSON list
of arguments, which the tests use).
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

from peekaudit.living import TARGETS, compare, make_baseline, summary_markdown

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "docs" / "results" / "living" / "baseline.json"
ISSUE_TITLE = "Living audit: a verdict changed"


def log(out: Path, msg: str) -> None:
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    with (out / "progress.log").open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def venv_python(venv: Path) -> Path:
    return venv / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def build_wheel(out: Path) -> Path:
    dist = out / "wheel"
    shutil.rmtree(dist, ignore_errors=True)
    subprocess.run(
        ["uv", "build", "--wheel", "--out-dir", str(dist), str(ROOT)],
        check=True,
        capture_output=True,
        text=True,
    )
    (wheel,) = dist.glob("peekaudit-*.whl")
    return wheel


# --no-config: this repository's [tool.uv] overrides (dependency pins lifted for its own
# lockfile) must not leak into the targets, which get what a user's install resolves to.
def install(target: str, wheel: Path, venv: Path, pins: dict[str, str]) -> None:
    t = TARGETS[target]
    if venv.exists():
        shutil.rmtree(venv)
    subprocess.run(
        ["uv", "venv", "--no-config", "--python", t.python, str(venv), "-q"],
        check=True,
        capture_output=True,
        text=True,
    )
    py = str(venv_python(venv))
    base = ["uv", "pip", "install", "-q", "--no-config", "--python", py]
    reqs = [pins.get(r, r) for r in t.requirements]
    if t.build_from_source:
        # Build only the pinned commit (no dependencies), then its dependencies as wheels.
        subprocess.run([*base, "--no-deps", *reqs], check=True, capture_output=True, text=True)
        reqs = json.loads(
            subprocess.run(
                [
                    py,
                    "-I",
                    "-c",
                    "import importlib.metadata as m, json; "
                    f"print(json.dumps(m.requires({t.dist!r}) or []))",
                ],
                check=True,
                capture_output=True,
                text=True,
            ).stdout
        )
    names = [r for r in t.requirements if "@" not in r]
    cmd = [
        *base,
        "--only-binary",
        ":all:",
        *(f"--refresh-package={n.split('<')[0]}" for n in names),
    ]
    subprocess.run([*cmd, str(wheel), *reqs], check=True, capture_output=True, text=True)


def run_target(target: str, out: Path, wheel: Path, pins: dict[str, str]) -> None:
    venv = out / "venvs" / target
    install(target, wheel, venv, pins)
    result = out / "results" / f"{target}.json"
    tmp = result.with_suffix(".partial")
    proc = subprocess.run(
        [
            str(venv_python(venv)),
            "-I",
            "-m",
            "peekaudit.living",
            "battery",
            target,
            "--work",
            str(out / "work" / target),
            "--result",
            str(tmp),
        ],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"battery exit {proc.returncode}: {proc.stderr.strip()[-800:]}")
    tmp.replace(result)
    log(out, f"{target}: {proc.stdout.strip()}")


def cmd_run(args: argparse.Namespace) -> int:
    out = Path(args.out)
    targets = args.targets.split(",") if args.targets else list(TARGETS)
    unknown = sorted(set(targets) - set(TARGETS))
    if unknown:
        print(f"error: unknown targets {unknown}; known: {sorted(TARGETS)}", file=sys.stderr)
        return 2
    pins = {}
    for spec in args.pin or []:
        name = spec.split("==")[0]
        if "==" not in spec or not any(name in TARGETS[t].requirements for t in targets):
            print(
                f"error: --pin {spec!r} is not NAME==VERSION of a target package", file=sys.stderr
            )
            return 2
        pins[name] = spec
    if args.fresh:
        shutil.rmtree(out, ignore_errors=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    config = {"targets": sorted(targets), "pins": sorted(pins.values())}
    cfg_path = out / "run.json"
    if cfg_path.exists() and json.loads(cfg_path.read_text("utf-8")) != config:
        print(
            f"error: {out} holds a run with other targets or pins; use --fresh or another --out",
            file=sys.stderr,
        )
        return 2
    cfg_path.write_text(json.dumps(config) + "\n", "utf-8")
    log(out, f"targets {','.join(targets)} pins {sorted(pins.values()) or 'none'}")
    failed: dict[str, str] = {}
    wheel: Path | None = None
    for t in targets:
        if (out / "results" / f"{t}.json").exists():
            log(out, f"{t}: result exists, skipped")
            continue
        try:
            wheel = wheel or build_wheel(out)
            t0 = time.time()
            run_target(t, out, wheel, pins)
            log(out, f"{t}: done in {time.time() - t0:.0f} s")
        except (subprocess.CalledProcessError, RuntimeError, OSError) as exc:
            detail = getattr(exc, "stderr", None) or str(exc)
            failed[t] = " ".join(str(detail).split())[-600:]
            log(out, f"{t}: FAILED {failed[t]}")
    results = {
        t: json.loads((out / "results" / f"{t}.json").read_text("utf-8"))
        for t in targets
        if (out / "results" / f"{t}.json").exists()
    }
    if args.update_baseline:
        if failed:
            print("error: not every target ran; baseline unchanged", file=sys.stderr)
            return 2
        BASELINE.parent.mkdir(parents=True, exist_ok=True)
        note = (
            "Verdicts of the living-audit battery (scripts/living.py). Regenerate with "
            "--update-baseline only after checking a change against the study."
        )
        BASELINE.write_text(json.dumps(make_baseline(results, note), indent=1) + "\n", "utf-8")
        log(out, f"baseline written: {len(results)} targets")
        return 0
    baseline = json.loads(BASELINE.read_text("utf-8")) if BASELINE.exists() else {}
    changes = compare(results, baseline)
    text = summary_markdown(results, baseline, changes, failed)
    (out / "summary.md").write_text(text, "utf-8")
    (out / "changes.json").write_text(
        json.dumps({"changes": [c.__dict__ for c in changes], "failed": failed}, indent=1),
        "utf-8",
    )
    print(text)
    if failed:
        return 2
    return 1 if changes else 0


def gh(*args: str) -> str:
    spec = os.environ.get("LIVING_GH", "gh")  # a command, or a JSON list of arguments
    base = json.loads(spec) if spec.startswith("[") else shlex.split(spec)
    proc = subprocess.run([*base, *args], check=True, capture_output=True, text=True)
    return proc.stdout


def cmd_notify(args: argparse.Namespace) -> int:
    out = Path(args.out)
    data = json.loads((out / "changes.json").read_text("utf-8"))
    if not data["changes"]:
        print("no verdict changed; no issue")
        return 0
    body = (out / "summary.md").read_text("utf-8")
    if args.run_url:
        body += f"\nRun and artifact: {args.run_url}\n"
    body_file = out / "issue.md"
    body_file.write_text(body, "utf-8")
    listed = json.loads(
        gh(
            "issue",
            "list",
            "--repo",
            args.repo,
            "--state",
            "open",
            "--search",
            f'"{ISSUE_TITLE}" in:title',
            "--json",
            "number,title",
        )
    )
    match = [i for i in listed if i["title"] == ISSUE_TITLE]
    if match:
        gh(
            "issue",
            "comment",
            str(match[0]["number"]),
            "--repo",
            args.repo,
            "--body-file",
            str(body_file),
        )
        print(f"commented on issue #{match[0]['number']}")
    else:
        gh(
            "issue",
            "create",
            "--repo",
            args.repo,
            "--title",
            ISSUE_TITLE,
            "--body-file",
            str(body_file),
        )
        print("opened an issue")
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["notify"]:
        n = argparse.ArgumentParser(prog="living.py notify")
        n.add_argument("out")
        n.add_argument("--repo", required=True)
        n.add_argument("--run-url", default="")
        return cmd_notify(n.parse_args(argv[1:]))
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", default=None, help="default tmp/living/<today>")
    p.add_argument("--targets", default="", help="comma-separated, default all")
    p.add_argument("--pin", action="append", help="NAME==VERSION instead of the newest")
    p.add_argument("--fresh", action="store_true", help="discard earlier results in --out")
    p.add_argument("--update-baseline", action="store_true")
    args = p.parse_args(argv)
    args.out = args.out or str(ROOT / "tmp" / "living" / time.strftime("%Y-%m-%d"))
    return cmd_run(args)


if __name__ == "__main__":
    sys.exit(main())
