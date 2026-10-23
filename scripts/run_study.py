"""Run a study with automatic resume-on-failure.

Usage: uv run python scripts/run_study.py [--config quick] [--out runs/quick] [--attempts 3]

Wraps ``peekaudit run`` in up to ``--attempts`` attempts. Each finished unit is checkpointed,
so a rerun only computes what is missing. Progress goes to ``<out>/progress.log``; this
wrapper appends its own attempt lines there too.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="quick")
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--attempts", type=int, default=3)
    args = ap.parse_args()
    out = Path(args.out or f"runs/{args.config}")
    out.mkdir(parents=True, exist_ok=True)
    log = out / "progress.log"
    cmd = [
        sys.executable, "-m", "peekaudit", "run", "--config", args.config,
        "--out", str(out), "--workers", str(args.workers),
    ]  # fmt: skip
    code = 1
    for attempt in range(1, args.attempts + 1):
        stamp = datetime.now(UTC).isoformat(timespec="seconds")
        with log.open("a", encoding="utf-8") as fh:
            fh.write(f"{stamp} wrapper attempt {attempt}/{args.attempts}: {' '.join(cmd)}\n")
        code = subprocess.call(cmd)
        if code == 0:
            return 0
        time.sleep(5)
    with log.open("a", encoding="utf-8") as fh:
        fh.write(
            f"{datetime.now(UTC).isoformat(timespec='seconds')} wrapper gave up (exit {code})\n"
        )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
