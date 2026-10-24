"""spotify.gst against like-for-like reference GSTs, on the A/A units of a finished run.

Usage: uv run python scripts/gst_reference_comparison.py [--run runs/full]
       [--out docs/results/full] [--workers 3] [--attempts 3] [--units 'bern_*']

Evaluates spotify.gst, ref.gst_pow2 (same rho = 2 spending, reference variances) and
ref.gst_obf on the run's own data (see ``peekaudit.gst_reference``). One checkpoint per unit
under ``<run>/gst_reference/``; a rerun skips finished units; failures are retried up to
``--attempts`` times. Writes ``gst_reference.{json,md}`` to ``--out``. The full run takes
about 10 minutes on 3 workers.
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback

from peekaudit.gst_reference import gst_reference_dir, run_gst_reference, write_summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", default="runs/full")
    ap.add_argument("--out", default="docs/results/full")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--attempts", type=int, default=3)
    ap.add_argument("--units", default=None, help="glob on unit ids")
    args = ap.parse_args()
    for attempt in range(1, args.attempts + 1):
        try:
            run_gst_reference(args.run, workers=args.workers, units=args.units)
            break
        except Exception:  # noqa: BLE001 - retried; checkpoints keep finished units
            log = gst_reference_dir(args.run) / "progress.log"
            with log.open("a", encoding="utf-8") as fh:
                fh.write(f"attempt {attempt}/{args.attempts} failed:\n{traceback.format_exc()}")
            if attempt == args.attempts:
                return 1
            time.sleep(5)
    summary = write_summary(args.run, args.out)
    print(f"{len(summary['cells'])} cells; run counts reproduced: {summary['all_match_run']}")
    return 0 if summary["all_match_run"] else 2


if __name__ == "__main__":
    sys.exit(main())
