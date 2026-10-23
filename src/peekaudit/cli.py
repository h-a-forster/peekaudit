"""Command-line interface.

Implemented: ``--version`` and ``versions`` (environment report: Python, core libraries and
the audited engines with their licences). Planned subcommands (``run``, ``status``,
``report``, ``conform``, ``verify``) are specified in docs/DESIGN.md, section "CLI".
"""

from __future__ import annotations

import argparse
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="peekaudit",
        description="Audit the sequential-testing guarantees of A/B testing engines.",
    )
    parser.add_argument("--version", action="version", version=f"peekaudit {__version__}")
    sub = parser.add_subparsers(dest="command")
    p_versions = sub.add_parser("versions", help="print Python, library and engine versions")
    p_versions.set_defaults(func=_cmd_versions)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        return 2
    return int(args.func(args))
