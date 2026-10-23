"""Method registry: stable names -> factories. Engine adapters are imported lazily."""

from __future__ import annotations

import importlib
from collections.abc import Callable

from peekaudit.methods.base import Method


def _ref(cls_name: str) -> Callable[[], Method]:
    def make() -> Method:
        mod = importlib.import_module("peekaudit.methods.reference")
        return getattr(mod, cls_name)()

    return make


def _adapter(module: str, factory: str, **kwargs: object) -> Callable[[], Method]:
    def make() -> Method:
        mod = importlib.import_module(f"peekaudit.methods.{module}")
        return getattr(mod, factory)(**kwargs)

    return make


#: The 11 methods of the study, in report order.
REGISTRY: dict[str, Callable[[], Method]] = {
    "gbstats.seq.rel": _adapter("gbstats_adapter", "GbstatsSequential", relative=True),
    "gbstats.seq.rel.tuned": _adapter(
        "gbstats_adapter", "GbstatsSequential", relative=True, tuned=True
    ),
    "gbstats.seq.abs": _adapter("gbstats_adapter", "GbstatsSequential", relative=False),
    "gbstats.ttest.rel": _adapter("gbstats_adapter", "GbstatsTTest"),
    "gbstats.bayes.ctw95": _adapter("gbstats_adapter", "GbstatsBayes"),
    "spotify.gst": _adapter("spotify_adapter", "SpotifyGST"),
    "spotify.z": _adapter("spotify_adapter", "SpotifyZ"),
    "ref.naive_z": _ref("NaiveZ"),
    "ref.msprt": _ref("MSPRT"),
    "ref.acs": _ref("ACS"),
    "ref.gst_obf": _ref("GstOBF"),
}

METHOD_NAMES: tuple[str, ...] = tuple(REGISTRY)


def get_method(spec: str) -> Method:
    """A registered name, or ``module:factory`` for a third-party method (factory() -> Method)."""
    if spec in REGISTRY:
        return REGISTRY[spec]()
    if ":" in spec:
        module, attr = spec.split(":", 1)
        obj = getattr(importlib.import_module(module), attr)
        method = obj() if callable(obj) and not isinstance(obj, Method) else obj
        if not isinstance(method, Method):
            raise TypeError(f"{spec} did not produce a peekaudit Method")
        return method
    raise KeyError(f"unknown method {spec!r}; known: {', '.join(METHOD_NAMES)}")


def all_methods(names: tuple[str, ...] | list[str] | None = None) -> list[Method]:
    return [get_method(n) for n in (names or METHOD_NAMES)]
