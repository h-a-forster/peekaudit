"""Audited engine adapters and reference methods behind one vectorised interface."""

from peekaudit.methods.base import (
    EvalContext,
    Method,
    MethodOutput,
    ScalarMethod,
    sample_moments,
    stopping_outcome,
)
from peekaudit.methods.registry import METHOD_NAMES, all_methods, get_method

__all__ = [
    "METHOD_NAMES",
    "EvalContext",
    "Method",
    "MethodOutput",
    "ScalarMethod",
    "all_methods",
    "get_method",
    "sample_moments",
    "stopping_outcome",
]
