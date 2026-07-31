"""Validation helpers for differential-expression ensembles."""

from __future__ import annotations

from collections.abc import Sequence
from numbers import Integral
from typing import TypeVar


SUPPORTED_DE_METHODS = ("limrots", "deqms", "proda", "limma", "rots")
DEFAULT_ENSEMBLE_METHODS = ("limrots", "deqms", "proda")
ResultValue = TypeVar("ResultValue")


def resolve_ensemble_methods(
    methods: Sequence[str] | None,
    min_k: int,
) -> tuple[str, ...]:
    """Return canonical ensemble members after validating the configuration."""
    selected = DEFAULT_ENSEMBLE_METHODS if methods is None else methods
    if isinstance(selected, (str, bytes)):
        raise ValueError("ensemble methods must be a sequence of member names")
    if not selected:
        raise ValueError("ensemble requires at least one member method")

    canonical: list[str] = []
    for method in selected:
        if not isinstance(method, str) or not method.strip():
            raise ValueError("ensemble member names must be non-empty strings")
        name = method.strip().lower()
        if name not in SUPPORTED_DE_METHODS:
            raise ValueError(
                f"unknown ensemble member '{method}'. Supported: {SUPPORTED_DE_METHODS}"
            )
        if name in canonical:
            raise ValueError(f"duplicate ensemble member '{name}'")
        canonical.append(name)

    if isinstance(min_k, bool) or not isinstance(min_k, Integral):
        raise ValueError(f"ensemble min_k must be an integer; got {min_k!r}")
    if not 1 <= min_k <= len(canonical):
        raise ValueError(
            "ensemble min_k must be between 1 and the number of member methods "
            f"({len(canonical)}); got {min_k!r}"
        )
    return tuple(canonical)


def validate_de_ensemble(
    *,
    enabled: bool,
    method: str,
    ensemble_methods: Sequence[str] | None,
    ensemble_min_k: int,
) -> tuple[str, ...] | None:
    """Validate ensemble options and return members when ensemble DE is enabled."""
    if not enabled:
        return None

    method_name = method.strip().lower()
    if method_name == "ensemble":
        return resolve_ensemble_methods(ensemble_methods, ensemble_min_k)
    if ensemble_methods is not None:
        raise ValueError("ensemble methods only apply when DE method is 'ensemble'")
    return None


def resolve_combined_results(
    results: dict[str, ResultValue],
    min_k: int,
) -> dict[str, ResultValue]:
    """Return consensus results keyed by validated canonical result labels."""
    if isinstance(min_k, bool) or not isinstance(min_k, Integral) or min_k < 1:
        raise ValueError(
            "ensemble min_k must be an integer greater than or equal to 1; "
            f"got {min_k!r}"
        )

    canonical: dict[str, ResultValue] = {}
    for label, result in results.items():
        if not isinstance(label, str) or not label.strip():
            raise ValueError("ensemble result labels must be non-empty strings")
        name = label.strip().lower()
        if name in canonical:
            raise ValueError(f"duplicate ensemble result key '{name}'")
        canonical[name] = result
    return canonical
