"""
Normalization implementations for the mokume package.

Only the IRS (Internal Reference Scaling) module survives in this self-contained
delegate copy: the rest of the normalization family plus the ``io.feature``
(duckdb) and SDRF / aggregation preprocessing re-exports are Rust-kernel-redundant
and are not vendored here. The IRS symbols are resolved lazily.
"""

# The IRS public symbols (the only ones the delegate commands use, via
# ``detect_condition_from_sdrf``).
_LAZY_NAMES = {
    "IRSNormalizer",
    "detect_pooled_from_sdrf",
    "detect_reference_by_column",
    "detect_reference_by_regex",
    "detect_plexes_from_sdrf",
    "detect_condition_from_sdrf",
}


def __getattr__(name):
    """Lazily import normalization symbols from the surviving ``irs`` module."""
    if name in _LAZY_NAMES:
        from mokume.normalization import irs

        return getattr(irs, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "IRSNormalizer",
    "detect_pooled_from_sdrf",
    "detect_reference_by_column",
    "detect_reference_by_regex",
    "detect_plexes_from_sdrf",
    "detect_condition_from_sdrf",
]
