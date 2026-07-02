"""
Input/Output utilities for the mokume package.

Only FASTA digestion survives in this self-contained delegate copy: the parquet
/ anndata IO and the duckdb-based feature backend were Rust-kernel-redundant and
have been removed. The FASTA symbols are resolved lazily so importing the
package does not force pyOpenMS to load until a symbol is actually accessed.
"""

_LAZY_NAMES = {
    "load_fasta",
    "digest_protein",
    "extract_fasta",
    "get_protein_molecular_weights",
}


def __getattr__(name):
    """Lazily import IO symbols from the surviving ``fasta`` module."""
    if name in _LAZY_NAMES:
        from mokume.io import fasta

        return getattr(fasta, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "load_fasta",
    "digest_protein",
    "extract_fasta",
    "get_protein_molecular_weights",
]
