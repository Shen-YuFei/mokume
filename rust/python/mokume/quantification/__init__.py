"""
Protein quantification methods for the mokume package.

Only the piBAQ peptide->protein path survives in this self-contained delegate
copy: the other quantification methods (Top3/TopN/MaxLFQ/DirectLFQ/AllPeptides/
Ratio/TMT/SpectralCount) and the ``get_quantification_method`` factory were
Rust-kernel-redundant and have been removed. The piBAQ symbols are resolved
lazily so importing the package does not force ``pibaq`` (pyOpenMS) to load until
a symbol is actually accessed.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .pibaq import compute_pibaq, normalize_pibaq

# The piBAQ public symbols and their defining submodule.
_LAZY_SUBMODULE = {
    "peptides_to_protein": "pibaq",
    "compute_pibaq": "pibaq",
    "normalize_pibaq": "pibaq",
    "extract_fasta": "pibaq",
    "ConcentrationWeightByProteomicRuler": "pibaq",
}


def __getattr__(name):
    """Lazily import quantification symbols from their defining submodule."""
    submodule = _LAZY_SUBMODULE.get(name)
    if submodule is not None:
        import importlib

        module = importlib.import_module(f"mokume.quantification.{submodule}")
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "peptides_to_protein",
    "compute_pibaq",
    "normalize_pibaq",
    "extract_fasta",
    "ConcentrationWeightByProteomicRuler",
]
