"""
Protein quantification methods for the mokume package.

This module provides implementations for various protein quantification
methods including iBAQ, Top3, TopN, MaxLFQ, DirectLFQ, and AllPeptides.

DirectLFQ is an optional dependency. Install with:
    pip install mokume[directlfq]
"""

from mokume.quantification.base import ProteinQuantificationMethod
from mokume.quantification.ibaq import (
    peptides_to_protein,
    normalize_ibaq,
    extract_fasta,
    ConcentrationWeightByProteomicRuler,
)
from mokume.quantification.top3 import Top3Quantification
from mokume.quantification.topn import TopNQuantification
from mokume.quantification.maxlfq import MaxLFQQuantification
from mokume.quantification.all_peptides import AllPeptidesQuantification
from mokume.quantification.ratio import RatioQuantification
from mokume.quantification.tmt_abundance import TMTAbundanceQuantification
from mokume.quantification.tmt_reporter import TMTReporterIntensityQuantification
from mokume.quantification.spectral_count import SpectralCountQuantification

# Lazy import for optional DirectLFQ
from mokume.quantification.directlfq import is_directlfq_available

__all__ = [
    # Base class
    "ProteinQuantificationMethod",
    # iBAQ
    "peptides_to_protein",
    "normalize_ibaq",
    "extract_fasta",
    "ConcentrationWeightByProteomicRuler",
    # Quantification methods
    "Top3Quantification",
    "TopNQuantification",
    "MaxLFQQuantification",
    "AllPeptidesQuantification",
    "RatioQuantification",
    "TMTAbundanceQuantification",
    "TMTReporterIntensityQuantification",
    "SpectralCountQuantification",
    # DirectLFQ (optional)
    "DirectLFQQuantification",
    "is_directlfq_available",
    # Factory function
    "get_quantification_method",
    "list_quantification_methods",
]


def __getattr__(name):
    """Lazy import for optional dependencies."""
    if name == "DirectLFQQuantification":
        from mokume.quantification.directlfq import DirectLFQQuantification

        return DirectLFQQuantification
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def get_quantification_method(method: str, **kwargs) -> ProteinQuantificationMethod:
    """
    Get a quantification method instance by name.

    Parameters
    ----------
    method : str
        Name of the quantification method. One of:
        'topN' (where N is any number, e.g., 'top3', 'top5', 'top10'),
        'maxlfq', 'directlfq', 'all', 'sum',
        'abd' / 'abundance', 'intensity' / 'reporter', 'ratio'.
    **kwargs
        Additional arguments passed to the quantification method constructor.

        For MaxLFQ:
            - min_peptides: int (default 2)
            - n_jobs: int (default -1, all cores)

        For TopN:
            - n: int (default 3, can also be parsed from method name)

        For DirectLFQ:
            - min_nonan: int (default 1)
            - num_cores: int (default None)
            - deactivate_normalization: bool (default False)

        For 'ratio' (TMT ratio quantification):
            - reference_samples: list[str] (required)
            - sample_to_plex: dict[str, str] (required)
            - fraction_merge_method: str (default 'mean')

        For 'abd' / 'abundance' and 'intensity' / 'reporter':
            No additional arguments.

    Returns
    -------
    ProteinQuantificationMethod
        An instance of the requested quantification method.

    Raises
    ------
    ValueError
        If the method name is not recognized.
    ImportError
        If DirectLFQ is requested but not installed.

    Examples
    --------
    >>> method = get_quantification_method("maxlfq", min_peptides=2, n_jobs=4)
    >>> result = method.quantify(peptide_df, ...)

    >>> # TopN with any N
    >>> method = get_quantification_method("top5")  # Top 5 peptides
    >>> method = get_quantification_method("top10")  # Top 10 peptides

    >>> # DirectLFQ requires optional install
    >>> method = get_quantification_method("directlfq", min_nonan=2)
    """
    import re

    method_lower = method.lower()

    # Handle topN methods (top3, top5, top10, etc.)
    if method_lower.startswith("top"):
        match = re.match(r"top(\d+)", method_lower)
        if match:
            n = int(match.group(1))
        else:
            n = kwargs.get("n", 3)
        return TopNQuantification(n=n)

    elif method_lower == "maxlfq":
        return MaxLFQQuantification(
            min_peptides=kwargs.get("min_peptides", 2),
            threads=kwargs.get("n_jobs", -1),
            verbose=kwargs.get("verbose", 0),
        )

    elif method_lower == "directlfq":
        from mokume.quantification.directlfq import DirectLFQQuantification

        return DirectLFQQuantification(
            min_nonan=kwargs.get("min_nonan", 1),
            num_cores=kwargs.get("num_cores", None),
            deactivate_normalization=kwargs.get("deactivate_normalization", False),
        )

    elif method_lower in ("all", "sum", "allpeptides"):
        return AllPeptidesQuantification()

    elif method_lower in ("abd", "abundance", "tmtabundance"):
        return TMTAbundanceQuantification()

    elif method_lower in ("intensity", "reporter", "tmtreporterintensity"):
        return TMTReporterIntensityQuantification()

    elif method_lower in ("spectralcount", "spectral_count", "count"):
        return SpectralCountQuantification()

    elif method_lower == "ratio":
        try:
            reference_samples = kwargs["reference_samples"]
            sample_to_plex = kwargs["sample_to_plex"]
        except KeyError as exc:
            raise ValueError(
                "RatioQuantification requires 'reference_samples' and "
                "'sample_to_plex' kwargs."
            ) from exc
        return RatioQuantification(
            reference_samples=reference_samples,
            sample_to_plex=sample_to_plex,
            fraction_merge_method=kwargs.get("fraction_merge_method", "mean"),
        )

    else:
        available = (
            "topN (e.g., top3, top5, top10), maxlfq, directlfq, all/sum, "
            "abd/abundance, intensity/reporter, ratio"
        )
        raise ValueError(
            f"Unknown quantification method: {method}. Available methods: {available}"
        )


def list_quantification_methods() -> dict:
    """
    List all available quantification methods.

    Returns
    -------
    dict
        Dictionary mapping method names to their availability status.
        For optional dependencies, shows whether they are installed.

    Examples
    --------
    >>> from mokume.quantification import list_quantification_methods
    >>> methods = list_quantification_methods()
    >>> print(methods)
    {'top3': True, 'topn': True, 'maxlfq': True, 'directlfq': False, 'sum': True}
    """
    return {
        "top3": True,
        "topn": True,
        "maxlfq": True,
        "directlfq": is_directlfq_available(),
        "sum": True,
        "abd": True,
        "intensity": True,
        "ratio": True,
        "spectral_count": True,
    }
