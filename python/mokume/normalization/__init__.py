"""
Normalization implementations for the mokume package.

This module provides implementations for feature-level, peptide-level,
protein-level, and hierarchical sample normalization.
"""

from mokume.normalization.feature import (
    normalize_runs,
    normalize_sample,
    normalize_replicates,
)
from mokume.normalization.peptide import (
    peptide_normalization,
)
from mokume.preprocessing.sdrf import analyse_sdrf
from mokume.preprocessing.aggregation import (
    remove_contaminants_entrapments_decoys,
    remove_protein_by_ids,
)
from mokume.io.feature import Feature
from mokume.normalization.hierarchical import (
    HierarchicalSampleNormalizer,
    HierarchicalIonAligner,
    DistanceMetric,
)
from mokume.normalization.irs import (
    IRSNormalizer,
    detect_pooled_from_sdrf,
    detect_reference_by_column,
    detect_reference_by_regex,
    detect_plexes_from_sdrf,
    detect_condition_from_sdrf,
)
from mokume.normalization.loess import (
    LOESSNormalizer,
    loess_normalize,
)
from mokume.normalization.protein import (
    median_center,
    quantile_normalize,
)
from mokume.normalization.rlr import (
    RlrNormalizer,
    rlr_normalize,
)
from mokume.normalization.distribution import (
    MeanCenterNormalizer,
    MedianCenterNormalizer,
    QuantileNormalizer,
)

__all__ = [
    # Feature normalization
    "normalize_runs",
    "normalize_sample",
    "normalize_replicates",
    # Peptide normalization
    "peptide_normalization",
    "analyse_sdrf",
    "remove_contaminants_entrapments_decoys",
    "remove_protein_by_ids",
    "Feature",
    # Hierarchical normalization (DirectLFQ-style, native mokume)
    "HierarchicalSampleNormalizer",
    "HierarchicalIonAligner",
    "DistanceMetric",
    # IRS normalization (multi-plex TMT)
    "IRSNormalizer",
    "detect_pooled_from_sdrf",
    "detect_reference_by_column",
    "detect_reference_by_regex",
    "detect_plexes_from_sdrf",
    "detect_condition_from_sdrf",
    # LOESS normalization
    "LOESSNormalizer",
    "loess_normalize",
    # Distribution-alignment normalization (quantile / median center / mean center)
    "quantile_normalize",
    "median_center",
    "QuantileNormalizer",
    "MedianCenterNormalizer",
    "MeanCenterNormalizer",
    # RLR normalization (NormalyzerDE)
    "RlrNormalizer",
    "rlr_normalize",
]
