"""
Post-processing utilities for the mokume package.

This module provides utilities for data reshaping, batch correction,
and multi-file combining.
"""

from mokume.postprocessing.reshape import (
    pivot_wider,
    pivot_longer,
    remove_samples_low_protein_number,
    remove_missing_values,
    describe_expression_metrics,
)
from mokume.postprocessing.batch_correction import (
    apply_batch_correction,
    compute_pca,
    get_batch_info_from_sample_names,
    remove_single_sample_batches,
    iterative_outlier_removal,
    is_batch_correction_available,
    is_inmoose_available,
    detect_batches,
    extract_covariates_from_sdrf,
)
from mokume.postprocessing.combiner import Combiner

__all__ = [
    # Reshape
    "pivot_wider",
    "pivot_longer",
    "remove_samples_low_protein_number",
    "remove_missing_values",
    "describe_expression_metrics",
    # Batch correction
    "apply_batch_correction",
    "compute_pca",
    "get_batch_info_from_sample_names",
    "remove_single_sample_batches",
    "iterative_outlier_removal",
    "is_batch_correction_available",
    "is_inmoose_available",
    "detect_batches",
    "extract_covariates_from_sdrf",
    # Combiner
    "Combiner",
]
