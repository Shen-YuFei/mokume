"""
Preprocessing utilities for the mokume package.

This module provides preprocessing filters, SDRF analysis, and data
aggregation tools for proteomics data analysis.
"""

from mokume.preprocessing.filters import (
    is_filter_config_available,
    load_filter_config,
    save_filter_config,
    generate_example_config,
    get_filter_pipeline,
    FilterPipeline,
    FilterResult,
    BaseFilter,
)
from mokume.preprocessing.sdrf import analyse_sdrf
from mokume.preprocessing.aggregation import (
    parse_uniprot_accession,
    get_canonical_peptide,
    remove_contaminants_entrapments_decoys,
    remove_protein_by_ids,
    reformat_quantms_feature_table_quant_labels,
    apply_initial_filtering,
    merge_fractions,
    get_peptidoform_normalize_intensities,
    sum_peptidoform_intensities,
)

__all__ = [
    # Filters
    "is_filter_config_available",
    "load_filter_config",
    "save_filter_config",
    "generate_example_config",
    "get_filter_pipeline",
    "FilterPipeline",
    "FilterResult",
    "BaseFilter",
    # SDRF analysis
    "analyse_sdrf",
    # Aggregation
    "parse_uniprot_accession",
    "get_canonical_peptide",
    "remove_contaminants_entrapments_decoys",
    "remove_protein_by_ids",
    "reformat_quantms_feature_table_quant_labels",
    "apply_initial_filtering",
    "merge_fractions",
    "get_peptidoform_normalize_intensities",
    "sum_peptidoform_intensities",
]
