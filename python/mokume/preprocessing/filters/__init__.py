"""
Preprocessing filters for the mokume package.

This module provides a comprehensive set of filters for quality control
and preprocessing of proteomics data at multiple levels:
- Intensity-based filtering
- Peptide-level filtering
- Protein-level filtering
- Run/Sample QC filtering
"""

from mokume.preprocessing.filters.base import BaseFilter, FilterResult
from mokume.preprocessing.filters.enums import (
    RazorPeptideHandling,
    ProteinGroupingStrategy,
    FilterLevel,
    FilterAction,
)
from mokume.preprocessing.filters.pipeline import FilterPipeline
from mokume.preprocessing.filters.io import (
    is_filter_config_available,
    load_filter_config,
    save_filter_config,
    generate_example_config,
)
from mokume.preprocessing.filters.factory import (
    create_intensity_filters,
    create_peptide_filters,
    create_protein_filters,
    create_run_qc_filters,
    get_filter_pipeline,
)

__all__ = [
    # Base classes
    "BaseFilter",
    "FilterResult",
    "FilterPipeline",
    # Enums
    "RazorPeptideHandling",
    "ProteinGroupingStrategy",
    "FilterLevel",
    "FilterAction",
    # I/O functions
    "is_filter_config_available",
    "load_filter_config",
    "save_filter_config",
    "generate_example_config",
    # Factory functions
    "create_intensity_filters",
    "create_peptide_filters",
    "create_protein_filters",
    "create_run_qc_filters",
    "get_filter_pipeline",
]
