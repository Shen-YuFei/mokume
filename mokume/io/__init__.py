"""
Input/Output utilities for the mokume package.

This module provides utilities for reading and writing various file formats
including Parquet, TSV/CSV, and FASTA files, as well as DuckDB-based
feature data access.
"""

from mokume.io.parquet import (
    create_anndata,
    combine_ibaq_tsv_files,
)
from mokume.io.fasta import (
    load_fasta,
    digest_protein,
    extract_fasta,
    get_protein_molecular_weights,
)
from mokume.io.feature import (
    Feature,
    SQLFilterBuilder,
)

__all__ = [
    "create_anndata",
    "combine_ibaq_tsv_files",
    "load_fasta",
    "digest_protein",
    "extract_fasta",
    "get_protein_molecular_weights",
    "Feature",
    "SQLFilterBuilder",
]
