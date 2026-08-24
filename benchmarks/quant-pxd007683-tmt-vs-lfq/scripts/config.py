"""
Configuration for PXD007683 benchmark.

This benchmark compares TMT vs LFQ quantification on the same biological samples
with known yeast spike-in ratios (ground truth for fold-change accuracy).

IMPORTANT: DirectLFQ Dependency
-------------------------------
When using DirectLFQ quantification, always use the external directlfq package
directly (pip install directlfq) rather than any fallback implementation.
This ensures reproducibility and uses the official Mann Lab algorithm.

Install with: pip install mokume-py[directlfq]
"""

from pathlib import Path

# =============================================================================
# Paths
# =============================================================================

BENCHMARK_DIR = Path(__file__).parent.parent
DATA_DIR = BENCHMARK_DIR / "data"
RESULTS_DIR = BENCHMARK_DIR / "results"
FIGURES_DIR = BENCHMARK_DIR / "figures"

# Use existing processed data from benchmarks-local
# These files contain peptide-level data ready for quantification
BENCHMARKS_LOCAL_DIR = Path(__file__).parent.parent.parent.parent / "benchmarks-local" / "PXD007683"
LOCAL_DATA_DIR = BENCHMARKS_LOCAL_DIR / "data" / "processed"

# Pre-quantified results (protein-level)
LOCAL_QUANTIFIED_DIR = BENCHMARKS_LOCAL_DIR / "results" / "quantification"

# Data URLs (QPX format from quantms) - for future reference
DATA_URLS = {
    "tmt": "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-benchmark/PXD007683-TMT/PXD007683-TMT.parquet",
    "lfq": "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-benchmark/PXD007683-LFQ/PXD007683-LFQ.parquet",
}

# FASTA file for iBAQ (human + yeast)
FASTA_URL = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-benchmark/PXD007683-LFQ/PXD007683-LFQ.fasta"
FASTA_PATH = DATA_DIR / "PXD007683.fasta"

# =============================================================================
# Experimental Design (Ground Truth)
# =============================================================================

# Sample to condition mapping (yeast spike-in percentage)
# Condition names based on yeast content: 10%, 5%, 3.3%
# Sample IDs in the existing data are just numbers: "1", "2", ..., "11"
SAMPLE_CONDITIONS = {
    # Samples 1-3: 10% yeast (3x relative to baseline)
    "1": "QY_10pct",
    "2": "QY_10pct",
    "3": "QY_10pct",
    # Samples 4-7: 5% yeast (1.5x relative to baseline)
    "4": "QY_5pct",
    "5": "QY_5pct",
    "6": "QY_5pct",
    "7": "QY_5pct",
    # Samples 8-11: 3.3% yeast (baseline)
    "8": "QY_3pct",
    "9": "QY_3pct",
    "10": "QY_3pct",
    "11": "QY_3pct",
}

# Expected fold-changes for YEAST proteins (ground truth)
# Based on spike-in ratios: 10% / 3.3% ≈ 3.0, 10% / 5% = 2.0, 5% / 3.3% ≈ 1.5
EXPECTED_FOLD_CHANGES = {
    ("QY_10pct", "QY_3pct"): 3.0,   # 10% vs 3.3%
    ("QY_10pct", "QY_5pct"): 2.0,   # 10% vs 5%
    ("QY_5pct", "QY_3pct"): 1.5,    # 5% vs 3.3%
}

# Species identification patterns
SPECIES_PATTERNS = {
    "human": ["HUMAN", "_HUMAN", "Homo sapiens"],
    "yeast": ["YEAST", "_YEAST", "Saccharomyces"],
}

# =============================================================================
# Methods to Test
# =============================================================================

# Quantification methods by technology
# NOTE: DirectLFQ should NOT be used for TMT data - it was designed for LFQ
# DirectLFQ assumptions (MS1 precursor intensities, missing values from stochastic
# sampling) do not apply to TMT (MS2/MS3 reporter ions, multiplexed samples)

QUANTIFICATION_METHODS_TMT = [
    "sum",          # Sum of all peptides (baseline)
    "top3",         # Top 3 peptides
    "top5",         # Top 5 peptides
    "top10",        # Top 10 peptides
    "maxlfq",       # MaxLFQ (peptide trace alignment)
    "ibaq",         # iBAQ (requires FASTA)
]

QUANTIFICATION_METHODS_LFQ = [
    "sum",          # Sum of all peptides (baseline)
    "top3",         # Top 3 peptides
    "top5",         # Top 5 peptides
    "top10",        # Top 10 peptides
    "maxlfq",       # MaxLFQ (peptide trace alignment)
    "directlfq",    # DirectLFQ (designed for LFQ data)
    "ibaq",         # iBAQ (requires FASTA)
]

# Combined list for backward compatibility
QUANTIFICATION_METHODS = list(set(QUANTIFICATION_METHODS_TMT + QUANTIFICATION_METHODS_LFQ))

# Sample-level normalization methods
NORMALIZATION_METHODS = [
    "none",              # No normalization
    "median",            # Median centering (log space)
    "quantile",          # Quantile normalization
    "total_intensity",   # Total intensity normalization
    "hierarchical",      # Hierarchical sample alignment
]

# Imputation methods
IMPUTATION_METHODS = [
    "none",      # No imputation (complete cases)
    "knn",       # K-nearest neighbors
    "min",       # Minimum value / 2 (MNAR assumption)
    "median",    # Median per protein
]

# =============================================================================
# Quality Thresholds
# =============================================================================

# Minimum peptides per protein
MIN_PEPTIDES = 2

# CV thresholds for quality assessment
CV_THRESHOLDS = {
    "excellent": 0.10,  # < 10%
    "good": 0.20,       # < 20%
    "acceptable": 0.30, # < 30%
    "poor": 0.50,       # < 50%
}

# =============================================================================
# Plotting Configuration
# =============================================================================

FIGURE_DPI = 150
FIGURE_FORMAT = "png"

# Color palettes
CONDITION_COLORS = {
    "QY_10pct": "#e41a1c",   # Red (highest yeast)
    "QY_5pct": "#377eb8",    # Blue (middle)
    "QY_3pct": "#4daf4a",    # Green (baseline)
}

TECHNOLOGY_COLORS = {
    "TMT": "#ff7f00",
    "LFQ": "#984ea3",
}

METHOD_COLORS = {
    "sum": "#a6cee3",
    "top3": "#1f78b4",
    "top5": "#b2df8a",
    "top10": "#33a02c",
    "maxlfq": "#fb9a99",
    "directlfq": "#e31a1c",
    "ibaq": "#fdbf6f",
}
