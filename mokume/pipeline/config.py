"""Pipeline configuration dataclasses.

This module provides nested configuration dataclasses for the
quantification pipeline, replacing the flat PipelineConfig.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class InputConfig:

    """Input file paths."""

    parquet: str
    sdrf: Optional[str] = None
    fasta_file: Optional[str] = None


@dataclass
class FilterConfig:

    """Peptide/protein filtering parameters."""

    min_aa: int = 7
    min_unique_peptides: int = 2
    remove_contaminants: bool = True


@dataclass
class NormalizationConfig:

    """Normalization method parameters."""

    run_method: str = "median"
    sample_method: str = "globalMedian"
    proteins_file: Optional[str] = None


@dataclass
class QuantificationConfig:

    """Quantification method parameters."""

    method: str = "maxlfq"
    ion_alignment: Optional[str] = None
    coverage_threshold: Optional[float] = None
    ratio_fraction_merge: str = "mean"
    # DirectLFQ-specific
    directlfq_num_cores: Optional[int] = None
    directlfq_min_nonan: int = 1
    directlfq_num_samples_quadratic: int = 50


@dataclass
class IRSConfig:

    """IRS (Internal Reference Scaling) normalization parameters."""

    enabled: bool = False
    reference_samples: Optional[list] = None
    sdrf_column: Optional[str] = None
    sdrf_values: Optional[list] = None
    reference_regex: str = "pool|powder|ref|reference|bridge"
    stat: str = "median"
    remove_reference: bool = False


@dataclass
class BatchCorrectionConfig:

    """Batch correction parameters."""

    enabled: bool = False
    method: str = "sample_prefix"
    column: Optional[str] = None
    covariates: Optional[list] = None
    parametric: bool = True
    mean_only: bool = False
    ref_batch: Optional[int] = None


@dataclass
class ImputationConfig:

    """Missing value imputation parameters."""

    enabled: bool = False
    method: str = "none"
    # MinProb parameters
    quantile: float = 0.01
    shift: float = 1.6
    scale: float = 0.3
    # KNN parameters
    n_neighbors: int = 5


@dataclass
class DEConfig:

    """Differential expression analysis parameters."""

    enabled: bool = False
    contrasts: Optional[list] = None
    method: str = "auto"
    log2fc_threshold: float = 0.5
    fdr_threshold: float = 0.05
    fdr_method: str = "bh"
    output: Optional[str] = None


@dataclass
class OutputConfig:

    """Output and export parameters."""

    export_peptides: Optional[str] = None
    export_ions: Optional[str] = None
    plot_dir: Optional[str] = None
    plot_volcano: bool = False
    plot_heatmap: bool = False
    plot_pca: bool = False
    highlight_genes: Optional[list] = None
    interactive_report: bool = False
    report_output: Optional[str] = None


@dataclass
class PipelineConfig:

    """Configuration for the quantification pipeline.

    Organizes settings into logical groups:
    - input: File paths (parquet, sdrf, fasta)
    - filtering: Peptide/protein filters
    - normalization: Run and sample normalization methods
    - quantification: Quant method and method-specific params
    - irs: Internal Reference Scaling for multi-plex TMT
    - batch: Batch correction (ComBat)
    - imputation: Missing value imputation (MinProb, MinDet, KNN)
    - de: Differential expression analysis
    - output: Export paths, plotting, reports
    """

    input: InputConfig
    filtering: FilterConfig = field(default_factory=FilterConfig)
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    quantification: QuantificationConfig = field(default_factory=QuantificationConfig)
    irs: IRSConfig = field(default_factory=IRSConfig)
    batch: BatchCorrectionConfig = field(default_factory=BatchCorrectionConfig)
    imputation: ImputationConfig = field(default_factory=ImputationConfig)
    de: DEConfig = field(default_factory=DEConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
