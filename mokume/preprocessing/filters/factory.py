"""
Factory functions for creating preprocessing filters.
"""

from typing import List

from mokume.model.filters import (
    PreprocessingFilterConfig,
    IntensityFilterConfig,
    PeptideFilterConfig,
    ProteinFilterConfig,
    RunQCFilterConfig,
)
from mokume.preprocessing.filters.base import BaseFilter
from mokume.preprocessing.filters.intensity import (
    MinIntensityFilter,
    CVThresholdFilter,
    ReplicateAgreementFilter,
    QuantileFilter,
)
from mokume.preprocessing.filters.peptide import (
    PeptideLengthFilter,
    ChargeStateFilter,
    ModificationFilter,
    MissedCleavageFilter,
    SearchScoreFilter,
    SequencePatternFilter,
    FDRFilter as PeptideFDRFilter,
)
from mokume.preprocessing.filters.protein import (
    ContaminantFilter,
    MinPeptideFilter,
    ProteinFDRFilter,
    CoverageFilter,
    RazorPeptideFilter,
)
from mokume.preprocessing.filters.run_qc import (
    RunIntensityFilter,
    MinFeaturesFilter,
    MissingRateFilter,
    SampleCorrelationFilter,
)
from mokume.preprocessing.filters.pipeline import FilterPipeline


def create_intensity_filters(config: IntensityFilterConfig) -> List[BaseFilter]:
    """
    Create intensity filters from configuration.

    Parameters
    ----------
    config : IntensityFilterConfig
        Intensity filter configuration.

    Returns
    -------
    List[BaseFilter]
        List of configured intensity filters.
    """
    filters = []

    if config.remove_zero_intensity or config.min_intensity > 0:
        min_val = max(
            config.min_intensity, 1e-10 if config.remove_zero_intensity else 0
        )
        filters.append(MinIntensityFilter(min_intensity=min_val))

    if config.cv_threshold is not None:
        filters.append(CVThresholdFilter(cv_threshold=config.cv_threshold))

    if config.min_replicate_agreement > 1:
        filters.append(
            ReplicateAgreementFilter(min_replicates=config.min_replicate_agreement)
        )

    if config.quantile_lower > 0 or config.quantile_upper < 1:
        filters.append(
            QuantileFilter(
                lower_quantile=config.quantile_lower,
                upper_quantile=config.quantile_upper,
            )
        )

    return filters


def create_peptide_filters(config: PeptideFilterConfig) -> List[BaseFilter]:
    """
    Create peptide filters from configuration.

    Parameters
    ----------
    config : PeptideFilterConfig
        Peptide filter configuration.

    Returns
    -------
    List[BaseFilter]
        List of configured peptide filters.
    """
    filters = []

    if config.min_peptide_length > 0 or config.max_peptide_length < 100:
        filters.append(
            PeptideLengthFilter(
                min_length=config.min_peptide_length,
                max_length=config.max_peptide_length,
            )
        )

    if config.allowed_charge_states:
        filters.append(ChargeStateFilter(allowed_charges=config.allowed_charge_states))

    if config.exclude_modifications:
        filters.append(
            ModificationFilter(exclude_modifications=config.exclude_modifications)
        )

    if config.max_missed_cleavages is not None:
        filters.append(
            MissedCleavageFilter(max_missed_cleavages=config.max_missed_cleavages)
        )

    if config.min_search_score is not None:
        filters.append(SearchScoreFilter(min_score=config.min_search_score))

    if config.exclude_sequence_patterns:
        filters.append(
            SequencePatternFilter(exclude_patterns=config.exclude_sequence_patterns)
        )

    # Peptide FDR filter - only add if column exists (checked at apply time)
    if config.fdr_threshold < 1.0:
        filters.append(PeptideFDRFilter(fdr_threshold=config.fdr_threshold))

    return filters


def create_protein_filters(config: ProteinFilterConfig) -> List[BaseFilter]:
    """
    Create protein filters from configuration.

    Parameters
    ----------
    config : ProteinFilterConfig
        Protein filter configuration.

    Returns
    -------
    List[BaseFilter]
        List of configured protein filters.
    """
    filters = []

    if config.remove_contaminants or config.remove_decoys:
        filters.append(
            ContaminantFilter(
                patterns=config.contaminant_patterns,
                remove_decoys=config.remove_decoys,
            )
        )

    if config.min_peptides > 0 or config.min_unique_peptides > 0:
        filters.append(
            MinPeptideFilter(
                min_peptides=config.min_peptides,
                min_unique_peptides=config.min_unique_peptides,
            )
        )

    if config.fdr_threshold < 1.0:
        filters.append(ProteinFDRFilter(fdr_threshold=config.fdr_threshold))

    if config.min_coverage > 0:
        filters.append(CoverageFilter(min_coverage=config.min_coverage))

    if config.razor_peptide_handling != "keep":
        filters.append(RazorPeptideFilter(handling=config.razor_peptide_handling))

    return filters


def create_run_qc_filters(config: RunQCFilterConfig) -> List[BaseFilter]:
    """
    Create run QC filters from configuration.

    Parameters
    ----------
    config : RunQCFilterConfig
        Run QC filter configuration.

    Returns
    -------
    List[BaseFilter]
        List of configured run QC filters.
    """
    filters = []

    if config.min_total_intensity > 0:
        filters.append(RunIntensityFilter(min_intensity=config.min_total_intensity))

    if config.min_identified_features > 0 or config.min_identified_proteins > 0:
        filters.append(
            MinFeaturesFilter(
                min_features=config.min_identified_features,
                min_proteins=config.min_identified_proteins,
            )
        )

    if config.max_missing_rate < 1.0:
        filters.append(MissingRateFilter(max_missing_rate=config.max_missing_rate))

    if config.min_sample_correlation is not None:
        filters.append(
            SampleCorrelationFilter(min_correlation=config.min_sample_correlation)
        )

    return filters


def get_filter_pipeline(config: PreprocessingFilterConfig) -> FilterPipeline:
    """
    Create a complete filter pipeline from configuration.

    The pipeline applies filters in the following order:
    1. Run/Sample QC filters (remove bad runs first)
    2. Intensity filters
    3. Peptide filters
    4. Protein filters

    Parameters
    ----------
    config : PreprocessingFilterConfig
        Complete preprocessing filter configuration.

    Returns
    -------
    FilterPipeline
        Configured filter pipeline ready to apply.
    """
    pipeline = FilterPipeline(name=config.name)

    if not config.enabled:
        return pipeline

    # Add filters in recommended order
    # 1. Run QC first (remove bad samples early)
    pipeline.add_filters(create_run_qc_filters(config.run_qc))

    # 2. Intensity filters
    pipeline.add_filters(create_intensity_filters(config.intensity))

    # 3. Peptide filters
    pipeline.add_filters(create_peptide_filters(config.peptide))

    # 4. Protein filters (after peptide filtering)
    pipeline.add_filters(create_protein_filters(config.protein))

    return pipeline
