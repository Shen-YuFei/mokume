"""Experiment executor for agentic analysis."""

import numpy as np
import pandas as pd

from mokume.agentic.state import CandidateConfig
from mokume.analysis.differential_expression import DifferentialExpression
from mokume.core.logger import get_logger
from mokume.imputation.censored import impute_censored

logger = get_logger("mokume.agentic.runner")


def _apply_imputation(
    protein_df: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    """Apply imputation to the protein matrix (log2 space)."""
    if method == "none":
        return protein_df

    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)

    # Ensure log2 space for imputation
    log2_matrix = matrix.copy()
    if log2_matrix.max().max() > 40:
        log2_matrix = np.log2(log2_matrix.replace(0, np.nan))

    imputed = impute_censored(log2_matrix, method=method)

    # Convert back if we transformed
    if matrix.max().max() > 40:
        imputed = 2**imputed

    return imputed.reset_index()


def run_experiment(
    config: CandidateConfig,
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    contrast: tuple[str, str],
    peptide_counts: pd.Series | None = None,
) -> pd.DataFrame:
    """Run a single DE experiment with the given configuration."""
    logger.info("Running experiment: %s", config.name)

    # 1. Apply imputation
    imputed_df = _apply_imputation(protein_df, config.imputation)

    # 2. Run DE
    de = DifferentialExpression(
        method=config.de_method,
        log2fc_threshold=config.log2fc_threshold,
        fdr_threshold=0.05,
        fdr_method=config.fdr_method,
        peptide_counts=peptide_counts,
    )

    result = de.run(imputed_df, sample_to_condition, contrast)
    logger.info(
        "Experiment %s complete: %d proteins tested",
        config.name,
        len(result),
    )
    return result
