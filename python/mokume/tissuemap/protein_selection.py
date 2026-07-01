"""
Protein filtering: NaN rate, contaminants, and PTM artifacts.
"""

from __future__ import annotations

import logging
import re

import numpy as np
import pandas as pd

from mokume.tissuemap.config import FilteringConfig

logger = logging.getLogger(__name__)


def filter_proteins(
    mat: pd.DataFrame,
    config: FilteringConfig,
) -> pd.DataFrame:
    """Filter proteins by NaN fraction and contaminant patterns.

    Parameters
    ----------
    mat : pd.DataFrame
        Proteins (rows) x samples (columns), log2-normalized values.
    config : FilteringConfig
        Filtering thresholds.

    Returns
    -------
    pd.DataFrame
        Filtered matrix (subset of rows).
    """
    n_before = mat.shape[0]

    # 1. NaN fraction filter
    nan_frac = mat.isna().mean(axis=1)
    keep_nan = nan_frac <= config.max_nan_frac
    n_nan_dropped = (~keep_nan).sum()

    # 2. Contaminant filter
    if config.remove_contaminants and config.contaminant_pattern:
        pattern = re.compile(config.contaminant_pattern)
        is_contam = mat.index.to_series().apply(lambda x: bool(pattern.search(str(x))))
        keep_clean = ~is_contam.values
        n_contam = is_contam.sum()
    else:
        keep_clean = np.ones(mat.shape[0], dtype=bool)
        n_contam = 0

    keep = keep_nan.values & keep_clean
    mat_filtered = mat.loc[keep]

    logger.info(
        "Protein selection: %d → %d (NaN>%.0f%%: -%d, contaminant: -%d)",
        n_before,
        mat_filtered.shape[0],
        config.max_nan_frac * 100,
        n_nan_dropped,
        n_contam,
    )
    return mat_filtered
