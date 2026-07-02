"""
pyCombat batch correction with NaN-safe handling.
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd
from combat.pycombat import pycombat

logger = logging.getLogger(__name__)


def _batches_valid(
    batch_labels: np.ndarray,
    min_batches: int = 2,
    min_per_batch: int = 2,
) -> bool:
    """Check that batch labels satisfy minimum requirements."""
    unique, counts = np.unique(batch_labels, return_counts=True)
    return bool(len(unique) >= min_batches and (counts >= min_per_batch).all())


def run_pycombat(
    log2_data: np.ndarray,
    sample_ids: np.ndarray,
    batch_labels: np.ndarray,
) -> np.ndarray:
    """Run pyCombat, handling NaN by temporary per-protein median fill.

    NaN positions are restored after correction so that downstream steps
    can distinguish truly missing values from corrected values.

    Parameters
    ----------
    log2_data : np.ndarray
        Samples (rows) x proteins (columns), log2 scale.
    sample_ids : np.ndarray
        Sample identifiers (used as DataFrame column names).
    batch_labels : np.ndarray
        Batch label per sample.

    Returns
    -------
    np.ndarray
        Batch-corrected matrix, same shape, NaN restored.
    """
    # pycombat expects genes (rows) x samples (columns)
    mat_t = log2_data.T.copy()
    nan_mask = np.isnan(mat_t)

    # Fill NaN with per-protein (row) median for combat stability
    row_medians = np.nanmedian(mat_t, axis=1, keepdims=True)
    row_medians = np.where(np.isnan(row_medians), 0, row_medians)
    mat_filled = np.where(nan_mask, row_medians, mat_t)

    df = pd.DataFrame(mat_filled, columns=sample_ids)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = pycombat(df, batch_labels)

    result_arr = result.values.T  # back to samples x proteins
    # Restore original NaN positions
    result_arr[np.isnan(log2_data)] = np.nan
    return result_arr


def batch_correct(
    log2_matrix: np.ndarray,
    sample_ids: np.ndarray,
    batch_labels: np.ndarray,
) -> np.ndarray:
    """Apply batch correction if batches are valid.

    Parameters
    ----------
    log2_matrix : np.ndarray
        Samples x proteins, log2 scale.
    sample_ids : np.ndarray
        Sample identifiers.
    batch_labels : np.ndarray
        Batch label per sample.

    Returns
    -------
    np.ndarray
        Corrected (or original) matrix.
    """
    if not _batches_valid(batch_labels):
        logger.warning(
            "Batch correction skipped: need >= 2 batches with >= 2 "
            "samples each (got %d unique batches)",
            len(np.unique(batch_labels)),
        )
        return log2_matrix

    logger.info(
        "Running pyCombat: %d samples, %d batches",
        len(sample_ids),
        len(np.unique(batch_labels)),
    )
    try:
        corrected = run_pycombat(log2_matrix, sample_ids, batch_labels)
        logger.info("Batch correction done")
        return corrected
    except (ValueError, TypeError, RuntimeError, ArithmeticError, KeyError) as exc:
        logger.error("pyCombat failed: %s — returning uncorrected", exc)
        return log2_matrix
