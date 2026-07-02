"""
Sequential KNN (SeqKNN) imputation.

Implements the sequential k-nearest neighbours imputation method of
Kim, Kim & Yi (2004) BMC Bioinformatics 5:160.

Algorithm
---------
1. Order rows (proteins) by missing-value count, ascending.
2. Sample the rows with **fewest** missing values into a complete set.
3. For each remaining row, impute its missing entries using the
   weighted mean of its ``k`` nearest neighbours from the **already
   complete** set.
4. After a row is imputed, append it to the complete set so it can
   serve as a neighbour for subsequent rows.

This sequential structure is the key difference from the standard
sklearn ``KNNImputer`` (a single pass over all rows simultaneously)
and matches the Bioconductor ``SeqKnn`` semantics.

References
----------
- Kim, K.-Y., Kim, B.-J., & Yi, G.-S. (2004). Reuse of imputed data
  in microarray analysis increases imputation efficiency.
  *BMC Bioinformatics*, 5, 160.
"""

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.imputation.seqknn")


def _weighted_neighbour_mean(
    target_row: np.ndarray,
    complete_rows: np.ndarray,
    k: int,
) -> np.ndarray:
    """Impute missing entries of ``target_row`` using k nearest complete rows."""
    missing_mask = np.isnan(target_row)
    observed_mask = ~missing_mask
    if not observed_mask.any() or len(complete_rows) == 0:
        return target_row

    # Distances on observed columns only (Euclidean in observed subspace)
    diffs = complete_rows[:, observed_mask] - target_row[observed_mask]
    distances = np.sqrt(np.sum(diffs * diffs, axis=1))
    # Avoid divide-by-zero for identical rows
    distances = np.where(distances == 0, 1e-9, distances)

    k_use = min(k, len(complete_rows))
    nearest_idx = np.argsort(distances)[:k_use]
    weights = 1.0 / distances[nearest_idx]
    neighbours = complete_rows[nearest_idx][:, missing_mask]
    weighted_sum = np.sum(neighbours * weights[:, None], axis=0)
    imputed = weighted_sum / weights.sum()

    result = target_row.copy()
    result[missing_mask] = imputed
    return result


def _seed_with_column_means(matrix: np.ndarray, sorted_matrix: np.ndarray) -> int:
    """Seed the first sorted row with column means so distance can start."""
    seed = sorted_matrix[0].copy()
    col_means = np.nanmean(matrix, axis=0)
    seed[np.isnan(seed)] = col_means[np.isnan(seed)]
    sorted_matrix[0] = seed
    logger.info("SeqKNN: no fully-observed rows; seeded with column means")
    return 1


def impute_seqknn(data: pd.DataFrame, k: int = 10) -> pd.DataFrame:
    """Sequentially impute missing values with weighted k-NN.

    ``data`` is a DataFrame with NaN-encoded missings; ``k`` is the
    neighbour count. Returns a same-shape DataFrame with NaNs filled.
    """
    if data.empty:
        return data.copy()
    if k < 1:
        raise ValueError("k must be >= 1")

    matrix = data.to_numpy(dtype=float, copy=True)
    missing_per_row = np.sum(np.isnan(matrix), axis=1)
    if missing_per_row.max() == 0:
        return data.copy()
    if (missing_per_row == matrix.shape[1]).all():
        logger.warning("SeqKNN: every row is fully missing; cannot impute")
        return data.copy()

    sort_order = np.argsort(missing_per_row, kind="stable")
    sorted_matrix = matrix[sort_order]
    sorted_missing = missing_per_row[sort_order]
    n_complete = int(np.sum(sorted_missing == 0))
    if n_complete == 0:
        n_complete = _seed_with_column_means(matrix, sorted_matrix)

    complete = sorted_matrix[:n_complete].copy()
    for idx in range(n_complete, len(sorted_matrix)):
        sorted_matrix[idx] = _weighted_neighbour_mean(sorted_matrix[idx], complete, k)
        complete = np.vstack([complete, sorted_matrix[idx][None, :]])

    restored = sorted_matrix[np.argsort(sort_order)]
    logger.info(
        "SeqKNN: imputed %d rows with k=%d neighbours",
        int(np.sum(missing_per_row > 0)),
        k,
    )
    return pd.DataFrame(restored, index=data.index, columns=data.columns)
