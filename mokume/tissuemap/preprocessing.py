"""
Log2 + median normalization and tissue label harmonization.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Explicit sub-region overrides.  None → exclude the sample.
_SUBREGION_MAP: dict[str, str | None] = {
    "colon; transverse colon": None,
    "esophagus; gastroesophageal junction": None,
    "minor salivary gland": "salivary gland",
    "heart ; atrial appendage": "heart",
    "heart; atrial appendage": "heart",
}

# Synonym mapping applied after canonicalization
_TISSUE_SYNONYMS: dict[str, str] = {
    "urinary bladder": "bladder",
    "thyroid gland": "thyroid",
    "skeletal muscle": "muscle",
    "gut": "intestine",
    "small intestine": "intestine",
    "tube": "fallopian tube",
}


def canonicalize_tissue(raw: str) -> str | None:
    """Return canonical tissue name or ``None`` to exclude the sample.

    1. Checks ``_SUBREGION_MAP`` for explicit overrides / exclusions.
    2. Falls back to the first token before ``";"``, stripped and lowercased.
    3. Applies ``_TISSUE_SYNONYMS``.
    """
    if not isinstance(raw, str) or not raw.strip():
        return None
    key = raw.strip().lower()
    if key in _SUBREGION_MAP:
        return _SUBREGION_MAP[key]
    canonical = key.split(";")[0].strip()
    return _TISSUE_SYNONYMS.get(canonical, canonical)


def harmonize_tissues(
    meta: pd.DataFrame, *, min_samples: int = 2
) -> pd.DataFrame:
    """Canonicalize tissue labels and drop rare tissues.

    Parameters
    ----------
    meta : pd.DataFrame
        Must contain a ``tissue`` column.
    min_samples : int
        Tissues with fewer samples are dropped.

    Returns
    -------
    pd.DataFrame
        Filtered copy of *meta* with updated ``tissue`` and an added
        ``tissue_original`` column preserving the raw label.
    """
    meta = meta.copy()
    meta["tissue_original"] = meta["tissue"].copy()
    meta["tissue"] = meta["tissue"].map(
        lambda t: canonicalize_tissue(str(t)) if pd.notna(t) else None
    )
    # Drop samples with excluded tissue (None)
    meta = meta[meta["tissue"].notna()]

    # Drop tissues with too few samples
    tc = meta["tissue"].value_counts()
    keep = tc[tc >= min_samples].index
    dropped = tc[tc < min_samples].index
    if len(dropped) > 0:
        logger.info(
            "Dropping %d tissues with <%d samples: %s",
            len(dropped),
            min_samples,
            list(dropped),
        )
    return meta[meta["tissue"].isin(keep)].copy()



def log2_median_normalize(mat: pd.DataFrame) -> pd.DataFrame:
    """Log2 transform + sample-wise median centering.

    Parameters
    ----------
    mat : pd.DataFrame
        Proteins (rows) x samples (columns), raw intensities.
        Zeros are treated as missing.

    Returns
    -------
    pd.DataFrame
        Normalized matrix (same shape), with NaN where input was 0 / NaN.
    """
    log_mat = np.log2(mat.replace(0, np.nan))
    medians = log_mat.median(axis=0)
    grand_median = medians.median()
    norm_mat = log_mat.sub(medians, axis=1).add(grand_median)
    n_nan = norm_mat.isna().sum().sum()
    total = norm_mat.size
    logger.info(
        "log2+median norm: %d proteins x %d samples, NaN=%.1f%%",
        norm_mat.shape[0],
        norm_mat.shape[1],
        100.0 * n_nan / total if total else 0,
    )
    return norm_mat
