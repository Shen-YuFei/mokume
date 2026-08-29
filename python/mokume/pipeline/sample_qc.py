"""Sample-level quality control on quantified protein matrices."""

from collections import defaultdict

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.pipeline.sample_qc")

_MIN_CORRELATION_OVERLAP = 3


def _is_reference_condition(condition: str) -> bool:
    condition = condition.lower()
    return "powder" in condition or "pool" in condition


def _condition_samples(
    sample_columns: list[str], sample_to_condition: dict[str, str]
) -> dict[str, list[str]]:
    missing = [sample for sample in sample_columns if sample not in sample_to_condition]
    if missing:
        raise ValueError(
            "Sample correlation filtering found no condition metadata for "
            + ", ".join(repr(sample) for sample in missing)
        )
    grouped: dict[str, list[str]] = defaultdict(list)
    for sample in sample_columns:
        condition = sample_to_condition[sample]
        if not _is_reference_condition(condition):
            grouped[condition].append(sample)
    if not grouped:
        raise ValueError(
            "Sample correlation filtering found no biological samples with "
            "condition metadata"
        )
    return dict(grouped)


def _log2_matrix(
    protein_df: pd.DataFrame, sample_columns: list[str], values_are_log2: bool
) -> pd.DataFrame:
    numeric = protein_df[sample_columns].apply(pd.to_numeric, errors="coerce")
    if values_are_log2:
        return numeric.where(np.isfinite(numeric))
    positive = numeric.where(np.isfinite(numeric) & (numeric > 0.0))
    return np.log2(positive)


def _validate_pair_correlations(
    log2_matrix: pd.DataFrame,
    correlations: pd.DataFrame,
    condition: str,
    sample: str,
    peers: list[str],
) -> None:
    for peer in peers:
        if pd.notna(correlations.at[sample, peer]):
            continue
        overlap = int(log2_matrix[[sample, peer]].notna().all(axis=1).sum())
        raise ValueError(
            f"Sample correlation between {sample!r} and {peer!r} in "
            f"condition {condition!r} is undefined: {overlap} "
            "pairwise-complete usable proteins "
            f"(minimum {_MIN_CORRELATION_OVERLAP})"
        )


def _condition_exclusions(
    log2_matrix: pd.DataFrame,
    condition: str,
    samples: list[str],
    threshold: float,
) -> list[str]:
    if len(samples) < 2:
        raise ValueError(
            "Sample correlation filtering requires at least two samples "
            f"in condition {condition!r}"
        )
    correlations = log2_matrix[samples].corr(
        method="pearson", min_periods=_MIN_CORRELATION_OVERLAP
    )
    excluded = []
    for sample in samples:
        peers = [peer for peer in samples if peer != sample]
        _validate_pair_correlations(log2_matrix, correlations, condition, sample, peers)
        mean_correlation = float(correlations.loc[sample, peers].mean())
        logger.info(
            "Sample correlation QC: sample=%s condition=%s "
            "mean_correlation=%.6f peers=%d threshold=%.6f",
            sample,
            condition,
            mean_correlation,
            len(peers),
            threshold,
        )
        if mean_correlation < threshold:
            excluded.append(sample)
    return excluded


def filter_samples_by_correlation(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    threshold: float,
    *,
    values_are_log2: bool = False,
) -> pd.DataFrame:
    """Drop samples below mean Pearson correlation to same-condition peers.

    Linear intensities are restricted to positive values and log2-transformed;
    already-log2 methods use all finite values directly. Every condition must
    contain at least two samples and each pair must share at least three
    proteins with non-constant values.
    """
    if not np.isfinite(threshold) or not -1.0 <= threshold <= 1.0:
        raise ValueError("min_sample_correlation must be between -1 and 1")
    protein_column = protein_df.columns[0]
    sample_columns = [
        column for column in protein_df.columns if column != protein_column
    ]
    grouped = _condition_samples(sample_columns, sample_to_condition)
    log2_matrix = _log2_matrix(protein_df, sample_columns, values_are_log2)
    excluded = []
    for condition, samples in grouped.items():
        excluded.extend(
            _condition_exclusions(log2_matrix, condition, samples, threshold)
        )

    logger.info(
        "Sample correlation filtering complete: evaluated_samples=%d "
        "excluded_samples=%d threshold=%.6f",
        sum(len(samples) for samples in grouped.values()),
        len(excluded),
        threshold,
    )
    return protein_df.drop(columns=excluded)
