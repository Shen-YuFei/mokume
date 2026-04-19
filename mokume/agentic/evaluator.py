"""Metric computation engine for agentic analysis."""

import numpy as np
import pandas as pd

from mokume.agentic.config import ScoreWeights
from mokume.agentic.state import CandidateConfig, EvaluationResult
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.evaluator")


def _direction_aware_tp_fp(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    expected_direction: str = "UP",
) -> tuple[int, int, int]:
    """Count TP/FP/FN using direction-aware significance."""
    sig = de_df[de_df["significance"] == expected_direction]
    protein_col = "protein" if "protein" in de_df.columns else de_df.columns[0]
    sig_proteins = set(sig[protein_col]) if protein_col in sig.columns else set()

    # Fallback: use index if protein column not found
    if not sig_proteins and not sig.empty:
        sig_proteins = set(sig.index)

    tp = len(sig_proteins & ground_truth)
    fp = len(sig_proteins - ground_truth)
    fn = len(ground_truth - sig_proteins)
    return tp, fp, fn


def _compute_auc(
    de_df: pd.DataFrame,
    ground_truth: set[str],
) -> float | None:
    """Compute ROC AUC for ground truth proteins."""
    try:
        from sklearn.metrics import roc_auc_score  # pylint: disable=import-outside-toplevel
    except ImportError:
        return None

    protein_col = "protein" if "protein" in de_df.columns else de_df.columns[0]
    if protein_col not in de_df.columns:
        return None

    labels = de_df[protein_col].isin(ground_truth).astype(int).values
    if labels.sum() == 0 or labels.sum() == len(labels):
        return None

    pvals = de_df["pvalue"].fillna(1.0).values
    scores = 1.0 - pvals  # higher score = more significant
    return float(roc_auc_score(labels, scores))


def _de_counts(de_df: pd.DataFrame) -> tuple[int, int]:
    """Count UP and DOWN significant proteins."""
    n_up = int((de_df.get("significance", pd.Series()) == "UP").sum())
    n_down = int((de_df.get("significance", pd.Series()) == "DOWN").sum())
    return n_up, n_down


def _median_cv_from_matrix(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
) -> float | None:
    """Compute median CV from protein matrix."""
    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)
    conditions = set(sample_to_condition.values())
    all_cvs = []
    for cond in conditions:
        cols = [s for s in matrix.columns if sample_to_condition.get(s) == cond]
        if len(cols) < 2:
            continue
        subset = matrix[cols]
        row_mean = subset.mean(axis=1)
        row_std = subset.std(axis=1, ddof=1)
        cv = row_std / row_mean.replace(0, np.nan)
        all_cvs.extend(cv.dropna().tolist())
    return float(np.median(all_cvs)) if all_cvs else None


def compute_score_ground_truth(
    result: EvaluationResult,
    max_tp: int,
    total_tested: int,
    weights: ScoreWeights,
) -> float:
    """Compute composite score for ground truth mode (Mode A)."""
    auc_term = weights.w_auc * (result.auc or 0.0)
    tp_term = weights.w_tp * ((result.tp or 0) / max(max_tp, 1))
    fp_term = weights.w_fp * ((result.fp or 0) / max(total_tested, 1))
    return auc_term + tp_term - fp_term


def compute_score_unsupervised(
    result: EvaluationResult,
    max_de: int,
    weights: ScoreWeights,
) -> float:
    """Compute composite score for unsupervised mode (Mode B)."""
    n_de = result.n_de_up + result.n_de_down
    de_term = weights.w_de * (n_de / max(max_de, 1))
    cv_term = weights.w_cv * (1.0 - min(result.median_cv or 1.0, 1.0))
    miss_term = weights.w_miss * result.missing_rate
    return de_term + cv_term - miss_term


def _fill_ground_truth(
    result: EvaluationResult,
    de_df: pd.DataFrame,
    ground_truth: set[str],
) -> None:
    """Populate ground truth metrics on an EvaluationResult."""
    tp, fp, fn = _direction_aware_tp_fp(de_df, ground_truth)
    result.tp = tp
    result.fp = fp
    result.fn = fn
    result.auc = _compute_auc(de_df, ground_truth)
    result.sensitivity = tp / max(tp + fn, 1)
    result.specificity = 1.0 - (fp / max(tp + fp + fn, 1))


def evaluate(
    config: CandidateConfig,
    de_df: pd.DataFrame,
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    ground_truth: set[str] | None = None,
) -> EvaluationResult:
    """Evaluate a single experiment result."""
    n_up, n_down = _de_counts(de_df)
    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)
    missing = float(matrix.isna().sum().sum() / matrix.size) if matrix.size else 0.0

    result = EvaluationResult(
        config_name=config.name, config=config.to_dict(),
        n_de_up=n_up, n_de_down=n_down,
        median_cv=_median_cv_from_matrix(protein_df, sample_to_condition),
        missing_rate=missing,
    )

    if ground_truth:
        _fill_ground_truth(result, de_df, ground_truth)

    logger.info(
        "Evaluated %s: UP=%d DOWN=%d TP=%s FP=%s AUC=%s",
        config.name, n_up, n_down, result.tp, result.fp, result.auc,
    )
    return result
