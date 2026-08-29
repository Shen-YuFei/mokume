"""Metric computation engine for agentic analysis."""

import math

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from mokume.agentic.adaptive_fdr import adjust_pvalues, estimate_pi0
from mokume.agentic.profiler import compute_median_cv
from mokume.agentic.state import (
    CandidateConfig,
    EvaluationResult,
    FdrCalibrationMetrics,
    GroundTruthMetrics,
)
from mokume.core.logger import get_logger

logger = get_logger("mokume.agentic.evaluator")

# Alpha grid for the recall-vs-empirical-FDR curve: 1% to 30% in 1% steps.
# The lower end covers the strict operating points people actually publish at
# (1-5%); the upper end is where a spike-in benchmark's recall has usually
# saturated, so extending it further adds no information.
EMP_FDR_ALPHA_GRID: np.ndarray = np.round(np.arange(0.01, 0.31, 0.01), 2)

# The adaptive-pi0 procedure used for the BH-vs-adaptive guardrail report.
ADAPTIVE_FDR_METHOD = "storey"


def _protein_col(de_df: pd.DataFrame) -> str | None:
    """Resolve the protein identifier column of a DE table."""
    if "protein" in de_df.columns:
        return "protein"
    return str(de_df.columns[0]) if len(de_df.columns) else None


def _direction_mask(de_df: pd.DataFrame, expected_direction: str | None) -> np.ndarray:
    """Boolean mask of rows moving in the expected direction.

    Falls back to "everything passes" when the table carries no ``log2FC``
    column or no direction is expected. Non-finite fold changes never pass.
    """
    n = len(de_df)
    if expected_direction not in ("UP", "DOWN") or "log2FC" not in de_df.columns:
        return np.ones(n, dtype=bool)
    fc = de_df["log2FC"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        return fc > 0 if expected_direction == "UP" else fc < 0


def _de_tp_fp(
    de_df: pd.DataFrame,
    ground_truth: set[str],
) -> tuple[int, int, int, int]:
    """Count TP/FP/FN/TN from two-sided differential calls."""
    sig = de_df[de_df["significance"].isin(("UP", "DOWN"))]
    protein_col = "protein" if "protein" in de_df.columns else de_df.columns[0]
    sig_proteins = set(sig[protein_col]) if protein_col in sig.columns else set()
    tested_proteins = set(de_df[protein_col]) if protein_col in de_df.columns else set()

    # Fallback: use index if protein column not found
    if not sig_proteins and not sig.empty:
        sig_proteins = set(sig.index)
    if not tested_proteins and not de_df.empty:
        tested_proteins = set(de_df.index)

    tested_truth = tested_proteins & ground_truth
    tp = len(sig_proteins & tested_truth)
    fp = len(sig_proteins - tested_truth)
    fn = len(tested_truth - sig_proteins)
    # Tested negatives (not in ground truth) that were not called significant.
    tn = len((tested_proteins - tested_truth) - sig_proteins)
    return tp, fp, fn, tn


def _truth_direction_diagnostics(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    expected_direction: str,
) -> tuple[int, int, float | None]:
    """Count correctly and oppositely directed calls among truth proteins."""
    protein_col = _protein_col(de_df)
    if protein_col is None or expected_direction not in {"UP", "DOWN"}:
        return 0, 0, None
    opposite_direction = "DOWN" if expected_direction == "UP" else "UP"
    truth = de_df[de_df[protein_col].isin(ground_truth)]
    correct = int(truth["significance"].eq(expected_direction).sum())
    incorrect = int(truth["significance"].eq(opposite_direction).sum())
    called = correct + incorrect
    accuracy = correct / called if called else None
    return correct, incorrect, accuracy


def _ranking_scores(de_df: pd.DataFrame, expected_direction: str | None) -> np.ndarray:
    """Significance ranking scores, optionally direction-aware.

    The base score is ``-log_pvalue`` when a stable log representation is
    available, otherwise ``-pvalue`` (higher = more significant). When a
    direction is expected, proteins moving the *wrong* way are demoted below
    every right-moving protein by offsetting their score, rather than being
    dropped: dropping rows would change the label set and therefore the AUC
    denominator, which would no longer be comparable across candidates. The
    offset only re-orders, so any rank-based metric (ROC AUC, partial AUC)
    stays well defined.
    """
    if "log_pvalue" in de_df.columns:
        log_pvalue = de_df["log_pvalue"].to_numpy(dtype=float)
        scores = np.where(np.isfinite(log_pvalue), -np.minimum(log_pvalue, 0.0), 0.0)
    else:
        scores = -de_df["pvalue"].fillna(1.0).to_numpy(dtype=float)
    if expected_direction not in ("UP", "DOWN"):
        return scores
    offset = float(np.ptp(scores)) + 1.0
    return np.where(
        _direction_mask(de_df, expected_direction),
        scores,
        scores - offset,
    )


def _auc_inputs(
    de_df: pd.DataFrame,
    ground_truth: set[str],
) -> tuple[np.ndarray, np.ndarray] | None:
    """Labels and ranking scores for ROC metrics, or None if not computable."""
    protein_col = _protein_col(de_df)
    if protein_col is None or "pvalue" not in de_df.columns:
        return None
    labels = de_df[protein_col].isin(ground_truth).astype(int).to_numpy()
    if labels.sum() == 0 or labels.sum() == len(labels):
        return None
    return labels, _ranking_scores(de_df, None)


def _compute_auc(
    de_df: pd.DataFrame,
    ground_truth: set[str],
) -> float | None:
    """Compute direction-blind ROC AUC for ground-truth proteins."""
    inputs = _auc_inputs(de_df, ground_truth)
    if inputs is None:
        return None
    labels, scores = inputs
    return float(roc_auc_score(labels, scores))


def _pauc(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    max_fpr: float = 0.1,
) -> float | None:
    """Partial ROC AUC over the low-false-positive region (FPR <= max_fpr)."""
    inputs = _auc_inputs(de_df, ground_truth)
    if inputs is None:
        return None
    labels, scores = inputs
    return float(roc_auc_score(labels, scores, max_fpr=max_fpr))


def _emp_fdr_ladder(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    expected_direction: str = "UP",
) -> tuple[np.ndarray, np.ndarray, int] | None:
    """Walk the p-value-ranked list and return (emp_fdr, tp_cum, n_gt).

    Proteins moving in the expected direction are ranked by p-value. The
    empirical FDR is evaluated only after a complete group of equal ranking
    scores, because no p-value threshold can select part of a tie. This keeps
    the result independent of input row order. ``n_gt`` is the number of
    ground-truth proteins actually tested. Returns None when there is nothing
    to score against.
    """
    protein_col = _protein_col(de_df)
    if protein_col is None or "pvalue" not in de_df.columns:
        return None
    n_gt = len(set(de_df[protein_col]) & ground_truth)
    if n_gt == 0:
        return None

    mask = (
        _direction_mask(de_df, expected_direction) & de_df["pvalue"].notna().to_numpy()
    )
    d = de_df.loc[mask]
    if d.empty:
        return np.empty(0), np.empty(0, dtype=int), n_gt
    scores = _ranking_scores(d, None)
    order = np.argsort(-scores, kind="stable")
    d = d.iloc[order]
    scores = scores[order]

    is_tp = d[protein_col].isin(ground_truth).to_numpy()
    tp_cum = np.cumsum(is_tp)
    depth = np.arange(1, len(d) + 1)
    group_ends = np.r_[scores[1:] != scores[:-1], True]
    return (
        ((depth - tp_cum) / depth)[group_ends],
        tp_cum[group_ends],
        n_gt,
    )


def _recall_from_ladder(
    ladder: tuple[np.ndarray, np.ndarray, int], alpha: float
) -> float:
    """Recall at the deepest cut whose empirical FDR stays <= alpha."""
    emp_fdr, tp_cum, n_gt = ladder
    ok = emp_fdr <= alpha
    max_tp = int(tp_cum[ok].max()) if ok.any() else 0
    return max_tp / n_gt


def _recall_at_emp_fdr(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    alpha: float = 0.05,
    expected_direction: str = "UP",
) -> float | None:
    """Recall at a controlled empirical FDR, using the background as the null.

    Returns the recall (over all detected ground-truth proteins) at the deepest
    cut whose empirical FDR stays <= alpha. This separates the recall the
    ranking can deliver at a calibrated error rate from the fixed operating
    point. See :func:`recall_at_emp_fdr_curve` for the whole curve.
    """
    ladder = _emp_fdr_ladder(de_df, ground_truth, expected_direction)
    return None if ladder is None else _recall_from_ladder(ladder, alpha)


def recall_at_emp_fdr_curve(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    alphas: np.ndarray | list[float] | None = None,
    expected_direction: str = "UP",
) -> list[tuple[float, float]] | None:
    """Recall-vs-empirical-FDR curve over an alpha grid.

    The scalar :func:`_recall_at_emp_fdr` only reports one operating point; the
    curve shows *where* two candidates actually differ — one can win at
    alpha=0.01 and lose at alpha=0.2. Defaults to :data:`EMP_FDR_ALPHA_GRID`.
    Returns ``[(alpha, recall), ...]``, or None when there is no ground truth
    among the tested proteins.
    """
    ladder = _emp_fdr_ladder(de_df, ground_truth, expected_direction)
    if ladder is None:
        return None
    grid = EMP_FDR_ALPHA_GRID if alphas is None else np.asarray(alphas, dtype=float)
    return [(float(a), _recall_from_ladder(ladder, float(a))) for a in grid]


def _adjust_pvalues(
    pvalues: np.ndarray, method: str, alpha: float
) -> tuple[np.ndarray, str, float]:
    """Adjust p-values by "bh" or "adaptive", returning (adjusted, method_used, pi0).

    Delegates to the agentic evaluation copy of the validated adaptive FDR code,
    already gates the adaptive procedure on pi0 reliability and reports the
    method it actually applied ("bh" on fallback). pi0 is estimated alongside
    purely so the report can show *why* the two branches differ.
    """
    if method != "adaptive":
        adjusted, method_used = adjust_pvalues(pvalues, "bh", alpha)
        return adjusted, method_used, 1.0

    finite = pvalues[np.isfinite(pvalues)]
    pi0 = estimate_pi0(finite) if finite.size else 1.0
    adjusted, method_used = adjust_pvalues(pvalues, ADAPTIVE_FDR_METHOD, alpha)
    return adjusted, method_used, pi0


def compare_fdr_calibration(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    alpha: float = 0.05,
    expected_direction: str = "UP",
) -> dict | None:
    """Report BH and adaptive-pi0 empirical FDR/recall side by side at one alpha.

    Both procedures are applied to the *same* p-values, so the pair shows how
    much extra recall the adaptive correction buys and whether its empirical FDR
    is still at or under the nominal alpha. This is a reporting guardrail, not a
    search axis: the candidate's own ``fdr_method`` is untouched.

    Returns a dict keyed "bh" and "adaptive", each with ``n_called`` (proteins
    called in the expected direction at adjusted p < alpha), ``emp_fdr``
    (fraction of those not in the ground truth), ``recall`` (over tested
    ground-truth proteins), ``method_used`` (the adaptive entry reads
    :data:`ADAPTIVE_FDR_METHOD`, or "bh" when the pi0 reliability gate forced a
    fallback) and ``pi0``. Returns None when no ground-truth protein was tested.
    """
    protein_col = _protein_col(de_df)
    if protein_col is None or "pvalue" not in de_df.columns:
        return None
    n_gt = len(set(de_df[protein_col]) & ground_truth)
    if n_gt == 0:
        return None

    pvalues = de_df["pvalue"].to_numpy(dtype=float)
    direction = _direction_mask(de_df, expected_direction)
    out: dict = {"alpha": float(alpha), "n_ground_truth_tested": n_gt}
    for method in ("bh", "adaptive"):
        adj, method_used, pi0 = _adjust_pvalues(pvalues, method, alpha)
        called = set(
            de_df.loc[direction & np.isfinite(adj) & (adj < alpha), protein_col]
        )
        n_called = len(called)
        out[method] = {
            "n_called": n_called,
            "emp_fdr": (len(called - ground_truth) / n_called) if n_called else 0.0,
            "recall": len(called & ground_truth) / n_gt,
            "method_used": method_used,
            "pi0": float(pi0),
        }
    return out


def de_signed_call_set(
    de_df: pd.DataFrame,
    expected_direction: str | None = None,
) -> set[tuple[str, str]]:
    """Return significant calls as ``(protein, direction)`` pairs.

    Direction is part of the identity: a protein called UP by one configuration
    and DOWN by another is disagreement, not replication. ``None`` includes
    both canonical directions.
    """
    protein_col = _protein_col(de_df)
    if protein_col is None or "significance" not in de_df.columns:
        return set()
    directions = (
        {expected_direction} if expected_direction in {"UP", "DOWN"} else {"UP", "DOWN"}
    )
    called = de_df.loc[
        de_df["significance"].isin(directions),
        [protein_col, "significance"],
    ]
    return {
        (str(protein), str(direction))
        for protein, direction in called.itertuples(index=False, name=None)
    }


def _method_sensitivity_rows(
    calls: dict[str, set[tuple[str, str]]],
) -> list[dict[str, int | float | str]]:
    """Build one support row for each signed call in the union."""
    candidate_names = list(calls)
    candidate_count = len(candidate_names)
    union = set().union(*calls.values()) if calls else set()
    rows = []
    for protein, direction in sorted(union):
        supporting = [
            name for name in candidate_names if (protein, direction) in calls[name]
        ]
        support_count = len(supporting)
        classification = "unassessed"
        if candidate_count >= 2:
            classification = (
                "shared" if support_count == candidate_count else "method_sensitive"
            )
        rows.append(
            {
                "protein": protein,
                "direction": direction,
                "support_count": support_count,
                "candidate_count": candidate_count,
                "support_fraction": support_count / candidate_count,
                "supporting_candidates": ";".join(supporting),
                "classification": classification,
            }
        )
    return rows


def _method_sensitivity_summary(
    table: pd.DataFrame,
    candidate_count: int,
) -> dict[str, int | bool | str]:
    """Summarize shared and method-sensitive signed calls."""
    return {
        "comparison_available": candidate_count >= 2,
        "candidate_count": candidate_count,
        "signed_call_union": len(table),
        "shared_signed_calls": int((table["classification"] == "shared").sum()),
        "method_sensitive_signed_calls": int(
            (table["classification"] == "method_sensitive").sum()
        ),
        "interpretation": (
            "Call sharing describes sensitivity to the tested methods and is not "
            "evidence of biological truth."
        ),
    }


def method_sensitivity(
    candidate_de: dict[str, pd.DataFrame],
    expected_direction: str | None = None,
) -> tuple[pd.DataFrame, dict[str, int | bool | str]]:
    """Describe shared and method-sensitive signed calls without ranking."""
    columns = [
        "protein",
        "direction",
        "support_count",
        "candidate_count",
        "support_fraction",
        "supporting_candidates",
        "classification",
    ]
    calls = {
        name: de_signed_call_set(table, expected_direction)
        for name, table in candidate_de.items()
    }
    candidate_count = len(calls)
    table = pd.DataFrame(_method_sensitivity_rows(calls), columns=columns)
    summary = _method_sensitivity_summary(table, candidate_count)
    return table, summary


def _de_counts(de_df: pd.DataFrame) -> tuple[int, int]:
    """Count UP and DOWN significant proteins."""
    n_up = int((de_df.get("significance", pd.Series()) == "UP").sum())
    n_down = int((de_df.get("significance", pd.Series()) == "DOWN").sum())
    return n_up, n_down


def normalized_mcc(tp: int, fp: int, fn: int, tn: int) -> float:
    """Matthews correlation coefficient rescaled to [0, 1].

    MCC is the only confusion-matrix summary that uses all four cells, which
    is what keeps it honest when the classes are as lopsided as a spike-in
    benchmark -- tens of true positives against thousands of background
    proteins. An undefined MCC (a whole row or column is zero) maps to 0.5,
    the no-information value.
    """
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if den == 0:
        return 0.5
    return ((tp * tn - fp * fn) / den + 1.0) / 2.0


def compute_score_ground_truth(
    result: EvaluationResult,
) -> float | None:
    """Return the absolute Score A compatibility summary.

    This is the arithmetic mean of pAUC at the 1%, 5%, and 10% false-positive
    cutoffs plus normalized MCC. Candidate selection is deliberately separate:
    the service ranks all candidates per metric, includes G-mean, and then
    averages those ranks.

    A candidate is scorable only when all four terms are finite. Averaging a
    partial subset would make candidates with different available evidence
    incomparable. Averaging is otherwise legitimate because all four terms
    already live on [0, 1].
    """
    metrics = result.truth_metrics
    terms = [
        metrics.pauc001,
        metrics.pauc005,
        metrics.pauc,
        metrics.nmcc,
    ]
    if any(term is None or not math.isfinite(term) for term in terms):
        return None
    return float(sum(terms) / len(terms))


def _ground_truth_metrics(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    alpha: float = 0.05,
    expected_direction: str = "UP",
) -> tuple[GroundTruthMetrics, FdrCalibrationMetrics]:
    """Build OpDEA-style Score A metrics and directional diagnostics."""
    tp, fp, fn, tn = _de_tp_fp(de_df, ground_truth)
    sensitivity = tp / max(tp + fn, 1)
    specificity = tn / max(tn + fp, 1)
    direction_correct, direction_incorrect, direction_accuracy = (
        _truth_direction_diagnostics(de_df, ground_truth, expected_direction)
    )
    metrics = GroundTruthMetrics(
        tp=tp,
        fp=fp,
        fn=fn,
        tn=tn,
        auc=_compute_auc(de_df, ground_truth),
        sensitivity=sensitivity,
        specificity=specificity,
        truth_direction_correct=direction_correct,
        truth_direction_incorrect=direction_incorrect,
        truth_direction_accuracy=direction_accuracy,
        recall_at_emp_fdr=_recall_at_emp_fdr(
            de_df, ground_truth, alpha, expected_direction
        ),
        pauc=_pauc(de_df, ground_truth),
        pauc001=_pauc(de_df, ground_truth, max_fpr=0.01),
        pauc005=_pauc(de_df, ground_truth, max_fpr=0.05),
        nmcc=normalized_mcc(tp, fp, fn, tn),
        gmean=math.sqrt(sensitivity * specificity),
        recall_emp_fdr_curve=tuple(
            recall_at_emp_fdr_curve(
                de_df,
                ground_truth,
                expected_direction=expected_direction,
            )
        ),
    )
    return metrics, _fdr_calibration(
        de_df,
        ground_truth,
        alpha,
        expected_direction,
    )


def _fdr_calibration(
    de_df: pd.DataFrame,
    ground_truth: set[str],
    alpha: float,
    expected_direction: str,
) -> FdrCalibrationMetrics:
    """Build the BH-vs-adaptive side-by-side fields."""
    comparison = compare_fdr_calibration(de_df, ground_truth, alpha, expected_direction)
    if comparison is None:
        return FdrCalibrationMetrics()
    bh, adaptive = comparison["bh"], comparison["adaptive"]
    return FdrCalibrationMetrics(
        n_called_bh=bh["n_called"],
        emp_fdr_bh=bh["emp_fdr"],
        recall_at_alpha_bh=bh["recall"],
        n_called_adaptive=adaptive["n_called"],
        emp_fdr_adaptive=adaptive["emp_fdr"],
        recall_at_alpha_adaptive=adaptive["recall"],
        adaptive_method_used=adaptive["method_used"],
        adaptive_pi0=adaptive["pi0"],
    )


def evaluate(
    config: CandidateConfig,
    de_df: pd.DataFrame,
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    truth: tuple[set[str], str] | None = None,
) -> EvaluationResult:
    """Evaluate a single experiment result.

    Ground-truth metrics are computed at the candidate's own alpha
    (``config.fdr_threshold``), not a hard-coded 0.05.
    """
    n_up, n_down = _de_counts(de_df)
    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)
    missing = float(matrix.isna().sum().sum() / matrix.size) if matrix.size else 0.0

    truth_metrics = GroundTruthMetrics()
    fdr_calibration = FdrCalibrationMetrics()
    if truth is not None:
        ground_truth, expected_direction = truth
        if expected_direction not in {"UP", "DOWN"}:
            raise ValueError(
                "expected_direction must be UP or DOWN when ground truth is supplied"
            )
        truth_metrics, fdr_calibration = _ground_truth_metrics(
            de_df,
            ground_truth,
            alpha=float(config.fdr_threshold),
            expected_direction=expected_direction,
        )

    result = EvaluationResult(
        config_name=config.name,
        config=config.to_dict(),
        truth_metrics=truth_metrics,
        fdr_calibration=fdr_calibration,
        n_de_up=n_up,
        n_de_down=n_down,
        median_cv=compute_median_cv(
            protein_df.set_index(protein_col), sample_to_condition
        ),
        missing_rate=missing,
    )

    logger.info(
        "Evaluated %s: UP=%d DOWN=%d TP=%s FP=%s AUC=%s",
        config.name,
        n_up,
        n_down,
        result.truth_metrics.tp,
        result.truth_metrics.fp,
        result.truth_metrics.auc,
    )
    return result
