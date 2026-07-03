"""Tests for mokume.agentic.evaluator."""

import pandas as pd

from mokume.agentic.config import ScoreWeights
from mokume.agentic.evaluator import (
    compute_score_ground_truth,
    compute_score_unsupervised,
    evaluate,
)
from mokume.agentic.state import CandidateConfig, EvaluationResult


def test_evaluate_no_ground_truth(synthetic_protein_df, sample_to_condition):
    """Evaluate without ground truth produces universal metrics."""
    cfg = CandidateConfig(name="test", de_method="deqms")
    # Fake DE result
    de_df = pd.DataFrame(
        {
            "protein": [f"P{i:05d}" for i in range(10)],
            "pvalue": [0.001] * 3 + [0.5] * 7,
            "adj_pvalue": [0.01] * 3 + [0.8] * 7,
            "log2fc": [2.0] * 3 + [0.1] * 7,
            "significance": ["UP"] * 3 + ["Not Sig"] * 7,
        }
    )
    result = evaluate(cfg, de_df, synthetic_protein_df, sample_to_condition)
    assert result.n_de_up == 3
    assert result.tp is None


def test_evaluate_with_ground_truth(
    synthetic_protein_df,
    sample_to_condition,
    ground_truth_proteins,
):
    """Evaluate with ground truth computes TP/FP/AUC."""
    cfg = CandidateConfig(name="test_gt", de_method="deqms")
    de_df = pd.DataFrame(
        {
            "protein": ["P00000", "P00001", "P00002", "P00010", "P00020"],
            "pvalue": [0.001, 0.002, 0.003, 0.01, 0.9],
            "adj_pvalue": [0.01, 0.02, 0.03, 0.05, 0.99],
            "log2fc": [3.0, 2.5, 2.0, 1.5, 0.1],
            "significance": ["UP", "UP", "UP", "UP", "Not Sig"],
        }
    )
    result = evaluate(
        cfg,
        de_df,
        synthetic_protein_df,
        sample_to_condition,
        ground_truth_proteins,
    )
    assert result.tp == 3  # P00000, P00001, P00002
    assert result.fp == 1  # P00010


def test_score_ground_truth():
    """Score formula for Mode A."""
    r = EvaluationResult(
        config_name="t",
        config={},
        tp=47,
        fp=1,
        auc=0.996,
    )
    score = compute_score_ground_truth(
        r, max_tp=48, total_tested=1300, weights=ScoreWeights()
    )
    assert score > 0.7


def test_score_unsupervised():
    """Score formula for Mode B."""
    r = EvaluationResult(
        config_name="t",
        config={},
        n_de_up=200,
        n_de_down=50,
        median_cv=0.15,
        missing_rate=0.10,
    )
    score = compute_score_unsupervised(r, max_de=300, weights=ScoreWeights())
    assert score > 0
