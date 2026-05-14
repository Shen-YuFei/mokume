"""Tests for mokume.agentic.optimizer (integration with mock)."""

from unittest.mock import patch

import pandas as pd

from mokume.agentic.config import AgenticConfig
from mokume.agentic.optimizer import RoundContext, optimize_contrast
from mokume.agentic.profiler import profile_data


def _mock_de_result(n=30):
    """Create a fake DE result DataFrame."""
    return pd.DataFrame(
        {
            "protein": [f"P{i:05d}" for i in range(n)],
            "pvalue": [0.001] * 5 + [0.5] * (n - 5),
            "adj_pvalue": [0.01] * 5 + [0.8] * (n - 5),
            "log2fc": [3.0] * 5 + [0.1] * (n - 5),
            "significance": ["UP"] * 5 + ["Not Sig"] * (n - 5),
        }
    )


def test_optimize_contrast_no_llm(synthetic_protein_df, sample_to_condition, tmp_path):
    """Full optimization loop with mocked runner, no LLM."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    config = AgenticConfig(
        use_llm=False,
        max_rounds=2,
        max_experiments=4,
        output_dir=str(tmp_path / "output"),
    )

    ctx = RoundContext(
        protein_df=synthetic_protein_df,
        sample_to_condition=sample_to_condition,
        contrast=("A", "B"),
        ground_truth={f"P{i:05d}" for i in range(5)},
        peptide_counts=None,
        config=config,
    )
    with patch(
        "mokume.agentic.optimizer.run_experiment", return_value=_mock_de_result()
    ):
        state = optimize_contrast(ctx, profile)

    assert state.total_experiments > 0
    assert state.best_config is not None
    assert len(state.rounds) >= 1
    assert len(state.audit_trail) > 0


def test_optimize_budget_respected(synthetic_protein_df, sample_to_condition, tmp_path):
    """Optimizer stops when budget is exhausted."""
    profile = profile_data(synthetic_protein_df, sample_to_condition)
    config = AgenticConfig(
        use_llm=False,
        max_rounds=10,
        max_experiments=3,
        output_dir=str(tmp_path / "budget"),
    )

    ctx = RoundContext(
        protein_df=synthetic_protein_df,
        sample_to_condition=sample_to_condition,
        contrast=("A", "B"),
        ground_truth=None,
        peptide_counts=None,
        config=config,
    )
    with patch(
        "mokume.agentic.optimizer.run_experiment", return_value=_mock_de_result()
    ):
        state = optimize_contrast(ctx, profile)

    assert state.total_experiments <= 3
