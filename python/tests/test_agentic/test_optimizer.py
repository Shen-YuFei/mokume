"""Tests for mokume.agentic.optimizer (integration with mock)."""

from unittest.mock import patch

import pandas as pd
import pytest

from mokume.agentic.config import AgenticConfig
from mokume.agentic.optimizer import (
    RoundContext,
    optimize_contrast,
    optimize_from_dataset,
)
from mokume.agentic.profiler import profile_data
from mokume.core.dataset import QpxDataset


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


def test_optimize_from_dataset_delegates(synthetic_protein_df, sample_to_condition):
    """optimize_from_dataset pulls the proteins level and delegates to optimize()."""
    config = AgenticConfig(use_llm=False)
    dataset = QpxDataset()
    dataset.proteins = synthetic_protein_df

    with patch(
        "mokume.agentic.optimizer.optimize", return_value={"A_vs_B": "state"}
    ) as mock_optimize:
        result = optimize_from_dataset(dataset, sample_to_condition, config)

    mock_optimize.assert_called_once()
    passed_df = mock_optimize.call_args.args[0]
    pd.testing.assert_frame_equal(passed_df, synthetic_protein_df)
    assert result == {"A_vs_B": "state"}


def test_optimize_from_dataset_without_proteins_raises(sample_to_condition):
    """A QpxDataset with no proteins level raises a clear error."""
    config = AgenticConfig(use_llm=False)
    with pytest.raises(ValueError, match="proteins"):
        optimize_from_dataset(QpxDataset(), sample_to_condition, config)
