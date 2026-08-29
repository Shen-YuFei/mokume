"""Regression tests for matrix-level sample correlation filtering."""

import pandas as pd
import pytest

from mokume.pipeline.sample_qc import filter_samples_by_correlation


def test_sample_correlation_filter_is_one_shot_and_preserves_values() -> None:
    proteins = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3", "P4"],
            "A1": [1.0, 2.0, 4.0, 8.0],
            "A2": [2.0, 4.0, 8.0, 16.0],
            "A3": [8.0, 4.0, 2.0, 1.0],
            "B1": [3.0, 6.0, 12.0, 24.0],
            "B2": [6.0, 12.0, 24.0, 48.0],
        }
    )
    conditions = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B"}

    filtered = filter_samples_by_correlation(proteins, conditions, threshold=-0.5)

    assert list(filtered.columns) == ["ProteinName", "A1", "A2", "B1", "B2"]
    pd.testing.assert_frame_equal(filtered, proteins.drop(columns="A3"))


def test_sample_correlation_filter_rejects_singleton_conditions() -> None:
    proteins = pd.DataFrame({"ProteinName": ["P1", "P2", "P3"], "A1": [1.0, 2.0, 3.0]})

    with pytest.raises(ValueError, match="at least two samples"):
        filter_samples_by_correlation(proteins, {"A1": "A"}, threshold=0.8)


def test_sample_correlation_filter_rejects_undefined_pairs() -> None:
    proteins = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3"],
            "A1": [1.0, None, None],
            "A2": [2.0, None, None],
        }
    )

    with pytest.raises(ValueError, match="1 pairwise-complete usable proteins"):
        filter_samples_by_correlation(proteins, {"A1": "A", "A2": "A"}, threshold=0.8)


def test_sample_correlation_filter_accepts_negative_log2_values() -> None:
    proteins = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3"],
            "A1": [-2.0, 0.0, 2.0],
            "A2": [-1.0, 1.0, 3.0],
        }
    )

    filtered = filter_samples_by_correlation(
        proteins,
        {"A1": "A", "A2": "A"},
        threshold=0.99,
        values_are_log2=True,
    )

    pd.testing.assert_frame_equal(filtered, proteins)
