"""Unit tests for ensemble DE, MBQN normalization, and SpectralCount quantification."""

import numpy as np
import pandas as pd
import pytest

from mokume.analysis.ensemble import combine_de_results
from mokume.normalization.mbqn import mbqn_normalize, MBQNNormalizer
from mokume.quantification.spectral_count import SpectralCountQuantification


# ---------- Fixtures ----------


@pytest.fixture
def mock_de_results():
    """Two mock DE result DataFrames with overlapping proteins."""
    df_a = pd.DataFrame(
        {
            "Protein": ["P1", "P2", "P3", "P4"],
            "log2FC": [2.0, -1.5, 0.3, 1.8],
            "pvalue": [0.001, 0.01, 0.5, 0.002],
            "adj_pvalue": [0.004, 0.02, 0.6, 0.006],
            "significance": ["UP", "DOWN", "Unchanged", "UP"],
        }
    )
    df_b = pd.DataFrame(
        {
            "Protein": ["P1", "P2", "P3", "P5"],
            "log2FC": [1.8, -1.2, 0.1, 2.5],
            "pvalue": [0.002, 0.02, 0.4, 0.001],
            "adj_pvalue": [0.006, 0.03, 0.5, 0.003],
            "significance": ["UP", "DOWN", "Unchanged", "UP"],
        }
    )
    return {"method_a": df_a, "method_b": df_b}


@pytest.fixture
def wide_matrix():
    """5 proteins x 4 samples log2 matrix."""
    rng = np.random.default_rng(42)
    return pd.DataFrame(
        rng.normal(20, 2, size=(5, 4)),
        index=[f"P{i}" for i in range(5)],
        columns=["S1", "S2", "S3", "S4"],
    )


@pytest.fixture
def peptide_df():
    """Simple peptide-level DataFrame for SpectralCount."""
    return pd.DataFrame(
        {
            "Protein": ["A", "A", "A", "B", "B", "A", "B"],
            "Peptide": ["p1", "p2", "p1", "p3", "p3", "p1", "p4"],
            "NormIntensity": [100, 200, 150, 300, 250, 180, 400],
            "Sample": ["S1", "S1", "S1", "S1", "S1", "S2", "S2"],
        }
    )


# ---------- Ensemble ----------


class TestCombineDeResults:
    def test_consensus_up(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p1 = result[result["Protein"] == "P1"].iloc[0]
        assert p1["n_methods_up"] == 2
        assert p1["significance"] == "UP"

    def test_consensus_down(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p2 = result[result["Protein"] == "P2"].iloc[0]
        assert p2["n_methods_down"] == 2

    def test_no_consensus_unchanged(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p3 = result[result["Protein"] == "P3"].iloc[0]
        assert p3["significance"] == "Unchanged"

    def test_single_method_below_k(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p5 = result[result["Protein"] == "P5"].iloc[0]
        assert p5["n_methods_up"] == 1

    def test_empty_input(self):
        result = combine_de_results({})
        assert result.empty

    def test_median_fc(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=1)
        p1 = result[result["Protein"] == "P1"].iloc[0]
        assert p1["log2FC"] == pytest.approx(np.median([2.0, 1.8]))


# ---------- MBQN ----------


class TestMBQN:
    def test_preserves_row_means(self, wide_matrix):
        row_means_before = wide_matrix.mean(axis=1)
        result = mbqn_normalize(wide_matrix)
        row_means_after = result.mean(axis=1)
        np.testing.assert_array_almost_equal(
            row_means_after.values, row_means_before.values, decimal=10
        )

    def test_shape_preserved(self, wide_matrix):
        result = mbqn_normalize(wide_matrix)
        assert result.shape == wide_matrix.shape

    def test_normalizer_class(self, wide_matrix):
        norm = MBQNNormalizer()
        result = norm.fit_transform(wide_matrix)
        assert result.shape == wide_matrix.shape


# ---------- SpectralCount ----------


class TestSpectralCount:
    def test_counts_spectra(self, peptide_df):
        sc = SpectralCountQuantification()
        result = sc.quantify(
            peptide_df,
            protein_column="Protein",
            peptide_column="Peptide",
            intensity_column="NormIntensity",
            sample_column="Sample",
        )
        a_s1 = result[(result["Protein"] == "A") & (result["Sample"] == "S1")][
            "Intensity"
        ].iloc[0]
        assert a_s1 == 3  # 3 PSMs for protein A in S1

    def test_name(self):
        assert SpectralCountQuantification().name == "SpectralCount"
