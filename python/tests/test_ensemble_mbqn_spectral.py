"""Unit tests for ensemble DE and count-based quantification."""

import numpy as np
import pandas as pd
import pytest

from mokume.analysis.ensemble import combine_de_results
from mokume.quantification.spectral_count import (
    PeptideCountQuantification,
    SpectralCountQuantification,
)


# ---------- Fixtures ----------


@pytest.fixture
def mock_de_results():
    """Two mock DE result DataFrames with overlapping proteins."""
    df_a = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3", "P4"],
            "log2FC": [2.0, -1.5, 0.3, 1.8],
            "pvalue": [0.001, 0.01, 0.5, 0.002],
            "adj_pvalue": [0.004, 0.02, 0.6, 0.006],
            "significance": ["UP", "DOWN", "Unchanged", "UP"],
        }
    )
    df_b = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3", "P5"],
            "log2FC": [1.8, -1.2, 0.1, 2.5],
            "pvalue": [0.002, 0.02, 0.4, 0.001],
            "adj_pvalue": [0.006, 0.03, 0.5, 0.003],
            "significance": ["UP", "DOWN", "Unchanged", "UP"],
        }
    )
    return {"method_a": df_a, "method_b": df_b}


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
        p1 = result[result["ProteinName"] == "P1"].iloc[0]
        assert p1["n_methods_up"] == 2
        assert p1["significance"] == "UP"

    def test_consensus_down(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p2 = result[result["ProteinName"] == "P2"].iloc[0]
        assert p2["n_methods_down"] == 2

    def test_no_consensus_unchanged(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p3 = result[result["ProteinName"] == "P3"].iloc[0]
        assert p3["significance"] == "Unchanged"

    def test_single_method_below_k(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=2)
        p5 = result[result["ProteinName"] == "P5"].iloc[0]
        assert p5["n_methods_up"] == 1

    def test_empty_input(self):
        result = combine_de_results({})
        assert result.empty

    def test_median_fc(self, mock_de_results):
        result = combine_de_results(mock_de_results, min_k=1)
        p1 = result[result["ProteinName"] == "P1"].iloc[0]
        assert p1["log2FC"] == pytest.approx(np.median([2.0, 1.8]))


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


class TestPeptideCount:
    def test_counts_distinct_peptides(self, peptide_df):
        result = PeptideCountQuantification().quantify(
            peptide_df,
            protein_column="Protein",
            peptide_column="Peptide",
            intensity_column="NormIntensity",
            sample_column="Sample",
        )
        a_s1 = result[(result["Protein"] == "A") & (result["Sample"] == "S1")][
            "Intensity"
        ].iloc[0]
        assert a_s1 == 2

    def test_name(self):
        assert PeptideCountQuantification().name == "PeptideCount"
