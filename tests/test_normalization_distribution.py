"""Unit tests for distribution-alignment normalizers.

Covers QuantileNormalizer, MedianCenterNormalizer, and MeanCenterNormalizer.
"""

import numpy as np
import pandas as pd
import pytest

from mokume.normalization import (
    MeanCenterNormalizer,
    MedianCenterNormalizer,
    QuantileNormalizer,
    median_center,
    quantile_normalize,
)
from mokume.model.normalization import PeptideNormalizationMethod


@pytest.fixture
def protein_matrix():
    """Wide protein x sample matrix in log2 space."""
    return pd.DataFrame(
        {
            "S1": [10.0, 12.0, 14.0, 16.0],
            "S2": [11.0, 13.0, 15.0, 17.0],
            "S3": [9.0, 11.0, 13.0, 15.0],
        },
        index=["P1", "P2", "P3", "P4"],
    )


class TestQuantileNormalizer:
    def test_fit_transform_aligns_distributions(self, protein_matrix):
        normalized = QuantileNormalizer().fit_transform(protein_matrix)
        # After quantile normalization, the sorted values per sample
        # should be identical across samples.
        sorted_per_sample = {
            col: sorted(normalized[col].dropna()) for col in normalized.columns
        }
        first_sample = sorted_per_sample["S1"]
        for sample, values in sorted_per_sample.items():
            np.testing.assert_allclose(values, first_sample, rtol=1e-9)

    def test_matches_underlying_function(self, protein_matrix):
        from_class = QuantileNormalizer().fit_transform(protein_matrix)
        from_function = quantile_normalize(protein_matrix)
        pd.testing.assert_frame_equal(from_class, from_function)

    def test_transform_without_fit_raises(self, protein_matrix):
        with pytest.raises(ValueError, match="not been fitted"):
            QuantileNormalizer().transform(protein_matrix)


class TestMedianCenterNormalizer:
    def test_each_sample_has_zero_median(self, protein_matrix):
        normalized = MedianCenterNormalizer().fit_transform(protein_matrix)
        # After centering, each sample column should have median 0.
        for col in normalized.columns:
            assert normalized[col].median() == pytest.approx(0.0)

    def test_matches_underlying_function(self, protein_matrix):
        from_class = MedianCenterNormalizer().fit_transform(protein_matrix)
        from_function = median_center(protein_matrix, axis=0)
        pd.testing.assert_frame_equal(from_class, from_function)


class TestMeanCenterNormalizer:
    def test_each_sample_has_zero_mean(self, protein_matrix):
        normalized = MeanCenterNormalizer().fit_transform(protein_matrix)
        for col in normalized.columns:
            assert normalized[col].mean() == pytest.approx(0.0)

    def test_relative_protein_order_preserved(self, protein_matrix):
        normalized = MeanCenterNormalizer().fit_transform(protein_matrix)
        # Centering preserves rank within each sample.
        for col in protein_matrix.columns:
            raw_rank = protein_matrix[col].rank()
            norm_rank = normalized[col].rank()
            pd.testing.assert_series_equal(raw_rank, norm_rank, check_names=False)


class TestPeptideNormalizationMethodEnum:
    @pytest.mark.parametrize(
        "member",
        [
            PeptideNormalizationMethod.Quantile,
            PeptideNormalizationMethod.MedianCenter,
            PeptideNormalizationMethod.MeanCenter,
            PeptideNormalizationMethod.Rlr,
        ],
    )
    def test_new_members_are_dataset_level(self, member):
        assert member.is_dataset_level is True

    def test_globalmedian_is_not_dataset_level(self):
        assert PeptideNormalizationMethod.GlobalMedian.is_dataset_level is False

    def test_from_str_recognizes_new_members(self):
        assert (
            PeptideNormalizationMethod.from_str("Quantile")
            == PeptideNormalizationMethod.Quantile
        )
        assert (
            PeptideNormalizationMethod.from_str("rlr") == PeptideNormalizationMethod.Rlr
        )
