"""
Tests for hierarchical sample normalization.

Tests the HierarchicalSampleNormalizer and HierarchicalIonAligner classes
which provide DirectLFQ-style normalization natively in mokume.
"""

import numpy as np
import pandas as pd
import pytest

from mokume.normalization.hierarchical import (
    HierarchicalSampleNormalizer,
    HierarchicalIonAligner,
    DistanceMetric,
)


class TestDistanceMetric:
    """Tests for DistanceMetric enum."""

    def test_from_str_median(self):
        assert DistanceMetric.from_str("median") == DistanceMetric.MEDIAN
        assert DistanceMetric.from_str("MEDIAN") == DistanceMetric.MEDIAN

    def test_from_str_variance(self):
        assert DistanceMetric.from_str("variance") == DistanceMetric.VARIANCE

    def test_from_str_overlap(self):
        assert DistanceMetric.from_str("overlap") == DistanceMetric.OVERLAP

    def test_from_str_none(self):
        assert DistanceMetric.from_str(None) == DistanceMetric.MEDIAN

    def test_from_str_invalid(self):
        with pytest.raises(ValueError, match="Unknown distance metric"):
            DistanceMetric.from_str("invalid")


class TestHierarchicalSampleNormalizer:
    """Tests for HierarchicalSampleNormalizer."""

    @pytest.fixture
    def simple_df(self):
        """Simple test DataFrame with systematic shifts."""
        # Use deterministic base values (not random) for predictable results
        n_features = 100
        base = np.linspace(15, 25, n_features)  # Base intensities from 15 to 25

        return pd.DataFrame({
            "sample1": base,
            "sample2": base + 0.5,  # Shifted up by 0.5
            "sample3": base - 0.3,  # Shifted down by 0.3
            "sample4": base + 1.0,  # Shifted up by 1.0
        })

    @pytest.fixture
    def df_with_nans(self):
        """DataFrame with missing values."""
        np.random.seed(42)
        n_features = 100
        base = np.random.randn(n_features) * 2 + 20

        df = pd.DataFrame({
            "sample1": base,
            "sample2": base + 0.5,
            "sample3": base - 0.3,
        })

        # Add some NaNs
        df.iloc[10:20, 0] = np.nan  # sample1
        df.iloc[30:40, 1] = np.nan  # sample2
        df.iloc[50:60, 2] = np.nan  # sample3

        return df

    @pytest.fixture
    def multiindex_df(self):
        """DataFrame with MultiIndex [protein, peptide]."""
        np.random.seed(42)
        data = []

        for protein in ["P1", "P2", "P3"]:
            for peptide in ["pep1", "pep2", "pep3"]:
                base = np.random.randn() * 2 + 20
                data.append({
                    "protein": protein,
                    "peptide": peptide,
                    "sample1": base,
                    "sample2": base + 0.5,
                    "sample3": base - 0.3,
                })

        df = pd.DataFrame(data)
        df = df.set_index(["protein", "peptide"])
        return df

    def test_fit_transform_aligns_samples(self, simple_df):
        """Test that normalization produces valid shift factors."""
        # Force quadratic optimization for better alignment
        normalizer = HierarchicalSampleNormalizer(num_samples_quadratic=100)
        _ = normalizer.fit_transform(simple_df)

        # Check that normalization factors are computed
        assert normalizer.normalization_factors_ is not None
        assert len(normalizer.normalization_factors_) == 4

        # Check that normalization factors approximately match expected shifts
        # sample1: base (shift = 0)
        # sample2: base + 0.5 (needs shift ≈ -0.5)
        # sample3: base - 0.3 (needs shift ≈ +0.3)
        # sample4: base + 1.0 (needs shift ≈ -1.0)
        factors = normalizer.normalization_factors_

        # One sample should be the reference (shift = 0 or close to it)
        # The relative differences between shifts should match the original offsets
        shifts = list(factors.values())
        shift_range = max(shifts) - min(shifts)

        # The shift range should be approximately equal to the original offset range
        # Original: max offset = 1.0, min offset = -0.3, range = 1.3
        assert 1.0 < shift_range < 1.6  # Allow some tolerance

    def test_fit_transform_preserves_shape(self, simple_df):
        """Test that normalization preserves DataFrame shape."""
        normalizer = HierarchicalSampleNormalizer()
        normalized = normalizer.fit_transform(simple_df)

        assert normalized.shape == simple_df.shape
        assert list(normalized.columns) == list(simple_df.columns)

    def test_fit_transform_with_nans(self, df_with_nans):
        """Test that normalization handles NaN values."""
        normalizer = HierarchicalSampleNormalizer(min_overlap=5)
        normalized = normalizer.fit_transform(df_with_nans)

        # NaNs should be preserved
        assert normalized.isna().sum().sum() == df_with_nans.isna().sum().sum()

    def test_fit_then_transform(self, simple_df):
        """Test separate fit and transform calls."""
        normalizer = HierarchicalSampleNormalizer()
        normalizer.fit(simple_df)

        assert normalizer.normalization_factors_ is not None
        assert len(normalizer.normalization_factors_) == len(simple_df.columns)

        normalized = normalizer.transform(simple_df)
        assert normalized.shape == simple_df.shape

    def test_transform_without_fit_raises(self, simple_df):
        """Test that transform without fit raises error."""
        normalizer = HierarchicalSampleNormalizer()

        with pytest.raises(ValueError, match="Must call fit"):
            normalizer.transform(simple_df)

    def test_selected_proteins(self, multiindex_df):
        """Test normalization on selected proteins only."""
        normalizer = HierarchicalSampleNormalizer(selected_proteins=["P1", "P2"])
        normalizer.fit(multiindex_df)

        # Should still produce factors for all samples
        assert len(normalizer.normalization_factors_) == 3

    def test_linear_vs_quadratic(self, simple_df):
        """Test that linear and quadratic both produce valid shift factors."""
        normalizer_linear = HierarchicalSampleNormalizer(num_samples_quadratic=1)
        normalizer_quad = HierarchicalSampleNormalizer(num_samples_quadratic=100)

        normalizer_linear.fit(simple_df)
        normalizer_quad.fit(simple_df)

        # Both should produce normalization factors
        assert normalizer_linear.normalization_factors_ is not None
        assert normalizer_quad.normalization_factors_ is not None

        # Both should have factors for all samples
        assert len(normalizer_linear.normalization_factors_) == 4
        assert len(normalizer_quad.normalization_factors_) == 4

        # The shift ranges should be similar (within a factor of 2)
        linear_range = (
            max(normalizer_linear.normalization_factors_.values()) -
            min(normalizer_linear.normalization_factors_.values())
        )
        quad_range = (
            max(normalizer_quad.normalization_factors_.values()) -
            min(normalizer_quad.normalization_factors_.values())
        )

        # Both should capture roughly the same shift range (1.0 - (-0.3) = 1.3)
        assert 0.5 < linear_range < 2.5
        assert 0.5 < quad_range < 2.5

    def test_distance_metric_variance(self, simple_df):
        """Test with variance distance metric."""
        normalizer = HierarchicalSampleNormalizer(
            distance_metric=DistanceMetric.VARIANCE
        )
        normalized = normalizer.fit_transform(simple_df)

        assert normalized.shape == simple_df.shape

    def test_distance_metric_overlap(self, simple_df):
        """Test with overlap distance metric."""
        normalizer = HierarchicalSampleNormalizer(
            distance_metric=DistanceMetric.OVERLAP
        )
        normalized = normalizer.fit_transform(simple_df)

        assert normalized.shape == simple_df.shape

    def test_distance_metric_string(self, simple_df):
        """Test distance metric specified as string."""
        normalizer = HierarchicalSampleNormalizer(distance_metric="median")
        normalized = normalizer.fit_transform(simple_df)

        assert normalized.shape == simple_df.shape

    def test_single_sample_returns_unchanged(self):
        """Test that single sample returns unchanged."""
        df = pd.DataFrame({"sample1": [10.0, 11.0, 12.0]})

        normalizer = HierarchicalSampleNormalizer()
        normalized = normalizer.fit_transform(df)

        # Values should be unchanged (shift is 0)
        assert normalizer.normalization_factors_["sample1"] == 0.0
        np.testing.assert_array_almost_equal(
            normalized["sample1"].values, df["sample1"].values
        )

    def test_two_samples(self):
        """Test normalization with just two samples."""
        # Use enough data points to exceed min_overlap
        n = 20
        df = pd.DataFrame({
            "sample1": np.linspace(10.0, 15.0, n),
            "sample2": np.linspace(10.5, 15.5, n),  # Shifted by 0.5
        })

        normalizer = HierarchicalSampleNormalizer(min_overlap=5)
        normalized = normalizer.fit_transform(df)

        # With two samples, the shift should align them
        # The medians should be equal (or very close) after normalization
        median_diff_after = abs(
            normalized["sample1"].median() - normalized["sample2"].median()
        )
        assert median_diff_after < 0.01  # Should be essentially zero

    def test_empty_after_protein_filter_raises(self, multiindex_df):
        """Test that empty DataFrame after filtering raises error."""
        normalizer = HierarchicalSampleNormalizer(
            selected_proteins=["NONEXISTENT"]
        )

        with pytest.raises(ValueError, match="No features remaining"):
            normalizer.fit(multiindex_df)


class TestHierarchicalIonAligner:
    """Tests for HierarchicalIonAligner."""

    @pytest.fixture
    def protein_df(self):
        """DataFrame representing ions for a single protein."""
        np.random.seed(42)
        return pd.DataFrame({
            "sample1": [10.0, 10.5, 11.0, 11.5],
            "sample2": [10.2, 10.7, np.nan, 11.7],
            "sample3": [9.8, 10.3, 10.8, 11.3],
        }, index=["ion1", "ion2", "ion3", "ion4"])

    def test_align_protein_ions(self, protein_df):
        """Test basic ion alignment."""
        aligner = HierarchicalIonAligner()
        aligned = aligner.align_protein_ions(protein_df)

        assert aligned.shape == protein_df.shape
        assert list(aligned.index) == list(protein_df.index)
        assert list(aligned.columns) == list(protein_df.columns)

    def test_single_ion_unchanged(self):
        """Test that single ion returns unchanged."""
        df = pd.DataFrame({
            "sample1": [10.0],
            "sample2": [10.5],
        }, index=["ion1"])

        aligner = HierarchicalIonAligner()
        aligned = aligner.align_protein_ions(df)

        pd.testing.assert_frame_equal(aligned, df)

    def test_align_all_proteins(self):
        """Test aligning all proteins in a dataset."""
        np.random.seed(42)

        # Create MultiIndex DataFrame
        index = pd.MultiIndex.from_tuples([
            ("P1", "ion1"), ("P1", "ion2"), ("P1", "ion3"),
            ("P2", "ion1"), ("P2", "ion2"),
        ], names=["protein", "ion"])

        df = pd.DataFrame({
            "sample1": [10.0, 10.5, 11.0, 20.0, 20.5],
            "sample2": [10.2, 10.7, 11.2, 20.2, 20.7],
            "sample3": [9.8, 10.3, 10.8, 19.8, 20.3],
        }, index=index)

        aligner = HierarchicalIonAligner()
        aligned = aligner.align_all_proteins(df)

        assert aligned.shape == df.shape
        assert list(aligned.index) == list(df.index)


class TestIntegration:
    """Integration tests combining sample and ion normalization."""

    def test_sample_then_ion_normalization(self):
        """Test pipeline: sample normalization then ion alignment."""
        np.random.seed(42)

        # Create test data
        index = pd.MultiIndex.from_tuples([
            ("P1", "ion1"), ("P1", "ion2"),
            ("P2", "ion1"), ("P2", "ion2"),
        ], names=["protein", "ion"])

        base = np.array([10.0, 10.5, 20.0, 20.5])
        df = pd.DataFrame({
            "sample1": base,
            "sample2": base + 0.5,
            "sample3": base - 0.3,
        }, index=index)

        # Step 1: Sample normalization
        sample_normalizer = HierarchicalSampleNormalizer()
        sample_normalized = sample_normalizer.fit_transform(df)

        # Step 2: Ion alignment per protein
        ion_aligner = HierarchicalIonAligner()
        final = ion_aligner.align_all_proteins(sample_normalized)

        assert final.shape == df.shape

    def test_reproducibility(self):
        """Test that normalization is deterministic."""
        np.random.seed(42)
        df = pd.DataFrame({
            "s1": np.random.randn(50) + 10,
            "s2": np.random.randn(50) + 10.5,
            "s3": np.random.randn(50) + 9.5,
        })

        normalizer1 = HierarchicalSampleNormalizer()
        normalizer2 = HierarchicalSampleNormalizer()

        result1 = normalizer1.fit_transform(df)
        result2 = normalizer2.fit_transform(df)

        pd.testing.assert_frame_equal(result1, result2)
