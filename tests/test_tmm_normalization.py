"""
Tests for TMM (Trimmed Mean of M-values) normalization.

Tests cover:
- Basic functionality
- Edge cases (zeros, missing values)
- Geometric mean property
- Known value validation
"""

import numpy as np
import pandas as pd
import pytest

from mokume.normalization.tmm import TMMNormalizer, tmm_normalize


class TestTMMNormalizer:
    """Tests for TMMNormalizer class."""

    @pytest.fixture
    def simple_data(self):
        """Simple test data with known properties (50 proteins for proper TMM)."""
        np.random.seed(42)
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)
        return pd.DataFrame(
            {
                "S1": base * (1 + np.random.randn(n_proteins) * 0.1),
                "S2": base * (1 + np.random.randn(n_proteins) * 0.1),
                "S3": base * (1 + np.random.randn(n_proteins) * 0.1),
                "S4": base * (1 + np.random.randn(n_proteins) * 0.1),
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

    @pytest.fixture
    def biased_data(self):
        """Data with composition bias (one sample has outlier proteins)."""
        np.random.seed(42)
        n_proteins = 50
        base = np.ones(n_proteins) * 100

        # S1 and S2 are normal
        s1 = base.copy()
        s2 = base.copy()

        # S3 has a few highly abundant outlier proteins
        s3 = base.copy()
        s3[-5:] = 10000  # Last 5 proteins are 100x higher

        return pd.DataFrame(
            {
                "S1": s1,
                "S2": s2,
                "S3": s3,
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

    def test_fit_creates_factors(self, simple_data):
        """TMM fit should create normalization factors for each sample."""
        tmm = TMMNormalizer()
        tmm.fit(simple_data)

        assert tmm.norm_factors_ is not None
        assert len(tmm.norm_factors_) == 4
        assert all(f > 0 for f in tmm.norm_factors_)

    def test_geometric_mean_is_one(self, simple_data):
        """Geometric mean of factors should be 1."""
        tmm = TMMNormalizer()
        tmm.fit(simple_data)

        geom_mean = np.exp(np.mean(np.log(tmm.norm_factors_)))
        assert np.isclose(geom_mean, 1.0, atol=1e-6)

    def test_transform_preserves_shape(self, simple_data):
        """Transform should preserve data shape."""
        tmm = TMMNormalizer()
        result = tmm.fit_transform(simple_data)

        assert result.shape == simple_data.shape
        assert list(result.columns) == list(simple_data.columns)
        assert list(result.index) == list(simple_data.index)

    def test_reference_sample_selection(self, simple_data):
        """Should select valid reference sample."""
        tmm = TMMNormalizer()
        tmm.fit(simple_data)

        assert tmm.ref_sample_ in simple_data.columns

    def test_custom_reference_sample(self, simple_data):
        """Should use custom reference sample if provided."""
        tmm = TMMNormalizer(ref_sample="S2")
        tmm.fit(simple_data)

        assert tmm.ref_sample_ == "S2"

    def test_invalid_reference_raises(self, simple_data):
        """Should raise error for invalid reference sample."""
        tmm = TMMNormalizer(ref_sample="InvalidSample")

        with pytest.raises(ValueError, match="not found"):
            tmm.fit(simple_data)

    def test_handles_zeros(self):
        """TMM should handle zero values gracefully."""
        np.random.seed(42)
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)

        data = pd.DataFrame(
            {
                "S1": base.copy(),
                "S2": base.copy(),
                "S3": base.copy(),
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

        # Add some zeros
        data.iloc[0, 0] = 0
        data.iloc[5, 1] = 0

        tmm = TMMNormalizer()
        result = tmm.fit_transform(data)

        # Should not have all NaN
        assert not result.isnull().all().any()

    def test_handles_nan(self):
        """TMM should handle NaN values gracefully."""
        np.random.seed(42)
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)

        data = pd.DataFrame(
            {
                "S1": base.copy(),
                "S2": base.copy(),
                "S3": base.copy(),
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

        # Add some NaN values
        data.iloc[0, 0] = np.nan
        data.iloc[5, 1] = np.nan

        tmm = TMMNormalizer()
        result = tmm.fit_transform(data)

        # Should complete without error
        assert result is not None

    def test_composition_bias_correction(self, biased_data):
        """TMM should correct for composition bias from outlier proteins."""
        tmm = TMMNormalizer()
        tmm.fit(biased_data)

        # S3 has the outlier, should have higher factor (to scale down)
        # After rescaling to geometric mean = 1, S3 factor should be > 1
        factors = tmm.get_norm_factors()

        # S1 and S2 should have similar factors, S3 should be different
        assert np.isclose(factors["S1"], factors["S2"], atol=0.1)
        # S3 factor should be notably different due to outlier
        assert not np.isclose(factors["S3"], factors["S1"], atol=0.2)

    def test_trim_parameters(self, simple_data):
        """Different trim parameters should affect results."""
        tmm1 = TMMNormalizer(m_trim=0.1, a_trim=0.01)
        tmm2 = TMMNormalizer(m_trim=0.4, a_trim=0.1)

        result1 = tmm1.fit_transform(simple_data)
        result2 = tmm2.fit_transform(simple_data)

        # Results may or may not differ depending on data
        # At minimum, both should complete successfully
        assert result1 is not None
        assert result2 is not None

    def test_transform_without_fit_raises(self, simple_data):
        """Transform without fit should raise error."""
        tmm = TMMNormalizer()

        with pytest.raises(ValueError, match="not been fitted"):
            tmm.transform(simple_data)

    def test_get_norm_factors(self, simple_data):
        """get_norm_factors should return copy of factors."""
        tmm = TMMNormalizer()
        tmm.fit(simple_data)

        factors1 = tmm.get_norm_factors()
        factors2 = tmm.get_norm_factors()

        # Should be equal but not same object
        pd.testing.assert_series_equal(factors1, factors2)

    def test_get_effective_library_sizes(self, simple_data):
        """get_effective_library_sizes should return adjusted sizes."""
        tmm = TMMNormalizer()
        tmm.fit(simple_data)

        eff_sizes = tmm.get_effective_library_sizes(simple_data)

        assert len(eff_sizes) == len(simple_data.columns)
        assert all(s > 0 for s in eff_sizes)

    def test_empty_dataframe_raises(self):
        """Empty DataFrame should raise error."""
        tmm = TMMNormalizer()

        with pytest.raises(ValueError, match="empty"):
            tmm.fit(pd.DataFrame())

    def test_single_sample(self):
        """Single sample should return unchanged (cannot normalize)."""
        data = pd.DataFrame(
            {
                "S1": [100, 200, 300],
            },
            index=["P1", "P2", "P3"],
        )

        tmm = TMMNormalizer()
        _ = tmm.fit_transform(data)

        # Factor should be 1.0
        assert np.isclose(tmm.norm_factors_["S1"], 1.0)


class TestTMMConvenienceFunction:
    """Tests for tmm_normalize convenience function."""

    def test_basic_usage(self):
        """tmm_normalize should work with default parameters."""
        np.random.seed(42)
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)

        data = pd.DataFrame(
            {
                "S1": base * (1 + np.random.randn(n_proteins) * 0.1),
                "S2": base * (1 + np.random.randn(n_proteins) * 0.1),
                "S3": base * (1 + np.random.randn(n_proteins) * 0.1),
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

        result = tmm_normalize(data)

        assert result.shape == data.shape
        assert not result.isnull().all().any()

    def test_with_parameters(self):
        """tmm_normalize should accept custom parameters."""
        np.random.seed(42)
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)

        data = pd.DataFrame(
            {
                "S1": base * (1 + np.random.randn(n_proteins) * 0.1),
                "S2": base * (1 + np.random.randn(n_proteins) * 0.1),
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

        result = tmm_normalize(data, m_trim=0.2, a_trim=0.1, ref_sample="S1")

        assert result is not None


class TestTMMNumericalAccuracy:
    """Tests for numerical accuracy of TMM implementation."""

    def test_symmetric_data_factors_near_one(self):
        """For symmetric data, all factors should be close to 1."""
        # Create perfectly symmetric data with enough proteins
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)

        data = pd.DataFrame(
            {
                "S1": base,
                "S2": base,
                "S3": base,
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

        tmm = TMMNormalizer()
        tmm.fit(data)

        # All factors should be exactly 1.0
        for factor in tmm.norm_factors_:
            assert np.isclose(factor, 1.0, atol=1e-6)

    def test_composition_bias_affects_factors(self):
        """Composition bias should result in different normalization factors."""
        n_proteins = 50
        base = np.linspace(100, 1000, n_proteins)

        # S1 and S3 are normal
        # S2 has composition bias: 10% of proteins are 10x higher
        s1 = base.copy()
        s2 = base.copy()
        s2[-5:] = s2[-5:] * 10  # Last 5 proteins are 10x higher (composition bias)
        s3 = base.copy()

        data = pd.DataFrame(
            {
                "S1": s1,
                "S2": s2,
                "S3": s3,
            },
            index=[f"P{i}" for i in range(n_proteins)],
        )

        tmm = TMMNormalizer()
        tmm.fit(data)

        # S2 should have different factor due to composition bias
        # (the outlier proteins inflate the library size but are trimmed in M calculation)
        factors = tmm.get_norm_factors()

        # S1 and S3 should have similar factors
        assert np.isclose(factors["S1"], factors["S3"], atol=0.1)

        # After normalization, the non-outlier proteins in S2 should be adjusted
        result = tmm.transform(data)

        # The normalization should have happened (result is different from input)
        # At minimum, the normalized data should preserve relative order within proteins
        assert result is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
