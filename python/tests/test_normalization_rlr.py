"""Unit tests for RLR (Robust Linear Regression) normalization."""

import numpy as np
import pandas as pd
import pytest

from mokume.normalization import RlrNormalizer, rlr_normalize


@pytest.fixture
def shifted_protein_matrix():
    """Synthetic protein x sample matrix where each sample is a known
    affine transform of a base profile in log2 space.

    Sample S0 is the reference (identity); S1 is shifted by +2 (log2);
    S2 is scaled by 1.5x slope. After RLR normalization the three
    samples should converge toward the same per-protein value.
    """
    rng = np.random.default_rng(42)
    base_log2 = pd.Series(rng.normal(loc=15, scale=2, size=20))
    base = 2**base_log2
    return pd.DataFrame(
        {
            "S0": base,
            "S1": base * 4.0,  # log2 shift +2
            "S2": 2 ** (base_log2 * 1.5),  # log2 slope 1.5
            "S3": base * 8.0,  # log2 shift +3
        }
    )


class TestRlrNormalizer:
    def test_fit_returns_self(self, shifted_protein_matrix):
        normalizer = RlrNormalizer()
        result = normalizer.fit(shifted_protein_matrix)
        assert result is normalizer

    def test_fit_records_per_sample_intercept_and_slope(self, shifted_protein_matrix):
        normalizer = RlrNormalizer().fit(shifted_protein_matrix)
        assert set(normalizer.intercepts_.index) == {"S0", "S1", "S2", "S3"}
        assert set(normalizer.slopes_.index) == {"S0", "S1", "S2", "S3"}
        assert normalizer.reference_ is not None

    def test_transform_aligns_shifted_samples(self, shifted_protein_matrix):
        normalized = RlrNormalizer().fit_transform(shifted_protein_matrix)
        log2_normalized = np.log2(normalized.replace(0, np.nan))
        # After RLR, residual standard deviation across samples should be
        # smaller than the raw input's per-protein std.
        raw_log2 = np.log2(shifted_protein_matrix.replace(0, np.nan))
        raw_std = raw_log2.std(axis=1).mean()
        norm_std = log2_normalized.std(axis=1).mean()
        assert norm_std < raw_std

    def test_transform_without_fit_raises(self, shifted_protein_matrix):
        with pytest.raises(ValueError, match="not been fitted"):
            RlrNormalizer().transform(shifted_protein_matrix)

    def test_empty_input_raises(self):
        with pytest.raises(ValueError, match="empty"):
            RlrNormalizer().fit(pd.DataFrame())

    def test_log_transform_pass_through(self, shifted_protein_matrix):
        log_input = np.log2(shifted_protein_matrix.replace(0, np.nan))
        normalized = RlrNormalizer(log_transform=True).fit_transform(log_input)
        # log_transform=True returns log space output, so values should
        # be roughly in log2 range (5-25), not raw intensity range.
        assert (normalized.dropna(how="all") < 30).all().all()

    def test_rlr_normalize_helper_matches_class(self, shifted_protein_matrix):
        result_helper = rlr_normalize(shifted_protein_matrix)
        result_class = RlrNormalizer().fit_transform(shifted_protein_matrix)
        pd.testing.assert_frame_equal(result_helper, result_class)
