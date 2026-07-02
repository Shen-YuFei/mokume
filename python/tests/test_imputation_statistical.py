"""Tests for statistical imputation methods (impSeq, impSeqRob, BPCA).

These methods are pure-Python reimplementations of the corresponding R
algorithms (rrcovNA::impSeq, rrcovNA::impSeqRob, pcaMethods::bpca).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mokume.imputation.impseq import impute_impseq
from mokume.imputation.impseqrob import impute_impseqrob
from mokume.imputation.bpca import impute_bpca
from mokume.imputation.censored import impute_censored


@pytest.fixture()
def imputation_matrix():
    """Create a synthetic protein matrix with controlled missingness."""
    rng = np.random.default_rng(42)
    data = rng.normal(23, 1, size=(30, 6))
    mask = rng.random(data.shape) < 0.08
    data[mask] = np.nan
    proteins = [f"Protein_{i:03d}" for i in range(30)]
    samples = [f"Sample_{i}" for i in range(6)]
    return pd.DataFrame(data, index=proteins, columns=samples)


class TestImpSeq:
    """Tests for the impSeq pure-Python reimplementation."""

    def test_fills_all_nan(self, imputation_matrix):
        result = impute_impseq(imputation_matrix)
        assert result.isna().sum().sum() == 0

    def test_preserves_shape(self, imputation_matrix):
        result = impute_impseq(imputation_matrix)
        assert result.shape == imputation_matrix.shape

    def test_preserves_index(self, imputation_matrix):
        result = impute_impseq(imputation_matrix)
        assert list(result.index) == list(imputation_matrix.index)
        assert list(result.columns) == list(imputation_matrix.columns)

    def test_observed_values_unchanged(self, imputation_matrix):
        result = impute_impseq(imputation_matrix)
        observed = imputation_matrix.notna()
        np.testing.assert_allclose(
            result.values[observed.values],
            imputation_matrix.values[observed.values],
            atol=1e-10,
        )

    def test_dispatcher_integration(self, imputation_matrix):
        result = impute_censored(imputation_matrix, method="impseq")
        assert result.isna().sum().sum() == 0
        assert result.shape == imputation_matrix.shape


class TestImpSeqRob:
    """Tests for the impSeqRob pure-Python reimplementation."""

    def test_fills_all_nan(self, imputation_matrix):
        result = impute_impseqrob(imputation_matrix)
        assert result.isna().sum().sum() == 0

    def test_preserves_shape(self, imputation_matrix):
        result = impute_impseqrob(imputation_matrix)
        assert result.shape == imputation_matrix.shape

    def test_alpha_parameter(self, imputation_matrix):
        r1 = impute_impseqrob(imputation_matrix, alpha=0.5)
        r2 = impute_impseqrob(imputation_matrix, alpha=0.9)
        assert r1.isna().sum().sum() == 0
        assert r2.isna().sum().sum() == 0

    def test_dispatcher_integration(self, imputation_matrix):
        result = impute_censored(imputation_matrix, method="impseqrob")
        assert result.isna().sum().sum() == 0

    def test_raises_above_max_variables_guard(self):
        from mokume.imputation.impseqrob import _MAX_VARIABLES_GUARD

        rng = np.random.default_rng(0)
        wide = pd.DataFrame(rng.normal(20, 1, size=(20, _MAX_VARIABLES_GUARD + 1)))
        wide.iloc[0, 0] = np.nan
        with pytest.raises(ValueError, match="max_variables"):
            impute_impseqrob(wide)

    def test_max_variables_none_disables_guard(self):
        rng = np.random.default_rng(0)
        wide = pd.DataFrame(rng.normal(20, 1, size=(10, 8)))
        wide.iloc[0, 0] = np.nan
        result = impute_impseqrob(wide, max_variables=None)
        assert result.isna().sum().sum() == 0


class TestBPCA:
    """Tests for the BPCA pure-Python reimplementation."""

    def test_fills_all_nan(self, imputation_matrix):
        result = impute_bpca(imputation_matrix)
        assert result.isna().sum().sum() == 0

    def test_preserves_shape(self, imputation_matrix):
        result = impute_bpca(imputation_matrix)
        assert result.shape == imputation_matrix.shape

    def test_preserves_index(self, imputation_matrix):
        result = impute_bpca(imputation_matrix)
        assert list(result.index) == list(imputation_matrix.index)
        assert list(result.columns) == list(imputation_matrix.columns)

    def test_npcs_parameter(self, imputation_matrix):
        result = impute_bpca(imputation_matrix, n_pcs=3)
        assert result.isna().sum().sum() == 0

    def test_dispatcher_integration(self, imputation_matrix):
        result = impute_censored(imputation_matrix, method="bpca")
        assert result.isna().sum().sum() == 0
