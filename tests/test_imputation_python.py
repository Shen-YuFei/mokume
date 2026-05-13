"""Unit tests for Python-native imputation methods (SeqKNN + missForest)."""

import numpy as np
import pandas as pd
import pytest

from mokume.imputation import (
    impute_censored,
    impute_gms,
    impute_mice,
    impute_missforest,
    impute_missing_values,
    impute_mle,
    impute_nbavg,
    impute_qrilc,
    impute_seqknn,
)


@pytest.fixture
def matrix_with_missing():
    """5 proteins x 4 samples matrix with a few NaN values seeded so
    sequential fill order is deterministic.
    """
    rng = np.random.default_rng(0)
    base = rng.normal(loc=20, scale=2, size=(5, 4))
    df = pd.DataFrame(
        base, index=[f"P{i}" for i in range(5)], columns=["S1", "S2", "S3", "S4"]
    )
    # Insert NaN: P0 row complete; P1 has 1 NaN; P2 has 2 NaN; P3 has 1 NaN; P4 complete
    df.iloc[1, 0] = np.nan
    df.iloc[2, 0] = np.nan
    df.iloc[2, 3] = np.nan
    df.iloc[3, 2] = np.nan
    return df


class TestSeqKNN:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_seqknn(matrix_with_missing, k=2)
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape

    def test_preserves_observed_values(self, matrix_with_missing):
        result = impute_seqknn(matrix_with_missing, k=2)
        observed_mask = matrix_with_missing.notna()
        np.testing.assert_array_almost_equal(
            result[observed_mask].dropna(how="all").values,
            matrix_with_missing[observed_mask].dropna(how="all").values,
        )

    def test_k_must_be_positive(self, matrix_with_missing):
        with pytest.raises(ValueError, match="k must be"):
            impute_seqknn(matrix_with_missing, k=0)

    def test_no_missing_returns_copy(self):
        df = pd.DataFrame(np.arange(12).reshape(3, 4).astype(float))
        result = impute_seqknn(df)
        pd.testing.assert_frame_equal(result, df)

    def test_empty_returns_empty(self):
        empty = pd.DataFrame()
        assert impute_seqknn(empty).empty

    def test_no_complete_rows_seeds_with_column_mean(self):
        # Every row has at least one NaN
        df = pd.DataFrame(
            {
                "S1": [1.0, 2.0, np.nan, 4.0],
                "S2": [np.nan, 3.0, 4.0, 5.0],
                "S3": [3.0, np.nan, 5.0, 6.0],
            }
        )
        result = impute_seqknn(df, k=2)
        assert not result.isna().any().any()


class TestMissForest:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_missforest(
            matrix_with_missing, n_estimators=20, max_iter=3, random_state=0
        )
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape

    def test_no_missing_returns_copy(self):
        df = pd.DataFrame(np.arange(12).reshape(3, 4).astype(float))
        result = impute_missforest(df)
        pd.testing.assert_frame_equal(result, df)

    def test_empty_returns_empty(self):
        assert impute_missforest(pd.DataFrame()).empty

    def test_reproducible_with_random_state(self, matrix_with_missing):
        a = impute_missforest(
            matrix_with_missing, n_estimators=20, max_iter=3, random_state=42
        )
        b = impute_missforest(
            matrix_with_missing, n_estimators=20, max_iter=3, random_state=42
        )
        pd.testing.assert_frame_equal(a, b)


class TestDispatchers:
    def test_censored_dispatcher_seqknn(self, matrix_with_missing):
        result = impute_censored(matrix_with_missing, method="seqknn", k=2)
        assert not result.isna().any().any()

    def test_censored_dispatcher_missforest(self, matrix_with_missing):
        result = impute_censored(
            matrix_with_missing,
            method="missforest",
            n_estimators=10,
            max_iter=2,
        )
        assert not result.isna().any().any()

    def test_censored_dispatcher_unknown_lists_new_methods(self):
        with pytest.raises(ValueError, match="seqknn"):
            impute_censored(pd.DataFrame({"a": [np.nan]}), method="not-a-method")

    def test_impute_missing_values_seqknn(self, matrix_with_missing):
        result = impute_missing_values(
            matrix_with_missing, method="seqknn", n_neighbors=2
        )
        assert isinstance(result, pd.DataFrame)
        assert not result.isna().any().any()

    def test_impute_missing_values_missforest(self, matrix_with_missing):
        result = impute_missing_values(matrix_with_missing, method="missforest")
        assert isinstance(result, pd.DataFrame)
        assert not result.isna().any().any()

    def test_censored_dispatcher_qrilc(self, matrix_with_missing):
        result = impute_censored(matrix_with_missing, method="qrilc")
        assert not result.isna().any().any()

    def test_censored_dispatcher_mice(self, matrix_with_missing):
        result = impute_censored(matrix_with_missing, method="mice")
        assert not result.isna().any().any()

    def test_censored_dispatcher_nbavg(self, matrix_with_missing):
        result = impute_censored(matrix_with_missing, method="nbavg")
        assert not result.isna().any().any()

    def test_censored_dispatcher_gms(self, matrix_with_missing):
        result = impute_censored(matrix_with_missing, method="gms")
        assert not result.isna().any().any()

    def test_censored_dispatcher_mle(self, matrix_with_missing):
        result = impute_censored(matrix_with_missing, method="mle")
        assert not result.isna().any().any()


class TestQRILC:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_qrilc(matrix_with_missing)
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape

    def test_preserves_observed_values(self, matrix_with_missing):
        result = impute_qrilc(matrix_with_missing)
        observed = matrix_with_missing.notna()
        np.testing.assert_array_almost_equal(
            result[observed].dropna(how="all").values,
            matrix_with_missing[observed].dropna(how="all").values,
        )

    def test_no_missing_returns_copy(self):
        df = pd.DataFrame(np.arange(12).reshape(3, 4).astype(float))
        result = impute_qrilc(df)
        pd.testing.assert_frame_equal(result, df)


class TestMICE:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_mice(matrix_with_missing, max_iter=3)
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape

    def test_reproducible(self, matrix_with_missing):
        a = impute_mice(matrix_with_missing, random_state=0, max_iter=3)
        b = impute_mice(matrix_with_missing, random_state=0, max_iter=3)
        pd.testing.assert_frame_equal(a, b)


class TestNbavg:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_nbavg(matrix_with_missing, k=2)
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape

    def test_preserves_observed_values(self, matrix_with_missing):
        result = impute_nbavg(matrix_with_missing, k=2)
        observed = matrix_with_missing.notna()
        np.testing.assert_array_almost_equal(
            result[observed].dropna(how="all").values,
            matrix_with_missing[observed].dropna(how="all").values,
        )


class TestGMS:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_gms(matrix_with_missing)
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape


class TestMLE:
    def test_fills_all_nans(self, matrix_with_missing):
        result = impute_mle(matrix_with_missing, max_iter=10)
        assert not result.isna().any().any()
        assert result.shape == matrix_with_missing.shape

    def test_preserves_observed_values(self, matrix_with_missing):
        result = impute_mle(matrix_with_missing)
        observed = matrix_with_missing.notna()
        np.testing.assert_array_almost_equal(
            result[observed].dropna(how="all").values,
            matrix_with_missing[observed].dropna(how="all").values,
            decimal=3,
        )
