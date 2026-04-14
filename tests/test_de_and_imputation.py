"""
Tests for new analysis, imputation, and normalization modules.

Covers: DifferentialExpression (limrots, deqms, proda, ihw),
        censored imputation (minprob, mindet, knn),
        LOESS normalization.
"""

import warnings

import numpy as np
import pandas as pd
import pytest

from mokume.analysis.differential_expression import (
    DifferentialExpression,
    _trigamma,
    _tetragamma,
)
from mokume.analysis.deqms import run_deqms
from mokume.analysis.limrots import run_limrots
from mokume.analysis.proda import run_proda, DropoutParams
from mokume.imputation.censored import (
    classify_missing,
    impute_minprob,
    impute_mindet,
    impute_censored,
)
from mokume.normalization.loess import LOESSNormalizer, loess_normalize


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_protein_matrix(n_proteins=50, n_samples_per_group=3, seed=42):
    """Create a synthetic protein intensity matrix with known DE proteins."""
    rng = np.random.default_rng(seed)

    n_de = 10  # first 10 proteins are DE
    n_total = n_samples_per_group * 2

    # Background ~ N(23, 1)
    data = rng.normal(23, 1, size=(n_proteins, n_total))

    # Inject fold change for DE proteins in group A
    data[:n_de, :n_samples_per_group] += 3.0  # log2FC ~ 3

    proteins = [f"Protein_{i:03d}" for i in range(n_proteins)]
    samples_a = [f"A_{i}" for i in range(n_samples_per_group)]
    samples_b = [f"B_{i}" for i in range(n_samples_per_group)]

    df = pd.DataFrame(data, index=proteins, columns=samples_a + samples_b)

    # Sprinkle some NaN (10% missing)
    mask = rng.random(df.shape) < 0.10
    df[mask] = np.nan

    return df, samples_a, samples_b, proteins[:n_de]


def _make_wide_df(matrix, samples_a, samples_b):
    """Convert index-based matrix to wide-format with ProteinName column."""
    wide = matrix.copy()
    wide.insert(0, "ProteinName", wide.index)
    wide = wide.reset_index(drop=True)
    s2c = {s: "A" for s in samples_a}
    s2c.update({s: "B" for s in samples_b})
    return wide, s2c


# ---------------------------------------------------------------------------
# _trigamma / _tetragamma
# ---------------------------------------------------------------------------

class TestMathHelpers:
    def test_trigamma_positive(self):
        assert _trigamma(1.0) > 0
        assert np.isfinite(_trigamma(1.0))

    def test_tetragamma_negative(self):
        assert _tetragamma(1.0) < 0
        assert np.isfinite(_tetragamma(1.0))


# ---------------------------------------------------------------------------
# DifferentialExpression — LimROTS
# ---------------------------------------------------------------------------

class TestDELimROTS:
    def test_limrots_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="limrots", log2fc_threshold=1.0, skip_log2=True, n_boot=20,
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_limrots_detects_de(self):
        mat, sa, sb, _de_proteins = _make_protein_matrix(n_samples_per_group=5)
        wide, s2c = _make_wide_df(mat, sa, sb)
        res_lr = DifferentialExpression(
            method="limrots", log2fc_threshold=1.0, skip_log2=True, n_boot=20,
        ).run(wide, s2c, ("A", "B"))
        sig = res_lr[res_lr["significance"] == "UP"]["ProteinName"].tolist()
        assert len(sig) > 0

    def test_run_limrots_direct(self):
        mat, sa, sb, _ = _make_protein_matrix()
        log2_mat = np.log2(mat.clip(lower=1))
        result = run_limrots(log2_mat, sa, sb, ("A", "B"), n_boot=10)
        assert len(result) > 0
        assert "pvalue" in result.columns


# ---------------------------------------------------------------------------
# DifferentialExpression — DEqMS
# ---------------------------------------------------------------------------

class TestDEDEqMS:
    def test_deqms_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="deqms", log2fc_threshold=1.0, skip_log2=True,
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_run_deqms_direct(self):
        mat, sa, sb, _ = _make_protein_matrix()
        log2_mat = np.log2(mat.clip(lower=1))
        result = run_deqms(log2_mat, sa, sb, ("A", "B"))
        assert len(result) > 0
        assert "pvalue" in result.columns

    def test_run_deqms_no_lowess_runtime_warning_without_peptide_counts(self):
        mat, sa, sb, _ = _make_protein_matrix()
        log2_mat = np.log2(mat.clip(lower=1))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = run_deqms(log2_mat, sa, sb, ("A", "B"))
        assert len(result) > 0
        runtime_warnings = [
            w for w in caught
            if issubclass(w.category, RuntimeWarning)
        ]
        assert runtime_warnings == []


# ---------------------------------------------------------------------------
# DifferentialExpression — proDA
# ---------------------------------------------------------------------------

class TestDEProDA:
    def test_proda_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(method="proda", log2fc_threshold=1.0, skip_log2=True)
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_run_proda_direct(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_proda(np.log2(mat.clip(lower=1)), sa, sb, "A", "B")
        assert len(result) > 0
        assert "pvalue" in result.columns

    def test_dropout_params_namedtuple(self):
        dp = DropoutParams(rho=np.array([1.0]), zeta=np.array([0.5]))
        assert dp.rho[0] == 1.0
        assert dp.zeta[0] == 0.5


# ---------------------------------------------------------------------------
# IHW correction
# ---------------------------------------------------------------------------

class TestIHW:
    def test_ihw_via_de(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="deqms", log2fc_threshold=1.0, fdr_method="ihw", skip_log2=True,
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_ihw_fallback_on_small_data(self):
        # Very few proteins → should fall back to BH
        mat, sa, sb, _ = _make_protein_matrix(n_proteins=5)
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="deqms", log2fc_threshold=1.0, fdr_method="ihw", skip_log2=True,
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0


# ---------------------------------------------------------------------------
# Censored imputation
# ---------------------------------------------------------------------------

class TestCensoredImputation:
    def test_minprob_fills_nan(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(
            rng.normal(23, 1, (20, 6)),
            columns=[f"S{i}" for i in range(6)],
        )
        data.iloc[0, 0] = np.nan
        data.iloc[5, 3] = np.nan
        result = impute_minprob(data)
        assert result.isna().sum().sum() == 0

    def test_mindet_fills_nan(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(
            rng.normal(23, 1, (20, 6)),
            columns=[f"S{i}" for i in range(6)],
        )
        data.iloc[0, 0] = np.nan
        result = impute_mindet(data)
        assert result.isna().sum().sum() == 0

    def test_impute_censored_dispatcher(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(rng.normal(23, 1, (20, 6)))
        data.iloc[0, 0] = np.nan
        for method in ("minprob", "mindet", "knn", "none"):
            result = impute_censored(data, method=method)
            if method != "none":
                assert result.isna().sum().sum() == 0

    def test_impute_censored_invalid_method(self):
        data = pd.DataFrame({"a": [1, 2, 3]})
        with pytest.raises(ValueError, match="Unknown imputation method"):
            impute_censored(data, method="invalid")

    def test_classify_missing(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(
            rng.normal(23, 1, (20, 6)),
            columns=[f"S{i}" for i in range(6)],
        )
        data.iloc[0, 0] = np.nan
        is_mnar = classify_missing(data)
        assert is_mnar.shape == data.shape
        assert is_mnar.dtypes.unique().tolist() == [np.dtype("bool")]


# ---------------------------------------------------------------------------
# LOESS normalization
# ---------------------------------------------------------------------------

class TestLOESS:
    def test_loess_normalize_shape(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(
            rng.normal(23, 1, (100, 6)),
            columns=[f"S{i}" for i in range(6)],
        )
        # Add systematic bias to one sample
        data["S0"] += np.linspace(0, 3, 100)
        result = loess_normalize(data)
        assert result.shape == data.shape
        assert result.isna().sum().sum() == 0

    def test_loess_normalizer_class(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(
            rng.normal(23, 1, (100, 6)),
            columns=[f"S{i}" for i in range(6)],
        )
        norm = LOESSNormalizer(frac=0.75, reference="median")
        result = norm.fit_transform(data)
        assert result.shape == data.shape

    def test_loess_reference_mean(self):
        rng = np.random.default_rng(42)
        data = pd.DataFrame(rng.normal(23, 1, (100, 6)))
        result = loess_normalize(data, reference="mean")
        assert result.shape == data.shape

    def test_loess_invalid_reference(self):
        data = pd.DataFrame({"a": [1, 2, 3]})
        with pytest.raises(ValueError, match="Unknown reference"):
            loess_normalize(data, reference="invalid")


# ---------------------------------------------------------------------------
# run_comparisons
# ---------------------------------------------------------------------------

class TestRunComparisons:
    def test_multiple_contrasts(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(method="deqms", log2fc_threshold=1.0, skip_log2=True)
        results = de.run_comparisons(wide, s2c, [("A", "B")])
        assert "A-B" in results
        assert len(results["A-B"]) > 0
