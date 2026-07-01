"""
Tests for analysis, imputation, and normalization modules.

Covers: DifferentialExpression (limrots, deqms, proda, ihw),
        censored imputation (minprob, mindet, knn),
        LOESS normalization.

All DE methods are pure-Python reimplementations; no rpy2/R required.
"""

import numpy as np
import pandas as pd
import pytest

from mokume.analysis.differential_expression import DifferentialExpression
from mokume.analysis.deqms import run_deqms
from mokume.analysis.limma import run_limma
from mokume.analysis.limrots import run_limrots
from mokume.analysis.proda import run_proda
from mokume.analysis.rots import run_rots
from mokume.analysis.ensemble import run_ensemble
from mokume.imputation.censored import (
    impute_minprob,
    impute_mindet,
    impute_censored,
)
from mokume.normalization.loess import LOESSNormalizer, loess_normalize


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_protein_matrix(n_proteins=100, n_samples_per_group=4, seed=42):
    """Create a synthetic protein intensity matrix with known DE proteins."""
    rng = np.random.default_rng(seed)

    n_de = 15  # first 15 proteins are DE
    n_total = n_samples_per_group * 2

    # Background ~ N(23, 1)
    data = rng.normal(23, 1, size=(n_proteins, n_total))

    # Inject fold change for DE proteins in group A
    data[:n_de, :n_samples_per_group] += 3.0  # log2FC ~ 3

    proteins = [f"Protein_{i:03d}" for i in range(n_proteins)]
    samples_a = [f"A_{i}" for i in range(n_samples_per_group)]
    samples_b = [f"B_{i}" for i in range(n_samples_per_group)]

    df = pd.DataFrame(data, index=proteins, columns=samples_a + samples_b)

    # Sprinkle some NaN (5% missing, light enough for R packages)
    mask = rng.random(df.shape) < 0.05
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
# DifferentialExpression — LimROTS
# ---------------------------------------------------------------------------


class TestDELimROTS:
    def test_limrots_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="limrots",
            log2fc_threshold=1.0,
            skip_log2=True,
            n_boot=20,
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_limrots_detects_de(self):
        mat, sa, sb, _de_proteins = _make_protein_matrix(n_samples_per_group=5)
        wide, s2c = _make_wide_df(mat, sa, sb)
        res_lr = DifferentialExpression(
            method="limrots",
            log2fc_threshold=1.0,
            skip_log2=True,
            n_boot=20,
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
            method="deqms",
            log2fc_threshold=1.0,
            skip_log2=True,
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

    def test_deqms_with_peptide_counts(self):
        rng = np.random.default_rng(42)
        mat, sa, sb, _ = _make_protein_matrix()
        log2_mat = np.log2(mat.clip(lower=1))
        counts = pd.Series(
            rng.integers(2, 10, len(log2_mat)),
            index=log2_mat.index,
        )
        result = run_deqms(log2_mat, sa, sb, ("A", "B"), peptide_counts=counts)
        assert len(result) > 0
        assert "peptide_count" in result.columns


# ---------------------------------------------------------------------------
# DifferentialExpression — proDA
# ---------------------------------------------------------------------------


class TestDEProDA:
    def test_proda_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="proda", log2fc_threshold=1.0, skip_log2=True
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_run_proda_direct(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_proda(np.log2(mat.clip(lower=1)), sa, sb, "A", "B")
        assert len(result) > 0
        assert "pvalue" in result.columns

    def test_proda_has_expected_columns(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_proda(np.log2(mat.clip(lower=1)), sa, sb, "A", "B")
        assert len(result) > 0
        for col in ("ProteinName", "log2FC", "pvalue", "adj_pvalue"):
            assert col in result.columns

    def test_proda_uses_moderated_df(self):
        # proDA must test against the per-protein empirical-Bayes-moderated df
        # (fit["df"]), not the naive constant n_samples - p; the reported df is
        # therefore per-protein and not all equal to that residual constant.
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_proda(np.log2(mat.clip(lower=1)), sa, sb, "A", "B")
        assert "df" in result.columns
        naive_df = (len(sa) + len(sb)) - 2  # p = intercept + group
        df_values = result["df"].values
        assert np.all(df_values > 0)
        assert not np.allclose(df_values, naive_df)


# ---------------------------------------------------------------------------
# IHW correction
# ---------------------------------------------------------------------------


class TestIHW:
    def test_ihw_differs_from_bh(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        res_bh = DifferentialExpression(
            method="deqms",
            log2fc_threshold=1.0,
            fdr_method="bh",
            skip_log2=True,
        ).run(wide, s2c, ("A", "B"))
        res_ihw = DifferentialExpression(
            method="deqms",
            log2fc_threshold=1.0,
            fdr_method="ihw",
            skip_log2=True,
        ).run(wide, s2c, ("A", "B"))
        assert len(res_ihw) > 0
        assert "adj_pvalue" in res_ihw.columns
        # IHW should produce different adjusted p-values than BH
        merged = res_bh.merge(res_ihw, on="ProteinName", suffixes=("_bh", "_ihw"))
        assert not np.allclose(
            merged["adj_pvalue_bh"].values,
            merged["adj_pvalue_ihw"].values,
            equal_nan=True,
        )

    def test_ihw_fallback_on_small_data(self):
        # Very few proteins → should fall back to BH
        mat, sa, sb, _ = _make_protein_matrix(n_proteins=5)
        wide, s2c = _make_wide_df(mat, sa, sb)
        res_bh = DifferentialExpression(
            method="deqms",
            log2fc_threshold=1.0,
            fdr_method="bh",
            skip_log2=True,
        ).run(wide, s2c, ("A", "B"))
        res_ihw = DifferentialExpression(
            method="deqms",
            log2fc_threshold=1.0,
            fdr_method="ihw",
            skip_log2=True,
        ).run(wide, s2c, ("A", "B"))
        assert len(res_ihw) > 0
        # With too few proteins, IHW falls back to BH → same results
        merged = res_bh.merge(res_ihw, on="ProteinName", suffixes=("_bh", "_ihw"))
        assert np.allclose(
            merged["adj_pvalue_bh"].values,
            merged["adj_pvalue_ihw"].values,
            equal_nan=True,
        )

    def test_limrots_fdr_unchanged_by_ihw(self):
        # LimROTS carries its own permutation-based FDR in ``adj_pvalue``;
        # requesting IHW must not overwrite it, so the adjusted p-values are
        # identical under bh and ihw (same seed → same permutation FDR).
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        res_bh = DifferentialExpression(
            method="limrots",
            log2fc_threshold=1.0,
            fdr_method="bh",
            skip_log2=True,
            n_boot=20,
        ).run(wide, s2c, ("A", "B"))
        res_ihw = DifferentialExpression(
            method="limrots",
            log2fc_threshold=1.0,
            fdr_method="ihw",
            skip_log2=True,
            n_boot=20,
        ).run(wide, s2c, ("A", "B"))
        merged = res_bh.merge(res_ihw, on="ProteinName", suffixes=("_bh", "_ihw"))
        assert np.allclose(
            merged["adj_pvalue_bh"].values,
            merged["adj_pvalue_ihw"].values,
            equal_nan=True,
        )

    def test_rots_fdr_unchanged_by_ihw(self):
        # Same guarantee for ROTS: its permutation FDR survives an IHW request.
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        res_bh = DifferentialExpression(
            method="rots",
            log2fc_threshold=1.0,
            fdr_method="bh",
            skip_log2=True,
            n_boot=20,
        ).run(wide, s2c, ("A", "B"))
        res_ihw = DifferentialExpression(
            method="rots",
            log2fc_threshold=1.0,
            fdr_method="ihw",
            skip_log2=True,
            n_boot=20,
        ).run(wide, s2c, ("A", "B"))
        merged = res_bh.merge(res_ihw, on="ProteinName", suffixes=("_bh", "_ihw"))
        assert np.allclose(
            merged["adj_pvalue_bh"].values,
            merged["adj_pvalue_ihw"].values,
            equal_nan=True,
        )


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
        de = DifferentialExpression(
            method="deqms", log2fc_threshold=1.0, skip_log2=True
        )
        results = de.run_comparisons(wide, s2c, [("A", "B")])
        assert "A-B" in results
        assert len(results["A-B"]) > 0


# ---------------------------------------------------------------------------
# DifferentialExpression — limma
# ---------------------------------------------------------------------------


class TestDELimma:
    def test_limma_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="limma", log2fc_threshold=1.0, skip_log2=True
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_run_limma_direct(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_limma(np.log2(mat.clip(lower=1)), sa, sb, "A", "B")
        assert len(result) > 0
        assert "pvalue" in result.columns

    def test_limma_has_expected_columns(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_limma(np.log2(mat.clip(lower=1)), sa, sb, "A", "B")
        for col in ("ProteinName", "log2FC", "pvalue", "adj_pvalue", "t_stat"):
            assert col in result.columns


# ---------------------------------------------------------------------------
# DifferentialExpression — ROTS
# ---------------------------------------------------------------------------


class TestDEROTS:
    def test_rots_returns_results(self):
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        de = DifferentialExpression(
            method="rots", log2fc_threshold=1.0, skip_log2=True, n_boot=20
        )
        result = de.run(wide, s2c, ("A", "B"))
        assert len(result) > 0
        assert "adj_pvalue" in result.columns

    def test_run_rots_direct(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_rots(np.log2(mat.clip(lower=1)), sa, sb, "A", "B", n_boot=20)
        assert len(result) > 0
        assert "pvalue" in result.columns

    def test_rots_has_expected_columns(self):
        mat, sa, sb, _ = _make_protein_matrix()
        result = run_rots(np.log2(mat.clip(lower=1)), sa, sb, "A", "B", n_boot=20)
        for col in ("ProteinName", "log2FC", "pvalue", "adj_pvalue", "d_stat"):
            assert col in result.columns


# ---------------------------------------------------------------------------
# DifferentialExpression — non-finite log2FC guard
# ---------------------------------------------------------------------------


class TestNonFiniteLog2FCGuard:
    """A protein imputed to a degenerate value in one group can yield an
    Inf/-Inf log2FC; it must be coerced to NaN and excluded from significance
    rather than passing the ``|log2FC| > threshold`` filter."""

    def _finalize(self, log2fc, pvalue):
        de = DifferentialExpression(
            method="limma", log2fc_threshold=0.5, fdr_method="bh"
        )
        df = pd.DataFrame(
            {
                "ProteinName": [f"P{i}" for i in range(len(log2fc))],
                "log2FC": log2fc,
                "pvalue": pvalue,
            }
        )
        return de._finalize_results(df)

    def test_inf_log2fc_coerced_to_nan(self):
        out = self._finalize([np.inf, -np.inf, np.nan, 3.0], [1e-6, 1e-6, 1e-6, 1e-6])
        out = out.set_index("ProteinName")
        assert np.isnan(out.loc["P0", "log2FC"])  # +inf -> NaN
        assert np.isnan(out.loc["P1", "log2FC"])  # -inf -> NaN
        assert np.isfinite(out.loc["P3", "log2FC"])  # finite untouched

    def test_inf_log2fc_not_counted_significant(self):
        out = self._finalize([np.inf, -np.inf, 3.0, -3.0], [1e-6, 1e-6, 1e-6, 1e-6])
        out = out.set_index("ProteinName")
        # non-finite folds stay Unchanged despite a significant p-value
        assert out.loc["P0", "significance"] == "Unchanged"
        assert out.loc["P1", "significance"] == "Unchanged"
        # genuine finite folds are still classified
        assert out.loc["P2", "significance"] == "UP"
        assert out.loc["P3", "significance"] == "DOWN"


# ---------------------------------------------------------------------------
# DifferentialExpression — Ensemble
# ---------------------------------------------------------------------------


class TestDEEnsemble:
    def test_ensemble_combines_without_protein_keyerror(self):
        # Regression: run_ensemble previously crashed with KeyError 'Protein'
        # because combine_de_results read df["Protein"] while every member
        # method emits "ProteinName". It must now run and return a combined
        # table keyed by ProteinName with the consensus columns.
        mat, sa, sb, _ = _make_protein_matrix()
        wide, s2c = _make_wide_df(mat, sa, sb)
        result = run_ensemble(
            wide,
            s2c,
            ("A", "B"),
            methods=("limma", "deqms"),
            min_k=1,
            log2fc_threshold=1.0,
        )
        assert not result.empty
        assert "ProteinName" in result.columns
        assert "Protein" not in result.columns
        for col in (
            "log2FC",
            "pvalue",
            "adj_pvalue",
            "significance",
            "n_methods_up",
            "n_methods_down",
            "methods_significant",
        ):
            assert col in result.columns
