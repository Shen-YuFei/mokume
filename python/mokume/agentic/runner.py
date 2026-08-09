"""Experiment executor for agentic analysis."""

import threading

import numpy as np
import pandas as pd

from mokume.agentic.state import CandidateConfig
from mokume.analysis.differential_expression import DifferentialExpression
from mokume.analysis.ensemble import run_ensemble
from mokume.analysis.ensemble_config import validate_de_ensemble
from mokume.core.logger import get_logger
from mokume.imputation.censored import impute_censored

logger = get_logger("mokume.agentic.runner")


class PreprocessCache:
    """Memoize (normalization, imputation) pairs within an optimization run.

    Many candidate configurations share the same (norm, imp) prefix and only
    differ in the DE method / thresholds. Caching the imputed matrix avoids
    re-running expensive normalization (LOESS / RLR) and imputation
    (KNN / impSeq / MinProb) work on every config. Thread-safe so it can back
    a joblib.Parallel(backend='threading') sweep without double-compute.
    """

    def __init__(self) -> None:
        self._store: dict[tuple[str, str], pd.DataFrame] = {}
        self._lock = threading.Lock()
        self.hits = 0
        self.misses = 0

    def get_or_compute(
        self,
        norm_method: str,
        imp_method: str,
        protein_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return imputed matrix for (norm, imp), computing on miss."""
        key = (norm_method, imp_method)
        with self._lock:
            cached = self._store.get(key)
            if cached is not None:
                self.hits += 1
                logger.debug(
                    "PreprocessCache hit for (%s, %s)", norm_method, imp_method
                )
                return cached.copy()
            self.misses += 1
        # Compute outside the lock so concurrent misses on distinct keys
        # can run in parallel; only the dict insert below is serialised.
        normed = _apply_normalization(protein_df, norm_method)
        imputed = _apply_imputation(normed, imp_method)
        with self._lock:
            # Another thread may have raced and populated the same key;
            # keep the first winner to keep results reproducible.
            self._store.setdefault(key, imputed)
            return self._store[key].copy()

    def stats(self) -> dict:
        """Return cache hit / miss counters for telemetry."""
        with self._lock:
            return {
                "hits": self.hits,
                "misses": self.misses,
                "unique_combos": len(self._store),
            }


def _apply_normalization(
    protein_df: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    """Apply normalization to the protein matrix."""
    if method == "none":
        return protein_df

    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)

    if method == "median":
        from mokume.normalization.distribution import MedianCenterNormalizer

        normalized = MedianCenterNormalizer().fit_transform(matrix)
    elif method == "quantile":
        from mokume.normalization.distribution import QuantileNormalizer

        normalized = QuantileNormalizer().fit_transform(matrix)
    elif method == "mean":
        from mokume.normalization.distribution import MeanCenterNormalizer

        normalized = MeanCenterNormalizer().fit_transform(matrix)
    elif method == "rlr":
        from mokume.normalization.rlr import rlr_normalize

        normalized = rlr_normalize(matrix)
    elif method == "loess":
        from mokume.normalization.loess import loess_normalize

        normalized = loess_normalize(matrix)
    else:
        logger.warning("Unknown normalization '%s', skipping", method)
        return protein_df

    return normalized.reset_index()


def _apply_imputation(
    protein_df: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    """Apply imputation to the protein matrix (log2 space)."""
    if method == "none":
        return protein_df

    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)

    # Ensure log2 space for imputation
    log2_matrix = matrix.copy()
    if log2_matrix.max().max() > 40:
        log2_matrix = np.log2(log2_matrix.replace(0, np.nan))

    imputed = impute_censored(log2_matrix, method=method)

    # Convert back if we transformed
    if matrix.max().max() > 40:
        imputed = 2**imputed

    return imputed.reset_index()


def run_experiment(
    config: CandidateConfig,
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    contrast: tuple[str, str],
    peptide_counts: pd.Series | None = None,
    cache: PreprocessCache | None = None,
) -> pd.DataFrame:
    """Run a single DE experiment with the given configuration.

    If ``cache`` is provided, (normalization, imputation) is memoized across
    configurations so the same prefix only runs once per round.
    """
    logger.info("Running experiment: %s", config.name)

    configured_ensemble = (
        None
        if config.ensemble.strip().lower() == "none"
        else config.ensemble.split(",")
    )
    ensemble_methods = validate_de_ensemble(
        enabled=True,
        method=config.de_method,
        ensemble_methods=configured_ensemble,
        ensemble_min_k=config.ensemble_k,
    )

    # 1+2. Normalization + imputation (cached when possible)
    if cache is not None:
        imputed_df = cache.get_or_compute(
            config.normalization, config.imputation, protein_df
        )
    else:
        normed_df = _apply_normalization(protein_df, config.normalization)
        imputed_df = _apply_imputation(normed_df, config.imputation)

    # 3. Run DE (single method or ensemble)
    if ensemble_methods is not None:
        result = run_ensemble(
            imputed_df,
            sample_to_condition,
            contrast,
            methods=ensemble_methods,
            min_k=config.ensemble_k,
            fdr_method=config.fdr_method,
            fdr_threshold=0.05,
            log2fc_threshold=config.log2fc_threshold,
            peptide_counts=peptide_counts,
        )
    else:
        de = DifferentialExpression(
            method=config.de_method,
            log2fc_threshold=config.log2fc_threshold,
            fdr_threshold=0.05,
            fdr_method=config.fdr_method,
            peptide_counts=peptide_counts,
        )
        result = de.run(imputed_df, sample_to_condition, contrast)
    logger.info(
        "Experiment %s complete: %d proteins tested",
        config.name,
        len(result),
    )
    return result
