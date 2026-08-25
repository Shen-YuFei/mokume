"""Rust-backed experiment executor for plugin recommendations."""

from __future__ import annotations

from dataclasses import dataclass
import threading
from typing import Any

import numpy as np
import pandas as pd

from mokume import differential_expression, impute_matrix, normalize_matrix
from mokume.agentic.contract import validate_config_values
from mokume.agentic.state import CandidateConfig
from mokume.core.logger import get_logger
from mokume.imputation.missforest import impute_missforest

logger = get_logger("mokume.agentic.runner")


@dataclass(frozen=True)
class ExperimentContext:
    """Immutable sample design and runtime settings for one candidate round."""

    sample_to_condition: dict[str, str]
    contrast: tuple[str, str]
    peptide_counts: pd.Series | None = None
    threads: int = 24


class PreprocessCache:
    """Memoize normalization and imputation pairs within one evaluation."""

    def __init__(self, threads: int = 1) -> None:
        self._store: dict[tuple[str, str], pd.DataFrame] = {}
        self._lock = threading.Lock()
        self._threads = threads
        self.hits = 0
        self.misses = 0

    def get_or_compute(
        self,
        norm_method: str,
        imp_method: str,
        protein_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return a preprocessed matrix, computing it once per method pair."""
        key = (norm_method, imp_method)
        with self._lock:
            cached = self._store.get(key)
            if cached is not None:
                self.hits += 1
                return cached.copy()
            self.misses += 1
        normalized = _apply_normalization(protein_df, norm_method, self._threads)
        imputed = _apply_imputation(normalized, imp_method, self._threads)
        with self._lock:
            self._store.setdefault(key, imputed)
            return self._store[key].copy()

    def stats(self) -> dict[str, int]:
        """Return cache counters for the audit artifact."""
        with self._lock:
            return {
                "hits": self.hits,
                "misses": self.misses,
                "unique_combos": len(self._store),
            }


def _apply_normalization(
    protein_df: pd.DataFrame,
    method: str,
    threads: int,
) -> pd.DataFrame:
    """Apply normalization through the Rust matrix API."""
    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)
    if np.isinf(matrix.to_numpy(dtype=float)).any():
        raise ValueError("Protein matrix contains an infinite intensity")
    normalized = normalize_matrix(
        _optional_values(matrix),
        method,
        [str(column) for column in matrix.columns],
        threads,
    )
    return _matrix_frame(protein_df, normalized)


def _apply_imputation(
    protein_df: pd.DataFrame,
    method: str,
    threads: int,
) -> pd.DataFrame:
    """Apply Rust imputation, using the documented missForest fallback."""
    protein_col = protein_df.columns[0]
    matrix = protein_df.set_index(protein_col)
    values = matrix.to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError("Protein matrix contains an infinite intensity")

    if method == "missforest":
        empty_cols = matrix.columns[matrix.isna().all()]
        workable = matrix.drop(columns=empty_cols)
        if workable.empty:
            return protein_df.copy()
        log2_matrix = np.log2(workable.where(workable > 0))
        imputed = impute_missforest(log2_matrix, n_jobs=threads)
        with np.errstate(over="raise", invalid="raise"):
            imputed = np.exp2(imputed)
        if len(empty_cols):
            imputed = imputed.join(matrix[empty_cols])[matrix.columns]
        return imputed.reset_index()

    imputed = impute_matrix(_optional_values(matrix), method, threads=threads)
    return _matrix_frame(protein_df, imputed)


def _optional_values(matrix: pd.DataFrame) -> list[list[float | None]]:
    """Encode non-finite matrix cells as missing for the PyO3 boundary."""
    return [
        [float(value) if np.isfinite(value) else None for value in row]
        for row in matrix.to_numpy(dtype=float)
    ]


def _matrix_frame(
    original: pd.DataFrame,
    values: list[list[float | None]],
) -> pd.DataFrame:
    """Restore a Rust matrix result to the input's wide-table shape."""
    protein_col = original.columns[0]
    matrix = pd.DataFrame(
        values,
        index=original[protein_col],
        columns=original.columns[1:],
        dtype=float,
    )
    matrix.index.name = protein_col
    return matrix.reset_index()


def run_experiment(
    config: CandidateConfig,
    protein_df: pd.DataFrame,
    context: ExperimentContext,
    cache: PreprocessCache | None = None,
) -> pd.DataFrame:
    """Run one validated preprocessing and differential-expression config."""
    validate_config_values(config.to_dict())
    ensemble_methods = (
        config.ensemble.split(",") if config.de_method == "ensemble" else None
    )
    if cache is None:
        normalized = _apply_normalization(
            protein_df, config.normalization, context.threads
        )
        processed = _apply_imputation(normalized, config.imputation, context.threads)
    else:
        processed = cache.get_or_compute(
            config.normalization,
            config.imputation,
            protein_df,
        )
    result = _run_rust_de(
        processed,
        config,
        context,
        ensemble_methods,
    )
    logger.info("Experiment %s complete: %d proteins tested", config.name, len(result))
    return result


def _run_rust_de(
    protein_df: pd.DataFrame,
    config: CandidateConfig,
    context: ExperimentContext,
    ensemble_methods: list[str] | None,
) -> pd.DataFrame:
    """Order one contrast and run the Rust matrix-level DE API."""
    group_a, group_b = context.contrast
    protein_col, samples_a, samples_b = _contrast_columns(protein_df, context)
    ordered = protein_df[[protein_col, *samples_a, *samples_b]]
    proteins = ordered[protein_col].astype(str).tolist()
    options: dict[str, Any] = {
        "fdr_threshold": config.fdr_threshold,
        "condition_a": group_a,
        "condition_b": group_b,
        "threads": context.threads,
    }
    if config.de_method not in {"rots", "limrots"}:
        options["fdr_method"] = config.fdr_method
    if config.de_method == "ensemble":
        options["ensemble_methods"] = ensemble_methods
        options["ensemble_min_k"] = config.ensemble_k
    if config.log2fc_threshold == "auto":
        options["effect_size_gate"] = "mixture"
    else:
        options["log2fc_threshold"] = float(config.log2fc_threshold)
    if context.peptide_counts is not None:
        options["peptide_counts"] = (
            context.peptide_counts.reindex(proteins).fillna(1.0).astype(float).tolist()
        )

    rows = differential_expression(
        proteins,
        _optional_values(ordered.iloc[:, 1:]),
        len(samples_a),
        len(samples_b),
        config.de_method,
        **options,
    )
    return pd.DataFrame(rows)


def _contrast_columns(
    protein_df: pd.DataFrame,
    context: ExperimentContext,
) -> tuple[Any, list[Any], list[Any]]:
    """Resolve and validate matrix columns for the requested contrast."""
    group_a, group_b = context.contrast
    protein_col = protein_df.columns[0]
    samples_a = [
        column
        for column in protein_df.columns[1:]
        if context.sample_to_condition.get(str(column)) == group_a
    ]
    samples_b = [
        column
        for column in protein_df.columns[1:]
        if context.sample_to_condition.get(str(column)) == group_b
    ]
    if len(samples_a) < 2 or len(samples_b) < 2:
        raise ValueError(
            f"Contrast {group_a!r} vs {group_b!r} requires at least two "
            "matrix columns per group"
        )
    return protein_col, samples_a, samples_b
