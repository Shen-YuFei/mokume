"""
Dimensionality reduction: protein selection → MNAR-aware imputation → PCA →
t-SNE / UMAP.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

try:
    import umap as _umap_module  # type: ignore[import-not-found]
except ImportError:  # pragma: no cover - optional dependency
    _umap_module = None

from mokume.imputation import (
    impute_bpca,
    impute_gms,
    impute_impseq,
    impute_impseqrob,
    impute_mice,
    impute_mindet,
    impute_minprob,
    impute_missing_values,
    impute_mle,
    impute_nbavg,
    impute_qrilc,
)
from mokume.tissuemap.config import EmbeddingConfig

logger = logging.getLogger(__name__)

_MIN_PROTEINS_FOR_PCA = 500

# Imputation methods that operate directly on a DataFrame and return one back.
# ``impute_missing_values`` already wraps "knn", "mean", "median", "constant",
# "most_frequent", "seqknn" and "missforest", so we only need to dispatch
# explicitly for the more specialised left-censored / iterative methods.
_DIRECT_IMPUTERS = {
    "mindet": impute_mindet,
    "minprob": impute_minprob,
    "qrilc": impute_qrilc,
    "bpca": impute_bpca,
    "gms": impute_gms,
    "mle": impute_mle,
    "mice": impute_mice,
    "nbavg": impute_nbavg,
    "impseq": impute_impseq,
    "impseqrob": impute_impseqrob,
}


def _resolve_nan_threshold(
    log2_matrix: np.ndarray,
    configured: float | None,
) -> float:
    """Determine the NaN fraction threshold for PCA protein selection.

    When *configured* is ``None`` (auto mode), start from 0.50 and relax
    upward in 0.05 steps until at least ``_MIN_PROTEINS_FOR_PCA`` proteins
    (or 30 % of total, whichever is smaller) are retained.
    """
    if configured is not None:
        return configured

    nan_frac = np.isnan(log2_matrix).mean(axis=0)
    target = min(_MIN_PROTEINS_FOR_PCA, int(log2_matrix.shape[1] * 0.3))
    target = max(target, 10)

    threshold = 0.50
    while threshold < 0.96:
        n_kept = int((nan_frac <= threshold).sum())
        if n_kept >= target:
            break
        threshold += 0.05

    logger.info(
        "Auto NaN threshold for PCA: %.0f%% (target >= %d proteins)",
        threshold * 100,
        target,
    )
    return threshold


def _select_low_nan_proteins(
    log2_matrix: np.ndarray,
    max_nan_frac: float,
) -> np.ndarray:
    """Return boolean mask of proteins with NaN fraction ≤ threshold."""
    nan_frac = np.isnan(log2_matrix).mean(axis=0)
    return nan_frac <= max_nan_frac


def _run_pca(
    adata: ad.AnnData,
    x_sub: np.ndarray,
    config: EmbeddingConfig,
) -> tuple[np.ndarray, float]:
    """Fit PCA and store results in adata. Return (pca_emb, var_explained)."""
    n_components = min(config.pca_components, x_sub.shape[0] - 1, x_sub.shape[1])
    pca = PCA(
        n_components=n_components,
        svd_solver="randomized",
        random_state=config.random_state,
    )
    pca_emb = pca.fit_transform(x_sub)
    adata.obsm["X_pca"] = pca_emb
    adata.uns["pca_variance_ratio"] = pca.explained_variance_ratio_.copy()
    var_explained = pca.explained_variance_ratio_.sum()
    logger.info(
        "PCA: %d components, variance explained: %.1f%%",
        n_components,
        var_explained * 100,
    )
    return pca_emb, var_explained


def _impute_for_pca(
    x_sub: np.ndarray,
    method: str,
    n_neighbors: int,
) -> np.ndarray:
    """Impute remaining NaN values in the protein-by-sample sub-matrix.

    The matrix is samples × proteins. We delegate to ``mokume.imputation`` so
    every CLI / config that already supports those methods stays consistent
    with the tissuemap embedding step.
    """
    if not np.isnan(x_sub).any():
        return x_sub.astype(np.float32)

    method = method.lower()
    df = pd.DataFrame(x_sub)

    try:
        if method in _DIRECT_IMPUTERS:
            imputed = _DIRECT_IMPUTERS[method](df)
        else:
            imputed = impute_missing_values(df, method=method, n_neighbors=n_neighbors)
    except (
        ValueError,
        RuntimeError,
        ArithmeticError,
        ImportError,
        KeyError,
        np.linalg.LinAlgError,
    ) as exc:
        # Imputers can fail on degenerate matrices (LinAlg), bad params (Value/
        # KeyError), missing optional deps (Import), or numerical edge cases
        # (Arithmetic/Runtime). AttributeError / TypeError are programmer bugs
        # and must NOT be swallowed here.
        logger.warning(
            "Imputation '%s' failed (%s); falling back to column median",
            method,
            exc,
        )
        col_medians = np.nanmedian(x_sub, axis=0, keepdims=True)
        col_medians = np.where(np.isnan(col_medians), 0.0, col_medians)
        return np.where(np.isnan(x_sub), col_medians, x_sub).astype(np.float32)

    out = np.asarray(imputed, dtype=np.float32)
    if np.isnan(out).any():
        # Some methods leave residual NaNs (e.g. all-NaN columns); patch with 0.
        out = np.where(np.isnan(out), 0.0, out)
    return out


def _run_tsne(
    adata: ad.AnnData,
    pca_emb: np.ndarray,
    config: EmbeddingConfig,
    n_jobs: int,
) -> None:
    """Fit t-SNE on PCA embedding and store in adata."""
    perplexity = min(config.tsne_perplexity, (adata.n_obs - 1) / 3.0)
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        learning_rate="auto",
        init="pca",
        random_state=config.random_state,
        max_iter=2000,
        n_jobs=n_jobs,
    )
    adata.obsm["X_tsne"] = tsne.fit_transform(pca_emb)
    logger.info("t-SNE done")


def _run_umap(
    adata: ad.AnnData,
    pca_emb: np.ndarray,
    config: EmbeddingConfig,
    n_jobs: int,
) -> None:
    """Fit UMAP on PCA embedding and store in adata.

    Writes the result to both ``X_umap`` (the canonical key) and ``X_tsne``
    (so downstream plotters that key off ``X_tsne`` keep working without
    changes). UMAP 1.0+ honours ``n_jobs`` only when ``random_state`` is
    ``None``; otherwise it must stay single-threaded for reproducibility.
    """
    if _umap_module is None:
        raise ImportError(
            "UMAP requested but 'umap-learn' is not installed. "
            "Install it with `pip install umap-learn`."
        )

    n_neighbors = min(config.umap_n_neighbors, max(2, adata.n_obs - 1))
    umap_jobs = 1 if config.random_state is not None else n_jobs
    reducer = _umap_module.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=config.umap_min_dist,
        random_state=config.random_state,
        n_jobs=umap_jobs,
    )
    embedding = reducer.fit_transform(pca_emb)
    adata.obsm["X_umap"] = embedding
    # Reuse the existing plot key so downstream code does not need changes.
    adata.obsm["X_tsne"] = embedding
    logger.info("UMAP done (n_neighbors=%d)", n_neighbors)


_EMBEDDERS = {
    "tsne": _run_tsne,
    "umap": _run_umap,
}


def embed(
    adata: ad.AnnData,
    config: EmbeddingConfig,
    n_jobs: int = 8,
) -> ad.AnnData:
    """Run protein selection → imputation → PCA → t-SNE / UMAP.

    Uses a low-NaN protein subset for PCA, then fills the remaining missing
    values with the configured imputation method (delegated to
    :mod:`mokume.imputation`). Defaults to ``mindet`` because proteomics
    missingness is mostly MNAR.

    Parameters
    ----------
    adata : ad.AnnData
        Must have ``layers["log2_corrected"]``.
    config : EmbeddingConfig
        Embedding parameters.

    Returns
    -------
    ad.AnnData
        Updated with ``obsm["X_pca"]``, ``obsm["X_tsne"]`` (always — UMAP
        results are also mirrored here), and ``uns["embedding_metrics"]``.
        When ``embedding_method == "umap"`` the result is also stored at
        ``obsm["X_umap"]``.
    """
    x_data = adata.layers["log2_corrected"].copy()
    nan_threshold = _resolve_nan_threshold(x_data, config.max_nan_frac_for_pca)

    keep_mask = _select_low_nan_proteins(x_data, nan_threshold)
    n_kept = keep_mask.sum()
    logger.info(
        "Embedding: using %d / %d proteins (NaN <= %.0f%%)",
        n_kept,
        x_data.shape[1],
        nan_threshold * 100,
    )
    if n_kept < 10:
        logger.warning("Too few proteins for embedding (%d), skipping", n_kept)
        return adata

    x_sub = x_data[:, keep_mask]
    logger.info("Imputation: %s", config.imputation_method)
    x_sub = _impute_for_pca(
        x_sub, config.imputation_method, config.imputation_n_neighbors
    )

    pca_emb, var_explained = _run_pca(adata, x_sub, config)

    method = (config.embedding_method or "tsne").lower()
    embedder = _EMBEDDERS.get(method)
    if embedder is None:
        logger.warning(
            "Unknown embedding_method '%s'; falling back to t-SNE",
            config.embedding_method,
        )
        method = "tsne"
        embedder = _EMBEDDERS["tsne"]
    embedder(adata, pca_emb, config, n_jobs)

    adata.uns["embedding_metrics"] = {
        "pca_var_explained": round(var_explained, 4),
        "n_proteins_used": int(n_kept),
        "imputation_method": config.imputation_method,
        "embedding_method": method,
    }
    return adata
