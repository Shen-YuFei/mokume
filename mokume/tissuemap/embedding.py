"""
Dimensionality reduction: PCA + t-SNE.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from mokume.tissuemap.config import EmbeddingConfig

logger = logging.getLogger(__name__)

_MIN_PROTEINS_FOR_PCA = 500


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
        threshold * 100, target,
    )
    return threshold


def _select_low_nan_proteins(
    log2_matrix: np.ndarray,
    max_nan_frac: float,
) -> np.ndarray:
    """Return boolean mask of proteins with NaN fraction ≤ threshold."""
    nan_frac = np.isnan(log2_matrix).mean(axis=0)
    return nan_frac <= max_nan_frac


def embed(
    adata: ad.AnnData,
    config: EmbeddingConfig,
    n_jobs: int = 8,
) -> ad.AnnData:
    """Run PCA + t-SNE on the corrected layer.

    Uses a low-NaN protein subset for PCA to avoid zero-fill artifacts.
    Remaining NaN in the subset are filled with 0 for PCA input only.

    Parameters
    ----------
    adata : ad.AnnData
        Must have ``layers["log2_corrected"]``.
    config : EmbeddingConfig
        Embedding parameters.

    Returns
    -------
    ad.AnnData
        Updated with ``obsm["X_pca"]``, ``obsm["X_tsne"]``,
        and ``uns["embedding_metrics"]``.
    """
    x_data = adata.layers["log2_corrected"].copy()

    # Resolve NaN threshold: auto or configured
    nan_threshold = _resolve_nan_threshold(x_data, config.max_nan_frac_for_pca)

    # Select low-NaN proteins for embedding
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
    # Per-protein median fill for PCA (preserves variance structure better than zero)
    col_medians = np.nanmedian(x_sub, axis=0, keepdims=True)
    col_medians = np.where(np.isnan(col_medians), 0.0, col_medians)
    x_sub = np.where(np.isnan(x_sub), col_medians, x_sub).astype(np.float32)

    # PCA
    n_components = min(config.pca_components, x_sub.shape[0] - 1, x_sub.shape[1])
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=config.random_state)
    pca_emb = pca.fit_transform(x_sub)
    adata.obsm["X_pca"] = pca_emb
    adata.uns["pca_variance_ratio"] = pca.explained_variance_ratio_.copy()
    var_explained = pca.explained_variance_ratio_.sum()
    logger.info(
        "PCA: %d components, variance explained: %.1f%%",
        n_components,
        var_explained * 100,
    )

    # t-SNE
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
    tsne_emb = tsne.fit_transform(pca_emb)
    adata.obsm["X_tsne"] = tsne_emb
    logger.info("t-SNE done")

    metrics = {
        "pca_var_explained": round(var_explained, 4),
        "n_proteins_used": int(n_kept),
    }
    adata.uns["embedding_metrics"] = metrics
    return adata
