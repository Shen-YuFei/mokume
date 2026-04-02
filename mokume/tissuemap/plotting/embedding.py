"""
Embedding diagnostic plots: PCA scree.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)

def plot_pca_scree(
    adata: ad.AnnData,
    out_dir: Path,
    *,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """Two-panel figure: A. scree bar chart, B. cumulative variance curve.

    Requires ``adata.uns["pca_variance_ratio"]`` (set by
    :func:`mokume.tissuemap.embedding.embed`).
    """
    var_ratio = adata.uns.get("pca_variance_ratio")
    if var_ratio is None:
        logger.warning("No pca_variance_ratio in adata.uns, skipping scree plot")
        return

    var_ratio = np.asarray(var_ratio, dtype=float)
    n_pcs = len(var_ratio)
    cumvar = np.cumsum(var_ratio)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Scree
    ax = axes[0]
    ax.bar(
        range(1, n_pcs + 1), var_ratio * 100,
        color="#5C6BC0", edgecolor="white", linewidth=0.5, alpha=0.85,
    )
    ax.set_xlabel("Principal Component", fontsize=10)
    ax.set_ylabel("Variance explained (%)", fontsize=10)
    ax.set_title("A. Scree plot", fontsize=12, fontweight="bold", loc="left")
    ax.set_xlim(0.5, n_pcs + 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Annotate top 3 PCs
    for i in range(min(3, n_pcs)):
        ax.text(
            i + 1, var_ratio[i] * 100 + 0.3, f"{var_ratio[i] * 100:.1f}%",
            ha="center", fontsize=8, color="#333",
        )

    # Panel B: Cumulative
    ax = axes[1]
    ax.plot(
        range(1, n_pcs + 1), cumvar * 100,
        "o-", color="#5C6BC0", markersize=4, linewidth=2,
    )
    ax.axhline(70, color="#E53935", linestyle="--", linewidth=1, alpha=0.7)
    ax.text(n_pcs * 0.8, 71.5, "70%", fontsize=9, color="#E53935")
    ax.axhline(50, color="#FFA726", linestyle="--", linewidth=1, alpha=0.7)
    ax.text(n_pcs * 0.8, 51.5, "50%", fontsize=9, color="#FFA726")
    # Mark key thresholds
    for n_pc, label in [(2, "PC1-2"), (10, "PC1-10"), (20, "PC1-20")]:
        if n_pc <= n_pcs:
            ax.plot(n_pc, cumvar[n_pc - 1] * 100, "o", color="#E53935", markersize=7, zorder=3)
            ax.annotate(
                f"{label}\n{cumvar[n_pc - 1] * 100:.1f}%",
                (n_pc, cumvar[n_pc - 1] * 100),
                textcoords="offset points",
                xytext=(10, -15 if n_pc > 5 else 10),
                fontsize=7.5, color="#333",
                arrowprops={"arrowstyle": "-", "color": "#999", "lw": 0.5},
            )
    # Mark last PC
    ax.plot(n_pcs, cumvar[-1] * 100, "o", color="#E53935", markersize=7, zorder=3)
    ax.annotate(
        f"PC1-{n_pcs}\n{cumvar[-1] * 100:.1f}%",
        (n_pcs, cumvar[-1] * 100),
        textcoords="offset points", xytext=(10, -15),
        fontsize=7.5, color="#333",
        arrowprops={"arrowstyle": "-", "color": "#999", "lw": 0.5},
    )
    ax.set_xlabel("Number of PCs", fontsize=10)
    ax.set_ylabel("Cumulative variance (%)", fontsize=10)
    ax.set_title(
        "B. Cumulative variance", fontsize=12, fontweight="bold", loc="left",
    )
    ax.set_xlim(0.5, n_pcs + 0.5)
    ax.set_ylim(0, min(cumvar[-1] * 100 + 15, 100))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_dir / "pca_scree.png", dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(out_dir / "pca_scree.pdf", bbox_inches="tight")
    plt.close()
    logger.info("Saved pca_scree.png")
