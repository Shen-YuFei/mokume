"""Embedding diagnostic plots: PCA scree and cumulative variance."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

from mokume.plotting._common import style_axes
from mokume.plotting.style import (
    SAMPLE_POINT_COLOR,
    SINGLE_SERIES_COLOR,
    TEXT_COLOR,
    WHITE,
)

logger = logging.getLogger(__name__)


def _plot_scree_panel(ax, var_ratio: np.ndarray, n_pcs: int) -> None:
    """Panel A: scree bar chart."""
    ax.bar(
        range(1, n_pcs + 1),
        var_ratio * 100,
        color=SINGLE_SERIES_COLOR,
        width=0.7,
        linewidth=0,
    )
    ax.set_xlabel("Principal Component", fontsize=10)
    ax.set_ylabel("Variance explained (%)", fontsize=10)
    ax.set_title("A. Scree plot", fontsize=12, fontweight="bold", loc="left")
    ax.set_xlim(0.5, n_pcs + 0.5)
    ax.margins(y=0.15)
    style_axes(ax)
    for i in range(min(3, n_pcs)):
        ax.text(
            i + 1.1,
            var_ratio[i] * 100 + 0.3,
            f"{var_ratio[i] * 100:.1f}%",
            ha="left",
            fontsize=8,
            color=TEXT_COLOR,
        )


def _plot_cumulative_panel(ax, cumvar: np.ndarray, n_pcs: int) -> None:
    """Panel B: cumulative variance curve with neutral threshold annotations."""
    ax.plot(
        range(1, n_pcs + 1),
        cumvar * 100,
        "o-",
        color=SINGLE_SERIES_COLOR,
        markersize=8,
        markeredgecolor=WHITE,
        markeredgewidth=2,
        linewidth=2,
    )
    for threshold in (50, 70):
        ax.axhline(threshold, color=SAMPLE_POINT_COLOR, linestyle="--", linewidth=1)
        ax.text(
            n_pcs * 0.8, threshold + 1.5, f"{threshold}%", fontsize=9, color=TEXT_COLOR
        )
    for n_pc in dict.fromkeys([2, 10, 20, n_pcs]):
        if n_pc > n_pcs:
            continue
        ax.annotate(
            f"PC1-{n_pc}\n{cumvar[n_pc - 1] * 100:.1f}%",
            (n_pc, cumvar[n_pc - 1] * 100),
            textcoords="offset points",
            xytext=(-8, -30) if n_pc == n_pcs else (8, -25),
            ha="right" if n_pc == n_pcs else "left",
            fontsize=8,
            color=TEXT_COLOR,
            arrowprops={"arrowstyle": "-", "color": SAMPLE_POINT_COLOR, "lw": 0.6},
        )
    ax.set_xlabel("Number of PCs", fontsize=10)
    ax.set_ylabel("Cumulative variance (%)", fontsize=10)
    ax.set_title("B. Cumulative variance", fontsize=12, fontweight="bold", loc="left")
    ax.set_xlim(0.5, n_pcs + 0.5)
    ax.set_ylim(0, min(cumvar[-1] * 100 + 15, 105))
    style_axes(ax)


def plot_pca_scree(
    adata: ad.AnnData,
    out_dir: Path,
    *,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """Plot PCA variance without changing the stored component statistics.

    Requires ``adata.uns["pca_variance_ratio"]``, set by the embedding stage.
    """
    var_ratio = adata.uns.get("pca_variance_ratio")
    if var_ratio is None:
        logger.warning("No pca_variance_ratio in adata.uns, skipping scree plot")
        return

    var_ratio = np.asarray(var_ratio, dtype=float)
    n_pcs = len(var_ratio)
    cumvar = np.cumsum(var_ratio)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    _plot_scree_panel(axes[0], var_ratio, n_pcs)
    _plot_cumulative_panel(axes[1], cumvar, n_pcs)
    fig.tight_layout()
    fig.savefig(out_dir / "pca_scree.png", dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(out_dir / "pca_scree.pdf", bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved pca_scree.png")
