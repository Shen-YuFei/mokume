"""
AdaTiSS tissue specificity distribution and per-tissue bar chart.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm as _norm

if TYPE_CHECKING:
    from mokume.tissuemap.tissue_specificity import GMMThresholds

logger = logging.getLogger(__name__)


def plot_ts_distribution(
    ts_df: pd.DataFrame,
    unique_tissues: list[str],
    out_dir: Path,
    *,
    ts_enriched: float = 2.5,
    ts_specific: float = 4.0,
    ts_housekeeping: float | None = None,
    gmm: GMMThresholds | None = None,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """Two-panel figure: TS score histogram + per-tissue specific proteins.

    Parameters
    ----------
    ts_df : pd.DataFrame
        Output of :func:`compute_ts_scores` with ``max_ts`` and
        ``enrichment_category`` columns.
    unique_tissues : list[str]
        Tissue names (column order).
    out_dir : Path
        Directory to save plots.
    gmm : GMMThresholds, optional
        If provided, overlay GMM component curves on Panel A.
    """
    if "max_ts" not in ts_df.columns:
        logger.warning("No max_ts in ts_df, skipping TS distribution plot")
        return

    max_ts = ts_df["max_ts"].dropna()
    categories = ts_df.get("enrichment_category", pd.Series(dtype=str))

    n_specific = (categories == "tissue-specific").sum()
    n_enriched = (categories.isin(["tissue-specific", "tissue-enriched"])).sum()

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: TS score distribution
    ax = axes[0]
    # Auto x-axis limit: 99th percentile, at least ts_specific * 2
    p99 = float(max_ts.quantile(0.99)) if len(max_ts) > 0 else ts_specific * 3
    x_upper = max(p99, ts_specific * 2)
    n_outliers = int((max_ts > x_upper).sum())

    ax.hist(
        max_ts.clip(upper=x_upper).values, bins=50,
        color="#5C6BC0", edgecolor="white", linewidth=0.5, alpha=0.6,
        density=True,
    )
    ax.axvline(
        ts_specific, color="#E53935", linestyle="--", linewidth=1.5,
        label=f"Tissue-specific TS≥{ts_specific:.2f} (n={n_specific})",
    )
    ax.axvline(
        ts_enriched, color="#FF9800", linestyle="--", linewidth=1.5,
        label=f"Tissue-enriched TS≥{ts_enriched:.2f} (n={n_enriched})",
    )
    if ts_housekeeping is not None:
        n_hk = int((categories == "house-keeping").sum())
        ax.axvline(
            ts_housekeeping, color="#43A047", linestyle="-.", linewidth=1.5,
            label=f"House-keeping |TS|<{ts_housekeeping:.2f} (n={n_hk})",
        )
    # Overlay GMM fit curves if available
    if gmm is not None:
        x_grid = np.linspace(0, x_upper, 500)
        p_bg = gmm.bg_weight * _norm.pdf(x_grid, gmm.bg_mean, gmm.bg_std)
        p_sp = gmm.sp_weight * _norm.pdf(x_grid, gmm.sp_mean, gmm.sp_std)
        ax.plot(x_grid, p_bg, "--", color="#1565C0", lw=1.5,
                label=f"Background (\u03bc={gmm.bg_mean:.2f})")
        ax.plot(x_grid, p_sp, "--", color="#C62828", lw=1.5,
                label=f"Specific (\u03bc={gmm.sp_mean:.2f})")
        ax.plot(x_grid, p_bg + p_sp, "-", color="#222", lw=2,
                label="GMM mixture", alpha=0.7)

    ax.set_xlim(left=max_ts.min() - 0.5, right=x_upper * 1.05)
    if n_outliers > 0:
        ax.text(
            0.97, 0.95, f"{n_outliers} proteins with TS > {x_upper:.1f}\n(clipped)",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="#777", style="italic",
        )
    ax.set_xlabel("AdaTiSS max TS score", fontsize=11)
    ax.set_ylabel("Density" if gmm is not None else "Number of proteins", fontsize=11)
    ax.set_title(
        "A. AdaTiSS tissue specificity score distribution",
        fontsize=12, fontweight="bold", loc="left",
    )
    ax.legend(fontsize=9, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: Number of tissue-specific proteins per tissue
    ax = axes[1]
    specific_mask = categories == "tissue-specific"
    if specific_mask.any():
        specific_per_tissue = (
            ts_df.loc[specific_mask, "max_tissue"]
            .value_counts()
            .reindex(unique_tissues)
            .fillna(0)
            .astype(int)
        )
    else:
        specific_per_tissue = pd.Series(0, index=unique_tissues)

    # Reverse for horizontal barplot (top-down reading)
    spt = specific_per_tissue.iloc[::-1]

    colors = ["#5C6BC0"] * len(spt)

    ax.barh(
        range(len(spt)), spt.values,
        color=colors, edgecolor="white", linewidth=0.5,
    )
    ax.set_yticks(range(len(spt)))
    ax.set_yticklabels(spt.index, fontsize=7)
    ax.set_xlabel(
        f"Number of tissue-specific proteins (TS ≥ {ts_specific:.2f})", fontsize=9,
    )
    ax.set_title(
        "B. Tissue-specific proteins per tissue",
        fontsize=12, fontweight="bold", loc="left",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Count labels
    for i, v in enumerate(spt.values):
        if v > 0:
            ax.text(v + 0.5, i, str(v), va="center", fontsize=6, color="#555")

    plt.tight_layout()
    fig.savefig(out_dir / "ts_distribution.png", dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(out_dir / "ts_distribution.pdf", bbox_inches="tight")
    plt.close()
    logger.info("Saved ts_distribution.png")


def plot_specific_per_tissue(
    ts_df: pd.DataFrame,
    unique_tissues: list[str],
    out_dir: Path,
    *,
    ts_specific: float = 4.0,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """Standalone per-tissue bar chart of tissue-specific protein counts.

    This is a simpler standalone version of Panel B above, useful when
    the combined figure is not needed.
    """
    categories = ts_df.get("enrichment_category", pd.Series(dtype=str))
    specific_mask = categories == "tissue-specific"

    if specific_mask.any():
        specific_per_tissue = (
            ts_df.loc[specific_mask, "max_tissue"]
            .value_counts()
            .reindex(unique_tissues)
            .fillna(0)
            .astype(int)
        )
    else:
        specific_per_tissue = pd.Series(0, index=unique_tissues)

    fig, ax = plt.subplots(figsize=(10, max(6, len(unique_tissues) * 0.3)))

    spt = specific_per_tissue.iloc[::-1]
    ax.barh(
        range(len(spt)), spt.values,
        color="#5C6BC0", edgecolor="white", linewidth=0.5,
    )
    ax.set_yticks(range(len(spt)))
    ax.set_yticklabels(spt.index, fontsize=8)
    ax.set_xlabel(
        f"Number of tissue-specific proteins (TS ≥ {ts_specific:.2f})", fontsize=10,
    )
    ax.set_title(
        "Tissue-specific proteins per tissue (AdaTiSS)",
        fontsize=13, fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, v in enumerate(spt.values):
        if v > 0:
            ax.text(v + 0.5, i, str(v), va="center", fontsize=7, color="#555")

    plt.tight_layout()
    fig.savefig(
        out_dir / "specific_per_tissue.png", dpi=dpi, bbox_inches="tight",
    )
    if save_pdf:
        fig.savefig(
            out_dir / "specific_per_tissue.pdf", bbox_inches="tight",
        )
    plt.close()
    logger.info("Saved specific_per_tissue.png")
