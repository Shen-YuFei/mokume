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


def _overlay_gmm_curves(ax, gmm, x_upper: float) -> None:
    """Draw GMM component density curves on an axes."""
    x_grid = np.linspace(0, x_upper, 500)
    p_bg = gmm.bg_weight * _norm.pdf(x_grid, gmm.bg_mean, gmm.bg_std)
    p_sp = gmm.sp_weight * _norm.pdf(x_grid, gmm.sp_mean, gmm.sp_std)
    ax.plot(
        x_grid,
        p_bg,
        "--",
        color="#1565C0",
        lw=1.5,
        label=f"Background (\u03bc={gmm.bg_mean:.2f})",
    )
    ax.plot(
        x_grid,
        p_sp,
        "--",
        color="#C62828",
        lw=1.5,
        label=f"Specific (\u03bc={gmm.sp_mean:.2f})",
    )
    ax.plot(
        x_grid, p_bg + p_sp, "-", color="#222", lw=2, label="GMM mixture", alpha=0.7
    )


def _plot_histogram_panel(
    ax,
    max_ts,
    categories,
    ts_enriched,
    ts_specific,
    ts_housekeeping,
    gmm,
) -> None:
    """Panel A: TS score histogram with threshold lines and optional GMM."""
    n_specific = int((categories == "tissue-specific").sum())
    n_enriched = int(categories.isin(["tissue-specific", "tissue-enriched"]).sum())

    p99 = float(max_ts.quantile(0.99)) if len(max_ts) > 0 else ts_specific * 3
    x_upper = max(p99, ts_specific * 2)
    n_outliers = int((max_ts > x_upper).sum())

    use_density = gmm is not None
    ax.hist(
        max_ts.clip(upper=x_upper).values,
        bins=50,
        color="#5C6BC0",
        edgecolor="white",
        linewidth=0.5,
        alpha=0.6,
        density=use_density,
    )
    ax.axvline(
        ts_specific,
        color="#E53935",
        linestyle="--",
        linewidth=1.5,
        label=f"Tissue-specific TS≥{ts_specific:.2f} (n={n_specific})",
    )
    ax.axvline(
        ts_enriched,
        color="#FF9800",
        linestyle="--",
        linewidth=1.5,
        label=f"Tissue-enriched TS≥{ts_enriched:.2f} (n={n_enriched})",
    )
    if ts_housekeeping is not None:
        n_hk = int((categories == "house-keeping").sum())
        ax.axvline(
            ts_housekeeping,
            color="#43A047",
            linestyle="-.",
            linewidth=1.5,
            label=f"House-keeping |TS|<{ts_housekeeping:.2f} (n={n_hk})",
        )
    if gmm is not None:
        _overlay_gmm_curves(ax, gmm, x_upper)

    left_bound = float(max_ts.min()) - 0.5 if len(max_ts) > 0 else 0.0
    ax.set_xlim(left=left_bound, right=x_upper * 1.05)
    if n_outliers > 0:
        ax.text(
            0.97,
            0.95,
            f"{n_outliers} proteins with TS > {x_upper:.1f}\n(clipped)",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color="#777",
            style="italic",
        )
    ax.set_xlabel("AdaTiSS max TS score", fontsize=11)
    ax.set_ylabel("Density" if use_density else "Number of proteins", fontsize=11)
    ax.set_title(
        "A. AdaTiSS tissue specificity score distribution",
        fontsize=12,
        fontweight="bold",
        loc="left",
    )
    ax.legend(fontsize=9, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _compute_specific_counts(
    ts_df: pd.DataFrame,
    unique_tissues: list[str],
) -> pd.Series:
    """Count tissue-specific proteins per tissue."""
    categories = ts_df.get("enrichment_category", pd.Series(dtype=str))
    specific_mask = categories == "tissue-specific"
    if specific_mask.any():
        return (
            ts_df.loc[specific_mask, "max_tissue"]
            .value_counts()
            .reindex(unique_tissues)
            .fillna(0)
            .astype(int)
        )
    return pd.Series(0, index=unique_tissues)


def _plot_bar_panel(
    ax,
    specific_per_tissue: pd.Series,
    ts_specific: float,
) -> None:
    """Panel B: horizontal bar chart of tissue-specific protein counts."""
    spt = specific_per_tissue.iloc[::-1]
    ax.barh(
        range(len(spt)),
        spt.values,
        color="#5C6BC0",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(range(len(spt)))
    ax.set_yticklabels(spt.index, fontsize=7)
    ax.set_xlabel(
        f"Number of tissue-specific proteins (TS ≥ {ts_specific:.2f})",
        fontsize=9,
    )
    ax.set_title(
        "B. Tissue-specific proteins per tissue",
        fontsize=12,
        fontweight="bold",
        loc="left",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for i, v in enumerate(spt.values):
        if v > 0:
            ax.text(v + 0.5, i, str(v), va="center", fontsize=6, color="#555")


def plot_ts_distribution(
    ts_df: pd.DataFrame,
    unique_tissues: list[str],
    out_dir: Path,
    *,
    ts_enriched: float = 2.5,
    ts_specific: float = 4.0,
    ts_housekeeping: float | None = None,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """Two-panel figure: TS score histogram + per-tissue specific proteins.

    Parameters
    ----------
    ts_df : pd.DataFrame
        Output of :func:`compute_ts_scores` with ``max_ts`` and
        ``enrichment_category`` columns.  GMM fit is read from
        ``ts_df.attrs["gmm"]`` if available.
    unique_tissues : list[str]
        Tissue names (column order).
    out_dir : Path
        Directory to save plots.
    """
    if "max_ts" not in ts_df.columns:
        logger.warning("No max_ts in ts_df, skipping TS distribution plot")
        return

    max_ts = ts_df["max_ts"].dropna()
    categories = ts_df.get("enrichment_category", pd.Series(dtype=str))
    gmm: GMMThresholds | None = ts_df.attrs.get("gmm")

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    _plot_histogram_panel(
        axes[0],
        max_ts,
        categories,
        ts_enriched,
        ts_specific,
        ts_housekeeping,
        gmm,
    )
    _plot_bar_panel(
        axes[1],
        _compute_specific_counts(ts_df, unique_tissues),
        ts_specific,
    )

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
    specific_per_tissue = _compute_specific_counts(ts_df, unique_tissues)

    fig, ax = plt.subplots(figsize=(10, max(6, len(unique_tissues) * 0.3)))

    spt = specific_per_tissue.iloc[::-1]
    ax.barh(
        range(len(spt)),
        spt.values,
        color="#5C6BC0",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(range(len(spt)))
    ax.set_yticklabels(spt.index, fontsize=8)
    ax.set_xlabel(
        f"Number of tissue-specific proteins (TS ≥ {ts_specific:.2f})",
        fontsize=10,
    )
    ax.set_title(
        "Tissue-specific proteins per tissue (AdaTiSS)",
        fontsize=13,
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, v in enumerate(spt.values):
        if v > 0:
            ax.text(v + 0.5, i, str(v), va="center", fontsize=7, color="#555")

    plt.tight_layout()
    fig.savefig(
        out_dir / "specific_per_tissue.png",
        dpi=dpi,
        bbox_inches="tight",
    )
    if save_pdf:
        fig.savefig(
            out_dir / "specific_per_tissue.pdf",
            bbox_inches="tight",
        )
    plt.close()
    logger.info("Saved specific_per_tissue.png")
