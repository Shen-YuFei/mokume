"""
Tissue marker analysis plots: heatmap, marker t-SNE, dotplot.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
import scanpy as sc

logger = logging.getLogger(__name__)

# Custom colormaps
_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "custom_heatmap", ["#2166AC", "#F7F7F7", "#B2182B"]
)
_EXPR_CMAP = LinearSegmentedColormap.from_list(
    "expr", ["#EEEEEE", "#FDD835", "#F4511E", "#B71C1C"]
)



def compute_markers(adata: ad.AnnData, min_group_size: int = 2) -> ad.AnnData:
    """Run Wilcoxon rank-sum test for tissue markers.

    Adds ``tissue_markers`` key to ``adata.uns``.
    Skips gracefully when too many tissues have < *min_group_size* samples
    (Wilcoxon requires ≥ 2 per group).
    """

    # Check group sizes — Wilcoxon needs ≥ 2 samples per group
    tissue_counts = adata.obs["tissue"].value_counts()
    valid = tissue_counts[tissue_counts >= min_group_size]
    if len(valid) < 2:
        n_single = (tissue_counts < min_group_size).sum()
        logger.warning(
            "Wilcoxon markers skipped: %d / %d tissues have < %d samples",
            n_single, len(tissue_counts), min_group_size,
        )
        return adata

    # Subset to tissues with enough samples
    if len(valid) < len(tissue_counts):
        dropped = tissue_counts[tissue_counts < min_group_size]
        logger.info(
            "Wilcoxon: using %d / %d tissues (dropped %d with < %d samples)",
            len(valid), len(tissue_counts), len(dropped), min_group_size,
        )
        mask = adata.obs["tissue"].isin(valid.index)
        adata_sub = adata[mask].copy()
    else:
        adata_sub = adata

    # Ensure X has no NaN for rank_genes_groups
    x_backup = adata_sub.X.copy()
    adata_sub.X = np.nan_to_num(adata_sub.X, nan=0.0).astype(np.float32)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        warnings.simplefilter("ignore", PerformanceWarning)
        sc.tl.rank_genes_groups(
            adata_sub, groupby="tissue", method="wilcoxon",
            key_added="tissue_markers", use_raw=False,
        )
    adata_sub.X = x_backup

    # Copy results back to original adata
    adata.uns["tissue_markers"] = adata_sub.uns["tissue_markers"]
    logger.info("Wilcoxon rank_genes_groups done")
    return adata



def _collect_top_markers(
    adata: ad.AnnData, tissue_order: list[str], n_top: int,
) -> list[tuple[str, str]]:
    """Collect top *n_top* unique markers per tissue from Wilcoxon results."""
    heat_proteins: list[tuple[str, str]] = []
    seen: set[str] = set()
    for tissue in tissue_order:
        try:
            result = sc.get.rank_genes_groups_df(
                adata, group=tissue, key="tissue_markers"
            )
        except (KeyError, ValueError):
            continue
        count = 0
        for _, row in result.iterrows():
            prot = row["names"]
            if prot not in seen and count < n_top:
                heat_proteins.append((prot, tissue))
                seen.add(prot)
                count += 1
    return heat_proteins


def _compute_zscore_matrix(
    tissue_means: pd.DataFrame,
    available_proteins: list[str],
    available_tissues: list[str],
) -> np.ndarray:
    """Compute row-wise z-scored matrix clipped to [-3, 3]."""
    heat_mat = tissue_means.loc[available_proteins, available_tissues].values.astype(float)
    row_mean = np.nanmean(heat_mat, axis=1, keepdims=True)
    row_std = np.nanstd(heat_mat, axis=1, keepdims=True)
    row_std[row_std == 0] = 1
    heat_z = (heat_mat - row_mean) / row_std
    heat_z = np.clip(heat_z, -3, 3)
    return np.nan_to_num(heat_z, nan=0.0)


def _draw_heatmap(
    ax, heat_z, available_proteins, available_tissues,
    protein_names, protein_tissues, tissue_colors, title,
) -> None:
    """Render the heatmap image with annotations."""
    im = ax.imshow(
        heat_z, aspect="auto", cmap=_HEATMAP_CMAP, vmin=-3, vmax=3,
        interpolation="nearest",
    )
    ax.set_yticks(range(len(available_proteins)))
    ax.set_yticklabels(available_proteins, fontsize=5)
    ax.set_xticks(range(len(available_tissues)))
    ax.set_xticklabels(available_tissues, rotation=60, ha="right", fontsize=7)

    for i, t in enumerate(available_tissues):
        color = tissue_colors.get(t, "#333")
        if isinstance(color, tuple):
            color = matplotlib.colors.rgb2hex(color[:3])
        ax.get_xticklabels()[i].set_color(color)
        ax.get_xticklabels()[i].set_fontweight("bold")

    for i, t in enumerate(available_tissues):
        c = tissue_colors.get(t, "#999")
        ax.add_patch(plt.Rectangle((i - 0.5, -1.5), 1, 1, color=c, clip_on=False))

    prev_tissue = protein_tissues[0] if protein_tissues else None
    for i, prot in enumerate(available_proteins):
        idx_in_heat = protein_names.index(prot)
        t = protein_tissues[idx_in_heat]
        if t != prev_tissue:
            ax.axhline(y=i - 0.5, color="#ddd", linewidth=0.5)
            prev_tissue = t

    cbar = plt.colorbar(im, ax=ax, shrink=0.3, aspect=20, pad=0.02)
    cbar.set_label("Z-score (mean tissue expression)", fontsize=9)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=20)


def _compute_tissue_means(
    adata: ad.AnnData, tissue_order: list[str],
) -> pd.DataFrame:
    """Per-tissue mean expression for each protein."""
    tissues = adata.obs["tissue"].values
    tissue_means = pd.DataFrame(index=adata.var.index)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for t in tissue_order:
            mask = tissues == t
            if mask.sum() > 0:
                tissue_means[t] = np.nanmean(adata.X[mask, :], axis=0)
    return tissue_means


def _prepare_heatmap_data(
    heat_proteins: list[tuple[str, str]],
    tissue_means: pd.DataFrame, tissue_order: list[str],
) -> tuple[list[str], list[str], list[str], list[str], np.ndarray] | None:
    """Filter markers to available proteins/tissues and compute z-scores.

    Returns None if no valid data remains.
    """
    protein_names, protein_tissues = zip(*heat_proteins)
    cols_set, idx_set = set(tissue_means.columns), set(tissue_means.index)
    available_tissues = [t for t in tissue_order if t in cols_set]
    available_proteins = [p for p in protein_names if p in idx_set]
    if not (available_proteins and available_tissues):
        return None
    heat_z = _compute_zscore_matrix(tissue_means, available_proteins, available_tissues)
    return list(protein_names), list(protein_tissues), available_tissues, available_proteins, heat_z


def plot_marker_heatmap(
    adata: ad.AnnData,
    tissue_order: list[str],
    out_dir: Path,
    *,
    n_top: int = 5,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """Top *n_top* Wilcoxon markers per tissue, z-scored heatmap."""
    if "tissue_markers" not in adata.uns:
        logger.warning("No tissue_markers in adata.uns, skipping heatmap")
        return

    tissue_colors = adata.uns.get("tissue_colors", {})
    tissue_means = _compute_tissue_means(adata, tissue_order)

    heat_proteins = _collect_top_markers(adata, tissue_order, n_top)
    if not heat_proteins:
        logger.warning("No markers found for heatmap")
        return

    prepared = _prepare_heatmap_data(heat_proteins, tissue_means, tissue_order)
    if prepared is None:
        return
    protein_names, protein_tissues, available_tissues, available_proteins, heat_z = prepared

    fig_height = max(12, len(available_proteins) * 0.14)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    heatmap_title = (
        f"Top {n_top} marker proteins per tissue (Wilcoxon, "
        f"{len(available_tissues)} tissues)\n"
        f"{adata.n_vars:,} proteins, {adata.n_obs} samples"
    )
    _draw_heatmap(
        ax, heat_z, available_proteins, available_tissues,
        protein_names, protein_tissues, tissue_colors, heatmap_title,
    )

    plt.tight_layout()
    fig.savefig(out_dir / "marker_heatmap.png", dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(out_dir / "marker_heatmap.pdf", bbox_inches="tight")
    plt.close()
    logger.info("Saved marker_heatmap.png")



def _select_showcase_markers(
    adata: ad.AnnData, tissue_order: list[str], n_showcase: int,
) -> list[tuple[str, str]]:
    """Pick one top marker per tissue (up to *n_showcase*)."""
    showcase: list[tuple[str, str]] = []
    seen_prots: set[str] = set()
    for t in tissue_order:
        try:
            result = sc.get.rank_genes_groups_df(
                adata, group=t, key="tissue_markers"
            )
        except (KeyError, ValueError):
            continue
        for _, row in result.iterrows():
            prot = row["names"]
            if prot not in seen_prots and prot in adata.var.index:
                showcase.append((prot, t))
                seen_prots.add(prot)
                break
        if len(showcase) >= n_showcase:
            break
    return showcase


def _render_marker_subplot(ax, tsne, adata, prot, tissue, var_index_map):
    """Render a single marker expression subplot."""
    prot_idx = var_index_map[prot]
    expr = np.nan_to_num(adata.X[:, prot_idx].copy(), nan=0.0)
    order = np.argsort(expr)
    vmin = float(np.percentile(expr, 5))
    vmax = float(np.percentile(expr, 95))
    scat = ax.scatter(
        tsne[order, 0], tsne[order, 1], c=expr[order],
        cmap=_EXPR_CMAP, s=18, alpha=0.85, edgecolors="none",
        vmin=vmin, vmax=vmax,
    )
    plt.colorbar(scat, ax=ax, shrink=0.6, aspect=15, pad=0.02)
    ts_val = ""
    if "max_ts" in adata.var.columns:
        ts_val = f" (TS={adata.var.loc[prot, 'max_ts']:.2f})"
    ax.set_title(f"{prot}\n{tissue} marker{ts_val}", fontsize=9, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def _create_subplot_grid(
    n_items: int, n_cols: int = 4,
) -> tuple[plt.Figure, list]:
    """Create a subplot grid and return (fig, flat_axes_list)."""
    n_rows = (n_items + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5.5 * n_rows))
    flat = np.asarray(axes).flatten().tolist()
    for ax in flat[n_items:]:
        ax.set_visible(False)
    return fig, flat


def plot_marker_tsne(
    adata: ad.AnnData,
    tissue_order: list[str],
    out_dir: Path,
    *,
    n_showcase: int = 8,
    dpi: int = 200,
    save_pdf: bool = True,
) -> None:
    """2x4 grid of t-SNE colored by showcase marker expression."""
    if "X_tsne" not in adata.obsm:
        logger.warning("Missing t-SNE, skipping marker_tsne")
        return
    if "tissue_markers" not in adata.uns:
        logger.warning("Missing markers, skipping marker_tsne")
        return

    showcase = _select_showcase_markers(adata, tissue_order, n_showcase)
    if not showcase:
        return

    tsne = adata.obsm["X_tsne"]
    fig, flat_axes = _create_subplot_grid(len(showcase))

    var_index_map = {name: i for i, name in enumerate(adata.var.index)}
    for idx, (prot, tissue) in enumerate(showcase):
        _render_marker_subplot(flat_axes[idx], tsne, adata, prot, tissue, var_index_map)

    fig.suptitle(
        "t-SNE — Top tissue marker expression",
        fontsize=14, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "marker_tsne.png", dpi=dpi, bbox_inches="tight")
    if save_pdf:
        fig.savefig(out_dir / "marker_tsne.pdf", bbox_inches="tight")
    plt.close()
    logger.info("Saved marker_tsne.png")



def _collect_dotplot_genes(
    adata: ad.AnnData, available: list[str], n_top: int,
) -> list[str]:
    """Collect top marker genes for dotplot (max 60)."""
    genes: list[str] = []
    seen: set[str] = set()
    max_total = n_top * len(available)
    for tissue in available:
        try:
            result = sc.get.rank_genes_groups_df(
                adata, group=tissue, key="tissue_markers"
            )
        except (KeyError, ValueError):
            continue
        count = 0
        for _, row in result.iterrows():
            g = row["names"]
            if g not in seen and g in adata.var.index:
                genes.append(g)
                seen.add(g)
                count += 1
                if count >= n_top:
                    break
            if len(genes) >= max_total:
                break
    return genes[:60]


def plot_dotplot(
    adata: ad.AnnData,
    tissue_order: list[str],
    out_dir: Path,
    *,
    n_top: int = 3,
    dpi: int = 150,
) -> None:
    """Scanpy dotplot of top markers per tissue."""
    if "tissue_markers" not in adata.uns:
        return

    available = [t for t in tissue_order if t in adata.obs["tissue"].values]
    dotplot_genes = _collect_dotplot_genes(adata, available, n_top)
    if not dotplot_genes:
        return

    x_backup = adata.X.copy()
    tissue_backup = adata.obs["tissue"].copy()
    try:
        adata.X = np.nan_to_num(adata.X, nan=0.0).astype(np.float32)
        adata.obs["tissue"] = pd.Categorical(
            adata.obs["tissue"], categories=available, ordered=True,
        )
        dp = sc.pl.dotplot(
            adata, var_names=dotplot_genes, groupby="tissue",
            standard_scale="var", show=False, return_fig=True,
            figsize=(20, 10),
        )
        dp.savefig(out_dir / "marker_dotplot.png", dpi=dpi, bbox_inches="tight")
        plt.close()
        logger.info("Saved marker_dotplot.png")
    except (ValueError, TypeError, RuntimeError, KeyError) as exc:
        logger.error("Dotplot failed: %s", exc)
    finally:
        adata.X = x_backup
        adata.obs["tissue"] = tissue_backup



def save_markers_csv(
    adata: ad.AnnData,
    tissue_order: list[str],
    out_dir: Path,
    *,
    n_top: int = 10,
) -> None:
    """Save top markers per tissue to CSV."""
    if "tissue_markers" not in adata.uns:
        return

    records: list[dict] = []
    for tissue in tissue_order:
        try:
            result = sc.get.rank_genes_groups_df(
                adata, group=tissue, key="tissue_markers"
            )
        except (KeyError, ValueError):
            continue
        top = result.head(n_top)
        for _, row in top.iterrows():
            records.append(
                {
                    "tissue": tissue,
                    "protein": row["names"],
                    "score": round(row["scores"], 2),
                    "logfoldchange": round(row["logfoldchanges"], 3),
                    "pval_adj": f"{row['pvals_adj']:.2e}",
                }
            )

    if records:
        df = pd.DataFrame(records)
        csv_path = out_dir / "tissue_markers_top10.csv"
        df.to_csv(csv_path, index=False)
        logger.info("Saved %s (%d markers)", csv_path.name, len(df))
