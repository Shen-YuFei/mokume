"""
Tissue atlas slide figure: t-SNE + dendrogram in one 16:9 panel.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.cluster.hierarchy import dendrogram as scipy_dendro
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

logger = logging.getLogger(__name__)

_SMALL_THRESHOLD = 8

_GROUP_PALETTE = {
    "Muscle/Cardiac": {"base": "#D32F2F", "cmap": matplotlib.colormaps["Reds"]},
    "Neural": {"base": "#1565C0", "cmap": matplotlib.colormaps["Blues"]},
    "Hepato-renal": {"base": "#8D6E63", "cmap": matplotlib.colormaps["copper"]},
    "Immune/Gut-associated": {"base": "#2E7D32", "cmap": matplotlib.colormaps["Greens"]},
    "Gonadal/Reproductive": {"base": "#7B1FA2", "cmap": matplotlib.colormaps["Purples"]},
    "Oropharyngeal": {"base": "#00838F", "cmap": matplotlib.colormaps["GnBu"]},
    "Pituitary": {"base": "#F57F17", "cmap": matplotlib.colormaps["YlOrBr"]},
    "Mucosal/Epithelial": {"base": "#E65100", "cmap": matplotlib.colormaps["Oranges"]},
    "Secretory/Glandular": {"base": "#558B2F", "cmap": matplotlib.colormaps["YlGn"]},
    "Stromal/Connective": {"base": "#78909C", "cmap": matplotlib.colormaps["cool"]},
}

_GROUP_HEURISTICS: list[tuple[set[str], str]] = [
    ({"heart", "muscle"}, "Muscle/Cardiac"),
    ({"brain", "eye", "spinal cord"}, "Neural"),
    ({"kidney", "liver"}, "Hepato-renal"),
    ({"intestine", "lymph node", "spleen"}, "Immune/Gut-associated"),
    ({"ovary", "testis", "uterus"}, "Gonadal/Reproductive"),
    ({"salivary gland", "tonsil"}, "Oropharyngeal"),
    ({"pituitary gland"}, "Pituitary"),
]


def _despine(ax: matplotlib.axes.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def _draw_convex_hull(ax, points, color, label=None, pad=2.5):
    if len(points) < 3:
        return
    try:
        hull = ConvexHull(points)
    except (ValueError, RuntimeError):
        return
    hull_pts = points[hull.vertices]
    hull_pts = np.vstack([hull_pts, hull_pts[0]])
    centroid = points.mean(axis=0)
    norms = np.linalg.norm(hull_pts - centroid, axis=1, keepdims=True).clip(1e-6)
    expanded = centroid + (hull_pts - centroid) * (1 + pad / norms)
    ax.plot(
        expanded[:, 0], expanded[:, 1],
        color=color, linestyle="--", linewidth=1.8, alpha=0.7, zorder=1,
    )
    if label:
        top_idx = np.argmax(expanded[:, 1])
        ax.text(
            expanded[top_idx, 0], expanded[top_idx, 1] + 1.5, label,
            fontsize=7, fontweight="bold", color=color, ha="center", va="bottom",
            bbox={
                "boxstyle": "round,pad=0.2", "facecolor": "white",
                "edgecolor": color, "alpha": 0.85, "linewidth": 0.8,
            },
        )


def _draw_ellipse(ax, points, color, label=None):
    if len(points) < 2:
        return
    cx, cy = points.mean(axis=0)
    sx = max(points[:, 0].std() * 2.5, 3)
    sy = max(points[:, 1].std() * 2.5, 3)
    ellipse = Ellipse(
        (cx, cy), width=sx * 2, height=sy * 2,
        fill=False, edgecolor=color, linestyle=":",
        linewidth=1.5, alpha=0.7, zorder=1,
    )
    ax.add_patch(ellipse)
    if label:
        ax.text(
            cx, cy + sy + 1.5, label,
            fontsize=6.5, fontweight="bold", color=color, ha="center", va="bottom",
            bbox={
                "boxstyle": "round,pad=0.15", "facecolor": "white",
                "edgecolor": color, "alpha": 0.8, "linewidth": 0.6,
            },
        )


def _name_cluster(
    members: list[str], cluster_id: int, used_names: set[str],
) -> str:
    """Assign a human-readable name to a cluster using heuristics."""
    member_set = set(members)
    # Exact subset match first
    for pattern, pname in _GROUP_HEURISTICS:
        if pattern <= member_set and pname not in used_names:
            return pname
    # Partial overlap fallback
    for pattern, pname in _GROUP_HEURISTICS:
        if pattern & member_set and pname not in used_names:
            return pname
    return f"Group {cluster_id}"


def _compute_tissue_groups(
    adata: ad.AnnData,
    max_clusters: int = 10,
) -> tuple[dict[str, list[str]], dict[str, str], np.ndarray]:
    """Compute data-driven tissue groups from hierarchical linkage.

    Returns
    -------
    proteomic_groups : dict[str, list[str]]
    tissue_to_group : dict[str, str]
    z_linkage : np.ndarray
        Ward linkage matrix (reused by the plotting function).
    """
    tissues = adata.obs["tissue"].values
    unique_tissues = sorted(np.unique(tissues))
    tissue_means = np.array(
        [adata.X[tissues == t].mean(axis=0) for t in unique_tissues]
    )

    tissue_means = np.nan_to_num(tissue_means, nan=0.0)
    dist = pdist(tissue_means, metric="correlation")
    finite_mask = np.isfinite(dist)
    if not finite_mask.all():
        dist[~finite_mask] = np.nanmax(dist[finite_mask]) if finite_mask.any() else 1.0

    z_linkage = linkage(dist, method="ward")
    n_groups = min(max_clusters, len(unique_tissues))
    cut_labels = fcluster(z_linkage, t=n_groups, criterion="maxclust")

    raw_groups: dict[int, list[str]] = {}
    for i, tissue in enumerate(unique_tissues):
        raw_groups.setdefault(cut_labels[i], []).append(tissue)

    proteomic_groups: dict[str, list[str]] = {}
    used_names: set[str] = set()
    for cl, members in raw_groups.items():
        name = _name_cluster(members, cl, used_names)
        used_names.add(name)
        proteomic_groups[name] = members

    tissue_to_group: dict[str, str] = {
        tissue: group_name
        for group_name, members in proteomic_groups.items()
        for tissue in members
    }
    return proteomic_groups, tissue_to_group, z_linkage


def _build_tissue_colors(
    proteomic_groups: dict[str, list[str]],
) -> dict[str, tuple]:
    """Assign per-tissue colors from group colormaps."""
    tissue_colors: dict[str, tuple] = {}
    for group_name, members in proteomic_groups.items():
        palette = _GROUP_PALETTE.get(
            group_name, {"base": "#999", "cmap": matplotlib.colormaps["Greys"]}
        )
        cmap = palette["cmap"]
        n = len(members)
        for i, t in enumerate(members):
            frac = 0.35 + 0.55 * (i / max(n - 1, 1))
            tissue_colors[t] = cmap(frac)
    return tissue_colors


def _draw_tsne_panel(
    ax, tsne_emb, tissues, tissue_order, tissue_colors, proteomic_groups,
) -> None:
    """Scatter t-SNE points and draw group hulls / ellipses."""
    for t in tissue_order:
        mask = tissues == t
        ax.scatter(
            tsne_emb[mask, 0], tsne_emb[mask, 1],
            c=[tissue_colors.get(t, "#999")], s=45, alpha=0.88,
            edgecolors="white", linewidths=0.4, zorder=2,
        )

    for group_name, members in proteomic_groups.items():
        present = [t for t in members if (tissues == t).sum() > 0]
        if not present:
            continue
        group_mask = np.isin(tissues, present)
        group_pts = tsne_emb[group_mask]
        color = _GROUP_PALETTE.get(group_name, {"base": "#999"})["base"]
        if group_mask.sum() >= _SMALL_THRESHOLD:
            _draw_convex_hull(ax, group_pts, color, label=group_name, pad=2.5)
        else:
            _draw_ellipse(ax, group_pts, color, label=group_name)


def _build_legend_elements(
    tissues, proteomic_groups, tissue_colors,
) -> list[Line2D]:
    """Build legend handles for tissue groups."""
    elements: list[Line2D] = []
    for group_name, members in proteomic_groups.items():
        elements.append(
            Line2D(
                [0], [0], marker="none", label=f"  {group_name}",
                color="none", markerfacecolor="none",
            )
        )
        for t in members:
            n = int((tissues == t).sum())
            if n > 0:
                elements.append(
                    Line2D(
                        [0], [0], marker="o", color="none",
                        markerfacecolor=tissue_colors.get(t, "#999"),
                        markeredgecolor="white", markeredgewidth=0.3,
                        markersize=8, label=f"    {t} ({n})",
                    )
                )
    return elements


def _draw_dendrogram_panel(
    ax, z_linkage, unique_tissues, tissue_to_group,
) -> None:
    """Render dendrogram and style labels by group color."""
    scipy_dendro(
        z_linkage, labels=unique_tissues, ax=ax,
        orientation="right", leaf_font_size=8,
        above_threshold_color="#BDBDBD", color_threshold=0,
    )
    for coll in ax.collections:
        coll.set_linewidth(0.8)
    for line in ax.lines:
        line.set_linewidth(0.8)
    for lbl in ax.get_yticklabels():
        t = lbl.get_text()
        grp = tissue_to_group.get(t, "")
        color = _GROUP_PALETTE.get(grp, {"base": "#333"})["base"]
        lbl.set_color(color)
        lbl.set_fontweight("bold")
        lbl.set_fontsize(8)
    ax.set_title("B  Tissue Dendrogram", fontsize=16, fontweight="bold", loc="left", pad=8)
    ax.set_xlabel("Ward distance (correlation)", fontsize=9, labelpad=5)
    ax.tick_params(axis="x", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)


def _save_atlas_figures(
    fig, ax_tsne, ax_dend, out_dir: Path, dpi: int, save_pdf: bool,
) -> None:
    """Save combined figure and individual sub-plots."""
    fig.savefig(
        out_dir / "slide_atlas_dendrogram.png",
        dpi=dpi, bbox_inches="tight", facecolor="white",
    )
    if save_pdf:
        fig.savefig(
            out_dir / "slide_atlas_dendrogram.pdf",
            bbox_inches="tight", facecolor="white",
        )
    for ax, name in [(ax_tsne, "tissue_atlas"), (ax_dend, "tissue_dendrogram")]:
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(
            fig.dpi_scale_trans.inverted(),
        )
        fig.savefig(
            out_dir / f"{name}.png",
            dpi=dpi, bbox_inches=extent, facecolor="white",
        )
        if save_pdf:
            fig.savefig(
                out_dir / f"{name}.pdf",
                bbox_inches=extent, facecolor="white",
            )
    plt.close()
    logger.info("Saved slide_atlas_dendrogram.png")


def plot_slide_atlas_dendrogram(
    adata: ad.AnnData,
    out_dir: Path,
    *,
    dpi: int = 300,
    save_pdf: bool = True,
) -> None:
    """Plot the combined tissue atlas and dendrogram figure."""
    if "X_tsne" not in adata.obsm:
        logger.warning("No t-SNE embedding, skipping slide figure")
        return

    tsne_emb = adata.obsm["X_tsne"]
    tissues = adata.obs["tissue"].values
    unique_tissues = sorted(np.unique(tissues))

    proteomic_groups, tissue_to_group, z_linkage = _compute_tissue_groups(adata)
    tissue_colors = _build_tissue_colors(proteomic_groups)

    adata.uns["proteomic_groups"] = proteomic_groups
    adata.uns["tissue_colors"] = {
        t: matplotlib.colors.rgb2hex(c[:3]) if isinstance(c, tuple) else c
        for t, c in tissue_colors.items()
    }

    tissue_order = [
        t for members in proteomic_groups.values()
        for t in members if (tissues == t).sum() > 0
    ]

    n_legend_entries = sum(
        1 + sum(1 for t in members if (tissues == t).sum() > 0)
        for members in proteomic_groups.values()
    )
    fig_width = max(20, 14 + n_legend_entries * 0.35)

    fig = plt.figure(figsize=(fig_width, 11))
    gs = GridSpec(
        1, 2, width_ratios=[1.4, 0.8], wspace=0.01,
        left=0.01, right=0.99, bottom=0.05, top=0.90,
    )

    ax_tsne = fig.add_subplot(gs[0])
    _draw_tsne_panel(ax_tsne, tsne_emb, tissues, tissue_order, tissue_colors, proteomic_groups)

    legend_elements = _build_legend_elements(tissues, proteomic_groups, tissue_colors)
    leg = ax_tsne.legend(
        handles=legend_elements, loc="center left", bbox_to_anchor=(1.01, 0.5),
        fontsize=7, frameon=True, framealpha=0.95, edgecolor="#ddd",
        handletextpad=0.5, labelspacing=0.2, borderpad=0.6,
        title="Proteomic group", title_fontsize=8,
    )
    leg.set_alignment("left")

    ax_tsne.set_title("A  Tissue Atlas", fontsize=16, fontweight="bold", loc="left", pad=8)
    ax_tsne.set_aspect("equal")
    _despine(ax_tsne)
    ax_tsne.set_xlabel("t-SNE 1", fontsize=10, color="#888", labelpad=5)
    ax_tsne.set_ylabel("t-SNE 2", fontsize=10, color="#888", labelpad=5)

    metrics = adata.uns.get("embedding_metrics", {})
    pca_var = metrics.get("pca_var_explained")
    pca_str = f"{pca_var:.3f}" if pca_var is not None else "N/A"
    stats = (
        f"{adata.n_vars:,} proteins | {adata.n_obs} samples"
        f" | {adata.obs['tissue'].nunique()} tissues"
        f" | PCA variance = {pca_str}"
    )
    ax_tsne.text(
        0.01, 0.01, stats, transform=ax_tsne.transAxes,
        fontsize=8.5, fontfamily="monospace", verticalalignment="bottom",
        bbox={
            "boxstyle": "round,pad=0.5", "facecolor": "#F5F5F5",
            "edgecolor": "#BDBDBD", "alpha": 0.95, "linewidth": 1,
        },
    )

    ax_dend = fig.add_subplot(gs[1])
    _draw_dendrogram_panel(ax_dend, z_linkage, unique_tissues, tissue_to_group)

    _save_atlas_figures(fig, ax_tsne, ax_dend, out_dir, dpi, save_pdf)
