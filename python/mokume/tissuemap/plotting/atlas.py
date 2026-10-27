"""Tissue atlas and dendrograms sharing proteomic-group color/shape codes."""

from __future__ import annotations

from collections import Counter
import logging
from pathlib import Path
from typing import NamedTuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox, ScaledTranslation
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.cluster.hierarchy import dendrogram as scipy_dendro
from scipy.spatial.distance import pdist

from mokume.plotting._common import panel_label, static_style, style_axes
from mokume.plotting.style import (
    CONDITION_COLORS,
    CONDITION_MARKERS,
    CONDITION_STYLE_INDICES,
    SAMPLE_POINT_COLOR,
    TEXT_COLOR,
    WHITE,
)

logger = logging.getLogger(__name__)


class _SaveOpts(NamedTuple):
    """Bundle the figure-output options shared by every atlas saver."""

    dpi: int
    save_pdf: bool


def _name_cluster(members: list[str], cluster_id: int, used_names: set[str]) -> str:
    """Name a cluster only when one organ family has a strict majority."""
    family, count = Counter(_organ_key(tissue) for tissue in members).most_common(1)[0]
    name = family.capitalize()
    if count > len(members) / 2 and name not in used_names:
        return name
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
    for cl, members in sorted(raw_groups.items()):
        name = _name_cluster(members, cl, used_names)
        used_names.add(name)
        proteomic_groups[name] = members

    tissue_to_group: dict[str, str] = {
        tissue: group_name
        for group_name, members in proteomic_groups.items()
        for tissue in members
    }
    return proteomic_groups, tissue_to_group, z_linkage


# Keyword order matters: specific organ families precede broader matches.
_ORGAN_ALIASES: list[tuple[str, str]] = [
    ("lymph node", "lymph node"),
    ("bone marrow", "bone marrow"),
    ("bone", "bone"),
    ("brain", "brain"),
    ("oral cavity", "oral cavity"),
    ("mouth", "oral cavity"),
    ("tonsil", "oral cavity"),
    ("lung", "lung"),
    ("pleura", "lung"),
    ("salivary gland", "salivary gland"),
    ("small intestine", "intestine"),
    ("large intestine", "intestine"),
    ("ileocecal", "intestine"),
    ("cecum", "intestine"),
    ("intestine", "intestine"),
    ("rectum", "rectum"),
    ("esophag", "esophagus"),
    ("ovar", "ovary"),
    ("fallopian", "uterus"),
    ("endometri", "uterus"),
    ("uter", "uterus"),
    ("cervix", "uterus"),
    ("vulva", "vulva"),
    ("vagina", "vagina"),
    ("renal pelvis", "kidney"),
    ("pararenal", "kidney"),
    ("kidney", "kidney"),
    ("liver", "liver"),
    ("pancrea", "pancreas"),
    ("colon", "colon"),
    ("pyloric", "stomach"),
    ("paragastric", "stomach"),
    ("gastric", "stomach"),
    ("stomach", "stomach"),
    ("breast", "breast"),
    ("prostate", "prostate"),
    ("thyroid", "thyroid"),
    ("bladder", "bladder"),
    ("muscle", "muscle"),
    ("skin", "skin"),
    ("epidermis", "skin"),
    ("hypoderm", "skin"),
    ("pharynx", "pharynx"),
    ("spleen", "spleen"),
    ("adrenal", "adrenal gland"),
    ("testis", "testis"),
    ("placenta", "placenta"),
    ("omentum", "omentum"),
    ("retroperiton", "peritoneum"),
    ("peritone", "peritoneum"),
    ("pericard", "heart"),
    ("heart", "heart"),
    ("cerebrospinal", "nervous system"),
    ("spinal", "nervous system"),
    ("nervous", "nervous system"),
    ("orbital", "eye"),
    ("eye", "eye"),
    ("synovi", "joint"),
    ("cartil", "joint"),
    ("head and neck", "head and neck"),
    ("neck", "head and neck"),
    ("larynx", "larynx"),
    ("blood", "blood"),
    ("ascites", "ascites"),
]


def _organ_key(tissue: str) -> str:
    """Map a tissue label to its organ family for majority-based cluster naming.

    Falls back to the text before the first comma (e.g. "colon, sigmoid" ->
    "colon") when no keyword matches.
    """
    low = tissue.lower()
    for keyword, canonical in _ORGAN_ALIASES:
        if keyword in low:
            return canonical
    return tissue.split(",")[0].strip()


def _group_styles(proteomic_groups: dict) -> dict[str, tuple[str, str]]:
    """Use circles for up to six groups, then unique color/shape pairs."""
    styles = {}
    for index, group in enumerate(proteomic_groups):
        color, marker = CONDITION_STYLE_INDICES[index]
        styles[group] = (
            CONDITION_COLORS[color],
            CONDITION_MARKERS[marker] if len(proteomic_groups) > 6 else "o",
        )
    return styles


def _label_position(size, center, placed, bounds, others):
    """Choose a nearby label box inside the axes, avoiding earlier labels."""
    width, height = size
    candidates = []
    for radius in (12, 28, 48, 72, 104, 144):
        for angle in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            point = center + (np.cos(angle), np.sin(angle)) * np.array(
                [radius + width / 2, radius + height / 2]
            )
            point = np.clip(
                point,
                [bounds.x0 + width / 2, bounds.y0 + height / 2],
                [bounds.x1 - width / 2, bounds.y1 - height / 2],
            )
            box = Bbox.from_bounds(
                point[0] - width / 2, point[1] - height / 2, width, height
            )
            distance = np.linalg.norm(point - center)
            score = sum(box.overlaps(other) for other in placed) * 1e6
            # A label nearer another group's anchor would read as that group's.
            score += 1e3 * (
                np.linalg.norm(others - point, axis=1).min(initial=np.inf) < distance
            )
            # Prefer above, then beside, so every label reads the same way.
            score += distance + 24 * (1 - np.sin(angle))
            candidates.append((score, point, box))
    _, point, box = min(candidates, key=lambda candidate: candidate[0])
    return point, box


def _label_groups(ax, anchors: dict) -> None:
    """Place neutral group labels near their anchors without crossing the axes."""
    labels = [
        ax.text(
            *center,
            group,
            color=TEXT_COLOR,
            ha="center",
            va="center",
            fontsize=8,
            zorder=4,
            bbox={"facecolor": WHITE, "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
        )
        for group, center in anchors.items()
    ]
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    bounds = ax.get_window_extent().padded(-6)
    # Treat every drawn point as an obstacle so labels never hide the data.
    points = [
        Bbox.from_bounds(*(xy - radius), 2 * radius, 2 * radius)
        for collection in ax.collections
        for xy, radius in zip(
            collection.get_offset_transform().transform(collection.get_offsets()),
            np.sqrt(np.resize(collection.get_sizes(), len(collection.get_offsets())))
            * ax.figure.dpi
            / 144,
        )
    ]
    centers = ax.transData.transform(np.array(list(anchors.values())))
    placed = []
    for index, text in enumerate(labels):
        extent = text.get_window_extent(renderer).padded(4)
        point, box = _label_position(
            (extent.width, extent.height),
            centers[index],
            placed + points,
            bounds,
            np.delete(centers, index, axis=0),
        )
        text.set_position(ax.transData.inverted().transform(point))
        placed.append(box)


def _draw_tsne_panel(ax, adata, proteomic_groups, group_styles) -> None:
    """Draw existing t-SNE coordinates with the proteomic-group encoding."""
    tsne = adata.obsm["X_tsne"]
    tissues = adata.obs["tissue"].values
    anchors = {}
    for group, members in proteomic_groups.items():
        points = tsne[np.isin(tissues, members)]
        color, marker = group_styles[group]
        ax.scatter(
            points[:, 0],
            points[:, 1],
            color=color,
            marker=marker,
            s=36,
            linewidths=0,
            zorder=3,
        )
        # Anchor on a real member: a split group's mean can sit in empty space.
        anchors[group] = points[
            np.argmin(np.linalg.norm(points - np.median(points, axis=0), axis=1))
        ]
    ax.set_title("Tissue atlas", loc="left", pad=12)
    ax.set(xlabel="t-SNE 1", ylabel="t-SNE 2")
    ax.margins(0.18)
    style_axes(ax)
    _label_groups(ax, anchors)


def _group_legend(ax, group_styles, **kwargs) -> None:
    """Show one neutral-text legend entry for each proteomic group."""
    if len(group_styles) < 2:
        return
    handles = [
        Line2D(
            [],
            [],
            marker=marker,
            linestyle="none",
            markersize=6,
            markerfacecolor=color,
            markeredgewidth=0,
            label=group,
        )
        for group, (color, marker) in group_styles.items()
    ]
    ax.legend(
        handles=handles,
        title="Proteomic group",
        labelcolor=TEXT_COLOR,
        frameon=False,
        **kwargs,
    )


def _draw_dendrogram_panel(
    ax, z_linkage, unique_tissues, tissue_to_group, group_styles
) -> None:
    """Use neutral leaf labels and the same group markers as the atlas."""
    scipy_dendro(
        z_linkage,
        labels=unique_tissues,
        ax=ax,
        orientation="right",
        leaf_font_size=8,
        above_threshold_color=SAMPLE_POINT_COLOR,
        color_threshold=0,
    )
    for collection in ax.collections:
        collection.set_linewidth(1)
    style_axes(ax)
    marker_transform = ax.get_yaxis_transform() + ScaledTranslation(
        -12 / 72, 0, ax.figure.dpi_scale_trans
    )
    for label in ax.get_yticklabels():
        label.set_color(TEXT_COLOR)
        color, marker = group_styles[tissue_to_group[label.get_text()]]
        ax.scatter(
            0,
            label.get_position()[1],
            transform=marker_transform,
            color=color,
            marker=marker,
            s=36,
            linewidths=0,
            clip_on=False,
            zorder=4,
        )
    ax.tick_params(axis="y", length=0, pad=25)
    ax.set_title("Tissue dendrogram", loc="left", pad=12)
    ax.set_xlabel("Ward distance (correlation)")


def _save_figure(fig, out_dir: Path, name: str, opts: _SaveOpts) -> None:
    """Save each atlas view with consistent output options."""
    fig.savefig(
        out_dir / f"{name}.png", dpi=opts.dpi, bbox_inches="tight", facecolor=WHITE
    )
    if opts.save_pdf:
        fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", facecolor=WHITE)
    plt.close(fig)
    logger.info("Saved %s.png", name)


def _save_individual_atlas(out_dir, adata, proteomic_groups, group_styles, opts):
    """Render a standalone atlas with a compact group legend."""
    fig, ax = plt.subplots(figsize=(4.5, 3.6), layout="constrained")
    _group_legend(ax, group_styles, loc="center left", bbox_to_anchor=(1.02, 0.5))
    _draw_tsne_panel(ax, adata, proteomic_groups, group_styles)
    _save_figure(fig, out_dir, "tissue_atlas", opts)


def _save_individual_dendrogram(out_dir, tree, group_styles, opts):
    """Render a standalone dendrogram with a group legend."""
    z_linkage, tissues, tissue_to_group = tree
    fig, ax = plt.subplots(
        figsize=(3.5, max(3, len(tissues) * 0.16)), layout="constrained"
    )
    _group_legend(ax, group_styles, loc="upper left", bbox_to_anchor=(1.02, 1))
    _draw_dendrogram_panel(ax, z_linkage, tissues, tissue_to_group, group_styles)
    _save_figure(fig, out_dir, "tissue_dendrogram", opts)


@static_style()
def plot_slide_atlas_dendrogram(
    adata: ad.AnnData,
    out_dir: Path,
    *,
    dpi: int = 300,
    save_pdf: bool = True,
) -> None:
    """Plot combined and standalone views of the same proteomic groups."""
    if "X_tsne" not in adata.obsm:
        logger.warning("No t-SNE embedding, skipping slide figure")
        return
    tissues = adata.obs["tissue"].values
    unique_tissues = sorted(np.unique(tissues))
    if len(unique_tissues) < 2:
        logger.warning(
            "Atlas dendrogram requires at least two tissues; found %d, "
            "skipping slide figure",
            len(unique_tissues),
        )
        return
    groups, tissue_to_group, z_linkage = _compute_tissue_groups(adata)
    styles = _group_styles(groups)
    adata.uns["proteomic_groups"] = groups
    adata.uns["tissue_colors"] = {
        tissue: styles[group][0] for tissue, group in tissue_to_group.items()
    }
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, max(3.6, len(unique_tissues) * 0.16)),
        gridspec_kw={"width_ratios": [1.3, 1]},
        layout="constrained",
    )
    _group_legend(
        axes[0], styles, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3
    )
    _draw_dendrogram_panel(axes[1], z_linkage, unique_tissues, tissue_to_group, styles)
    _draw_tsne_panel(axes[0], adata, groups, styles)
    for ax, letter in zip(axes, "ab"):
        panel_label(ax, letter)
    opts = _SaveOpts(dpi, save_pdf)
    _save_figure(fig, out_dir, "slide_atlas_dendrogram", opts)
    _save_individual_atlas(out_dir, adata, groups, styles, opts)
    _save_individual_dendrogram(
        out_dir, (z_linkage, unique_tissues, tissue_to_group), styles, opts
    )
