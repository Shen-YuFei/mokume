"""
Plotting functions for differential expression analysis results.

Provides volcano plots, heatmaps, and PCA visualizations colored by condition.
"""

from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.transforms import Bbox
from sklearn.decomposition import PCA

from mokume.plotting._common import static_style, style_axes as _style_axes
from mokume.plotting.style import (
    CONDITION_COLORS,
    CONDITION_MARKERS,
    CONDITION_STYLE_INDICES,
    MISSING_COLOR,
    SIGNIFICANCE_COLORS,
    TEXT_COLOR,
)


@static_style()
def plot_volcano(
    de_results: pd.DataFrame,
    log2fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
    highlight_genes: Optional[list[str]] = None,
    title: str = "Volcano Plot",
    output_file: Optional[str] = None,
    figsize: tuple = (3.5, 3.2),
) -> plt.Figure:
    """
    Create a volcano plot from differential expression results.

    Parameters
    ----------
    de_results : pd.DataFrame
        DE results with columns: ProteinName, log2FC, pvalue, adj_pvalue, significance.
    log2fc_threshold : float
        Threshold for fold change significance lines.
    fdr_threshold : float
        Threshold for FDR significance line.
    highlight_genes : list[str], optional
        Protein names to label on the plot.
    title : str
        Plot title.
    output_file : str, optional
        Path to save the figure.
    figsize : tuple
        Figure dimensions.

    Returns
    -------
    matplotlib.figure.Figure
        The volcano plot figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    df = de_results.copy()
    df["neg_log10_p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    # Color mapping
    color_map = SIGNIFICANCE_COLORS
    colors = df["significance"].map(color_map)

    # One mark per protein: small, translucent and unringed so overlaps read as density.
    ax.scatter(
        df["log2FC"],
        df["neg_log10_p"],
        c=colors,
        alpha=0.6,
        s=8,
        linewidths=0,
        rasterized=True,
    )

    # Significance thresholds
    cutoff = df.loc[df["adj_pvalue"] < fdr_threshold, "pvalue"].max()
    if pd.notna(cutoff):
        ax.axhline(
            y=-np.log10(max(cutoff, 1e-300)),
            color="grey",
            linestyle="--",
            linewidth=0.8,
        )
    ax.axvline(x=log2fc_threshold, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(x=-log2fc_threshold, color="grey", linestyle="--", linewidth=0.8)

    # Label highlighted genes
    if highlight_genes:
        for gene in highlight_genes:
            matches = df[df["ProteinName"].str.contains(gene, case=False, na=False)]
            ax.scatter(
                matches["log2FC"],
                matches["neg_log10_p"],
                c=matches["significance"].map(color_map),
                s=36,
                edgecolors=TEXT_COLOR,
                linewidths=0.6,
                zorder=3,
            )
            for _, row in matches.iterrows():
                ax.annotate(
                    row["ProteinName"],
                    xy=(row["log2FC"], row["neg_log10_p"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
                )

    # Legend
    from matplotlib.patches import Patch

    counts = df["significance"].value_counts()
    legend_elements = [
        Patch(facecolor=color, label=f"{label} ({counts[label]})")
        for label, color in color_map.items()
        if counts.get(label, 0)
    ]
    if len(legend_elements) > 1:
        # A single row under the axes leaves the plot its full width.
        ax.legend(
            handles=legend_elements,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.2),
            ncol=len(legend_elements),
            frameon=False,
            handlelength=1,
            columnspacing=1,
            labelcolor=TEXT_COLOR,
        )
    elif legend_elements:
        title += f" ({legend_elements[0].get_label()})"

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)
    _style_axes(ax)
    fig.text(
        0.5,
        0.01,
        "Colors: adjusted p-value and |log2FC| cutoffs",
        ha="center",
    )

    # Reserve the bottom band for the caption so it clears the x-axis label.
    plt.tight_layout(rect=(0, 0.03, 1, 1))

    if output_file:
        plt.savefig(output_file, bbox_inches="tight")
        plt.close(fig)

    return fig


@static_style()
def plot_heatmap(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    proteins: Optional[list[str]] = None,
    title: str = "Protein Heatmap",
    output_file: Optional[str] = None,
    figsize: tuple = (6, 9),
    top_n: int = 50,
) -> plt.Figure:
    """
    Create a clustered heatmap of protein intensities.

    Parameters
    ----------
    protein_df : pd.DataFrame
        Wide-format protein matrix.
    sample_to_condition : dict
        Sample-to-condition mapping for annotation.
    proteins : list[str], optional
        Specific proteins to include. If None, uses top_n most variable.
    title : str
        Plot title.
    output_file : str, optional
        Path to save the figure.
    figsize : tuple
        Figure dimensions.
    top_n : int
        Number of most variable proteins to show (if proteins is None).

    Returns
    -------
    matplotlib.figure.Figure
        The heatmap figure.
    """
    protein_col = protein_df.columns[0]
    sample_cols = [c for c in protein_df.columns if c != protein_col]

    # Filter to samples with conditions
    valid_samples = [s for s in sample_cols if s in sample_to_condition]
    if not valid_samples:
        raise ValueError(
            "Heatmap: no protein matrix columns matched the condition mapping. "
            f"Protein columns (first 5): {sample_cols[:5]}, "
            f"Condition keys (first 5): {list(sample_to_condition.keys())[:5]}."
        )
    intensity = protein_df.set_index(protein_col)[valid_samples]

    # Log2 transform
    intensity = intensity.replace(0, np.nan)
    intensity = np.log2(intensity)

    # Select proteins
    if proteins:
        # Match by substring
        mask = intensity.index.to_series().apply(
            lambda x: any(p in str(x) for p in proteins)
        )
        intensity = intensity[mask]
    else:
        # Top N most variable
        variance = intensity.var(axis=1).dropna()
        top_proteins = variance.nlargest(top_n).index
        intensity = intensity.loc[top_proteins]

    if intensity.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No proteins to display", ha="center", va="center")
        if output_file:
            plt.savefig(output_file, bbox_inches="tight")
            plt.close(fig)
        return fig

    # Row-scale (z-score per protein)
    # Drop proteins with zero std (constant) or all-NaN before z-scoring
    row_std = intensity.std(axis=1)
    intensity = intensity[row_std.notna() & (row_std > 0)]
    intensity_z = intensity.subtract(intensity.mean(axis=1), axis=0).divide(
        intensity.std(axis=1), axis=0
    )
    # Fill remaining per-sample NaN with 0 (neutral in z-score) for clustering,
    # but mask those cells so missing values are not drawn as average intensity.
    missing = intensity_z.isna()
    intensity_z = intensity_z.fillna(0)

    if intensity_z.shape[0] < 2 or intensity_z.shape[1] < 2:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5,
            0.5,
            f"Not enough data for heatmap clustering\n"
            f"({intensity_z.shape[0]} proteins × {intensity_z.shape[1]} samples after filtering)",
            ha="center",
            va="center",
        )
        ax.set_title(title)
        if output_file:
            plt.savefig(output_file, bbox_inches="tight")
            plt.close(fig)
        return fig

    # Create condition color annotation
    conditions = [sample_to_condition.get(s, "Unknown") for s in valid_samples]
    condition_colors = _condition_palette(conditions)
    col_colors = pd.Series(
        [condition_colors[c] for c in conditions], index=valid_samples, name="Condition"
    )

    g = sns.clustermap(
        intensity_z,
        cmap="RdBu_r",
        center=0,
        mask=missing,
        col_colors=col_colors,
        figsize=figsize,
        dendrogram_ratio=0.1,
        cbar_pos=(0.02, 0.87, 0.03, 0.1),
        yticklabels=True,
        xticklabels=True,
    )
    g.ax_heatmap.set_facecolor(MISSING_COLOR)
    g.ax_heatmap.tick_params(colors=TEXT_COLOR)
    # Protein rows are dense; 7 pt keeps one label per row legible.
    g.ax_heatmap.tick_params(axis="y", labelsize=7)
    g.fig.suptitle(title, y=1.02, color=TEXT_COLOR)

    # Beside the column dendrogram the legend clears the row labels.
    _add_heatmap_legend(g.ax_col_dendrogram, condition_colors, missing.to_numpy().any())

    if output_file:
        g.savefig(output_file, bbox_inches="tight")
        plt.close(g.fig)

    return g.fig


def _add_heatmap_legend(ax, condition_colors, has_missing):
    """Add the condition legend, plus a Missing entry when cells are masked."""
    legend_patches = [
        plt.matplotlib.patches.Patch(color=color, label=cond)
        for cond, color in condition_colors.items()
    ]
    if has_missing:
        legend_patches.append(
            plt.matplotlib.patches.Patch(color=MISSING_COLOR, label="Missing")
        )
    ax.legend(
        handles=legend_patches,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Condition",
        labelcolor=TEXT_COLOR,
    ).set_visible(len(legend_patches) > 1)


def _condition_palette(conditions):
    """Assign colors from the nonrepeating color/shape combinations."""
    return {
        condition: CONDITION_COLORS[CONDITION_STYLE_INDICES[index][0]]
        for index, condition in enumerate(sorted(set(conditions)))
    }


def _label_pca_samples(ax, pca_df):
    """Label outliers first; drop labels that would cover another label or sample.

    Call after the final layout, since overlaps are measured in display space.
    """
    if len(pca_df) > 30:
        return
    y = pca_df["PC2"] if "PC2" in pca_df else pca_df["PC1"]
    points = ax.transData.transform(np.column_stack([pca_df["PC1"], y]))
    radius = np.sqrt(ax.collections[0].get_sizes()[0]) * ax.figure.dpi / 144
    dots = [Bbox.from_bounds(*(xy - radius), 2 * radius, 2 * radius) for xy in points]
    renderer = ax.figure.canvas.get_renderer()
    placed = []
    for index in np.argsort(-np.linalg.norm(points - points.mean(axis=0), axis=1)):
        label = ax.annotate(
            pca_df.index[index],
            (pca_df["PC1"].iloc[index], y.iloc[index]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )
        box = label.get_window_extent(renderer)
        if box.x1 > ax.get_window_extent().x1:
            # Near the right edge the label goes left so it stays clear of the legend.
            label.set_horizontalalignment("right")
            label.xyann = (-4, 4)
            box = label.get_window_extent(renderer)
        # A label never hides another sample's point; its own point sits clear of it.
        others = placed + dots[:index] + dots[index + 1 :]
        if any(box.overlaps(other) for other in others):
            label.remove()
        else:
            placed.append(box)


@static_style()
def plot_pca_conditions(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    title: str = "PCA by Condition",
    output_file: Optional[str] = None,
    figsize: tuple = (3.5, 3),
) -> plt.Figure:
    """
    Create a PCA scatter plot colored by experimental condition.

    Parameters
    ----------
    protein_df : pd.DataFrame
        Wide-format protein matrix.
    sample_to_condition : dict
        Sample-to-condition mapping for coloring.
    title : str
        Plot title.
    output_file : str, optional
        Path to save the figure.
    figsize : tuple
        Figure dimensions.

    Returns
    -------
    matplotlib.figure.Figure
        The PCA scatter plot figure.
    """
    protein_col = protein_df.columns[0]
    sample_cols = [c for c in protein_df.columns if c != protein_col]

    # Filter samples with conditions
    valid_samples = [s for s in sample_cols if s in sample_to_condition]
    if not valid_samples:
        raise ValueError(
            "PCA: no protein matrix columns matched the condition mapping. "
            f"Protein columns (first 5): {sample_cols[:5]}, "
            f"Condition keys (first 5): {list(sample_to_condition.keys())[:5]}. "
            "Check that detect_condition_from_sdrf maps run file names correctly."
        )

    intensity = protein_df.set_index(protein_col)[valid_samples]

    # Log2 transform
    intensity = intensity.replace(0, np.nan)
    intensity = np.log2(intensity)

    # Transpose: samples as rows, proteins as columns
    data = intensity.T.dropna(axis=1)

    if data.shape[1] < 2:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "Not enough features for PCA", ha="center", va="center")
        if output_file:
            plt.savefig(output_file, bbox_inches="tight")
            plt.close(fig)
        return fig

    # Compute PCA
    n_components = min(2, data.shape[0], data.shape[1])
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data)

    pca_df = pd.DataFrame(
        pca_result[:, :2],
        columns=["PC1", "PC2"] if n_components >= 2 else ["PC1"],
        index=data.index,
    )
    pca_df["Condition"] = [sample_to_condition[s] for s in pca_df.index]

    fig, ax = plt.subplots(figsize=figsize)
    conditions = sorted(pca_df["Condition"].unique())
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2" if "PC2" in pca_df.columns else "PC1",
        hue="Condition",
        style="Condition",
        markers={
            condition: CONDITION_MARKERS[CONDITION_STYLE_INDICES[index][1]]
            for index, condition in enumerate(conditions)
        },
        palette=_condition_palette(pca_df["Condition"]),
        s=40,
        edgecolor="#ffffff",
        linewidth=0.6,
        ax=ax,
    )

    var_explained = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    if n_components >= 2:
        ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
    ax.set_title(f"{title}\nn = {data.shape[1]:,} proteins, complete cases")
    ax.legend(
        loc="center left", bbox_to_anchor=(1.05, 0.5), labelcolor=TEXT_COLOR
    ).set_visible(len(conditions) > 1)
    _style_axes(ax)

    plt.tight_layout()
    _label_pca_samples(ax, pca_df)

    if output_file:
        plt.savefig(output_file, bbox_inches="tight")
        plt.close(fig)

    return fig
