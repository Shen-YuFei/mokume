"""
Plotting functions for differential expression analysis results.

Provides volcano plots, heatmaps, and PCA visualizations colored by condition.
"""

from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from sklearn.decomposition import PCA

matplotlib.use("Agg")


def plot_volcano(
    de_results: pd.DataFrame,
    log2fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
    highlight_genes: Optional[list[str]] = None,
    title: str = "Volcano Plot",
    output_file: Optional[str] = None,
    figsize: tuple = (10, 8),
) -> plt.Figure:
    """
    Create a volcano plot from differential expression results.

    Parameters
    ----------
    de_results : pd.DataFrame
        DE results with columns: ProteinName, log2FC, adj_pvalue, significance.
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
    df["-log10(adj_pvalue)"] = -np.log10(df["adj_pvalue"].clip(lower=1e-300))

    # Color mapping
    color_map = {
        "UP": "#E74C3C",
        "DOWN": "#3498DB",
        "Unchanged": "#BDC3C7",
        "NotTested": "#666666",
    }
    colors = df["significance"].map(color_map)

    ax.scatter(
        df["log2FC"],
        df["-log10(adj_pvalue)"],
        c=colors,
        alpha=0.6,
        s=20,
        edgecolors="none",
    )

    # Significance thresholds
    ax.axhline(y=-np.log10(fdr_threshold), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(x=log2fc_threshold, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(x=-log2fc_threshold, color="grey", linestyle="--", linewidth=0.8)

    # Label highlighted genes
    if highlight_genes:
        for gene in highlight_genes:
            matches = df[df["ProteinName"].str.contains(gene, case=False, na=False)]
            for _, row in matches.iterrows():
                ax.annotate(
                    row["ProteinName"],
                    xy=(row["log2FC"], row["-log10(adj_pvalue)"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
                )

    # Legend
    from matplotlib.patches import Patch

    n_up = (df["significance"] == "UP").sum()
    n_down = (df["significance"] == "DOWN").sum()
    n_ns = (df["significance"] == "Unchanged").sum()
    n_not_tested = (df["significance"] == "NotTested").sum()
    legend_elements = [
        Patch(facecolor=color_map["UP"], label=f"UP ({n_up})"),
        Patch(facecolor=color_map["DOWN"], label=f"DOWN ({n_down})"),
        Patch(facecolor=color_map["Unchanged"], label=f"Unchanged ({n_ns})"),
    ]
    if n_not_tested:
        legend_elements.append(
            Patch(
                facecolor=color_map["NotTested"],
                label=f"NotTested ({n_not_tested})",
            )
        )
    ax.legend(handles=legend_elements, loc="upper right")

    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax.set_title(title, fontsize=14)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)

    return fig


def plot_heatmap(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    proteins: Optional[list[str]] = None,
    title: str = "Protein Heatmap",
    output_file: Optional[str] = None,
    figsize: tuple = (12, 10),
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
            plt.savefig(output_file, dpi=150, bbox_inches="tight")
            plt.close(fig)
        return fig

    # Row-scale (z-score per protein)
    # Drop proteins with zero std (constant) or all-NaN before z-scoring
    row_std = intensity.std(axis=1)
    intensity = intensity[row_std.notna() & (row_std > 0)]
    intensity_z = intensity.subtract(intensity.mean(axis=1), axis=0).divide(
        intensity.std(axis=1), axis=0
    )
    # Fill remaining per-sample NaN with 0 (neutral in z-score) for clustering
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
            plt.savefig(output_file, dpi=150, bbox_inches="tight")
            plt.close(fig)
        return fig

    # Create condition color annotation
    conditions = [sample_to_condition.get(s, "Unknown") for s in valid_samples]
    unique_conditions = sorted(set(conditions))
    palette = sns.color_palette("Set2", len(unique_conditions))
    condition_colors = {c: palette[i] for i, c in enumerate(unique_conditions)}
    col_colors = pd.Series(
        [condition_colors[c] for c in conditions], index=valid_samples, name="Condition"
    )

    g = sns.clustermap(
        intensity_z,
        cmap="RdBu_r",
        center=0,
        col_colors=col_colors,
        figsize=figsize,
        dendrogram_ratio=0.1,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        yticklabels=True,
        xticklabels=True,
    )
    g.fig.suptitle(title, y=1.02, fontsize=14)

    # Add condition legend
    legend_patches = [
        plt.matplotlib.patches.Patch(color=condition_colors[c], label=c)
        for c in unique_conditions
    ]
    g.ax_heatmap.legend(
        handles=legend_patches,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Condition",
    )

    if output_file:
        g.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(g.fig)

    return g.fig


def plot_pca_conditions(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    title: str = "PCA by Condition",
    output_file: Optional[str] = None,
    figsize: tuple = (8, 6),
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
            plt.savefig(output_file, dpi=150, bbox_inches="tight")
            plt.close(fig)
        return fig

    # Fill remaining NaN with column median
    data = data.fillna(data.median())

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
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2" if "PC2" in pca_df.columns else "PC1",
        hue="Condition",
        palette="Set2",
        s=100,
        ax=ax,
    )

    var_explained = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)", fontsize=12)
    if n_components >= 2:
        ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)

    return fig
