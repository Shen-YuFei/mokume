"""Locally styled PCA and t-SNE plots with shared condition color/shape codes."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from mokume.plotting._common import static_style, style_axes
from mokume.plotting.style import (
    CONDITION_COLORS,
    CONDITION_MARKERS,
    CONDITION_STYLE_INDICES,
    SAMPLE_POINT_COLOR,
    SINGLE_SERIES_COLOR,
    WHITE,
)


def _embedding_plot(frame, axes, group_column, options):
    """Plot existing coordinates without modifying their calculation."""
    fig, ax = plt.subplots(figsize=options.get("figsize", (8, 6)))
    groups = sorted(frame[group_column].unique())
    if len(groups) > len(CONDITION_STYLE_INDICES):
        # Beyond the color/shape combinations groups are indistinguishable.
        ax.scatter(
            frame[axes[0]],
            frame[axes[1]],
            color=SAMPLE_POINT_COLOR,
            s=64,
            linewidths=2,
            edgecolors=WHITE,
        )
        groups = []
    for index, group in enumerate(groups):
        color, shape = CONDITION_STYLE_INDICES[index]
        subset = frame[frame[group_column] == group]
        ax.scatter(
            subset[axes[0]],
            subset[axes[1]],
            label=group,
            color=CONDITION_COLORS[color],
            s=64,
            linewidths=2,
            edgecolors=WHITE,
            marker=CONDITION_MARKERS[shape] if len(groups) > 3 else "o",
        )
    ax.set(xlabel=axes[0], ylabel=axes[1], title=options["title"])
    if len(groups) > 1:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(options["output"], bbox_inches="tight")
    plt.close(fig)


@static_style()
def plot_pca(
    df_pca: pd.DataFrame,
    output_file: str,
    x_col: str = "PC1",
    y_col: str = "PC2",
    hue_col: str = "batch",
    palette: str = "Set2",
    title: str = "PCA plot",
    figsize: tuple = (8, 6),
) -> None:
    """
    Plot a PCA scatter plot and save it to a file.

    Groups use the shared condition color/shape combinations.

    Parameters
    ----------
    df_pca : pd.DataFrame
        DataFrame containing PCA components and metadata.
    output_file : str
        Path to save the output figure.
    x_col : str, optional
        Column name for x-axis (default: "PC1").
    y_col : str, optional
        Column name for y-axis (default: "PC2").
    hue_col : str, optional
        Column name for color grouping (default: "batch").
    palette : str, optional
        Ignored; kept for backward compatibility.
    title : str, optional
        Plot title (default: "PCA plot").
    figsize : tuple, optional
        Figure size in inches (default: (8, 6)).
    """
    del palette
    _embedding_plot(
        df_pca,
        (x_col, y_col),
        hue_col,
        {"output": output_file, "title": title, "figsize": figsize},
    )


@static_style()
def compute_pca_with_plot(df: pd.DataFrame, n_components: int = 5) -> pd.DataFrame:
    """
    Compute principal components and display a variance explained plot.

    This function performs PCA on the input dataframe and generates
    a plot showing the cumulative variance explained by each component.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe where rows are samples and columns are features.
    n_components : int, optional
        Number of principal components to compute (default: 5).

    Returns
    -------
    pd.DataFrame
        DataFrame with principal component values, indexed by the input index.
    """
    pca = PCA(n_components=n_components)
    pca.fit(df)
    df_pca = pca.transform(df)
    df_pca = pd.DataFrame(
        df_pca, index=df.index, columns=[f"PC{i}" for i in range(1, n_components + 1)]
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    xi = np.arange(1, n_components + 1, step=1)
    ax.plot(
        xi,
        np.cumsum(pca.explained_variance_ratio_),
        marker="o",
        linewidth=2,
        color=SINGLE_SERIES_COLOR,
        markersize=8,
        markeredgecolor=WHITE,
        markeredgewidth=2,
    )
    ax.set(
        ylim=(0, 1.1),
        xlabel="Number of components",
        ylabel="Cumulative variance",
        title="PCA explained variance",
        xticks=xi,
    )
    ax.axhline(0.95, color=SAMPLE_POINT_COLOR, linestyle="--", linewidth=2)
    ax.text(0.02, 0.9, "95% reference", transform=ax.transAxes)
    style_axes(ax)
    fig.tight_layout()
    plt.show()
    return df_pca


@static_style()
def plot_tsne(
    df: pd.DataFrame, x_col: str, y_col: str, hue_col: str, file_name: str
) -> None:
    """
    Generate and save a t-SNE scatter plot from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing t-SNE coordinates and metadata.
    x_col : str
        Column name for x-axis (typically "tSNE1").
    y_col : str
        Column name for y-axis (typically "tSNE2").
    hue_col : str
        Column name for color grouping.
    file_name : str
        Output file path for the saved plot.
    """
    _embedding_plot(
        df,
        (x_col, y_col),
        hue_col,
        {"title": f"t-SNE of {len(df)} samples", "output": file_name},
    )
