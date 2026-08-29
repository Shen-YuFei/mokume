"""
PCA and t-SNE plotting functions for the mokume package.

This module provides functions for visualizing PCA and t-SNE results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


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
        Color palette to use (default: "Set2").
    title : str, optional
        Plot title (default: "PCA plot").
    figsize : tuple, optional
        Figure size in inches (default: (8, 6)).
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(x=x_col, y=y_col, hue=hue_col, data=df_pca, palette=palette, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


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

    plt.rcParams["figure.figsize"] = (12, 6)
    _, ax = plt.subplots()
    xi = np.arange(1, n_components + 1, step=1)
    y = np.cumsum(pca.explained_variance_ratio_)

    plt.ylim(0.0, 1.1)
    plt.plot(xi, y, marker="o", linestyle="--", color="b")
    plt.xlabel("Number of Components")
    plt.xticks(np.arange(0, n_components, step=1))
    plt.ylabel("Cumulative variance (%)")
    plt.title("The number of components needed to explain variance")

    plt.axhline(y=0.95, color="r", linestyle="-")
    plt.text(0.5, 0.85, "95% cut-off threshold", color="red", fontsize=16)
    ax.grid(axis="x")
    plt.show()

    return df_pca


def plot_tsne(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    hue_col: str,
    file_name: str,
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
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    sns.scatterplot(
        x=x_col, y=y_col, hue=hue_col, data=df, ax=ax, markers=["o", "+", "x"]
    )
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"{x_col} vs {y_col} with {hue_col} information")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=8)
    plt.subplots_adjust(right=0.8)
    plt.savefig(file_name)
    plt.close(fig)
