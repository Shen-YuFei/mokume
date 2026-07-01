"""
Distribution plotting functions for the mokume package.

This module provides functions for plotting data distributions,
box plots, and violin plots for QC reports.
"""

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.figure import Figure


def plot_distributions(
    dataset: pd.DataFrame,
    field: str,
    class_field: str,
    title: str = "",
    log2: bool = True,
    width: float = 10,
) -> Figure:
    """
    Print the quantile plot for the dataset.

    Parameters
    ----------
    dataset : pd.DataFrame
        DataFrame to plot.
    field : str
        Field that would be used in the dataframe to plot the quantile.
    class_field : str
        Field to group the quantile into classes.
    title : str, optional
        Title of the box plot.
    log2 : bool, optional
        Log the intensity values.
    width : float, optional
        Size of the plot.

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    pd.set_option("mode.chained_assignment", None)
    normalize = dataset[[field, class_field]].reset_index(drop=True)
    if log2:
        normalize[field] = np.log2(normalize[field])
    normalize.dropna(subset=[field], inplace=True)
    plt.figure(dpi=500, figsize=(width, 8))
    fig = sns.kdeplot(
        data=normalize, x=field, hue=class_field, palette="Paired", linewidth=2
    )
    sns.despine(ax=fig, top=True, right=True)
    plt.title(title)
    pd.set_option("mode.chained_assignment", "warn")

    return plt.gcf()


def plot_box_plot(
    dataset: pd.DataFrame,
    field: str,
    class_field: str,
    log2: bool = False,
    width: float = 10,
    rotation: int = 30,
    title: str = "",
    violin: bool = False,
) -> Figure:
    """
    Plot a box plot of two values field and classes field.

    Parameters
    ----------
    dataset : pd.DataFrame
        Dataframe with peptide intensities.
    field : str
        Intensity field.
    class_field : str
        Class to group the peptides.
    log2 : bool, optional
        Transform peptide intensities to log scale.
    width : float, optional
        Size of the plot.
    rotation : int, optional
        Rotation of the x-axis.
    title : str, optional
        Title of the box plot.
    violin : bool, optional
        Also add violin on top of box plot.

    Returns
    -------
    Figure
        Matplotlib figure object.
    """
    pd.set_option("mode.chained_assignment", None)
    normalized = dataset[[field, class_field]]
    np.seterr(divide="ignore")
    plt.figure(figsize=(width, 14))
    if log2:
        normalized[field] = np.log2(normalized[field])

    if violin:
        chart = sns.violinplot(
            x=class_field,
            y=field,
            data=normalized,
            boxprops=dict(alpha=0.3),
            palette="muted",
        )
    else:
        chart = sns.boxplot(
            x=class_field,
            y=field,
            data=normalized,
            boxprops=dict(alpha=0.3),
            palette="muted",
        )

    chart.set(title=title)
    chart.set_xticklabels(chart.get_xticklabels(), rotation=rotation, ha="right")
    pd.set_option("mode.chained_assignment", "warn")

    return plt.gcf()
