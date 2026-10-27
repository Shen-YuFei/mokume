"""Abundance ranks and multi-contrast summaries without extra dependencies."""

from collections import Counter

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patheffects

from mokume.plotting._common import static_style, style_axes
from mokume.plotting.style import (
    CONDITION_COLORS,
    CONDITION_MARKERS,
    INACTIVE_COLOR,
    SAMPLE_POINT_COLOR,
    SINGLE_SERIES_COLOR,
    TEXT_COLOR,
    WHITE,
)


@static_style()
def plot_abundance_rank(log2_values, highlight_genes=None):
    """Rank proteins by observed mean log2 abundance."""
    abundance = log2_values.mean(axis=1).dropna().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    ranks = np.arange(1, len(abundance) + 1)
    ax.plot(ranks, abundance, color=SINGLE_SERIES_COLOR, linewidth=2)
    for protein in highlight_genes or []:
        if protein in abundance.index:
            rank = abundance.index.get_loc(protein) + 1
            ax.scatter(
                rank,
                abundance[protein],
                color=SINGLE_SERIES_COLOR,
                s=36,
                edgecolors=TEXT_COLOR,
                linewidths=0.6,
                zorder=3,
            )
            ax.annotate(
                protein,
                (rank, abundance[protein]),
                xytext=(5, 5),
                textcoords="offset points",
                color=TEXT_COLOR,
            )
    ax.set(xlabel="Protein rank", ylabel="Mean log2 intensity", title="Abundance rank")
    style_axes(ax)
    fig.tight_layout()
    return fig


def upset_counts(sets):
    """Count exclusive membership patterns, as in a standard UpSet plot."""
    return Counter(
        tuple(name for name, members in sets.items() if protein in members)
        for protein in set().union(*sets.values())
    )


@static_style()
def plot_upset(sets, output_file):
    """Show intersection sizes, set sizes and the membership matrix."""
    counts = upset_counts(sets)
    intersections = sorted(counts, key=lambda key: (-counts[key], key))
    names = list(sets)
    fig = plt.figure(figsize=(max(4.5, len(intersections) * 0.32), 3.6))
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 3], height_ratios=[2, 1])
    top, left, matrix = (
        fig.add_subplot(grid[position]) for position in ((0, 1), (1, 0), (1, 1))
    )
    bars = top.bar(
        range(len(intersections)),
        [counts[key] for key in intersections],
        color=SINGLE_SERIES_COLOR,
        width=0.7,
        linewidth=0,
    )
    top.bar_label(bars, color=TEXT_COLOR)
    top.margins(y=0.12)
    left.barh(
        range(len(names)),
        [len(sets[name]) for name in names],
        color=SINGLE_SERIES_COLOR,
        height=0.7,
        linewidth=0,
    )
    _intersection_matrix(matrix, names, intersections)
    top.set(
        ylabel="Intersection size",
        title="Significant proteins across contrasts",
        xticks=[],
    )
    left.set(xlabel="Set size", yticks=range(len(names)), yticklabels=names)
    matrix.set(xticks=[], yticks=range(len(names)), yticklabels=[])
    matrix.tick_params(left=False)
    left.invert_xaxis()
    for ax in (top, left, matrix):
        style_axes(ax)
    matrix.set_ylim(left.get_ylim())
    top.set_xlim(matrix.get_xlim())
    fig.tight_layout()
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


def _intersection_matrix(ax, names, intersections):
    """Draw active membership points over a light inactive grid."""
    for x, intersection in enumerate(intersections):
        enabled = [index for index, name in enumerate(names) if name in intersection]
        ax.scatter(
            [x] * len(names),
            range(len(names)),
            s=64,
            color=INACTIVE_COLOR,
            edgecolors=WHITE,
            linewidths=2,
        )
        ax.plot([x] * len(enabled), enabled, color=SAMPLE_POINT_COLOR, linewidth=2)
        ax.scatter(
            [x] * len(enabled),
            enabled,
            s=64,
            color=SAMPLE_POINT_COLOR,
            edgecolors=WHITE,
            linewidths=2,
            zorder=3,
        )


def _fold_change_series(ax, tables, protein, color, marker):
    """Use solid and hollow markers for significant and nonsignificant values."""
    rows = [
        table.set_index("ProteinName").reindex([protein]).iloc[0]
        for table in tables.values()
    ]
    values = [row["log2FC"] for row in rows]
    significant = [row["significance"] in ("UP", "DOWN") for row in rows]
    ax.plot(range(len(rows)), values, color=color, linewidth=2, label=protein)
    for index, value in enumerate(values):
        effects = (
            []
            if significant[index]
            else [
                patheffects.Stroke(linewidth=6, foreground=WHITE),
                patheffects.Normal(),
            ]
        )
        ax.scatter(
            index,
            value,
            s=64,
            marker=marker,
            linewidths=2,
            edgecolors=WHITE if significant[index] else color,
            facecolors=color if significant[index] else WHITE,
            path_effects=effects,
            zorder=3,
        )


@static_style()
def plot_fold_changes(tables, proteins, output_file):
    """Compare fold changes, switching to small multiples above six proteins."""
    multiple = len(proteins) > 6
    columns, rows = (3, int(np.ceil(len(proteins) / 3))) if multiple else (1, 1)
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(7.2, rows * 2.4) if multiple else (3.5, 3),
        squeeze=False,
        sharey=True,
    )
    for index, protein in enumerate(proteins):
        ax = axes.flat[index] if multiple else axes[0, 0]
        color = SINGLE_SERIES_COLOR if multiple else CONDITION_COLORS[index]
        marker = "o" if multiple else CONDITION_MARKERS[index]
        _fold_change_series(ax, tables, protein, color, marker)
        ax.set_title(protein if multiple else "log2FC across contrasts")
    rotation = 30 if max(map(len, tables)) > 8 else 0
    for ax in axes.flat:
        ax.axhline(0, color=SAMPLE_POINT_COLOR, linewidth=0.6)
        ax.set_xticks(
            range(len(tables)),
            list(tables),
            rotation=rotation,
            ha="right" if rotation else "center",
        )
        style_axes(ax)
    fig.supylabel("log2 fold change", color=TEXT_COLOR)
    if multiple:
        for ax in list(axes.flat)[len(proteins) :]:
            ax.set_visible(False)
    elif len(proteins) > 1:
        axes[0, 0].legend(
            labelcolor=TEXT_COLOR, loc="center left", bbox_to_anchor=(1.02, 0.5)
        )
    fig.text(
        0.5,
        0.01,
        "Solid: significant; hollow: nonsignificant. No confidence "
        "intervals:\nDE tables contain no standard errors.",
        ha="center",
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)
