"""Peptide/protein profiles using SDRF conditions and log2 intensities."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from mokume.plotting._common import (
    finite_list,
    log2_matrix,
    ordered_samples,
    read_frame,
    static_style,
    style_axes,
)
from mokume.plotting.style import (
    PEPTIDE_COLOR,
    SINGLE_SERIES_COLOR,
    TEXT_COLOR,
)


def peptide_matrix(source, sample_to_condition):
    """Read four- or six-column exports; conditions always come from SDRF."""
    frame = read_frame(source)
    frame = frame[frame["SampleID"].isin(sample_to_condition)].copy()
    frame["Condition"] = frame["SampleID"].map(sample_to_condition)
    values = frame.pivot_table(
        index=["ProteinName", "PeptideCanonical"],
        columns="SampleID",
        values="NormIntensity",
        aggfunc="sum",
    )
    return np.log2(values.where(values > 0)).reindex(
        columns=ordered_samples(values.columns, sample_to_condition)
    )


def peptide_payload(source, sample_to_condition):
    """Encode peptide traces with two decimals and null missing values."""
    values = peptide_matrix(source, sample_to_condition)
    return {
        str(protein): [finite_list(row) for row in group.to_numpy()]
        for protein, group in values.groupby(level="ProteinName", sort=False)
    }, list(values.columns)


def condition_groups(samples, conditions):
    """Return group boundaries and (center, condition) labels for ordered samples."""
    starts = [
        index
        for index in range(len(samples))
        if index == 0 or conditions[samples[index]] != conditions[samples[index - 1]]
    ]
    ends = starts[1:] + [len(samples)]
    labels = [
        [(start + end - 1) / 2, conditions[samples[start]]]
        for start, end in zip(starts, ends)
    ]
    return [start - 0.5 for start in starts[1:]], labels


@static_style()
def plot_profile(protein_df, peptides, sample_to_condition, protein, output_file):
    """Plot gray peptide trajectories behind the protein summary."""
    values = peptide_matrix(peptides, sample_to_condition)
    proteins = log2_matrix(protein_df)
    samples = ordered_samples(
        [sample for sample in proteins if sample in sample_to_condition],
        sample_to_condition,
    )
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    has_peptides = protein in values.index.get_level_values("ProteinName")
    if has_peptides:
        ax.plot(
            values.loc[protein].reindex(columns=samples).T.to_numpy(),
            color=PEPTIDE_COLOR,
            linewidth=1,
        )
    ax.plot(
        proteins.loc[protein, samples],
        color=SINGLE_SERIES_COLOR,
        linewidth=2,
        marker="o",
        markersize=6,
        markeredgewidth=0,
    )
    boundaries, labels = condition_groups(samples, sample_to_condition)
    for boundary in boundaries:
        ax.axvline(boundary, color=PEPTIDE_COLOR, linewidth=0.6)
    ax.set_xticks(range(len(samples)), samples, rotation=60, ha="right")
    ax.set(xlabel="Sample", ylabel="log2 intensity", title=protein)
    conditions = ax.secondary_xaxis("top")
    conditions.set_xticks(*zip(*labels))
    conditions.tick_params(length=0, colors=TEXT_COLOR)
    conditions.spines["top"].set_visible(False)
    if has_peptides:
        ax.legend(
            handles=[
                Line2D([], [], color=PEPTIDE_COLOR, label="Peptides", lw=1),
                Line2D([], [], color=SINGLE_SERIES_COLOR, label="Protein", lw=2),
            ]
        )
    style_axes(ax)
    fig.tight_layout()
    _hide_crowded_labels(fig, ax, conditions, boundaries, len(samples))
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


def _hide_crowded_labels(fig, ax, conditions, boundaries, count):
    """Drop sample names that would overlap and condition names wider than their block."""

    def step():
        return np.diff(ax.transData.transform([(0, 0), (1, 0)])[:, 0])[0]

    size = ax.xaxis.get_ticklabels()[0].get_fontsize() * fig.dpi / 72
    if step() < 1.2 * size / np.sin(np.radians(60)):
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        fig.tight_layout()
    fig.canvas.draw()
    widths = np.diff([-0.5, *boundaries, count - 0.5]) * step()
    names = conditions.get_xticklabels()
    conditions.set_xticks(
        conditions.get_xticks(),
        [
            name.get_text() if name.get_window_extent().width <= width else ""
            for name, width in zip(names, widths)
        ],
    )
