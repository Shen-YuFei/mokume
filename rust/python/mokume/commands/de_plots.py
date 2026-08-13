#!/usr/bin/env python3
"""DE-plotting command (volcano / heatmap / PCA) for the mokume wheel.

The Rust ``features2proteins`` pipeline owns the numbers: it writes the protein
intensity matrix (``--output``) and one differential-expression (DE) result CSV
per contrast (``--de-output`` / the path the CLI synthesizes inside the plot
directory). This script is the thin "copy-py" glue that reproduces the volcano /
heatmap / PCA plots the Python ``features_to_proteins`` path emits via
``mokume/pipeline/stages.py`` ``generate_plots`` (lines 2061-2085) and its
helpers, calling the same shared plotting functions
``mokume.plotting.differential_expression.plot_volcano`` / ``plot_heatmap`` /
``plot_pca_conditions``.

Design (matches ``peptides2protein_qc.py``: Rust writes the tables, Python
only draws): we never recompute DE or the protein matrix here. We read the CSVs
Rust produced and call the shared mokume plotting helpers, so the cells in the
plots are the cells in the Rust CSVs. No mokume algorithm is duplicated -- a
system ``pip install mokume-rs[plotting]`` provides the helpers.

argv contract (runnable via ``python -m mokume.commands.de_plots``):

    python -m mokume.commands.de_plots \
        --protein-matrix <proteins.csv> \
        --plot-dir <plots/> \
        [--sdrf <sdrf.tsv>] \
        [--volcano] [--heatmap] [--pca] \
        [--log2fc-threshold 0.5] [--fdr-threshold 0.05] \
        [--highlight-genes GENE1,GENE2] \
        [--irs-remove-reference] \
        [--contrast <key> <condA> <condB> <de_csv>] ...

Each ``--contrast`` flag carries one contrast: its output-file key (Python's
``f"{condA}-{condB}"`` or the single-contrast ``--de-output`` stem), the two
condition labels, and the path to the DE result CSV Rust wrote. Volcano and
heatmap plots are per-contrast; the PCA plot is one figure over all conditions
and needs no DE table.

Exit code 0 on success; non-zero (with a message on stderr) otherwise.
"""

from __future__ import annotations

import argparse
import sys


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="de_plots.py",
        description="Render features2proteins DE plots from Rust-written CSVs.",
    )
    parser.add_argument(
        "--protein-matrix",
        required=True,
        help="Protein intensity matrix CSV written by the Rust pipeline.",
    )
    parser.add_argument(
        "--plot-dir",
        required=True,
        help="Output directory for the PNG plots (created if absent).",
    )
    parser.add_argument(
        "--sdrf",
        default=None,
        help="SDRF file; required for heatmap (condition annotation) and PCA.",
    )
    parser.add_argument("--volcano", action="store_true", help="Render volcano plots.")
    parser.add_argument("--heatmap", action="store_true", help="Render DE heatmaps.")
    parser.add_argument(
        "--pca", action="store_true", help="Render the PCA-by-condition plot."
    )
    parser.add_argument(
        "--log2fc-threshold",
        type=float,
        default=0.5,
        help="|log2FC| significance threshold (matches --de-log2fc).",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR significance threshold (matches --de-fdr).",
    )
    parser.add_argument(
        "--highlight-genes",
        default=None,
        help="Comma-separated protein names to label on volcano plots.",
    )
    parser.add_argument(
        "--irs-remove-reference",
        action="store_true",
        help=(
            "Restrict the condition mapping to matrix columns "
            "(mirrors stages._plot_sample_conditions when IRS removed reference)."
        ),
    )
    parser.add_argument(
        "--contrast",
        nargs=4,
        action="append",
        default=[],
        metavar=("KEY", "COND_A", "COND_B", "DE_CSV"),
        help="One contrast: output key, condition A, condition B, DE result CSV.",
    )
    return parser.parse_args(argv)


def _plot_sample_conditions(sample_to_condition, protein_df, irs_remove_reference):
    """Mirror ``stages.MokumePostprocessing._plot_sample_conditions``.

    Without ``--irs-remove-reference`` the full SDRF mapping is returned; with
    it, only conditions for samples that survive as matrix columns are kept.
    """
    if not irs_remove_reference:
        return sample_to_condition
    protein_col = protein_df.columns[0]
    available = [column for column in protein_df.columns if column != protein_col]
    return {
        sample: condition
        for sample, condition in sample_to_condition.items()
        if sample in available
    }


def _significant_de_proteins(de_df, log2fc_threshold, fdr_threshold):
    """Mirror ``stages._significant_de_proteins`` (adj_pvalue + |log2FC| gate)."""
    return de_df[
        (de_df["adj_pvalue"] < fdr_threshold)
        & (de_df["log2FC"].abs() > log2fc_threshold)
    ]


def _generate_volcano_plots(pd, contrasts, plot_dir, args):
    """Mirror ``stages._generate_volcano_plots``: one volcano per contrast."""
    from mokume.plotting.differential_expression import plot_volcano

    highlight = (
        [gene.strip() for gene in args.highlight_genes.split(",") if gene.strip()]
        if args.highlight_genes
        else None
    )
    for key, _cond_a, _cond_b, de_csv in contrasts:
        de_df = pd.read_csv(de_csv)
        output_file = str(plot_dir / "volcano_{0}.png".format(key))
        plot_volcano(
            de_df,
            log2fc_threshold=args.log2fc_threshold,
            fdr_threshold=args.fdr_threshold,
            highlight_genes=highlight,
            title="Volcano Plot: {0}".format(key),
            output_file=output_file,
        )
        print("Volcano plot saved to {0}".format(output_file))


def _generate_heatmap_plots(
    pd, contrasts, protein_df, sample_to_condition, plot_dir, args
):
    """Mirror ``stages._generate_heatmap_plots``: per-contrast DE heatmap."""
    from mokume.plotting.differential_expression import plot_heatmap

    protein_col = protein_df.columns[0]
    for key, cond_a, cond_b, de_csv in contrasts:
        de_df = pd.read_csv(de_csv)
        sig = _significant_de_proteins(de_df, args.log2fc_threshold, args.fdr_threshold)
        if sig.empty:
            print("Heatmap skipped for {0}: no significant proteins".format(key))
            continue

        sig_proteins = (
            sig.sort_values(
                by="log2FC", key=lambda series: series.abs(), ascending=False
            )["ProteinName"]
            .head(50)
            .tolist()
        )
        contrast_samples = [
            sample
            for sample in protein_df.columns
            if sample != protein_col
            and sample_to_condition.get(sample) in (cond_a, cond_b)
        ]
        contrast_mapping = {
            sample: condition
            for sample, condition in sample_to_condition.items()
            if condition in (cond_a, cond_b)
        }
        output_file = str(plot_dir / "heatmap_{0}.png".format(key))
        plot_heatmap(
            protein_df[[protein_col] + contrast_samples],
            contrast_mapping,
            proteins=sig_proteins,
            title="DE Heatmap: {0} (top {1}/{2} sig)".format(
                key, len(sig_proteins), len(sig)
            ),
            output_file=output_file,
        )
        print("Heatmap saved to {0}".format(output_file))


def _generate_pca_plot(protein_df, sample_to_condition, plot_dir):
    """Mirror ``stages._generate_pca_plot``: one PCA figure over all conditions."""
    from mokume.plotting.differential_expression import plot_pca_conditions

    output_file = str(plot_dir / "pca_conditions.png")
    plot_pca_conditions(
        protein_df,
        sample_to_condition,
        title="PCA by Condition",
        output_file=output_file,
    )
    print("PCA plot saved to {0}".format(output_file))


def main(argv=None):
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    from pathlib import Path

    # Import pandas / mokume.plotting lazily so a clear message is printed when
    # the plotting extra is missing, rather than an opaque ImportError trace.
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "DE plots aborted: pandas is not installed ({0}). "
            "Install with: pip install mokume-rs[plotting]".format(exc)
        )

    try:
        from mokume.plotting import is_plotting_available
    except ImportError as exc:
        raise SystemExit(
            "DE plots aborted: mokume plotting dependencies are not installed "
            "({0}). Install with: pip install mokume-rs[plotting]".format(exc)
        )

    if not is_plotting_available():
        raise SystemExit(
            "DE plots aborted: mokume plotting dependencies (matplotlib, seaborn) "
            "are not installed. Install with: pip install mokume-rs[plotting]"
        )

    plot_dir = Path(args.plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    protein_df = pd.read_csv(args.protein_matrix)

    # Sample -> condition mapping (heatmap annotation + PCA coloring) mirrors
    # ``stages``: derive it from the SDRF, then apply the IRS-remove-reference
    # restriction exactly as ``_plot_sample_conditions`` does.
    sample_to_condition = {}
    if args.sdrf:
        from mokume.normalization.irs import detect_condition_from_sdrf

        sample_to_condition = detect_condition_from_sdrf(args.sdrf)
    sample_to_condition = _plot_sample_conditions(
        sample_to_condition, protein_df, args.irs_remove_reference
    )

    if args.volcano and args.contrast:
        _generate_volcano_plots(pd, args.contrast, plot_dir, args)
    if args.heatmap and args.contrast and sample_to_condition:
        _generate_heatmap_plots(
            pd, args.contrast, protein_df, sample_to_condition, plot_dir, args
        )
    if args.pca and sample_to_condition:
        _generate_pca_plot(protein_df, sample_to_condition, plot_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
