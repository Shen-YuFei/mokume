#!/usr/bin/env python3
"""Interactive HTML DE-report command for the mokume wheel.

The Rust ``features2proteins`` pipeline owns the numbers: it writes the protein
intensity matrix (``--output``) and one differential-expression (DE) result CSV
per contrast (``--de-output`` / the path the CLI synthesizes). This script is the
thin "copy-py" glue that reproduces the standalone interactive HTML report the
Python ``features_to_proteins`` path emits via ``mokume/pipeline/stages.py``
``generate_interactive_report`` (lines 2087-2129), calling the same shared
``mokume.reports.interactive.generate_de_report``.

Design (matches ``peptides2protein_qc.py`` / ``de_plots.py``: Rust writes
the tables, Python only draws): we never recompute DE or the protein matrix here.
We read the CSVs Rust produced and call the shared mokume report builder, so the
cells in the report are the cells in the Rust CSVs. No mokume algorithm is
duplicated -- a system ``pip install mokume[reports]`` provides the builder
(plotly is loaded from a CDN by the generated HTML; the Python builder itself
needs only pandas + numpy, but the availability gate matches Python and checks
for plotly).

argv contract (runnable via ``python -m mokume.commands.interactive_report``):

    python -m mokume.commands.interactive_report \
        --protein-matrix <proteins.csv> \
        --sdrf <sdrf.tsv> \
        [--report-output <report.html>] \
        [--plot-dir <plots/>] \
        [--log2fc-threshold 0.5] [--fdr-threshold 0.05] \
        [--highlight-genes GENE1,GENE2] \
        [--contrast <key> <condA> <condB> <de_csv>] ...

Each ``--contrast`` flag carries one contrast: its output-file key, the two
condition labels, and the path to the DE result CSV Rust wrote. The report is
per-contrast. The output path mirrors ``stages.generate_interactive_report``:
with ``--report-output`` and a single contrast the path is used verbatim; with
multiple contrasts the contrast key is inserted before the extension; without
``--report-output`` the report lands in ``<plot-dir or .>/report_<key>.html``.

Exit code 0 on success; non-zero (with a message on stderr) otherwise.
"""

from __future__ import annotations

import argparse
import sys


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="interactive_report.py",
        description="Render the features2proteins interactive HTML report.",
    )
    parser.add_argument(
        "--protein-matrix",
        required=True,
        help="Protein intensity matrix CSV written by the Rust pipeline.",
    )
    parser.add_argument(
        "--sdrf",
        required=True,
        help="SDRF file; the report needs the sample-to-condition mapping.",
    )
    parser.add_argument(
        "--report-output",
        default=None,
        help="Output HTML path (default: <plot-dir or .>/report_<key>.html).",
    )
    parser.add_argument(
        "--plot-dir",
        default=None,
        help="Plot directory used to place the report when --report-output is unset.",
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
        help="Comma-separated protein names to highlight on the volcano plot.",
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


def _resolve_output_html(report_output, plot_dir, key, n_contrasts):
    """Mirror the output-path logic in ``stages.generate_interactive_report``."""
    from pathlib import Path

    if report_output:
        if n_contrasts == 1:
            return report_output
        base = report_output.rsplit(".", 1)
        return "{0}_{1}.html".format(base[0], key)
    plot_dir = plot_dir or "."
    return str(Path(plot_dir) / "report_{0}.html".format(key))


def _validate_args(args):
    if args.report_output and args.plot_dir:
        raise SystemExit(
            "Interactive report aborted: choose --report-output or --plot-dir, not both"
        )
    if args.log2fc_threshold < 0:
        raise SystemExit(
            "Interactive report aborted: --log2fc-threshold must be non-negative"
        )
    if not 0 <= args.fdr_threshold <= 1:
        raise SystemExit(
            "Interactive report aborted: --fdr-threshold must be between 0 and 1"
        )


def _require_sample_conditions(sample_to_condition):
    if not sample_to_condition:
        raise SystemExit(
            "Interactive report aborted: the SDRF yielded no sample conditions"
        )


def _validate_contrast_conditions(sample_to_condition, cond_a, cond_b):
    missing_conditions = [
        condition
        for condition in (cond_a, cond_b)
        if condition not in set(sample_to_condition.values())
    ]
    if missing_conditions:
        raise SystemExit(
            "Interactive report aborted: contrast conditions absent from SDRF: "
            + ", ".join(missing_conditions)
        )


def main(argv=None):
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    _validate_args(args)

    # Import pandas / mokume.reports lazily so a clear message is printed when
    # the reports extra is missing, rather than an opaque ImportError trace.
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "Interactive report aborted: pandas is not installed ({0}). "
            "Install with: pip install mokume[reports]".format(exc)
        )

    try:
        from mokume.reports import is_interactive_available
    except ImportError as exc:
        raise SystemExit(
            "Interactive report aborted: mokume reports dependencies are not "
            "installed ({0}). Install with: pip install mokume[reports]".format(exc)
        )

    if not is_interactive_available():
        raise SystemExit(
            "Interactive report aborted: report dependencies (plotly) are not "
            "installed. Install with: pip install mokume[reports]"
        )

    if not args.contrast:
        raise SystemExit(
            "Interactive report aborted: no contrasts supplied. The report is "
            "generated per differential-expression contrast."
        )

    from mokume.reports.interactive import generate_de_report
    from mokume.normalization.irs import detect_condition_from_sdrf

    protein_df = pd.read_csv(args.protein_matrix)
    sample_to_condition = detect_condition_from_sdrf(args.sdrf)
    _require_sample_conditions(sample_to_condition)
    highlight = (
        [gene.strip() for gene in args.highlight_genes.split(",") if gene.strip()]
        if args.highlight_genes
        else None
    )

    n_contrasts = len(args.contrast)
    for key, cond_a, cond_b, de_csv in args.contrast:
        _validate_contrast_conditions(sample_to_condition, cond_a, cond_b)
        de_df = pd.read_csv(de_csv, float_precision="round_trip")
        output_html = _resolve_output_html(
            args.report_output, args.plot_dir, key, n_contrasts
        )
        generate_de_report(
            de_results=de_df,
            protein_df=protein_df,
            sample_to_condition=sample_to_condition,
            output_html=output_html,
            title="DE Report: {0} ({1} vs {2})".format(key, cond_a, cond_b),
            highlight_genes=highlight,
            log2fc_threshold=args.log2fc_threshold,
            fdr_threshold=args.fdr_threshold,
        )
        print("Interactive report saved to {0}".format(output_html))

    return 0


if __name__ == "__main__":
    sys.exit(main())
