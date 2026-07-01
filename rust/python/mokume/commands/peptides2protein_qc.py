#!/usr/bin/env python3
"""iBAQ QC-report command (density + box plots) for the mokume wheel.

The Rust ``peptides2protein`` iBAQ path computes the protein table itself and
writes it as a tab-separated file. This script is the thin "copy-py" glue that
reproduces the QC report that the Python ``peptides_to_protein(verbose=True)``
path emits (``mokume/quantification/ibaq.py`` lines 969-1058): density (KDE) and
box plots, log2-scaled, one pair per quantification column, written into a single
PDF.

Design (blueprint section 2.5, option ii -- "Rust writes the TSV, Python only
draws"): the Rust kernel stays the single source of the numbers, so the cells in
the QC plots are byte-for-byte the cells in the Rust TSV. We do **not** recompute
iBAQ here; we only read the table Rust produced and call the shared plotting
helpers ``mokume.plotting.plot_distributions`` / ``plot_box_plot`` (exactly
the functions the Python verbose block uses). No mokume algorithm code lives in
Rust -- the first-class ``mokume.plotting`` helpers draw the figures (their
third-party deps come from the ``plotting`` extra).

argv contract (runnable via ``python -m mokume.commands.peptides2protein_qc``):

    python -m mokume.commands.peptides2protein_qc \
        --protein-table <proteins.tsv> \
        --qc-report <QCprofile.pdf> \
        [--plot-column Ibaq|IbaqPpb] \
        [--tpa] [--ruler]

Exit code 0 on success; non-zero (with a message on stderr) otherwise.
"""

import argparse
import sys


# Column headers, mirroring ``mokume.core.constants`` so the command reads the
# same output schema the Rust kernel writes (kept inline to avoid importing the
# mokume extension just for string constants -- the plotting import below is the
# only hard mokume need).
SAMPLE_ID = "SampleID"
IBAQ = "Ibaq"
IBAQ_PPB = "IbaqPpb"
TPA = "TPA"
COPYNUMBER = "CopyNumber"
CONCENTRATION_NM = "Concentration[nM]"


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Render the peptides2protein iBAQ QC report from a protein table."
    )
    parser.add_argument(
        "--protein-table",
        required=True,
        help="Tab-separated protein table produced by the Rust iBAQ path.",
    )
    parser.add_argument(
        "--qc-report",
        required=True,
        help="Destination PDF for the QC images.",
    )
    parser.add_argument(
        "--plot-column",
        default=IBAQ,
        help=(
            "Primary quantification column to plot. The Python verbose path uses "
            "IbaqPpb when --normalize is set, otherwise Ibaq."
        ),
    )
    parser.add_argument(
        "--tpa",
        action="store_true",
        help="Also plot the TPA distribution (matches the Python verbose block).",
    )
    parser.add_argument(
        "--ruler",
        action="store_true",
        help=(
            "Also plot the CopyNumber and Concentration[nM] distributions "
            "(matches the Python verbose block)."
        ),
    )
    return parser.parse_args(argv)


def _require_column(frame, column, table_path):
    if column not in frame.columns:
        raise SystemExit(
            "QC report aborted: column '{0}' is not present in '{1}'. "
            "Available columns: {2}".format(
                column, table_path, ", ".join(map(str, frame.columns))
            )
        )


def main(argv=None):
    args = _parse_args(sys.argv[1:] if argv is None else argv)

    # Import pandas / mokume.plotting lazily so a clear message is printed when
    # the plotting extra is missing, rather than an opaque ImportError trace.
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - environment guard
        raise SystemExit(
            "QC report aborted: pandas is not installed ({0}). "
            "Install with: pip install mokume[plotting]".format(exc)
        )

    try:
        from mokume.plotting import PdfPages, plot_box_plot, plot_distributions
    except ImportError as exc:
        raise SystemExit(
            "QC report aborted: mokume plotting dependencies are not "
            "installed ({0}). Install with: pip install mokume[plotting]".format(exc)
        )

    res = pd.read_csv(args.protein_table, sep="\t")
    _require_column(res, SAMPLE_ID, args.protein_table)
    _require_column(res, args.plot_column, args.protein_table)

    # ``plot_width = len(set(SampleID)) * 0.5 + 10`` -- identical to
    # ``mokume/quantification/ibaq.py`` line 978.
    plot_width = len(set(res[SAMPLE_ID])) * 0.5 + 10

    pdf = PdfPages(args.qc_report)
    try:
        # Primary iBAQ column: density + box, exactly as ibaq.py lines 980-998.
        density = plot_distributions(
            res,
            args.plot_column,
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="{0} Distribution".format(args.plot_column),
        )
        box = plot_box_plot(
            res,
            args.plot_column,
            SAMPLE_ID,
            log2=True,
            width=plot_width,
            title="{0} Distribution".format(args.plot_column),
            violin=False,
        )
        pdf.savefig(density, bbox_inches="tight")
        pdf.savefig(box, bbox_inches="tight")

        if args.tpa:
            _require_column(res, TPA, args.protein_table)
            density_tpa = plot_distributions(
                res,
                TPA,
                SAMPLE_ID,
                log2=True,
                width=plot_width,
                title="TPA Distribution",
            )
            box_tpa = plot_box_plot(
                res,
                TPA,
                SAMPLE_ID,
                log2=True,
                width=plot_width,
                title="{0} Distribution".format(TPA),
                violin=False,
            )
            pdf.savefig(density_tpa, bbox_inches="tight")
            pdf.savefig(box_tpa, bbox_inches="tight")

        if args.ruler:
            for column in (COPYNUMBER, CONCENTRATION_NM):
                _require_column(res, column, args.protein_table)
                density_col = plot_distributions(
                    res,
                    column,
                    SAMPLE_ID,
                    width=plot_width,
                    log2=True,
                    title="{0} Distribution".format(column),
                )
                box_col = plot_box_plot(
                    res,
                    column,
                    SAMPLE_ID,
                    width=plot_width,
                    log2=True,
                    title="{0} Distribution".format(column),
                    violin=False,
                )
                pdf.savefig(density_col, bbox_inches="tight")
                pdf.savefig(box_col, bbox_inches="tight")
    finally:
        pdf.close()

    print(
        "QC report written to {0} ({1} sample(s))".format(
            args.qc_report, len(set(res[SAMPLE_ID]))
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
