"""Shared CLI arguments for DE plots and interactive reports."""


def add_protein_matrix_argument(parser):
    """Add the protein-matrix input shared by downstream visualizations."""
    parser.add_argument(
        "-p",
        "--protein-matrix",
        metavar="<FILE>",
        required=True,
        help="CSV input.",
    )


def add_de_result_arguments(parser):
    """Add thresholds, highlights, and contrast tables shared by DE outputs."""
    parser.add_argument(
        "--log2fc",
        metavar="<VALUE>",
        type=float,
        default=0.5,
        help="|log2FC| threshold.",
    )
    parser.add_argument(
        "--fdr",
        metavar="<FRACTION>",
        type=float,
        default=0.05,
        help="Significance threshold.",
    )
    parser.add_argument(
        "--highlight-protein",
        metavar="<PROTEIN>",
        action="append",
        default=[],
        help="May be repeated.",
    )
    parser.add_argument(
        "--contrast",
        nargs=4,
        action="append",
        default=[],
        metavar=("<KEY>", "<GROUP_A>", "<GROUP_B>", "<DE_FILE>"),
        help="May be repeated.",
    )
