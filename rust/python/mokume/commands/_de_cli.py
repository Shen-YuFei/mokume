"""Shared CLI arguments for DE plots and interactive reports."""


def add_protein_matrix_argument(parser):
    """Add the protein-matrix input shared by downstream visualizations."""
    parser.add_argument(
        "-p",
        "--protein-matrix",
        required=True,
        help="Protein intensity matrix CSV written by the Rust pipeline.",
    )


def add_de_result_arguments(parser):
    """Add thresholds, highlights, and contrast tables shared by DE outputs."""
    parser.add_argument(
        "--log2fc",
        type=float,
        default=0.5,
        help="|log2FC| significance threshold (matches --de-log2fc).",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        help="FDR significance threshold (matches --de-fdr).",
    )
    parser.add_argument(
        "--highlight-protein",
        action="append",
        default=[],
        help="Protein name to highlight; repeat for multiple proteins.",
    )
    parser.add_argument(
        "--contrast",
        nargs=4,
        action="append",
        default=[],
        metavar=("KEY", "COND_A", "COND_B", "DE_CSV"),
        help="One contrast: output key, condition A, condition B, DE result CSV.",
    )
