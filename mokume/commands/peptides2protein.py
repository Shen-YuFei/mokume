"""
CLI command for computing protein quantification values.
"""

import logging

import click
import pandas as pd

from mokume.quantification.ibaq import peptides_to_protein
from mokume.quantification import (
    get_quantification_method,
    is_directlfq_available,
)
from mokume.model.organism import OrganismDescription
from mokume.core.constants import (
    PROTEIN_NAME,
    SAMPLE_ID,
    CONDITION,
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    is_parquet,
)

logger = logging.getLogger(__name__)


# Base methods always available
QUANTIFICATION_METHODS = ["ibaq", "top3", "topn", "maxlfq", "sum", "directlfq"]


def get_available_methods():
    """Get list of available quantification methods based on installed packages."""
    methods = ["ibaq", "top3", "topn", "maxlfq", "sum"]
    if is_directlfq_available():
        methods.append("directlfq")
    return methods


@click.command("peptides2protein", short_help="Compute protein quantification values")
@click.option(
    "-f",
    "--fasta",
    help="Protein database to compute IBAQ values (required for ibaq method)",
    type=click.Path(exists=True),
)
@click.option(
    "-p",
    "--peptides",
    help="Peptide identifications with intensities following the peptide intensity output",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "--method",
    help="Quantification method to use (directlfq requires: pip install mokume[directlfq])",
    type=click.Choice(QUANTIFICATION_METHODS, case_sensitive=False),
    default="ibaq",
)
@click.option(
    "-e",
    "--enzyme",
    help="Enzyme used during the analysis of the dataset (default: Trypsin)",
    default="Trypsin",
)
@click.option(
    "-n",
    "--normalize",
    help="Normalize quantification values",
    is_flag=True,
)
@click.option(
    "--min_aa", help="Minimum number of amino acids to consider a peptide", default=7
)
@click.option(
    "--max_aa", help="Maximum number of amino acids to consider a peptide", default=30
)
@click.option(
    "-t", "--tpa", help="Whether calculate TPA (iBAQ method only)", is_flag=True
)
@click.option(
    "-r",
    "--ruler",
    help="Whether to use ProteomicRuler (iBAQ method only)",
    is_flag=True,
)
@click.option("-i", "--ploidy", help="Ploidy number (default: 2)", default=2)
@click.option(
    "-m",
    "--organism",
    help="Organism source of the data (default: human)",
    type=click.Choice(
        sorted(map(str.lower, OrganismDescription.registered_organisms())),
        case_sensitive=False,
    ),
    default="human",
)
@click.option(
    "-c",
    "--cpc",
    help="Cellular protein concentration(g/L) (default: 200)",
    default=200,
)
@click.option(
    "-o", "--output", help="Output file with the proteins and quantification values"
)
@click.option(
    "--verbose",
    help="Print additional information about the distributions of the intensities",
    is_flag=True,
)
@click.option(
    "--qc_report",
    help="PDF file to store multiple QC images (iBAQ method only)",
    default="QCprofile.pdf",
)
@click.option(
    "--topn_n",
    help="Number of top peptides to use for TopN method (default: 3)",
    default=3,
    type=int,
)
@click.option(
    "--threads",
    help="Number of parallel threads for MaxLFQ (-1 for all cores, default: -1)",
    default=-1,
    type=int,
)
@click.option(
    "--min_nonan",
    help="Minimum non-NaN ion intensities per protein for DirectLFQ (default: 1)",
    default=1,
    type=int,
)
@click.option(
    "--families",
    "families_yaml",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help=(
        "Optional YAML file declaring explicit protein family overrides "
        "for piBAQ (paralog-aware iBAQ; schema: families: [{name, members}])."
    ),
)
@click.option(
    "--min-shared",
    "min_shared",
    type=int,
    default=2,
    show_default=True,
    help=(
        "Minimum number of distinct peptides two proteins must share to be "
        "automatically grouped into the same iBAQ family."
    ),
)
@click.option(
    "--min-anchors",
    "min_anchors",
    type=int,
    default=1,
    show_default=True,
    help=(
        "Minimum proteotypic ('anchor') peptides a family member needs for the "
        "family to stay on the per-protein branch. The family rolls up to a "
        "single iBAQ only when NO member reaches this threshold (iBAQ only)."
    ),
)
@click.option(
    "--high-anchor-threshold",
    "high_anchor_threshold",
    type=int,
    default=3,
    show_default=True,
    help=(
        "Minimum anchor count (weakest member) for a family to be labelled "
        "EvidenceLevel='high' (iBAQ only)."
    ),
)
@click.pass_context
def peptides2protein(
    click_context,
    fasta: str,
    peptides: str,
    method: str,
    enzyme: str,
    normalize: bool,
    min_aa: int,
    max_aa: int,
    tpa: bool,
    ruler: bool,
    organism: str,
    ploidy: int,
    cpc: float,
    output: str,
    verbose: bool,
    qc_report: str,
    topn_n: int,
    threads: int,
    min_nonan: int,
    families_yaml: str,
    min_shared: int,
    min_anchors: int,
    high_anchor_threshold: int,
) -> None:
    """
    Compute protein quantification values from peptide intensity data.

    This command processes peptide identifications and computes protein
    quantification values using various methods:

    \b
    - ibaq: Intensity-Based Absolute Quantification (default, requires FASTA)
    - top3: Average of the 3 most intense peptides
    - topn: Average of the N most intense peptides (use --topn_n)
    - maxlfq: MaxLFQ delayed normalization algorithm (parallelized)
    - sum: Sum of all peptide intensities
    - directlfq: DirectLFQ intensity traces (requires: pip install mokume[directlfq])

    For the iBAQ method, a FASTA file is required. Other methods can work
    without a FASTA file.

    \b
    Examples:
        # Using iBAQ (requires FASTA)
        mokume peptides2protein --method ibaq -f proteome.fasta -p peptides.csv -o proteins.tsv

        # Using MaxLFQ with 4 threads
        mokume peptides2protein --method maxlfq --threads 4 -p peptides.csv -o proteins.tsv

        # Using DirectLFQ (requires optional install)
        mokume peptides2protein --method directlfq -p peptides.csv -o proteins.tsv
    """
    method_lower = method.lower()

    # Check DirectLFQ availability
    if method_lower == "directlfq" and not is_directlfq_available():
        raise click.UsageError(
            "DirectLFQ is not installed. Install with: pip install mokume[directlfq]"
        )

    if method_lower == "ibaq":
        # Use the existing iBAQ implementation
        if not fasta:
            raise click.UsageError("The --fasta option is required for the iBAQ method")

        peptides_to_protein(
            fasta=fasta,
            peptides=peptides,
            enzyme=enzyme,
            normalize=normalize,
            min_aa=min_aa,
            max_aa=max_aa,
            tpa=tpa,
            ruler=ruler,
            ploidy=ploidy,
            cpc=cpc,
            organism=organism,
            output=output,
            verbose=verbose,
            qc_report=qc_report,
            families_yaml=families_yaml,
            min_shared=min_shared,
            min_anchors=min_anchors,
            high_anchor_threshold=high_anchor_threshold,
        )
    else:
        # Use the generic quantification methods
        click.echo(f"Using {method} quantification method")

        # Load peptide data
        if is_parquet(peptides):
            peptide_df = pd.read_parquet(peptides)
        else:
            peptide_df = pd.read_csv(peptides)

        click.echo(f"Loaded {len(peptide_df)} peptide measurements")

        # Get the quantification method with appropriate parameters
        if method_lower == "topn":
            quant_method = get_quantification_method(method, n=topn_n)
        elif method_lower == "maxlfq":
            quant_method = get_quantification_method(
                method, threads=threads, min_peptides=2
            )
        elif method_lower == "directlfq":
            quant_method = get_quantification_method(method, min_nonan=min_nonan)
        else:
            quant_method = get_quantification_method(method)

        # Determine column names (try to auto-detect)
        protein_col = (
            PROTEIN_NAME if PROTEIN_NAME in peptide_df.columns else "ProteinName"
        )
        sample_col = SAMPLE_ID if SAMPLE_ID in peptide_df.columns else "SampleID"
        intensity_col = (
            NORM_INTENSITY if NORM_INTENSITY in peptide_df.columns else "NormIntensity"
        )
        peptide_col = (
            PEPTIDE_CANONICAL
            if PEPTIDE_CANONICAL in peptide_df.columns
            else "PeptideSequence"
        )

        # Check for required columns
        for col, name in [
            (protein_col, "protein"),
            (sample_col, "sample"),
            (intensity_col, "intensity"),
        ]:
            if col not in peptide_df.columns:
                raise click.UsageError(
                    f"Could not find {name} column '{col}' in peptide file"
                )

        # Run quantification
        click.echo(f"Quantifying {peptide_df[protein_col].nunique()} proteins...")
        result_df = quant_method.quantify(
            peptide_df,
            protein_column=protein_col,
            peptide_column=peptide_col,
            intensity_column=intensity_col,
            sample_column=sample_col,
        )

        # Add condition if available
        if CONDITION in peptide_df.columns:
            condition_map = (
                peptide_df[[sample_col, CONDITION]]
                .drop_duplicates()
                .set_index(sample_col)[CONDITION]
            )
            result_df[CONDITION] = result_df[sample_col].map(condition_map)

        # Normalize if requested
        if normalize:
            # Find the intensity column (last column that contains 'Intensity')
            intensity_cols = [c for c in result_df.columns if "Intensity" in c]
            if intensity_cols:
                intensity_col_out = intensity_cols[-1]
                result_df[f"{intensity_col_out}Norm"] = result_df.groupby(sample_col)[
                    intensity_col_out
                ].transform(lambda x: x / x.sum())

        # Save output
        if output:
            if output.endswith(".parquet"):
                result_df.to_parquet(output, index=False)
            else:
                result_df.to_csv(output, sep="\t", index=False)
            click.echo(f"Results saved to {output}")
        else:
            click.echo(result_df.to_string())
