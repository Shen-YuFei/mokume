"""
CLI command for batch effect correction.
"""

import logging
import re
from pathlib import Path
from typing import Union

import click
import pandas as pd

from mokume.io.parquet import create_anndata, combine_ibaq_tsv_files
from mokume.core.constants import (
    SAMPLE_ID_REGEX,
    SAMPLE_ID,
    PROTEIN_NAME,
    IBAQ,
    IBAQ_BEC,
)
from mokume.postprocessing.reshape import pivot_wider, pivot_longer
from mokume.postprocessing.batch_correction import apply_batch_correction


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def is_valid_sample_id(
    samples: Union[str, list, pd.Series], sample_id_pattern: str = SAMPLE_ID_REGEX
) -> bool:
    """Validate sample IDs against a specified pattern."""
    sample_pattern = re.compile(sample_id_pattern)

    if isinstance(samples, str):
        samples = [samples]
    elif isinstance(samples, pd.Series):
        samples = samples.tolist()

    invalid_samples = [
        sample for sample in samples if not sample_pattern.fullmatch(sample)
    ]

    if invalid_samples:
        logger.error("The following sample IDs are invalid:")
        for invalid_sample in invalid_samples:
            logger.error(" - %s", invalid_sample)
        return False
    return True


def get_batch_id_from_sample_names(samples: list) -> list:
    """Extract batch IDs from a list of sample names."""
    batch_ids = []
    for sample in samples:
        parts = sample.split("-")
        if not parts or not parts[0]:
            raise ValueError(f"Invalid sample name format: {sample}")
        batch_id = parts[0]
        if not re.match(r"^[A-Za-z0-9]+$", batch_id):
            raise ValueError(f"Invalid batch ID format: {batch_id}")
        batch_ids.append(batch_id)
    return pd.factorize(batch_ids)[0]


def run_batch_correction(
    folder: str,
    pattern: str,
    comment: str,
    sep: str,
    output: str,
    sample_id_column: str = SAMPLE_ID,
    protein_id_column: str = PROTEIN_NAME,
    ibaq_raw_column: str = IBAQ,
    ibaq_corrected_column: str = IBAQ_BEC,
    export_anndata: bool = False,
) -> pd.DataFrame:
    """Run batch correction on iBAQ data from TSV files."""
    logger.info("Loading iBAQ data from TSV files in folder '%s'", folder)

    try:
        df_ibaq = combine_ibaq_tsv_files(
            folder, pattern=pattern, comment=comment, sep=sep
        )
    except Exception as e:
        raise ValueError(f"Failed to load input files: {str(e)}") from e

    df_wide = pivot_wider(
        df_ibaq,
        row_name=protein_id_column,
        col_name=sample_id_column,
        values=ibaq_raw_column,
        fillna=True,
    )

    if not is_valid_sample_id(df_wide.columns, SAMPLE_ID_REGEX):
        raise ValueError("Invalid sample IDs found in the data.")

    batch_ids = get_batch_id_from_sample_names(df_wide.columns)

    logger.info("Applying batch correction to iBAQ values")
    df_corrected = apply_batch_correction(df_wide, list(batch_ids), kwargs={})

    df_corrected = df_corrected.reset_index()
    df_corrected_long = pivot_longer(
        df_corrected,
        row_name=protein_id_column,
        col_name=sample_id_column,
        values=ibaq_corrected_column,
    )

    df_ibaq = df_ibaq.merge(
        df_corrected_long, how="left", on=[sample_id_column, protein_id_column]
    )

    if output:
        try:
            df_ibaq.to_csv(output, sep=sep, index=False)
        except Exception as e:
            raise ValueError(f"Failed to save output file: {str(e)}") from e

    if export_anndata:
        logger.info("Exporting raw and corrected iBAQ values to an AnnData object")
        output_path = Path(output)
        if not output_path.exists():
            raise FileNotFoundError(f"Output file {output} does not exist!")
        adata = create_anndata(
            df_ibaq,
            obs_col=sample_id_column,
            var_col=protein_id_column,
            value_col=ibaq_raw_column,
            layer_cols=[ibaq_corrected_column],
        )
        adata_filename = output_path.with_suffix(".h5ad")
        try:
            adata.write(adata_filename)
        except Exception as e:
            raise ValueError(f"Failed to write AnnData object: {e}") from e

    logger.info("Batch correction completed...")
    return df_ibaq


@click.command("correct-batches", short_help="Batch effect correction for iBAQ values.")
@click.option(
    "-f",
    "--folder",
    help="Folder that contains all TSV files with raw iBAQ values",
    required=True,
    default=None,
)
@click.option(
    "-p",
    "--pattern",
    help="Pattern for the TSV files with raw iBAQ values",
    required=True,
    default="*ibaq.tsv",
)
@click.option(
    "--comment",
    help="Comment character for the TSV files",
    required=False,
    default="#",
)
@click.option("--sep", help="Separator for the TSV files", required=False, default="\t")
@click.option(
    "-o",
    "--output",
    help="Output file name for the combined iBAQ corrected values",
    required=True,
)
@click.option(
    "-sid",
    "--sample_id_column",
    help="Sample ID column name",
    required=False,
    default=SAMPLE_ID,
)
@click.option(
    "-pid",
    "--protein_id_column",
    help="Protein ID column name",
    required=False,
    default=PROTEIN_NAME,
)
@click.option(
    "-ibaq",
    "--ibaq_raw_column",
    help="Name of the raw iBAQ column",
    required=False,
    default=IBAQ,
)
@click.option(
    "--ibaq_corrected_column",
    help="Name for the corrected iBAQ column",
    required=False,
    default=IBAQ_BEC,
)
@click.option(
    "--export_anndata",
    help="Export the raw and corrected iBAQ values to an AnnData object",
    is_flag=True,
)
@click.pass_context
def correct_batches(
    _ctx,
    folder: str,
    pattern: str,
    comment: str,
    sep: str,
    output: str,
    sample_id_column: str,
    protein_id_column: str,
    ibaq_raw_column: str,
    ibaq_corrected_column: str,
    export_anndata: bool,
):
    """Correcting batch effects in iBAQ data."""
    run_batch_correction(
        folder=folder,
        pattern=pattern,
        comment=comment,
        sep=sep,
        output=output,
        sample_id_column=sample_id_column,
        protein_id_column=protein_id_column,
        ibaq_raw_column=ibaq_raw_column,
        ibaq_corrected_column=ibaq_corrected_column,
        export_anndata=export_anndata,
    )
