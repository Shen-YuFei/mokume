import logging
from pathlib import Path

import pandas as pd

from mokume.io.parquet import combine_pibaq_tsv_files, create_anndata
from mokume.core.constants import (
    SAMPLE_ID,
    PROTEIN_NAME,
    PIBAQ,
    PIBAQ_NORMALIZED,
    PIBAQ_LOG,
)

TESTS_DIR = Path(__file__).parent


def test_combine_pibaq_tsv_files():
    """
    Test functions for combining iBAQ TSV files and creating AnnData objects.

    Functions:
    - test_combine_pibaq_tsv_files: Tests the combination of multiple TSV files
      into a single DataFrame and verifies the shape of the resulting DataFrame.
    - test_create_anndata: Tests the creation of an AnnData object from a DataFrame
      with specified observation and variable columns, additional layers, and metadata.
    """
    ibaq_dir = TESTS_DIR / "ibaq-raw-hela"
    files_pattern = "*ibaq.tsv"
    df_pibaq = combine_pibaq_tsv_files(
        dir_path=str(ibaq_dir), pattern=files_pattern, sep="\t"
    )
    logging.info(df_pibaq.head())
    if df_pibaq.shape != (83725, 14):
        raise AssertionError(f"Expected shape (83725, 14), got {df_pibaq.shape}")


def test_create_anndata():
    """
    Test functions for combining iBAQ TSV files and creating AnnData objects.

    Functions:
    - test_combine_pibaq_tsv_files: Tests the combination of multiple TSV files
      into a single DataFrame and verifies the shape of the resulting DataFrame.
    - test_create_anndata: Tests the creation of an AnnData object from a DataFrame
      with specified observation and variable columns, additional layers, and metadata.
    """
    df = pd.read_csv(TESTS_DIR / "ibaq-raw-hela/PXD000396.ibaq.tsv", sep="\t")
    df = df.rename(
        columns={"Ibaq": PIBAQ, "IbaqNorm": PIBAQ_NORMALIZED, "IbaqLog": PIBAQ_LOG}
    )
    obs_col = SAMPLE_ID
    var_col = PROTEIN_NAME
    value_col = PIBAQ
    layers = [PIBAQ_NORMALIZED, PIBAQ_LOG]
    adata = create_anndata(
        df=df,
        obs_col=obs_col,
        var_col=var_col,
        value_col=value_col,
        layer_cols=layers,
        obs_metadata_cols=["Condition"],
        var_metadata_cols=[],
    )
    logging.info(adata)
    if adata.shape != (12, 3096):
        raise AssertionError(f"Expected shape (12, 3096), got {adata.shape}")
    if adata.layers[PIBAQ_NORMALIZED].shape != (12, 3096):
        raise AssertionError(
            f"Expected PIBAQ_NORMALIZED shape (12, 3096), got {adata.layers[PIBAQ_NORMALIZED].shape}"
        )
    if adata.layers[PIBAQ_LOG].shape != (12, 3096):
        raise AssertionError(
            f"Expected PIBAQ_LOG shape (12, 3096), got {adata.layers[PIBAQ_LOG].shape}"
        )
    if "HeLa" not in adata.obs["Condition"].values:
        raise AssertionError("'HeLa' not found in adata.obs['Condition'].values")
