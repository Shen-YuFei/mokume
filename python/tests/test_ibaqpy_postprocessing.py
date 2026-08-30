from pathlib import Path
import pandas as pd

from mokume.core.constants import PIBAQ_NORMALIZED, SAMPLE_ID
from mokume.postprocessing.reshape import (
    remove_samples_low_protein_number,
    remove_missing_values,
    describe_expression_metrics,
)
import logging

TESTS_DIR = Path(__file__).parent


def _load_pibaq_fixture() -> pd.DataFrame:
    """Load the historical iBAQ fixture under the canonical piBAQ schema."""
    return pd.read_csv(
        TESTS_DIR / "example/PXD017834-example-ibaq.tsv", sep="\t"
    ).rename(columns={"IbaqNorm": PIBAQ_NORMALIZED})


def test_remove_samples_low_protein_number():
    """
    Test functions for post-processing iBAQ data.

    These tests validate the functionality of the following operations:
    - Removing samples with a low number of proteins.
    - Removing samples with a high percentage of missing values.
    - Describing expression metrics across samples.

    Each test reads a sample iBAQ dataset, applies the respective function,
    and logs the number of samples before and after processing.
    """
    pibaq_df = _load_pibaq_fixture()
    number_samples = len(pibaq_df[SAMPLE_ID].unique())
    logging.info("The number of samples in the dataframe {}".format(number_samples))

    new_pibaq = remove_samples_low_protein_number(pibaq_df, min_protein_num=286)

    number_samples = len(new_pibaq[SAMPLE_ID].unique())
    logging.info(
        "The number of samples with number of proteins higher than 286 is {}".format(
            number_samples
        )
    )


def test_remove_missing_values():
    """
    Test functions for post-processing iBAQ data.

    These tests validate the functionality of the following operations:
    - Removing samples with a low number of proteins.
    - Removing samples with a high percentage of missing values.
    - Describing expression metrics across samples.

    Each test reads a sample iBAQ dataset, applies the respective function,
    and logs the number of samples before and after processing.
    """
    pibaq_df = _load_pibaq_fixture()
    number_samples = len(pibaq_df[SAMPLE_ID].unique())
    logging.info("The number of samples in the dataframe {}".format(number_samples))
    new_pibaq = remove_missing_values(
        pibaq_df, missingness_percentage=1, expression_column=PIBAQ_NORMALIZED
    )
    number_samples = len(new_pibaq[SAMPLE_ID].unique())
    logging.info(
        "The number of samples with less than 1% of missing values is {}".format(
            number_samples
        )
    )


def test_describe_expression_metrics():
    """
    Test functions for post-processing iBAQ data.

    These tests validate the functionality of the following operations:
    - Removing samples with a low number of proteins.
    - Removing samples with a high percentage of missing values.
    - Describing expression metrics across samples.

    Each test reads a sample iBAQ dataset, applies the respective function,
    and logs the number of samples before and after processing.
    """
    pibaq_df = _load_pibaq_fixture()

    metrics = describe_expression_metrics(pibaq_df)
    logging.info(metrics)
