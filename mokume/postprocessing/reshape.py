"""
Data reshaping utilities for the mokume package.
"""

import logging
from typing import Union

import pandas as pd

from mokume.core.constants import (
    IBAQ,
    IBAQ_NORMALIZED,
    IBAQ_PPB,
    IBAQ_LOG,
    TPA,
    COPYNUMBER,
    PROTEIN_NAME,
    SAMPLE_ID,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def remove_samples_low_protein_number(
    ibaq_df: pd.DataFrame, min_protein_num: int
) -> pd.DataFrame:
    """Remove samples with a low number of unique proteins."""
    protein_num = ibaq_df.groupby(SAMPLE_ID)[PROTEIN_NAME].nunique()
    samples_to_keep = protein_num[protein_num >= min_protein_num].index
    samples_to_remove = protein_num[protein_num < min_protein_num].index
    logger.info(
        "The number of samples with number of proteins lower than {} is {}".format(
            min_protein_num, len(samples_to_remove)
        )
    )
    ibaq_df = ibaq_df[ibaq_df["SampleID"].isin(samples_to_keep)]
    return ibaq_df


def remove_missing_values(
    ibaq_df: pd.DataFrame,
    missingness_percentage: float = 30,
    expression_column: str = IBAQ,
) -> pd.DataFrame:
    """Remove samples based on missing values in the expression column."""
    if not isinstance(ibaq_df, pd.DataFrame):
        raise ValueError("The input ibaq_df must be a pandas DataFrame.")
    if expression_column not in ibaq_df.columns:
        raise ValueError(
            f"The expression column '{expression_column}' is not in the DataFrame."
        )

    initial_sample_count = ibaq_df["SampleID"].nunique()
    logger.info(f"Initial number of samples: {initial_sample_count}")

    pivot_df = ibaq_df.pivot_table(
        index=PROTEIN_NAME, columns=SAMPLE_ID, values=expression_column
    )
    non_missing_samples = pivot_df.columns[pivot_df.notna().any(axis=0)]
    missingness = pivot_df[non_missing_samples].isna().sum() / len(pivot_df) * 100
    valid_samples = missingness[missingness <= missingness_percentage].index
    filtered_df = ibaq_df[ibaq_df[SAMPLE_ID].isin(valid_samples)]

    final_sample_count = filtered_df[SAMPLE_ID].nunique()
    logger.info(f"Final number of samples: {final_sample_count}")
    logger.info(
        f"Number of samples removed: {initial_sample_count - final_sample_count}"
    )

    return filtered_df


def describe_expression_metrics(ibaq_df: pd.DataFrame) -> pd.DataFrame:
    """Generate descriptive statistics for expression metrics."""
    possible_expression_values = [
        IBAQ,
        IBAQ_NORMALIZED,
        IBAQ_LOG,
        IBAQ_PPB,
        TPA,
        COPYNUMBER,
    ]
    expression_columns = [
        col for col in ibaq_df.columns if col in possible_expression_values
    ]
    metrics = ibaq_df.groupby(SAMPLE_ID)[expression_columns].describe()
    return metrics


def pivot_wider(
    df: pd.DataFrame,
    row_name: str,
    col_name: str,
    values: str,
    fillna: Union[int, float, bool] = False,
) -> pd.DataFrame:
    """Create a matrix from a DataFrame given the row, column, and value columns."""
    missing_columns = {row_name, col_name, values} - set(df.columns)
    if missing_columns:
        raise ValueError(f"Columns {missing_columns} not found in the DataFrame.")

    duplicates = df.groupby([row_name, col_name]).size()
    if (duplicates > 1).any():
        raise ValueError(
            f"Found duplicate combinations of {row_name} and {col_name}. "
            "Use an aggregation function to handle duplicates."
        )

    matrix = df.pivot_table(
        index=row_name, columns=col_name, values=values, aggfunc="first"
    )

    if fillna is True:
        matrix = matrix.fillna(0)
    elif fillna not in [None, False]:
        matrix = matrix.fillna(fillna)

    return matrix


def pivot_longer(
    df: pd.DataFrame, row_name: str, col_name: str, values: str
) -> pd.DataFrame:
    """Transforms a wide-format DataFrame into a long-format DataFrame."""
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame")
    if row_name not in df.columns:
        raise ValueError(f"Row name '{row_name}' not found in DataFrame")

    matrix_reset = df.reset_index()
    long_df = pd.melt(
        matrix_reset, id_vars=[row_name], var_name=col_name, value_name=values
    )

    if long_df[values].isna().any():
        logging.warning(
            f"Found {long_df[values].isna().sum()} missing values in the result"
        )

    return long_df
