"""
Parquet and file I/O utilities for the mokume package.
"""

import glob
import logging
import warnings
from typing import List, Optional, TYPE_CHECKING

import pandas as pd

from mokume.postprocessing.reshape import pivot_wider

if TYPE_CHECKING:
    import anndata as an

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def create_anndata(
    df: pd.DataFrame,
    obs_col: str,
    var_col: str,
    value_col: str,
    layer_cols: Optional[List[str]] = None,
    obs_metadata_cols: Optional[List[str]] = None,
    var_metadata_cols: Optional[List[str]] = None,
) -> "an.AnnData":
    """
    Create an AnnData object from a long-format DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input data in long format.
    obs_col : str
        Column name representing observation IDs.
    var_col : str
        Column name representing variable IDs.
    value_col : str
        Column name representing the main data values.
    layer_cols : Optional[List[str]]
        List of column names to add as additional layers.
    obs_metadata_cols : Optional[List[str]]
        List of column names to add as observation metadata.
    var_metadata_cols : Optional[List[str]]
        List of column names to add as variable metadata.

    Returns
    -------
    anndata.AnnData
        The constructed AnnData object.
    """
    try:
        import anndata as an
    except ImportError as exc:
        raise ImportError(
            "anndata is required for AnnData/.h5ad export. "
            "Install it with: pip install mokume-py[h5ad]"
        ) from exc

    if df.empty:
        raise ValueError("Cannot create AnnData object from empty DataFrame")

    required_cols = [obs_col, var_col, value_col]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"The following required columns are missing: {missing}")

    df_matrix = pivot_wider(
        df, row_name=obs_col, col_name=var_col, values=value_col, fillna=True
    )
    if df_matrix.empty:
        raise ValueError("Pivot operation resulted in an empty DataFrame")
    if df_matrix.shape[0] == 0 or df_matrix.shape[1] == 0:
        raise ValueError("Pivot operation resulted in a DataFrame with zero dimensions")

    adata = an.AnnData(
        X=df_matrix.to_numpy(),
        obs=df_matrix.index.to_frame(),
        var=df_matrix.columns.to_frame(),
    )

    def add_metadata(
        metadata_df: pd.DataFrame, key: str, cols: List[str]
    ) -> pd.DataFrame:
        for col in cols:
            if col not in df.columns:
                warnings.warn(f"Column '{col}' not found. Skipping.")
                continue
            mapping = df[[key, col]].drop_duplicates().set_index(key)[col]
            metadata_df[col] = metadata_df.index.map(mapping)
        return metadata_df

    if obs_metadata_cols:
        adata.obs = add_metadata(adata.obs, obs_col, obs_metadata_cols)

    if var_metadata_cols:
        adata.var = add_metadata(adata.var, var_col, var_metadata_cols)

    if layer_cols:
        for layer_col in layer_cols:
            if layer_col not in df.columns:
                warnings.warn(f"Layer column '{layer_col}' not found. Skipping.")
                continue
            df_layer = pivot_wider(
                df, row_name=obs_col, col_name=var_col, values=layer_col, fillna=True
            )
            adata.layers[layer_col] = df_layer.to_numpy()

    logger.info("Created AnnData object:\n %s", adata)
    return adata


def combine_pibaq_tsv_files(
    dir_path: str, pattern: str = "*", comment: str = "#", sep: str = "\t"
) -> pd.DataFrame:
    """
    Combine multiple TSV files from a directory into a single DataFrame.

    Parameters
    ----------
    dir_path : str
        Directory path containing the TSV files.
    pattern : str, optional
        Pattern to match files in the directory. Default is '*'.
    comment : str, optional
        Character indicating the start of a comment line. Default is '#'.
    sep : str, optional
        Delimiter to use for reading the TSV files. Default is '\t'.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing data from all TSV files.
    """
    file_paths = glob.glob(f"{dir_path}/{pattern}")

    if not file_paths:
        raise FileNotFoundError(
            f"No files found in the directory '{dir_path}' matching the pattern '{pattern}'."
        )

    dataframes = []
    first_schema = None

    for file_path in file_paths:
        try:
            df = pd.read_csv(file_path, sep=sep, comment=comment)

            if first_schema is None:
                first_schema = set(df.columns)
            elif set(df.columns) != first_schema:
                raise ValueError(
                    f"Schema mismatch in file '{file_path}'. "
                    f"Expected columns: {sorted(first_schema)}, "
                    f"got: {sorted(df.columns)}"
                )

            dataframes.append(df)
        except Exception as e:
            raise ValueError(f"Error reading file '{file_path}': {str(e)}") from e

    combined_df = pd.concat(dataframes, ignore_index=True)
    return combined_df
