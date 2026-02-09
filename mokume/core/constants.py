"""
Constants and common utilities for the mokume package.

This module defines column names, mappings, and utility functions used
throughout the package for data processing and analysis.
"""

import os

import pandas as pd


# Column name constants
PROTEIN_NAME = "ProteinName"
PEPTIDE_SEQUENCE = "PeptideSequence"
PEPTIDE_CANONICAL = "PeptideCanonical"
PEPTIDE_CHARGE = "PrecursorCharge"
CHANNEL = "Channel"
CONDITION = "Condition"
BIOREPLICATE = "BioReplicate"
TECHREPLICATE = "TechReplicate"
RUN = "Run"
FRACTION = "Fraction"
INTENSITY = "Intensity"
NORM_INTENSITY = "NormIntensity"
REFERENCE = "Reference"
SAMPLE_ID = "SampleID"
SAMPLE_ID_REGEX = r"^[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*$"
SEARCH_ENGINE = "searchScore"
IBAQ = "Ibaq"
IBAQ_NORMALIZED = "IbaqNorm"
IBAQ_LOG = "IbaqLog"
IBAQ_BEC = "IbaqBec"
IBAQ_PPB = "IbaqPpb"
TPA = "TPA"
MOLECULARWEIGHT = "MolecularWeight"
COPYNUMBER = "CopyNumber"
CONCENTRATION_NM = "Concentration[nM]"
WEIGHT_NG = "Weight[ng]"
MOLES_NMOL = "Moles[nmol]"
GLOBALMEDIAN = "globalMedian"
CONDITIONMEDIAN = "conditionMedian"

# Aggregation level constants
AGGREGATION_LEVEL_SAMPLE = "sample"
AGGREGATION_LEVEL_RUN = "run"


# Parquet column names (QPX compatible)
PARQUET_COLUMNS = [
    "pg_accessions",
    "peptidoform",
    "sequence",
    "precursor_charge",
    "channel",
    "condition",
    "biological_replicate",
    "run",
    "fraction",
    "intensity",
    "reference_file_name",
    "sample_accession",
]


# Mapping from parquet column names to internal column names
parquet_map = {
    "pg_accessions": PROTEIN_NAME,
    "peptidoform": PEPTIDE_SEQUENCE,
    "sequence": PEPTIDE_CANONICAL,
    "precursor_charge": PEPTIDE_CHARGE,
    "channel": CHANNEL,
    "condition": CONDITION,
    "biological_replicate": BIOREPLICATE,
    "run": RUN,
    "fraction": FRACTION,
    "intensity": INTENSITY,
    "reference_file_name": REFERENCE,
    "sample_accession": SAMPLE_ID,
}


def get_accession(identifier: str) -> str:
    """
    Get protein accession from the identifier (e.g. sp|P12345|PROT_NAME).

    Parameters
    ----------
    identifier : str
        Protein identifier.

    Returns
    -------
    str
        Protein accession.
    """
    identifier_lst = identifier.split("|")
    if len(identifier_lst) == 1:
        return identifier_lst[0]
    else:
        return identifier_lst[1]


# Functions needed by Combiner
def load_sdrf(sdrf_path: str) -> pd.DataFrame:
    """
    Load SDRF TSV as a dataframe.

    Parameters
    ----------
    sdrf_path : str
        Path to SDRF TSV.

    Returns
    -------
    pd.DataFrame
        Loaded SDRF data.
    """
    if not os.path.exists(sdrf_path):
        raise FileNotFoundError(f"{sdrf_path} does not exist!")
    sdrf_df = pd.read_csv(sdrf_path, sep="\t")
    sdrf_df.columns = [col.lower() for col in sdrf_df.columns]
    return sdrf_df


def load_feature(feature_path: str) -> pd.DataFrame:
    """
    Load feature file as a dataframe.

    Parameters
    ----------
    feature_path : str
        Path to feature file.

    Returns
    -------
    pd.DataFrame
        Loaded feature data.

    Raises
    ------
    ValueError
        If the provided file's suffix is not supported, either "parquet" or "csv".
    """
    suffix = os.path.splitext(feature_path)[1][1:]
    if suffix == "parquet":
        return pd.read_parquet(feature_path)
    elif suffix == "csv":
        return pd.read_csv(feature_path)
    else:
        raise ValueError(
            f"{suffix} is not allowed as input, please provide msstats_in or feature parquet."
        )


def is_parquet(path: str) -> bool:
    """
    Check if a file is in Parquet format.

    This function attempts to open the specified file and read its header
    to determine if it matches the Parquet file signature.

    Parameters
    ----------
    path : str
        The file path to check.

    Returns
    -------
    bool
        True if the file is a Parquet file, False otherwise.
    """
    try:
        with open(path, "rb") as fh:
            header = fh.read(4)
        return header == b"PAR1"
    except IOError:
        return False
