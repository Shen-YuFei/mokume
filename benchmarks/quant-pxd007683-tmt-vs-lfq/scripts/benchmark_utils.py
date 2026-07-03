"""
Shared utility functions for PXD007683 benchmark scripts.

Provides standalone implementations of normalization, imputation, filtering,
and FASTA parsing that were previously part of mokume's internal Python API.
These functions operate on pandas DataFrames / numpy arrays directly.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# Constants
# =============================================================================

ENZYME_CLEAVAGE_RULES = {
    "Trypsin": r"(?<=[KR])(?!P)",
    "Trypsin/P": r"(?<=[KR])",
    "Lys-C": r"(?<=K)",
    "Lys-N": r"(?=K)",
    "Arg-C": r"(?<=R)",
    "Asp-N": r"(?=[BD])",
    "Chymotrypsin": r"(?<=[FYWL])(?!P)",
    "Glu-C": r"(?<=[DE])",
}

MIN_AA = 7
MAX_AA = 50

# =============================================================================
# FASTA Parsing
# =============================================================================

def parse_fasta(fasta_path: str) -> dict:
    """
    Parse a FASTA file and return {accession: sequence}.

    Handles UniProt-style headers: >sp|P12345|NAME_SPECIES ...
    as well as generic headers: >ACCESSION description...

    Parameters
    ----------
    fasta_path : str or Path
        Path to the FASTA file.

    Returns
    -------
    dict
        Mapping from protein accession to amino acid sequence.
    """
    proteins = {}
    current_acc = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_acc is not None:
                    proteins[current_acc] = "".join(current_seq)
                header = line[1:].split()[0]
                # Extract accession from UniProt-style header
                parts = header.split("|")
                if len(parts) >= 2:
                    current_acc = parts[1]
                else:
                    current_acc = header
                current_seq = []
            else:
                current_seq.append(line)

    if current_acc is not None:
        proteins[current_acc] = "".join(current_seq)

    return proteins

def count_theoretical_peptides(
    sequence: str,
    enzyme: str = "Trypsin",
    min_aa: int = MIN_AA,
    max_aa: int = MAX_AA,
) -> int:
    """
    Count the number of theoretical tryptic peptides for a protein sequence.

    Parameters
    ----------
    sequence : str
        Amino acid sequence.
    enzyme : str
        Enzyme name (must be in ENZYME_CLEAVAGE_RULES).
    min_aa : int
        Minimum peptide length.
    max_aa : int
        Maximum peptide length.

    Returns
    -------
    int
        Number of theoretical peptides within the length range.
    """
    pattern = ENZYME_CLEAVAGE_RULES.get(enzyme)
    if pattern is None:
        raise ValueError(f"Unknown enzyme: {enzyme}")

    peptides = re.split(pattern, sequence)
    count = sum(1 for p in peptides if min_aa <= len(p) <= max_aa)
    return count

def get_accession(protein_name: str) -> str:
    """
    Extract UniProt accession from a protein name string.

    Handles formats like:
    - sp|P12345|NAME_SPECIES
    - tr|A0A0B4J2F0|A0A0B4J2F0_HUMAN
    - P12345

    Parameters
    ----------
    protein_name : str
        Protein name or header.

    Returns
    -------
    str or None
        Extracted accession, or None if not parseable.
    """
    parts = protein_name.split("|")
    if len(parts) >= 2:
        return parts[1]
    name = protein_name.strip()
    if name and not name.startswith(">"):
        return name
    return None

# =============================================================================
# Normalization
# =============================================================================

def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Quantile normalization: rank-based, forces identical distributions.

    Algorithm:
    1. Sort each column independently.
    2. Compute row means of the sorted matrix.
    3. Reassign values by original rank order.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame (proteins x samples).

    Returns
    -------
    pd.DataFrame
        Quantile-normalized DataFrame.
    """
    rank_mean = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    result = df.rank(method="min").stack().astype(int).map(rank_mean).unstack()
    return result

def median_center(df: pd.DataFrame, axis: int = 0) -> pd.DataFrame:
    """
    Median centering: subtract the column (or row) median.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame.
    axis : int
        0 for column-wise centering, 1 for row-wise.

    Returns
    -------
    pd.DataFrame
        Median-centered DataFrame.
    """
    return df.subtract(df.median(axis=axis), axis=1 - axis)

def total_intensity_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Total intensity normalization: scale each sample to the median total.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame (proteins x samples).

    Returns
    -------
    pd.DataFrame
        Normalized DataFrame.
    """
    col_sums = df.sum()
    target = col_sums.median()
    return df * (target / col_sums)

def global_median_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Global median normalization: scale each sample's median to the global median.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame (proteins x samples).

    Returns
    -------
    pd.DataFrame
        Normalized DataFrame.
    """
    global_median = df.median().median()
    sample_medians = df.median()
    factors = global_median / sample_medians
    return df * factors

def hierarchical_normalize(df_log: pd.DataFrame, min_overlap: int = 10) -> pd.DataFrame:
    """
    Hierarchical sample normalization.

    For each sample pair with sufficient overlap, compute median log-ratio
    to a reference sample (the one with the most non-NaN values), then scale.

    Parameters
    ----------
    df_log : pd.DataFrame
        Log-transformed numeric DataFrame (proteins x samples).
    min_overlap : int
        Minimum number of shared non-NaN values between sample pairs.

    Returns
    -------
    pd.DataFrame
        Normalized log-transformed DataFrame.
    """
    result = df_log.copy()
    # Choose reference as sample with most non-NaN values
    non_nan_counts = result.notna().sum()
    ref_col = non_nan_counts.idxmax()

    for col in result.columns:
        if col == ref_col:
            continue
        # Find overlapping proteins
        mask = result[ref_col].notna() & result[col].notna()
        if mask.sum() < min_overlap:
            continue
        # Median log-ratio
        shift = (result.loc[mask, ref_col] - result.loc[mask, col]).median()
        result[col] = result[col] + shift

    return result

def impute_values(df: pd.DataFrame, method: str, **kwargs) -> pd.DataFrame:
    """
    Impute missing values in a numeric DataFrame.

    Supported methods:
    - "mean": column-wise mean imputation
    - "median": column-wise median imputation
    - "knn": K-Nearest Neighbors imputation
    - "constant": fill with a constant value (pass fill_value=...)
    - "missforest": iterative random forest imputation (via sklearn IterativeImputer)

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame with NaN for missing values.
    method : str
        Imputation method name.
    **kwargs
        Additional arguments (e.g., n_neighbors for knn, fill_value for constant).

    Returns
    -------
    pd.DataFrame
        Imputed DataFrame with no NaN values.
    """
    if method in ("mean", "median"):
        from sklearn.impute import SimpleImputer
        imputer = SimpleImputer(strategy=method)
        result = pd.DataFrame(
            imputer.fit_transform(df),
            index=df.index,
            columns=df.columns,
        )
        return result

    elif method == "knn":
        from sklearn.impute import KNNImputer
        n_neighbors = kwargs.get("n_neighbors", 5)
        imputer = KNNImputer(n_neighbors=n_neighbors)
        result = pd.DataFrame(
            imputer.fit_transform(df),
            index=df.index,
            columns=df.columns,
        )
        return result

    elif method == "constant":
        fill_value = kwargs.get("fill_value", 0.0)
        return df.fillna(fill_value)

    elif method == "missforest":
        from sklearn.experimental import enable_iterative_imputer  # noqa: F401
        from sklearn.impute import IterativeImputer
        from sklearn.ensemble import RandomForestRegressor
        max_iter = kwargs.get("max_iter", 10)
        imputer = IterativeImputer(
            estimator=RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1),
            max_iter=max_iter,
            random_state=42,
        )
        result = pd.DataFrame(
            imputer.fit_transform(df),
            index=df.index,
            columns=df.columns,
        )
        return result

    else:
        raise ValueError(f"Unknown imputation method: {method}")

# =============================================================================
# Preprocessing Filters
# =============================================================================

def filter_peptide_length(
    df: pd.DataFrame,
    min_length: int = 7,
    max_length: int = 50,
    sequence_column: str = "PeptideSequence",
) -> pd.DataFrame:
    """
    Filter peptides by sequence length.

    Parameters
    ----------
    df : pd.DataFrame
        Peptide-level DataFrame.
    min_length : int
        Minimum peptide length (inclusive).
    max_length : int
        Maximum peptide length (inclusive).
    sequence_column : str
        Name of the column containing peptide sequences.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    lengths = df[sequence_column].str.len()
    mask = (lengths >= min_length) & (lengths <= max_length)
    return df[mask].copy()

def filter_missed_cleavages(
    df: pd.DataFrame,
    max_missed_cleavages: int = 2,
    sequence_column: str = "PeptideSequence",
    enzyme: str = "trypsin",
) -> pd.DataFrame:
    """
    Filter peptides by number of missed cleavage sites.

    For trypsin: counts K or R not followed by P (internal to the peptide).

    Parameters
    ----------
    df : pd.DataFrame
        Peptide-level DataFrame.
    max_missed_cleavages : int
        Maximum allowed missed cleavages.
    sequence_column : str
        Name of the column containing peptide sequences.
    enzyme : str
        Enzyme name (currently only "trypsin" is supported).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    def count_missed(seq: str) -> int:
        if enzyme.lower() in ("trypsin",):
            # Count K or R not followed by P, excluding the last residue
            internal = seq[:-1] if len(seq) > 1 else ""
            count = 0
            for i, aa in enumerate(internal):
                if aa in ("K", "R"):
                    # Check if followed by P
                    if i + 1 < len(internal) and internal[i + 1] == "P":
                        continue
                    count += 1
            return count
        else:
            return 0

    mc_counts = df[sequence_column].apply(count_missed)
    return df[mc_counts <= max_missed_cleavages].copy()

def filter_min_peptides(
    df: pd.DataFrame,
    min_unique_peptides: int = 2,
    protein_column: str = "ProteinName",
    peptide_column: str = "PeptideSequence",
) -> pd.DataFrame:
    """
    Filter proteins by minimum number of unique peptides.

    Parameters
    ----------
    df : pd.DataFrame
        Peptide-level DataFrame.
    min_unique_peptides : int
        Minimum number of unique peptides per protein.
    protein_column : str
        Name of the protein column.
    peptide_column : str
        Name of the peptide sequence column.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    peptide_counts = df.groupby(protein_column)[peptide_column].nunique()
    passing_proteins = peptide_counts[peptide_counts >= min_unique_peptides].index
    return df[df[protein_column].isin(passing_proteins)].copy()

def filter_min_cv(
    df_wide: pd.DataFrame,
    max_cv: float,
    condition_map: dict,
) -> pd.DataFrame:
    """
    Filter proteins by maximum CV within conditions.

    Parameters
    ----------
    df_wide : pd.DataFrame
        Wide-format DataFrame (proteins x samples).
    max_cv : float
        Maximum allowed CV.
    condition_map : dict
        Mapping from sample names to condition labels.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    keep = pd.Series(True, index=df_wide.index)

    conditions = set(condition_map.values())
    for cond in conditions:
        cols = [c for c in df_wide.columns if condition_map.get(str(c)) == cond]
        if len(cols) < 2:
            continue
        subset = df_wide[cols]
        cv = subset.std(axis=1) / subset.mean(axis=1)
        cv = cv.replace([np.inf, -np.inf], np.nan).fillna(1.0)
        keep = keep & (cv <= max_cv)

    return df_wide[keep].copy()

def filter_min_intensity(
    df: pd.DataFrame,
    min_intensity: float,
    intensity_column: str = "Intensity",
) -> pd.DataFrame:
    """
    Filter rows by minimum intensity value.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with an intensity column.
    min_intensity : float
        Minimum intensity threshold.
    intensity_column : str
        Name of the intensity column.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    return df[df[intensity_column] >= min_intensity].copy()

def filter_charge_state(
    df: pd.DataFrame,
    allowed_charges: set,
    charge_column: str = "Charge",
) -> pd.DataFrame:
    """
    Filter rows by allowed charge states.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a charge column.
    allowed_charges : set
        Set of allowed charge state integers.
    charge_column : str
        Name of the charge column.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    return df[df[charge_column].isin(allowed_charges)].copy()

# =============================================================================
# Quantification Helpers
# =============================================================================

def pivot_wider(
    df: pd.DataFrame,
    row_name: str,
    col_name: str,
    values: str,
    fillna: bool = False,
) -> pd.DataFrame:
    """
    Pivot a long-format DataFrame to wide format.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame.
    row_name : str
        Column to use as row index.
    col_name : str
        Column to use as column headers.
    values : str
        Column containing values.
    fillna : bool
        Whether to fill NaN with 0.

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame.
    """
    result = df.pivot_table(index=row_name, columns=col_name, values=values, aggfunc="first")
    if fillna:
        result = result.fillna(0)
    return result
