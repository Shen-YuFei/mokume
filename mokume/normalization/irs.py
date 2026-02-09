"""
Internal Reference Scaling (IRS) normalization for multi-plex TMT experiments.

IRS corrects between-plex batch effects by scaling protein intensities so that
reference channels (pooled samples present across all plexes) have equal
intensity across plexes.

References
----------
Plubell et al. (2017). Extended Multiplexing of TMT Labeling Reveals Age and
High Fat Diet Specific Proteome Changes in Mouse Epididymal Adipose Tissue.
Molecular & Cellular Proteomics, 16(5), 873-890.
"""

import re
from typing import Optional

import numpy as np
import pandas as pd

from mokume.core.constants import load_sdrf
from mokume.core.logger import get_logger

logger = get_logger("mokume.normalization.irs")


class IRSNormalizer:
    """
    Internal Reference Scaling normalizer for multi-plex TMT data.

    For each protein, IRS computes the reference intensity per plex (using
    reference channels), then scales all samples in each plex so that
    reference intensities are equalized across plexes.

    Parameters
    ----------
    reference_samples : list[str]
        Sample names (source names) that are reference channels.
    stat : str
        Statistic to compute reference intensity: "median" or "mean".
    """

    def __init__(self, reference_samples: list[str], stat: str = "median"):
        self.reference_samples = reference_samples
        self.stat = stat
        self.scaling_factors_ = None

    def fit(self, protein_df: pd.DataFrame, sample_to_plex: dict[str, str]):
        """
        Compute per-protein, per-plex IRS scaling factors.

        Parameters
        ----------
        protein_df : pd.DataFrame
            Wide-format protein matrix with a protein column and sample columns.
        sample_to_plex : dict
            Mapping from sample name to plex identifier.
        """
        protein_col = protein_df.columns[0]
        sample_cols = [c for c in protein_df.columns if c != protein_col]

        # Build plex->reference_samples mapping
        plexes = sorted(set(sample_to_plex.values()))
        plex_ref_samples = {}
        for plex in plexes:
            refs = [s for s in self.reference_samples if sample_to_plex.get(s) == plex]
            if not refs:
                raise ValueError(
                    f"No reference samples found in plex '{plex}'. "
                    f"Reference samples: {self.reference_samples}, "
                    f"Plex samples: {[s for s, p in sample_to_plex.items() if p == plex]}"
                )
            plex_ref_samples[plex] = refs

        logger.info(
            f"IRS fit: {len(plexes)} plexes, "
            f"{len(self.reference_samples)} reference samples, stat={self.stat}"
        )
        for plex, refs in plex_ref_samples.items():
            logger.info(f"  Plex '{plex}': reference channels = {refs}")

        # Compute reference intensity per protein per plex
        intensity_matrix = protein_df.set_index(protein_col)[sample_cols]
        plex_ref_intensity = {}

        for plex, refs in plex_ref_samples.items():
            ref_cols = [r for r in refs if r in intensity_matrix.columns]
            if not ref_cols:
                raise ValueError(f"Reference samples {refs} not found in protein matrix columns")
            ref_data = intensity_matrix[ref_cols]
            if self.stat == "median":
                plex_ref_intensity[plex] = ref_data.median(axis=1)
            else:
                plex_ref_intensity[plex] = ref_data.mean(axis=1)

        # Compute global reference intensity (geometric mean across plexes)
        ref_df = pd.DataFrame(plex_ref_intensity)
        # Replace 0 and negatives with NaN for geometric mean
        ref_df = ref_df.replace(0, np.nan)
        log_refs = np.log(ref_df)
        global_ref = np.exp(log_refs.mean(axis=1))

        # Compute scaling factors: global_ref / plex_ref for each plex
        scaling_factors = {}
        for plex in plexes:
            scaling_factors[plex] = global_ref / plex_ref_intensity[plex]

        self.scaling_factors_ = scaling_factors
        self._sample_to_plex = sample_to_plex
        self._protein_col = protein_col

        return self

    def transform(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply IRS scaling to protein intensities.

        Parameters
        ----------
        protein_df : pd.DataFrame
            Wide-format protein matrix.

        Returns
        -------
        pd.DataFrame
            IRS-normalized protein matrix.
        """
        if self.scaling_factors_ is None:
            raise RuntimeError("Must call fit() before transform()")

        protein_col = self._protein_col
        result = protein_df.copy()
        result = result.set_index(protein_col)

        for sample, plex in self._sample_to_plex.items():
            if sample in result.columns and plex in self.scaling_factors_:
                result[sample] = result[sample] * self.scaling_factors_[plex]

        return result.reset_index()

    def fit_transform(
        self, protein_df: pd.DataFrame, sample_to_plex: dict[str, str]
    ) -> pd.DataFrame:
        """Fit and transform in one step."""
        return self.fit(protein_df, sample_to_plex).transform(protein_df)


# ---------------------------------------------------------------------------
# SDRF helper functions for reference channel and plex detection
# ---------------------------------------------------------------------------


def detect_pooled_from_sdrf(sdrf_path: str) -> Optional[list[str]]:
    """
    Detect reference samples from SDRF ``characteristics[pooled sample]`` column.

    Per the SDRF-Proteomics specification, samples with value ``pooled`` or
    ``SN=sample1;SN=sample2;...`` are reference channels.

    Parameters
    ----------
    sdrf_path : str
        Path to the SDRF TSV file.

    Returns
    -------
    list[str] or None
        List of reference sample source names, or None if column not present.
    """
    sdrf = load_sdrf(sdrf_path)

    pooled_col = "characteristics[pooled sample]"
    if pooled_col not in sdrf.columns:
        return None

    # Pooled if value is "pooled" or starts with "SN=" (not "not pooled")
    lowered = sdrf[pooled_col].str.lower()
    mask = (lowered == "pooled") | (lowered.str.startswith("sn=", na=False))
    ref_samples = sdrf.loc[mask, "source name"].unique().tolist()

    if ref_samples:
        logger.info(
            f"Detected {len(ref_samples)} pooled reference samples from SDRF "
            f"'{pooled_col}' column: {ref_samples}"
        )
        return ref_samples

    return None


def detect_reference_by_column(
    sdrf_path: str, column: str, values: list[str]
) -> list[str]:
    """
    Detect reference samples by matching values in a specific SDRF column.

    Parameters
    ----------
    sdrf_path : str
        Path to the SDRF TSV file.
    column : str
        SDRF column name (e.g., ``"factor value[disease]"``).
    values : list[str]
        Values that indicate a reference sample.

    Returns
    -------
    list[str]
        List of reference sample source names.
    """
    sdrf = load_sdrf(sdrf_path)
    column_lower = column.lower()

    if column_lower not in sdrf.columns:
        raise ValueError(f"Column '{column}' not found in SDRF. Available: {list(sdrf.columns)}")

    values_lower = [v.lower() for v in values]
    mask = sdrf[column_lower].str.lower().isin(values_lower)
    ref_samples = sdrf.loc[mask, "source name"].unique().tolist()

    logger.info(
        f"Detected {len(ref_samples)} reference samples from column "
        f"'{column}' with values {values}: {ref_samples}"
    )
    return ref_samples


def detect_reference_by_regex(
    sdrf_path: str, regex: str = r"pool|powder|ref|reference|bridge"
) -> list[str]:
    """
    Detect reference samples by regex scan across all factor/characteristic columns.

    Parameters
    ----------
    sdrf_path : str
        Path to the SDRF TSV file.
    regex : str
        Regex pattern to match reference sample values.

    Returns
    -------
    list[str]
        List of reference sample source names.
    """
    sdrf = load_sdrf(sdrf_path)

    # Scan all factor value[*] and characteristics[*] columns
    scan_cols = [
        c for c in sdrf.columns
        if c.startswith("factor value[") or c.startswith("characteristics[")
    ]

    pattern = re.compile(regex, re.IGNORECASE)
    ref_mask = pd.Series(False, index=sdrf.index)

    for col in scan_cols:
        col_match = sdrf[col].astype(str).str.contains(pattern, na=False)
        ref_mask = ref_mask | col_match

    ref_samples = sdrf.loc[ref_mask, "source name"].unique().tolist()

    logger.info(
        f"Detected {len(ref_samples)} reference samples by regex '{regex}' "
        f"across {len(scan_cols)} columns: {ref_samples}"
    )
    return ref_samples


def detect_plexes_from_sdrf(sdrf_path: str) -> dict[str, str]:
    """
    Detect plex membership from SDRF ``source name`` prefix convention.

    Parses quantms-style source names (e.g., ``p1_1`` -> plex ``p1``,
    ``p2_1`` -> plex ``p2``).

    Parameters
    ----------
    sdrf_path : str
        Path to the SDRF TSV file.

    Returns
    -------
    dict[str, str]
        Mapping from sample name (source name) to plex identifier.
    """
    sdrf = load_sdrf(sdrf_path)

    sample_to_plex = {}
    for source_name in sdrf["source name"].unique():
        # Parse prefix: p1_1 -> p1, p2_10 -> p2
        match = re.match(r"^(p\d+)_\d+$", source_name)
        if match:
            sample_to_plex[source_name] = match.group(1)
        else:
            # Fallback: use everything before the last underscore+digit
            parts = source_name.rsplit("_", 1)
            if len(parts) == 2 and parts[1].isdigit():
                sample_to_plex[source_name] = parts[0]
            else:
                sample_to_plex[source_name] = "plex1"
                logger.warning(
                    f"Could not parse plex from source name '{source_name}', "
                    "assigning to 'plex1'"
                )

    plexes = sorted(set(sample_to_plex.values()))
    logger.info(f"Detected {len(plexes)} plexes: {plexes}")
    for plex in plexes:
        members = [s for s, p in sample_to_plex.items() if p == plex]
        logger.info(f"  Plex '{plex}': {sorted(members)}")

    return sample_to_plex


def detect_condition_from_sdrf(sdrf_path: str) -> dict[str, str]:
    """
    Detect sample-to-condition mapping from SDRF factor value columns.

    Looks for the first ``factor value[*]`` column and maps each source name
    to its condition.

    Parameters
    ----------
    sdrf_path : str
        Path to the SDRF TSV file.

    Returns
    -------
    dict[str, str]
        Mapping from sample name (source name) to condition.
    """
    sdrf = load_sdrf(sdrf_path)

    # Find factor value columns
    factor_cols = [c for c in sdrf.columns if c.startswith("factor value[")]
    if not factor_cols:
        raise ValueError("No 'factor value[*]' columns found in SDRF")

    factor_col = factor_cols[0]
    sample_to_condition = {}
    for _, row in sdrf.drop_duplicates(subset=["source name"]).iterrows():
        sample_to_condition[row["source name"]] = str(row[factor_col])

    logger.info(f"Condition mapping from '{factor_col}': {len(sample_to_condition)} samples")
    return sample_to_condition
