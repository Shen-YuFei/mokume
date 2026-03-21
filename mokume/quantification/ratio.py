"""
Ratio-based quantification for TMT experiments.

Implements the Proteome Sciences (PS) post-processing protocol where
quantification operates in **ratio space**: PSM intensity / reference
intensity -> log2, with median aggregation up to protein level.

This is a self-contained pipeline (like DirectLFQ) that does NOT go
through the standard _load_and_process_peptides() -> _quantify() flow.

References
----------
Proteome Sciences post-processing protocol for TMT data.
"""

import numpy as np
import pandas as pd
import duckdb

from mokume.core.constants import PROTEIN_NAME, PEPTIDE_CANONICAL, PEPTIDE_CHARGE, SAMPLE_ID
from mokume.core.logger import get_logger
from mokume.io.feature import SQLFilterBuilder

logger = get_logger("mokume.quantification.ratio")


class RatioQuantification:
    """
    Ratio-based TMT quantification following the PS protocol.

    Computes log2(sample / reference) per PSM per plex, then aggregates
    PSM -> peptide -> protein using median. The output is a wide protein
    x sample matrix of log2 ratios (already in log2 space).

    Parameters
    ----------
    reference_samples : list[str]
        Sample names (source names) that are reference channels.
    sample_to_plex : dict[str, str]
        Mapping from sample name to plex identifier.
    fraction_merge_method : str
        How to merge fractions: "mean" (PS protocol) or "max" (mokume default).
    """

    def __init__(
        self,
        reference_samples: list[str],
        sample_to_plex: dict[str, str],
        fraction_merge_method: str = "mean",
    ):
        self.reference_samples = reference_samples
        self.sample_to_plex = sample_to_plex
        self.fraction_merge_method = fraction_merge_method

    def quantify(self, psm_long_df: pd.DataFrame) -> pd.DataFrame:
        """
        Full ratio quantification pipeline.

        Parameters
        ----------
        psm_long_df : pd.DataFrame
            Long-format PSM data with columns: ProteinName, PeptideSequence
            (or PeptideCanonical), PrecursorCharge, SampleID, Fraction, Intensity.

        Returns
        -------
        pd.DataFrame
            Wide protein x sample matrix of log2 ratios, with ProteinName
            as the first column.
        """
        # Normalize column names
        df = psm_long_df.copy()
        if "PeptideSequence" in df.columns and PEPTIDE_CANONICAL not in df.columns:
            df.rename(columns={"PeptideSequence": PEPTIDE_CANONICAL}, inplace=True)
        if "Fraction" not in df.columns:
            df["Fraction"] = "1"

        n_psms_before = len(df)
        logger.info(f"Ratio quantification: {n_psms_before} PSM rows")

        # 1. Average fractions
        df = self._average_fractions(df)
        logger.info(f"After fraction averaging: {len(df)} rows")

        # 2. Compute reference intensity per plex per PSM group
        ref_df = self._compute_reference_intensity(df)

        # 3. Compute log2 ratios
        df = self._compute_log2_ratios(df, ref_df)
        logger.info(f"After log2 ratio computation: {len(df)} rows")

        # 4. Median aggregate PSM -> peptide (by sequence, ignoring charge)
        df = self._aggregate_psm_to_peptide(df)
        logger.info(f"After PSM->peptide aggregation: {len(df)} rows")

        # 5. Median aggregate peptide -> protein
        df = self._aggregate_peptide_to_protein(df)
        logger.info(f"After peptide->protein aggregation: {len(df)} rows")

        # 6. Pivot to wide format
        wide = df.pivot_table(
            index=PROTEIN_NAME,
            columns=SAMPLE_ID,
            values="log2ratio",
            aggfunc="median",
        )

        n_proteins = len(wide)
        n_samples = len(wide.columns)
        logger.info(f"Ratio quantification complete: {n_proteins} proteins, {n_samples} samples")

        return wide.reset_index()

    def _average_fractions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Average (or max) intensities across fractions per PSM per sample.

        Matches R's rowMeans(na.rm=TRUE) for the PS protocol.
        """
        group_cols = [PROTEIN_NAME, PEPTIDE_CANONICAL, PEPTIDE_CHARGE, SAMPLE_ID]

        if self.fraction_merge_method == "max":
            result = df.groupby(group_cols, as_index=False)["Intensity"].max()
        else:
            # mean — matches R's rowMeans(na.rm=TRUE)
            result = df.groupby(group_cols, as_index=False)["Intensity"].mean()

        return result

    def _compute_reference_intensity(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute reference intensity per plex per PSM group.

        For each plex, reference intensity = mean of reference channel
        intensities per (ProteinName, PeptideCanonical, PrecursorCharge).
        """
        # Assign plex to each row
        df["_plex"] = df[SAMPLE_ID].map(self.sample_to_plex)

        # Filter to reference samples only
        ref_mask = df[SAMPLE_ID].isin(self.reference_samples)
        ref_data = df[ref_mask]

        if len(ref_data) == 0:
            raise ValueError(
                f"No reference sample data found. Reference samples: {self.reference_samples}, "
                f"Available samples: {sorted(df[SAMPLE_ID].unique())}"
            )

        # Mean of reference channels per plex per PSM group
        group_cols = [PROTEIN_NAME, PEPTIDE_CANONICAL, PEPTIDE_CHARGE, "_plex"]
        ref_intensity = (
            ref_data.groupby(group_cols, as_index=False)["Intensity"]
            .mean()
            .rename(columns={"Intensity": "_ref_intensity"})
        )

        plexes = sorted(ref_intensity["_plex"].unique())
        for plex in plexes:
            plex_refs = [s for s in self.reference_samples if self.sample_to_plex.get(s) == plex]
            logger.info(f"  Plex '{plex}': reference channels = {plex_refs}")

        return ref_intensity

    def _compute_log2_ratios(self, df: pd.DataFrame, ref_df: pd.DataFrame) -> pd.DataFrame:
        """Compute log2(sample_intensity / reference_intensity) per PSM per plex."""
        df["_plex"] = df[SAMPLE_ID].map(self.sample_to_plex)

        # Remove reference samples from the output (log2(ref/ref) ~ 0)
        non_ref_mask = ~df[SAMPLE_ID].isin(self.reference_samples)
        df = df[non_ref_mask]

        # Merge reference intensity
        merge_cols = [PROTEIN_NAME, PEPTIDE_CANONICAL, PEPTIDE_CHARGE, "_plex"]
        df = df.merge(ref_df, on=merge_cols, how="left")

        # Replace zeros with NaN before division
        df["Intensity"] = df["Intensity"].replace(0, np.nan)
        df["_ref_intensity"] = df["_ref_intensity"].replace(0, np.nan)

        # Compute ratio and log2
        ratio = df["Intensity"] / df["_ref_intensity"]
        df["log2ratio"] = np.log2(ratio)

        # Replace inf/-inf with NaN
        df["log2ratio"] = df["log2ratio"].replace([np.inf, -np.inf], np.nan)

        # Drop helper columns
        df = df.drop(columns=["_plex", "_ref_intensity", "Intensity"])

        # Drop rows with NaN log2ratio
        df = df.dropna(subset=["log2ratio"])

        return df

    def _aggregate_psm_to_peptide(self, df: pd.DataFrame) -> pd.DataFrame:
        """Median aggregate PSM -> peptide (by sequence, ignoring charge)."""
        group_cols = [PROTEIN_NAME, PEPTIDE_CANONICAL, SAMPLE_ID]
        result = df.groupby(group_cols, as_index=False)["log2ratio"].median()
        return result

    def _aggregate_peptide_to_protein(self, df: pd.DataFrame) -> pd.DataFrame:
        """Median aggregate peptide -> protein."""
        group_cols = [PROTEIN_NAME, SAMPLE_ID]
        result = df.groupby(group_cols, as_index=False)["log2ratio"].median()
        return result


def apply_coverage_filter(
    protein_df: pd.DataFrame,
    sample_to_condition: dict[str, str],
    threshold: float = 0.65,
) -> pd.DataFrame:
    """
    Remove proteins with insufficient coverage in any condition.

    A protein passes if, in ALL conditions, the fraction of non-missing
    (and non-zero) values >= threshold.

    Matches R logic:
        all(sapply(gene_by_group, \\(gv) sum(!is.na(gv) & gv != 0) / length(gv) >= 0.65))

    Parameters
    ----------
    protein_df : pd.DataFrame
        Wide-format protein matrix (ProteinName | sample1 | sample2 | ...).
    sample_to_condition : dict
        Mapping from sample name to condition.
    threshold : float
        Minimum fraction of non-missing values required per condition.

    Returns
    -------
    pd.DataFrame
        Filtered protein matrix.
    """
    protein_col = protein_df.columns[0]
    sample_cols = [c for c in protein_df.columns if c != protein_col]

    # Group sample columns by condition
    condition_samples = {}
    for sample in sample_cols:
        cond = sample_to_condition.get(sample)
        if cond is not None:
            condition_samples.setdefault(cond, []).append(sample)

    if not condition_samples:
        logger.warning("No condition mapping for samples, skipping coverage filter")
        return protein_df

    intensity_matrix = protein_df.set_index(protein_col)

    # For each protein, check coverage in each condition
    keep_mask = pd.Series(True, index=intensity_matrix.index)

    for cond, samples in condition_samples.items():
        cond_data = intensity_matrix[samples]
        # Count non-missing and non-zero values
        valid_count = ((cond_data.notna()) & (cond_data != 0)).sum(axis=1)
        coverage = valid_count / len(samples)
        keep_mask = keep_mask & (coverage >= threshold)

    n_before = len(protein_df)
    result = protein_df[keep_mask.values]
    n_after = len(result)

    logger.info(
        f"Coverage filter ({threshold:.0%}): {n_before} -> {n_after} proteins "
        f"({n_before - n_after} removed)"
    )

    return result.reset_index(drop=True)


def load_psm_data(
    parquet_path: str,
    sdrf_path: str,
    min_aa: int = 7,
    min_unique_peptides: int = 2,
    remove_contaminants: bool = True,
) -> pd.DataFrame:
    """
    Load PSM data from QPX parquet, unnest intensities, join SDRF metadata.

    Uses DuckDB directly following the existing Feature class pattern.

    Parameters
    ----------
    parquet_path : str
        Path to the QPX parquet file.
    sdrf_path : str
        Path to the SDRF TSV file.
    min_aa : int
        Minimum peptide sequence length.
    min_unique_peptides : int
        Minimum unique peptides per protein.
    remove_contaminants : bool
        Whether to remove contaminants and decoys.

    Returns
    -------
    pd.DataFrame
        Long-format PSM data with columns: ProteinName, PeptideCanonical,
        PrecursorCharge, SampleID, Fraction, Intensity.
    """
    # Build SQL filters (where_clause built after is_decoy detection below)
    filter_builder = SQLFilterBuilder(
        remove_contaminants=remove_contaminants,
        min_peptide_length=min_aa,
        require_unique=True,
    )

    # Load SDRF for fraction info
    sdrf_df = pd.read_csv(sdrf_path, sep="\t")
    sdrf_df.columns = [c.lower() for c in sdrf_df.columns]

    # Build fraction mapping if available
    has_fraction = "comment[fraction identifier]" in sdrf_df.columns

    import re as _re

    def _strip_raw_ext(name: str) -> str:
        return _re.sub(
            r'\.(raw|mzML|mzml|d|wiff|RAW)$',
            '', str(name).replace('\\', '/').split('/')[-1],
        )

    conn = duckdb.connect()
    try:
        # Detect QPX format using parameterized read_parquet
        cols = [
            r[0] for r in conn.execute(
                "SELECT column_name FROM (DESCRIBE SELECT * FROM read_parquet(?))",
                [parquet_path],
            ).fetchall()
        ]
        is_new_qpx = "charge" in cols or "run_file_name" in cols

        # Set has_is_decoy before building WHERE clause so DECOY filter is optimal
        if "is_decoy" in cols:
            filter_builder.has_is_decoy = True
        where_clause, where_params = filter_builder.build_where_clause()

        # Detect if pg_accessions is list<struct{accession,...}> (new QPX)
        pg_is_struct = False
        if "pg_accessions" in cols:
            try:
                type_str = conn.execute(
                    "SELECT typeof(pg_accessions) FROM read_parquet(?) LIMIT 1",
                    [parquet_path],
                ).fetchone()[0].lower()
                pg_is_struct = "struct" in type_str
            except Exception as exc:
                logger.debug("Could not detect pg_accessions type: %s", exc)
        pg_col = (
            "list_transform(pg_accessions, x -> x.accession) as pg_accessions"
            if pg_is_struct
            else "pg_accessions"
        )

        # Predefined query templates (no user-controlled data)
        _QUERY_NEW_QPX = "".join([
            "SELECT ", pg_col, ", sequence,",
            " charge as precursor_charge,",
            " run_file_name as run_file_name,",
            " unnest.label as label,",
            " unnest.intensity as intensity",
            " FROM read_parquet(?) AS parquet_raw, UNNEST(intensities) as unnest",
            " WHERE unnest.intensity IS NOT NULL AND ",
        ])
        _QUERY_OLD_QPX = "".join([
            "SELECT ", pg_col, ", sequence,",
            " precursor_charge as precursor_charge,",
            " unnest.sample_accession as sample_accession,",
            " reference_file_name as run_file_name,",
            " unnest.channel as label,",
            " unnest.intensity as intensity",
            " FROM read_parquet(?) AS parquet_raw, UNNEST(intensities) as unnest",
            " WHERE unnest.intensity IS NOT NULL AND ",
        ])

        base_query = _QUERY_NEW_QPX if is_new_qpx else _QUERY_OLD_QPX
        query = "".join((base_query, where_clause))

        df = conn.execute(query, [parquet_path] + where_params).df()
    finally:
        conn.close()

    if len(df) == 0:
        raise ValueError("No PSM data after filtering. Check parquet file and filters.")

    # Build SDRF mapping: (run_file, label) -> source name, fraction
    sdrf_run_file = sdrf_df['comment[data file]'].apply(_strip_raw_ext) if 'comment[data file]' in sdrf_df.columns else sdrf_df.get('source name', pd.Series())
    sdrf_label = sdrf_df.get('comment[label]', pd.Series(dtype=str))
    sdrf_source = sdrf_df['source name']
    sdrf_fraction = sdrf_df.get('comment[fraction identifier]', pd.Series('1', index=sdrf_df.index))

    if is_new_qpx:
        # Map (run_file_name, label) -> source name via SDRF
        sdrf_map = pd.DataFrame({
            'run_file_name': sdrf_run_file.values,
            'label': sdrf_label.values if len(sdrf_label) > 0 else '',
            'source_name': sdrf_source.values,
            'fraction': sdrf_fraction.values,
        })
        df = df.merge(
            sdrf_map, on=['run_file_name', 'label'], how='left',
        )
        df['sample_accession'] = df['source_name'].fillna(df['run_file_name'])
        df['Fraction'] = df['fraction'].fillna('1')
    else:
        # Legacy: sample_accession already extracted from unnest
        if has_fraction:
            sdrf_fraction_map = {}
            for _, row in sdrf_df.iterrows():
                for col in ["comment[data file]", "comment[spectrum file]"]:
                    if col in sdrf_df.columns and pd.notna(row.get(col)):
                        sdrf_fraction_map[row[col]] = str(row["comment[fraction identifier]"])
                        break
            df["Fraction"] = df["run_file_name"].map(sdrf_fraction_map).fillna("1")
        else:
            df["Fraction"] = "1"

    # Parse protein accessions (take first accession, extract UniProt ID)
    first_acc = df["pg_accessions"].str[0].fillna("")
    df[PROTEIN_NAME] = np.where(
        first_acc.str.contains("|", regex=False),
        first_acc.str.split("|").str[1],
        first_acc,
    )

    # Rename columns to internal names
    df[PEPTIDE_CANONICAL] = df["sequence"]
    df[PEPTIDE_CHARGE] = df["precursor_charge"]
    df[SAMPLE_ID] = df["sample_accession"]
    df["Intensity"] = df["intensity"]

    # Keep only needed columns
    df = df[[PROTEIN_NAME, PEPTIDE_CANONICAL, PEPTIDE_CHARGE, SAMPLE_ID, "Fraction", "Intensity"]]

    # Filter multi-protein groups (keep only single-protein entries)
    df = df[~df[PROTEIN_NAME].str.contains(";")]

    # Filter by min unique peptides per protein
    peptide_counts = df.groupby(PROTEIN_NAME)[PEPTIDE_CANONICAL].nunique()
    valid_proteins = peptide_counts[peptide_counts >= min_unique_peptides].index
    df = df[df[PROTEIN_NAME].isin(valid_proteins)]

    n_proteins = df[PROTEIN_NAME].nunique()
    n_peptides = df[PEPTIDE_CANONICAL].nunique()
    n_samples = df[SAMPLE_ID].nunique()
    logger.info(
        f"Loaded PSM data: {len(df)} rows, {n_proteins} proteins, "
        f"{n_peptides} peptides, {n_samples} samples"
    )

    return df
