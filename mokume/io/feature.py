"""
Feature data access via DuckDB for proteomics parquet files.

This module provides the Feature class for reading and querying quantms.io/qpx
wide-format parquet files using DuckDB, and the SQLFilterBuilder for constructing
SQL-level filters for normalization pre-computations.
"""

import os
from dataclasses import dataclass, field
from typing import Iterator, Optional

import pandas as pd
import numpy as np
import duckdb

from mokume.model.labeling import QuantificationCategory, IsobaricLabel
from mokume.core.constants import load_sdrf
from mokume.core.logger import get_logger

logger = get_logger("mokume.io.feature")


@dataclass
class SQLFilterBuilder:
    """Builds SQL WHERE clauses for filtering parquet data at query time.

    This class is used to ensure that normalization factors (median maps, peptide
    frequencies, IRS scaling) are computed on filtered data that excludes contaminants,
    decoys, and other artifacts.

    Attributes
    ----------
    remove_contaminants : bool
        Whether to exclude rows containing contaminant patterns in pg_accessions.
    contaminant_patterns : list[str]
        List of substring patterns to identify contaminants (e.g., ["CONTAMINANT", "DECOY"]).
    min_intensity : float
        Minimum intensity threshold (0.0 means only filter intensity > 0).
    min_peptide_length : int
        Minimum peptide sequence length.
    require_unique : bool
        Whether to require unique peptides only (unique = 1).
    """

    remove_contaminants: bool = True
    contaminant_patterns: list[str] = field(
        default_factory=lambda: ["CONTAMINANT", "ENTRAP", "DECOY"]
    )
    min_intensity: float = 0.0
    min_peptide_length: int = 7
    require_unique: bool = True

    def build_where_clause(self) -> str:
        """Build SQL WHERE clause string for DuckDB queries.

        Returns
        -------
        str
            A SQL WHERE clause (without the WHERE keyword) that can be used
            in DuckDB queries to filter the parquet data.
        """
        conditions = []

        # Always filter intensity > 0
        conditions.append("intensity > 0")

        # Min intensity threshold
        if self.min_intensity > 0:
            conditions.append(f"intensity >= {self.min_intensity}")

        # Peptide length filter
        if self.min_peptide_length > 0:
            conditions.append(f'LENGTH("sequence") >= {self.min_peptide_length}')

        # Unique peptides only
        if self.require_unique:
            conditions.append('"unique" = 1')

        # Contaminant/decoy filter - cast pg_accessions array to text for LIKE matching
        if self.remove_contaminants and self.contaminant_patterns:
            pattern_conditions = []
            for pattern in self.contaminant_patterns:
                # Escape any SQL special characters in the pattern
                safe_pattern = pattern.replace("'", "''")
                pattern_conditions.append(f"pg_accessions::text NOT LIKE '%{safe_pattern}%'")
            conditions.append(f"({' AND '.join(pattern_conditions)})")

        return " AND ".join(conditions) if conditions else "1=1"


class Feature:
    """
    Represents a feature in a proteomics dataset, providing methods for data manipulation
    and analysis using a DuckDB database connection to a Parquet file.

    This class expects the quantms.io/qpx wide format where intensities are stored
    in a nested array. Supports both new QPX format intensities[{label, intensity}, ...]
    and legacy format intensities[{sample_accession, channel, intensity}, ...]

    Attributes
    ----------
    filter_builder : Optional[SQLFilterBuilder]
        If provided, this filter builder will be used to apply SQL-level filtering
        when computing normalization factors (median maps, peptide frequencies, etc.).
        This ensures normalization is computed on clean data without contaminants/decoys.
    """

    labels: Optional[list[str]]
    label: Optional[QuantificationCategory]
    choice: Optional[IsobaricLabel]
    technical_repetitions: Optional[int]
    filter_builder: Optional[SQLFilterBuilder]

    def __init__(
        self, database_path: str, filter_builder: Optional[SQLFilterBuilder] = None
    ):
        if not os.path.exists(database_path):
            raise FileNotFoundError(f"the file {database_path} does not exist.")

        self.parquet_db = duckdb.connect()

        safe_path = database_path.replace("'", "''")
        self.parquet_db.execute(
            "CREATE VIEW parquet_db_raw AS SELECT * FROM parquet_scan('{}')".format(safe_path)
        )

        self._detect_qpx_format()
        self._create_unnest_view()

        self.samples = self.get_unique_samples()
        self.filter_builder = filter_builder

    def _detect_qpx_format(self) -> None:
        """Detect whether the parquet uses new or legacy QPX schema."""
        cols = [
            r[0]
            for r in self.parquet_db.execute(
                "SELECT column_name FROM (DESCRIBE parquet_db_raw)"
            ).fetchall()
        ]
        self._is_new_qpx = "charge" in cols or "run_file_name" in cols
        self._charge_col = "charge" if self._is_new_qpx else "precursor_charge"
        self._run_col = "run_file_name" if self._is_new_qpx else "reference_file_name"

    def _create_unnest_view(self) -> None:
        """Create the long-format DuckDB view by unnesting intensities."""
        if self._is_new_qpx:
            unnest_sql = (
                "run_file_name as sample_accession,\n"
                "                    unnest.label as channel,\n"
                "                    unnest.intensity"
            )
            sa_default = "run_file_name"
        else:
            unnest_sql = (
                "unnest.sample_accession,\n"
                "                    unnest.channel,\n"
                "                    unnest.intensity"
            )
            sa_default = "unnest.sample_accession"

        charge_col, run_col = self._charge_col, self._run_col
        self.parquet_db.execute(f"""
            CREATE VIEW parquet_db AS
            SELECT
                sequence,
                peptidoform,
                pg_accessions,
                {charge_col} as charge,
                {run_col} as run_file_name,
                "unique",
                {unnest_sql},
                -- Defaults (can be enriched with SDRF later)
                {run_col} as run,
                {sa_default} as condition,
                1 as biological_replicate,
                '1' as fraction,
                split_part({sa_default}, '_', 1) as mixture
            FROM parquet_db_raw, UNNEST(intensities) as unnest
            WHERE unnest.intensity IS NOT NULL AND unnest.intensity > 0
        """)

    def enrich_with_sdrf(self, sdrf_path: str) -> None:
        """Enrich parquet data with SDRF metadata (condition, biological_replicate, etc.).

        Parameters
        ----------
        sdrf_path : str
            Path to the SDRF file containing sample metadata.
        """
        import re as _re

        sdrf_df = load_sdrf(sdrf_path)

        # Find the condition column (try factor value first, then characteristics)
        condition_col = None
        for col in sdrf_df.columns:
            if "factor value" in col:
                condition_col = col
                break
        if condition_col is None:
            for col in sdrf_df.columns:
                if "organism part" in col and "characteristics" in col:
                    condition_col = col
                    break

        # Strip common raw-file extensions so names match QPX run_file_name
        def _strip_raw_ext(name: str) -> str:
            return _re.sub(
                r"\.(raw|mzML|mzml|d|wiff|RAW)$",
                "",
                str(name).replace("\\", "/").split("/")[-1],
            )

        # Prepare SDRF mapping with both join keys
        sdrf_mapping = pd.DataFrame(
            {
                "sdrf_run_file": sdrf_df["comment[data file]"].apply(_strip_raw_ext),
                "sdrf_label": sdrf_df.get("comment[label]", ""),
                "sdrf_sample_accession": sdrf_df["source name"],
                "sdrf_condition": (
                    sdrf_df[condition_col] if condition_col else sdrf_df["source name"]
                ),
                "sdrf_biological_replicate": sdrf_df.get(
                    "characteristics[biological replicate]", 1
                ),
                "sdrf_fraction": sdrf_df.get("comment[fraction identifier]", "1"),
            }
        )

        self.parquet_db.register("sdrf_mapping", sdrf_mapping)

        # Build format-aware UNNEST SQL
        charge_col = self._charge_col
        run_col = self._run_col

        if self._is_new_qpx:
            unnest_cols = "unnest.label as channel,\n                unnest.intensity"
            extra_cols = ""
            join_clause = (
                "ON p.run_file_name = s.sdrf_run_file\n"
                "                AND (p.channel = s.sdrf_label\n"
                "                     OR s.sdrf_label IS NULL\n"
                "                     OR s.sdrf_label = '')"
            )
            sa_fallback = "p.run_file_name"
        else:
            unnest_cols = "unnest.channel,\n                unnest.intensity"
            extra_cols = ",\n                unnest.sample_accession as _legacy_sa"
            join_clause = "ON p._legacy_sa = s.sdrf_sample_accession"
            sa_fallback = "p._legacy_sa"

        # Create intermediate view for unnested data
        self.parquet_db.execute(f"""
            CREATE OR REPLACE VIEW parquet_db_unnested AS
            SELECT
                sequence,
                peptidoform,
                pg_accessions,
                {charge_col} as charge,
                {run_col} as run_file_name,
                "unique",
                {unnest_cols},
                {run_col} as run{extra_cols}
            FROM parquet_db_raw, UNNEST(intensities) as unnest
            WHERE unnest.intensity IS NOT NULL AND unnest.intensity > 0
        """)

        # Recreate main view with SDRF data joined
        self.parquet_db.execute("DROP VIEW IF EXISTS parquet_db")
        self.parquet_db.execute(f"""
            CREATE VIEW parquet_db AS
            SELECT
                p.sequence,
                p.peptidoform,
                p.pg_accessions,
                p.charge,
                p.run_file_name,
                p."unique",
                COALESCE(s.sdrf_sample_accession, {sa_fallback}) as sample_accession,
                p.channel,
                p.intensity,
                p.run,
                COALESCE(s.sdrf_condition, {sa_fallback}) as condition,
                COALESCE(CAST(s.sdrf_biological_replicate AS INTEGER), 1) as biological_replicate,
                COALESCE(CAST(s.sdrf_fraction AS VARCHAR), '1') as fraction,
                split_part(COALESCE(s.sdrf_sample_accession, {sa_fallback}), '_', 1) as mixture
            FROM parquet_db_unnested p
            LEFT JOIN sdrf_mapping s
                {join_clause}
        """)

        logger.info("Enriched parquet data with SDRF metadata from %s", sdrf_path)

    @staticmethod
    def standardize_df(df: pd.DataFrame) -> pd.DataFrame:
        """Standardizes column names in the given DataFrame."""
        return df.rename({"protein_accessions": "pg_accessions"}, axis=1)

    @property
    def experimental_inference(
        self,
    ) -> tuple[int, QuantificationCategory, list[str], Optional[IsobaricLabel]]:
        """Infers experimental details from the dataset."""
        self.labels = self.get_unique_labels()
        self.label, self.choice = QuantificationCategory.classify(self.labels)
        self.technical_repetitions = self.get_unique_tec_reps()
        return len(self.technical_repetitions), self.label, self.samples, self.choice

    def get_low_frequency_peptides(self, percentage: float = 0.2) -> tuple:
        """Identifies peptides that occur with low frequency across samples.

        If a filter_builder is set on this Feature instance, it will be used to
        exclude contaminants, decoys, and other artifacts from the frequency
        calculation.

        Parameters
        ----------
        percentage : float
            The frequency threshold. Peptides appearing in less than this
            fraction of samples are considered low frequency. Default is 0.2 (20%).

        Returns
        -------
        tuple
            A tuple of (protein_accession, sequence) pairs for low frequency peptides.
        """
        where_clause = self.filter_builder.build_where_clause() if self.filter_builder else "1=1"

        f_table = self.parquet_db.sql(f"""
            SELECT "sequence", "pg_accessions", COUNT(DISTINCT sample_accession) as "count"
            FROM parquet_db
            WHERE {where_clause}
            GROUP BY "sequence", "pg_accessions"
            """).df()
        f_table.dropna(subset=["pg_accessions"], inplace=True)
        try:
            f_table["pg_accessions"] = f_table["pg_accessions"].apply(lambda x: x[0].split("|")[1])
        except IndexError:
            f_table["pg_accessions"] = f_table["pg_accessions"].apply(lambda x: x[0])
        except Exception as e:
            raise ValueError(
                "Some errors occurred when parsing pg_accessions column in feature parquet!"
            ) from e
        f_table.set_index(["sequence", "pg_accessions"], inplace=True)
        f_table.drop(
            f_table[f_table["count"] >= (percentage * len(self.samples))].index,
            inplace=True,
        )
        f_table.reset_index(inplace=True)
        return tuple(zip(f_table["pg_accessions"], f_table["sequence"]))

    @staticmethod
    def csv2parquet(csv):
        """Converts a CSV file to a Parquet file using DuckDB."""
        parquet_path = os.path.splitext(csv)[0] + ".parquet"
        duckdb.read_csv(csv).to_parquet(parquet_path)

    def get_report_from_database(self, samples: list, columns: list = None):
        """Retrieves a standardized report from the database for specified samples."""
        cols = ",".join(columns) if columns is not None else "*"
        database = self.parquet_db.sql(
            """SELECT {} FROM parquet_db WHERE sample_accession IN {}""".format(
                cols, tuple(samples)
            )
        )
        report = database.df()
        return Feature.standardize_df(report)

    def iter_samples(
        self, sample_num: int = 20, columns: list = None
    ) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """Iterates over samples in batches."""
        ref_list = [
            self.samples[i : i + sample_num] for i in range(0, len(self.samples), sample_num)
        ]
        for refs in ref_list:
            batch_df = self.get_report_from_database(refs, columns)
            yield refs, batch_df

    def get_unique_samples(self) -> list[str]:
        """Retrieves a list of unique sample accessions from the Parquet database."""
        unique = self.parquet_db.sql("SELECT DISTINCT sample_accession FROM parquet_db").df()
        return unique["sample_accession"].tolist()

    def get_unique_labels(self) -> list[str]:
        """Retrieves a list of unique channel labels from the Parquet database."""
        unique = self.parquet_db.sql("SELECT DISTINCT channel FROM parquet_db").df()
        return unique["channel"].tolist()

    def get_unique_tec_reps(self) -> list[int]:
        """Retrieves a list of unique technical repetition identifiers.

        Attempts to extract technical replicate numbers from run names in order:
        1. If run name contains '_', try to extract the last part as an integer
        2. If run name is numeric, use it directly
        3. Fall back to sequential integers (1, 2, 3, ...) based on unique runs
        """
        unique = self.parquet_db.sql("SELECT DISTINCT run FROM parquet_db").df()

        try:
            # Try to extract last part after underscore as tech rep
            if unique["run"].str.contains("_").all():
                # Get the last part after splitting by underscore
                last_parts = unique["run"].str.split("_").str.get(-1)
                unique["run"] = last_parts.astype("int")
            else:
                # Try to convert directly to int
                unique["run"] = unique["run"].astype("int")
        except (ValueError, TypeError):
            # Fall back to sequential integers
            unique["run"] = list(range(1, len(unique) + 1))

        return unique["run"].tolist()

    def get_median_map(self) -> dict[str, float]:
        """Computes a median intensity map for samples.

        If a filter_builder is set on this Feature instance, it will be used to
        exclude contaminants, decoys, and other artifacts from the median
        calculation. This ensures normalization factors are computed on clean data.

        Returns
        -------
        dict[str, float]
            A dictionary mapping sample accessions to their normalization factors
            (sample median / global median).
        """
        where_clause = self.filter_builder.build_where_clause() if self.filter_builder else "1=1"

        # Use SQL aggregation with filtering for efficiency
        result = self.parquet_db.sql(f"""
            SELECT sample_accession, MEDIAN(intensity) as median_intensity
            FROM parquet_db
            WHERE {where_clause}
            GROUP BY sample_accession
            """).df()

        med_map = dict(zip(result["sample_accession"], result["median_intensity"]))
        global_med = np.median(list(med_map.values()))

        for sample, med in med_map.items():
            med_map[sample] = med / global_med

        return med_map

    def get_report_condition_from_database(self, cons: list, columns: list = None) -> pd.DataFrame:
        """Retrieves a standardized report from the database for specified conditions."""
        cols = ",".join(columns) if columns is not None else "*"
        database = self.parquet_db.sql(
            f"""SELECT {cols} FROM parquet_db WHERE condition IN {tuple(cons)}"""
        )
        report = database.df()
        return Feature.standardize_df(report)

    def iter_conditions(
        self, conditions: int = 10, columns: list = None
    ) -> Iterator[tuple[list[str], pd.DataFrame]]:
        """Iterates over experimental conditions in batches."""
        condition_list = self.get_unique_conditions()
        ref_list = [
            condition_list[i : i + conditions] for i in range(0, len(condition_list), conditions)
        ]
        for refs in ref_list:
            batch_df = self.get_report_condition_from_database(refs, columns)
            yield refs, batch_df

    def get_unique_conditions(self) -> list[str]:
        """Retrieves a list of unique experimental conditions from the Parquet database."""
        unique = self.parquet_db.sql("SELECT DISTINCT condition FROM parquet_db").df()
        return unique["condition"].tolist()

    def get_median_map_to_condition(self) -> dict[str, dict[str, float]]:
        """Computes a median intensity map for each experimental condition.

        If a filter_builder is set on this Feature instance, it will be used to
        exclude contaminants, decoys, and other artifacts from the median
        calculation. This ensures normalization factors are computed on clean data.

        Returns
        -------
        dict[str, dict[str, float]]
            A nested dictionary mapping conditions to sample normalization factors.
            For each condition, samples are normalized to the condition mean.
        """
        where_clause = self.filter_builder.build_where_clause() if self.filter_builder else "1=1"

        # Use SQL aggregation with filtering for efficiency
        result = self.parquet_db.sql(f"""
            SELECT condition, sample_accession, MEDIAN(intensity) as median_intensity
            FROM parquet_db
            WHERE {where_clause}
            GROUP BY condition, sample_accession
            """).df()

        med_map = {}
        for condition in result["condition"].unique():
            cond_data = result[result["condition"] == condition]
            meds = pd.Series(
                cond_data["median_intensity"].values,
                index=cond_data["sample_accession"].values,
            )
            meds = meds / meds.mean()
            med_map[condition] = meds.to_dict()

        return med_map

    def get_irs_scaling_factors(
        self,
        irs_channel: str,
        irs_stat: str = "median",
        irs_scope: str = "global",
    ) -> dict[int, float]:
        """Compute IRS (Internal Reference Scaling) factors with filtering applied.

        If a filter_builder is set on this Feature instance, it will be used to
        exclude contaminants, decoys, and other artifacts from the IRS calculation.

        Parameters
        ----------
        irs_channel : str
            The channel label to use as internal reference (e.g., "126").
        irs_stat : str
            Statistic to use for computing reference values: "median" or "mean".
        irs_scope : str
            Scope of normalization: "global", "by_mixture", or "two_stage".

        Returns
        -------
        dict[int, float]
            Dictionary mapping technical replicate indices to scaling factors.
        """
        stat_fn = "median" if (irs_stat or "").lower() == "median" else "avg"

        # Build filter conditions for contaminants only (not unique peptide requirement)
        # since IRS uses specific channel which may have different characteristics
        filter_conditions = ["intensity > 0"]

        if self.filter_builder and self.filter_builder.remove_contaminants:
            for pattern in self.filter_builder.contaminant_patterns:
                safe_pattern = pattern.replace("'", "''")
                filter_conditions.append(f"pg_accessions::text NOT LIKE '%{safe_pattern}%'")

        if self.filter_builder and self.filter_builder.min_intensity > 0:
            filter_conditions.append(f"intensity >= {self.filter_builder.min_intensity}")

        # Add channel filter
        filter_conditions.append(f"channel = '{irs_channel}'")
        where_clause = " AND ".join(filter_conditions)

        irs_df = self.parquet_db.sql(f"""
            SELECT run, {stat_fn}(intensity) as irs_value, mixture, techreplicate as techrep_guess
            FROM (
                SELECT *,
                       CASE WHEN position('_' in run) > 0 THEN CAST(split_part(run, '_', 2) AS INTEGER)
                            ELSE CAST(run AS INTEGER) END AS techreplicate
                FROM parquet_db
                WHERE {where_clause}
            )
            GROUP BY run, mixture, techrep_guess
            """).df()

        irs_scale_by_techrep: dict[int, float] = {}

        if len(irs_df.index) > 0:
            irs_df = irs_df[irs_df["irs_value"] > 0]

            if irs_scope.lower() == "by_mixture":
                transform_fn = "median" if stat_fn == "median" else "mean"
                irs_df["mixture_center"] = irs_df.groupby("mixture")["irs_value"].transform(
                    transform_fn
                )
                irs_df["scale"] = irs_df["mixture_center"] / irs_df["irs_value"]
            elif irs_scope.lower() == "two_stage":
                transform_fn = "median" if stat_fn == "median" else "mean"
                irs_df["mixture_center"] = irs_df.groupby("mixture")["irs_value"].transform(
                    transform_fn
                )
                irs_df["scale_stage1"] = irs_df["mixture_center"] / irs_df["irs_value"]
                mixture_center_df = irs_df[["mixture", "mixture_center"]].drop_duplicates()
                if stat_fn == "median":
                    global_center = mixture_center_df["mixture_center"].median()
                else:
                    global_center = mixture_center_df["mixture_center"].mean()
                mixture_center_df["scale_stage2"] = (
                    global_center / mixture_center_df["mixture_center"]
                )
                irs_df = irs_df.merge(
                    mixture_center_df, on="mixture", how="left", suffixes=("", "_mix")
                )
                irs_df["scale"] = irs_df["scale_stage1"] * irs_df["scale_stage2"]
            else:
                # Global scope
                if stat_fn == "median":
                    global_center = irs_df["irs_value"].median()
                else:
                    global_center = irs_df["irs_value"].mean()
                irs_df["scale"] = global_center / irs_df["irs_value"]

            irs_scale_by_techrep = dict(
                zip(irs_df["techrep_guess"].tolist(), irs_df["scale"].tolist())
            )

        return irs_scale_by_techrep
