"""
Pipeline stage classes for the quantification pipeline.

Each stage handles a distinct phase of the proteomics quantification workflow:
- LoadingStage: Data loading and filtering
- NormalizationStage: Run-level and sample-level normalization
- QuantificationStage: Protein quantification methods
- PostprocessingStage: Batch correction, DE, plotting, reports
"""

import os
import re
from collections import Counter

import numpy as np
import pandas as pd
import pyarrow as pa
from pathlib import Path
from typing import Optional

from mokume.io.feature import Feature, SQLFilterBuilder
from mokume.preprocessing.sdrf import analyse_sdrf
from mokume.preprocessing.aggregation import (
    remove_contaminants_entrapments_decoys,
    apply_initial_filtering,
    get_peptidoform_normalize_intensities,
    sum_peptidoform_intensities,
    reformat_quantms_feature_table_quant_labels,
)
from mokume.normalization.hierarchical import HierarchicalSampleNormalizer
from mokume.model.normalization import (
    FeatureNormalizationMethod,
    PeptideNormalizationMethod,
)
from mokume.model.labeling import QuantificationCategory
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
    INTENSITY,
    PARQUET_COLUMNS,
    AGGREGATION_LEVEL_SAMPLE,
    CONDITION,
)
from mokume.core.logger import get_logger
from mokume.postprocessing.batch_correction import (
    is_batch_correction_available,
    detect_batches,
    extract_covariates_from_sdrf,
    apply_batch_correction,
)
from mokume.model.batch_correction import BatchDetectionMethod
from mokume.pipeline.config import PipelineConfig

logger = get_logger("mokume.pipeline")

# Per-(sample, tech_rep) replicate metric expressions for SQL-first run
# normalization. Each value must produce a scalar metric per group so that the
# cross-replicate factor ``intensity * sample_avg / run_metric`` matches the
# legacy ``normalize_sample`` / ``normalize_replicates`` contract from
# :mod:`mokume.model.normalization`. ``max_min`` is intentionally absent because
# its replicate function returns a per-row Series rather than a scalar metric;
# workflows using it fall back to the legacy pandas loop.
_SQL_RUN_METRIC_EXPR: dict[str, str] = {
    "median": 'MEDIAN("Intensity")',
    "mean": 'AVG("Intensity")',
    "max": 'MAX("Intensity")',
    "global": 'SUM("Intensity")',
    "iqr": (
        '(QUANTILE_CONT("Intensity", 0.75) + QUANTILE_CONT("Intensity", 0.25)) / 2'
    ),
}
_SQL_RUN_NORM_SUPPORTED: frozenset[str] = frozenset({"none", ""}) | frozenset(
    _SQL_RUN_METRIC_EXPR
)

_SQLFIRST_LOAD_QUERY_PREFIX = """
        WITH base AS (
            SELECT
                array_to_string(
                    list_transform(
                        pg_accessions,
                        x -> CASE
                            WHEN length(string_split(x, '|')) = 3
                                THEN string_split(x, '|')[2]
                            ELSE x
                        END
                    ),
                    ';'
                ) AS "ProteinName",
                sequence       AS "PeptideCanonical",
                peptidoform    AS "PeptideSequence",
                charge         AS "PrecursorCharge",
                intensity      AS "Intensity",
                COALESCE(condition, 'Empty') AS "Condition",
                COALESCE(CAST(biological_replicate AS INTEGER), 1) AS "BioReplicate",
                sample_accession AS "SampleID",
                run            AS "Run"
            FROM parquet_db
            WHERE ("""

_SQLFIRST_LOAD_QUERY_AFTER_WHERE = """)
              AND pg_accessions IS NOT NULL
              AND length(sequence) >= ?
              AND (condition IS NULL OR condition != 'Empty')
              AND run IS NOT NULL
        ),
        no_contam AS (
            SELECT * FROM base
            WHERE "ProteinName" NOT LIKE '%CONTAMINANT%'
              AND "ProteinName" NOT LIKE '%ENTRAP%'
              AND "ProteinName" NOT LIKE '%DECOY%'
        ),
        with_n_unique AS (
            SELECT *,
                COUNT(DISTINCT "PeptideCanonical") OVER (
                    PARTITION BY "ProteinName", "SampleID"
                ) AS _n_unique_pep
            FROM no_contam
        ),
        min_unique_kept AS (
            SELECT * FROM with_n_unique WHERE _n_unique_pep >= ?
        ),
        peptidoform_max AS (
            SELECT * FROM (
                SELECT """

_SQLFIRST_LOAD_QUERY_AFTER_SELECT = """,
                    ROW_NUMBER() OVER (
                        PARTITION BY mu."PeptideSequence", mu."PrecursorCharge",
                                     mu."SampleID", mu."Condition", mu."BioReplicate"
                        ORDER BY """

_SQLFIRST_LOAD_QUERY_AFTER_ORDERBY = """ DESC, mu."Run"
                    ) AS _rn
                FROM min_unique_kept mu
                """

_SQLFIRST_LOAD_QUERY_SUFFIX = """
            )
            WHERE _rn = 1
        )
        SELECT
            "ProteinName",
            "PeptideCanonical",
            "SampleID",
            "BioReplicate",
            "Condition",
            SUM("Intensity") AS "NormIntensity"
        FROM peptidoform_max
        GROUP BY "ProteinName", "PeptideCanonical", "SampleID",
                 "BioReplicate", "Condition"
        """

_MAXLFQ_DIRECTLFQ_QUERY_SUFFIX = """
            )
            WHERE _rn = 1
        )
        SELECT
            "ProteinName" AS protein,
            "PeptideCanonical" AS ion,
            "SampleID" AS sample_accession,
            SUM("Intensity") AS intensity
        FROM peptidoform_max
        GROUP BY "ProteinName", "PeptideCanonical", "SampleID"
        """

_RUN_SCALE_QUERY_PREFIX = """
        WITH base AS (
            SELECT
                array_to_string(
                    list_transform(
                        pg_accessions,
                        x -> CASE
                            WHEN length(string_split(x, '|')) = 3
                                THEN string_split(x, '|')[2]
                            ELSE x
                        END
                    ),
                    ';'
                ) AS "ProteinName",
                sequence         AS "PeptideCanonical",
                intensity        AS "Intensity",
                sample_accession AS "SampleID",
                run              AS "Run"
            FROM parquet_db
            WHERE ("""

_RUN_SCALE_QUERY_AFTER_WHERE = """)
              AND pg_accessions IS NOT NULL
              AND length(sequence) >= ?
              AND (condition IS NULL OR condition != 'Empty')
              AND run IS NOT NULL
        ),
        no_contam AS (
            SELECT * FROM base
            WHERE "ProteinName" NOT LIKE '%CONTAMINANT%'
              AND "ProteinName" NOT LIKE '%ENTRAP%'
              AND "ProteinName" NOT LIKE '%DECOY%'
        ),
        with_n_unique AS (
            SELECT *,
                COUNT(DISTINCT "PeptideCanonical") OVER (
                    PARTITION BY "ProteinName", "SampleID"
                ) AS _n_unique_pep
            FROM no_contam
        ),
        kept AS (
            SELECT "SampleID", "Run", "Intensity"
            FROM with_n_unique WHERE _n_unique_pep >= ?
        ),
        sample_run_meta AS (
            SELECT
                "SampleID",
                BOOL_AND("Run" LIKE '%\\_%' ESCAPE '\\') AS all_have_underscore,
                BOOL_AND(
                    TRY_CAST(split_part("Run", '_', -1) AS INTEGER) IS NOT NULL
                ) AS all_parse_last,
                BOOL_AND(TRY_CAST("Run" AS INTEGER) IS NOT NULL) AS all_parse_full
            FROM kept
            GROUP BY "SampleID"
        ),
        with_tech_rep AS (
            SELECT
                m."SampleID",
                m."Run",
                m."Intensity",
                CASE
                    WHEN srm.all_have_underscore AND srm.all_parse_last
                        THEN CAST(split_part(m."Run", '_', -1) AS INTEGER)
                    WHEN srm.all_parse_full
                        THEN CAST(m."Run" AS INTEGER)
                    ELSE
                        DENSE_RANK() OVER (
                            PARTITION BY m."SampleID" ORDER BY m."Run"
                        )
                END AS tech_rep
            FROM kept m
            JOIN sample_run_meta srm USING ("SampleID")
        ),
        run_metric AS (
            SELECT "SampleID", tech_rep, """

_RUN_SCALE_QUERY_SUFFIX = """ AS rm
            FROM with_tech_rep
            GROUP BY "SampleID", tech_rep
        ),
        sample_avg AS (
            SELECT "SampleID", AVG(rm) AS sa, COUNT(*) AS n_tech_reps
            FROM run_metric
            GROUP BY "SampleID"
        ),
        run_to_tech_rep AS (
            SELECT DISTINCT "SampleID", "Run", tech_rep
            FROM with_tech_rep
        )
        SELECT
            r."SampleID" AS "SampleID",
            r."Run"         AS "Run",
            CASE
                WHEN sa.n_tech_reps > 1
                    THEN sa.sa / NULLIF(rm.rm, 0)
                ELSE CAST(1.0 AS DOUBLE)
            END AS scale_factor
        FROM run_to_tech_rep r
        JOIN run_metric rm USING ("SampleID", tech_rep)
        JOIN sample_avg  sa USING ("SampleID")
        """


def _build_sqlfirst_load_query(
    where_clause: str,
    scale_select: str,
    scale_orderby: str,
    scale_join_clause: str,
    query_suffix: str = _SQLFIRST_LOAD_QUERY_SUFFIX,
) -> str:
    """Build the SQL-first peptide loading query from trusted fragments."""
    return "".join(
        [
            _SQLFIRST_LOAD_QUERY_PREFIX,
            where_clause,
            _SQLFIRST_LOAD_QUERY_AFTER_WHERE,
            scale_select,
            _SQLFIRST_LOAD_QUERY_AFTER_SELECT,
            scale_orderby,
            _SQLFIRST_LOAD_QUERY_AFTER_ORDERBY,
            scale_join_clause,
            query_suffix,
        ]
    )


def _build_run_scale_query(where_clause: str, metric_expr: str) -> str:
    """Build the run normalization scale-factor query from trusted fragments."""
    return "".join(
        [
            _RUN_SCALE_QUERY_PREFIX,
            where_clause,
            _RUN_SCALE_QUERY_AFTER_WHERE,
            metric_expr,
            _RUN_SCALE_QUERY_SUFFIX,
        ]
    )


class LoadingStage:
    """Loads and filters data from parquet + SDRF."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def load_for_mokume(self) -> pd.DataFrame:
        """Load data and apply normalization for mokume quantification methods.

        Dispatcher that routes to a SQL-first path when the workflow allows
        pushing the entire filter+aggregation pipeline into a single DuckDB query
        (LFQ label, no per-run normalization needed); otherwise falls back to the
        legacy per-sample pandas loop. Both paths produce equivalent peptide-level
        DataFrames (validated on PXD003539: rows + keys identical, NormIntensity
        agreeing to relative 1e-7, the precision limit of the float32 intensity
        column stored in parquet).
        """
        keep_shared_peptides = self.config.quantification.method.lower() == "ibaq"
        filter_builder = SQLFilterBuilder(
            remove_contaminants=self.config.filtering.remove_contaminants,
            min_peptide_length=self.config.filtering.min_aa,
            require_unique=not keep_shared_peptides,
        )

        feature = Feature(
            self.config.input.parquet,
            filter_builder=filter_builder,
            duckdb_memory=self.config.runtime.duckdb_memory,
            duckdb_threads=self.config.runtime.duckdb_threads,
        )

        if self.config.input.sdrf:
            feature.enrich_with_sdrf(self.config.input.sdrf)
            technical_repetitions, label, _sample_names, choice = analyse_sdrf(
                self.config.input.sdrf
            )
        else:
            (
                technical_repetitions,
                label,
                _sample_names,
                choice,
            ) = feature.experimental_inference

        # Build normalization config snapshot (shared between both load paths).
        run_norm = self.config.normalization.run_method.lower()
        sample_norm = self.config.normalization.sample_method.lower()
        sample_norm_method = PeptideNormalizationMethod.from_str(sample_norm)

        med_map: dict = {}
        if sample_norm == "globalmedian":
            med_map = feature.get_median_map()
        elif sample_norm == "conditionmedian":
            med_map = feature.get_median_map_to_condition()

        # SQL-first path is eligible when (a) the label is LFQ (TMT/iTRAQ need
        # channel-aware reformatting that the SQL prototype does not implement),
        # and (b) the requested run normalization method has a SQL translation
        # registered in :data:`_SQL_RUN_METRIC_EXPR`. Workflows using ``max_min``
        # (whose replicate function returns a per-row Series) fall back to the
        # legacy pandas loop.
        sqlfirst_eligible = (
            label == QuantificationCategory.LFQ and run_norm in _SQL_RUN_NORM_SUPPORTED
        )

        if sqlfirst_eligible:
            logger.info(
                "Loading with SQL-first path (label=%s, tech_reps=%d, run_norm=%s)",
                label,
                technical_repetitions,
                run_norm,
            )
            combined_df = self._load_for_mokume_sqlfirst(
                feature,
                sample_norm,
                sample_norm_method,
                med_map,
                run_norm=run_norm,
                technical_repetitions=technical_repetitions,
            )
        else:
            logger.info(
                "Loading with legacy per-sample pandas path "
                "(label=%s, tech_reps=%d, run_norm=%s)",
                label,
                technical_repetitions,
                run_norm,
            )
            combined_df = self._load_for_mokume_legacy_pandas(
                feature,
                label,
                choice,
                technical_repetitions,
                sample_norm,
                sample_norm_method,
                med_map,
            )

        return self._apply_dataset_normalization(combined_df, sample_norm_method)

    def _load_for_mokume_sqlfirst(
        self,
        feature: Feature,
        sample_norm: str,
        sample_norm_method: PeptideNormalizationMethod,
        med_map: dict,
        run_norm: str = "none",
        technical_repetitions: int = 1,
    ) -> pd.DataFrame:
        """Push the full filter+aggregation pipeline into a single DuckDB query.

        Translates the steps that the legacy pandas loop performs per sample
        (parse_uniprot_accession, length / condition / run filters,
        contaminant/decoy removal, per-sample min-unique-peptides filter,
        optional per-(sample, tech_rep) run normalization, peptidoform-max
        keep, sum aggregation at the sample level) into one SQL query against
        the SDRF-enriched ``parquet_db`` view. Only the small post-aggregation
        peptide-level table lands in pandas, where per-sample normalization is
        applied if needed.

        Run normalization is enabled when both ``run_norm`` resolves to a
        registered SQL metric and ``technical_repetitions > 1`` (matching the
        global guard in :meth:`FeatureNormalizationMethod.normalize_runs`).
        """
        where_clause, where_params = feature.filter_builder.build_where_clause()
        min_aa = self.config.filtering.min_aa
        min_unique = (
            0
            if self.config.quantification.method.lower() == "ibaq"
            else self.config.filtering.min_unique_peptides
        )

        run_norm_active = run_norm in _SQL_RUN_METRIC_EXPR and technical_repetitions > 1

        # Run normalization is applied via a small per-(SampleID, Run)
        # scale-factor lookup table registered with DuckDB and joined into the
        # main query, rather than materialising a full-sized scaled intensity
        # CTE (which on PXD030304-scale inputs forced ~95 GB of disk spill).
        scale_df = None
        scale_join_clause = ""
        scale_select = "mu.*"
        scale_orderby = 'mu."Intensity"'
        if run_norm_active:
            scale_df = self._compute_run_scale_map(feature, run_norm)
            feature.parquet_db.register("run_scale_map", scale_df)
            scale_join_clause = (
                "LEFT JOIN run_scale_map rsm "
                'ON rsm."SampleID" = mu."SampleID" '
                'AND rsm."Run" = mu."Run"'
            )
            scale_expr = '(mu."Intensity" * COALESCE(rsm.scale_factor, 1.0))'
            scale_select = f'mu.* REPLACE ({scale_expr} AS "Intensity")'
            scale_orderby = scale_expr

        # NOTE: this UniProt mid-section extraction (`sp|UNIPROT|name` -> `UNIPROT`)
        # is duplicated in two more places. Keep all three in sync when the
        # protein-accession parsing rules change:
        #   - _compute_run_scale_map (SQL, this file)
        #   - convert_to_directlfq_format (polars, this file)
        query = _build_sqlfirst_load_query(
            where_clause,
            scale_select,
            scale_orderby,
            scale_join_clause,
        )
        query_params = [*where_params, min_aa, min_unique]
        try:
            combined_df = feature.parquet_db.execute(query, query_params).df()
        finally:
            if scale_df is not None:
                feature.parquet_db.unregister("run_scale_map")

        # Per-sample normalization (skip dataset-level methods that operate on
        # the full table later in _apply_dataset_normalization).
        if not sample_norm_method.is_dataset_level and sample_norm != "none":
            chunks = []
            for sample_id, sample_df in combined_df.groupby(
                SAMPLE_ID, observed=True, sort=False
            ):
                chunks.append(
                    sample_norm_method(sample_df.copy(), str(sample_id), med_map)
                )
            combined_df = pd.concat(chunks, ignore_index=True)

        # Match legacy column dtypes (apply_initial_filtering casts these to
        # Categorical in the pandas path; downstream consumers may rely on it).
        combined_df[CONDITION] = pd.Categorical(combined_df[CONDITION])
        combined_df[SAMPLE_ID] = pd.Categorical(combined_df[SAMPLE_ID])
        return combined_df

    def _compute_run_scale_map(self, feature: Feature, run_norm: str) -> pd.DataFrame:
        """Compute the per-``(SampleID, Run)`` multiplicative scale factor
        applied by :meth:`FeatureNormalizationMethod.normalize_runs`.

        Returns a small DataFrame (at most ``n_runs`` rows) with columns
        ``SampleID``, ``Run`` and ``scale_factor``. The TechReplicate ladder
        mirrors :func:`apply_initial_filtering` precisely: per sample, prefer
        the integer trailing run-name segment after ``_``; fall back to casting
        the full run name; fall back further to a per-sample DENSE_RANK over
        unique run names. ``scale_factor`` is ``sample_avg_metric / run_metric``
        for samples with more than one TechReplicate, and ``1.0`` otherwise
        (single-tech-rep samples pass through unchanged), so the caller can
        unconditionally multiply Intensity by ``COALESCE(scale_factor, 1.0)``.

        Compared to the previous in-query CTE chain this aggregates straight
        from the filtered base to a small ``(SampleID, tech_rep) -> metric``
        table and never materialises a full-sized scaled-intensity CTE, which
        is what drove the >95 GB peak resident set on PXD030304.
        """
        where_clause, where_params = feature.filter_builder.build_where_clause()
        min_aa = self.config.filtering.min_aa
        min_unique = (
            0
            if self.config.quantification.method.lower() == "ibaq"
            else self.config.filtering.min_unique_peptides
        )
        metric_expr = _SQL_RUN_METRIC_EXPR[run_norm]

        # NOTE: we project only the columns needed for the run-norm aggregate
        # (SampleID, Run, sequence-for-min-unique gating, Intensity) so DuckDB
        # never has to carry the peptidoform/charge/condition columns through
        # the window-function pipeline used to enforce min_unique_peptides.
        #
        # The UniProt mid-section extraction below mirrors _load_for_mokume_sqlfirst
        # and convert_to_directlfq_format (polars). Keep all three in sync.
        sql = _build_run_scale_query(where_clause, metric_expr)
        sql_params = [*where_params, min_aa, min_unique]
        return feature.parquet_db.execute(sql, sql_params).df()

    def _load_for_mokume_legacy_pandas(
        self,
        feature: Feature,
        label,
        choice,
        technical_repetitions: int,
        sample_norm: str,
        sample_norm_method: PeptideNormalizationMethod,
        med_map: dict,
    ) -> pd.DataFrame:
        """Original per-sample pandas loop, kept verbatim for fallback paths.

        Used when the SQL-first dispatcher cannot apply (TMT/iTRAQ labels or
        an active per-run normalization stage that needs row-level intensities
        before peptidoform-max).
        """
        all_peptides = []

        for samples, df in feature.iter_samples():
            df.dropna(subset=["pg_accessions"], inplace=True)

            for sample in samples:
                dataset_df = df[df["sample_accession"] == sample].copy()
                if self.config.quantification.method.lower() != "ibaq":
                    dataset_df = dataset_df[dataset_df["unique"] == 1]
                dataset_df = dataset_df[PARQUET_COLUMNS]

                dataset_df = reformat_quantms_feature_table_quant_labels(
                    dataset_df, label, choice
                )
                dataset_df = apply_initial_filtering(
                    dataset_df, self.config.filtering.min_aa, AGGREGATION_LEVEL_SAMPLE
                )

                # Filter by min unique peptides
                if self.config.quantification.method.lower() != "ibaq":
                    dataset_df = dataset_df.groupby(PROTEIN_NAME).filter(
                        lambda x: (
                            len(set(x[PEPTIDE_CANONICAL]))
                            >= self.config.filtering.min_unique_peptides
                        )
                    )

                if self.config.filtering.remove_contaminants:
                    dataset_df = remove_contaminants_entrapments_decoys(dataset_df)

                dataset_df.rename(columns={INTENSITY: NORM_INTENSITY}, inplace=True)

                # Apply run normalization
                run_norm = self.config.normalization.run_method.lower()
                if run_norm not in ("none", "") and technical_repetitions > 1:
                    run_method = FeatureNormalizationMethod.from_str(run_norm)
                    dataset_df = run_method(dataset_df, technical_repetitions)

                dataset_df = get_peptidoform_normalize_intensities(dataset_df)
                dataset_df = sum_peptidoform_intensities(
                    dataset_df, AGGREGATION_LEVEL_SAMPLE
                )

                # Apply per-sample normalization (skip dataset-level methods)
                if not sample_norm_method.is_dataset_level and sample_norm != "none":
                    dataset_df = sample_norm_method(dataset_df, sample, med_map)

                all_peptides.append(dataset_df)

        return pd.concat(all_peptides, ignore_index=True)

    def _apply_dataset_normalization(
        self,
        combined_df: pd.DataFrame,
        sample_norm_method: PeptideNormalizationMethod,
    ) -> pd.DataFrame:
        """Apply dataset-level normalization (operates on the combined table).

        Dispatches to the matching ``NormalizationStage`` method when the
        configured peptide normalization is dataset-level (Hierarchical,
        Quantile, MedianCenter, MeanCenter, Rlr, Loess, TMM).
        """
        if sample_norm_method == PeptideNormalizationMethod.Hierarchical:
            combined_df = NormalizationStage(self.config).apply_hierarchical(
                combined_df
            )
        elif sample_norm_method == PeptideNormalizationMethod.Quantile:
            combined_df = NormalizationStage(self.config).apply_quantile(combined_df)
        elif sample_norm_method == PeptideNormalizationMethod.MedianCenter:
            combined_df = NormalizationStage(self.config).apply_median_center(
                combined_df
            )
        elif sample_norm_method == PeptideNormalizationMethod.MeanCenter:
            combined_df = NormalizationStage(self.config).apply_mean_center(combined_df)
        elif sample_norm_method == PeptideNormalizationMethod.Rlr:
            combined_df = NormalizationStage(self.config).apply_rlr(combined_df)
        elif sample_norm_method == PeptideNormalizationMethod.Loess:
            combined_df = NormalizationStage(self.config).apply_loess(combined_df)
        elif sample_norm_method == PeptideNormalizationMethod.TMM:
            combined_df = NormalizationStage(self.config).apply_tmm(combined_df)

        return combined_df

    def load_for_maxlfq_directlfq(self) -> pa.Table:
        """Load MaxLFQ input through SQL-first Arrow without pandas long table."""
        filter_builder = SQLFilterBuilder(
            remove_contaminants=self.config.filtering.remove_contaminants,
            min_peptide_length=self.config.filtering.min_aa,
            require_unique=True,
        )
        feature = Feature(
            self.config.input.parquet,
            filter_builder=filter_builder,
            duckdb_memory=self.config.runtime.duckdb_memory,
            duckdb_threads=self.config.runtime.duckdb_threads,
        )

        try:
            if self.config.input.sdrf:
                feature.enrich_with_sdrf(self.config.input.sdrf)
                technical_repetitions, label, _sample_names, _choice = analyse_sdrf(
                    self.config.input.sdrf
                )
            else:
                technical_repetitions, label, _sample_names, _choice = (
                    feature.experimental_inference
                )

            run_norm = self.config.normalization.run_method.lower()
            if label != QuantificationCategory.LFQ:
                raise ValueError(
                    f"MaxLFQ DirectLFQ-streaming path requires LFQ input; got {label}."
                )
            if run_norm not in _SQL_RUN_NORM_SUPPORTED:
                raise ValueError(
                    "MaxLFQ DirectLFQ-streaming path requires SQL-supported "
                    f"run normalization; got {run_norm}."
                )

            logger.info(
                "Loading MaxLFQ DirectLFQ input with SQL-first Arrow path "
                "(label=%s, tech_reps=%d, run_norm=%s)",
                label,
                technical_repetitions,
                run_norm,
            )

            where_clause, where_params = filter_builder.build_where_clause()
            min_aa = self.config.filtering.min_aa
            min_unique = self.config.filtering.min_unique_peptides
            run_norm_active = (
                run_norm in _SQL_RUN_METRIC_EXPR and technical_repetitions > 1
            )

            scale_df = None
            scale_join_clause = ""
            scale_select = "mu.*"
            scale_orderby = 'mu."Intensity"'
            if run_norm_active:
                scale_df = self._compute_run_scale_map(feature, run_norm)
                feature.parquet_db.register("run_scale_map", scale_df)
                scale_join_clause = (
                    "LEFT JOIN run_scale_map rsm "
                    'ON rsm."SampleID" = mu."SampleID" '
                    'AND rsm."Run" = mu."Run"'
                )
                scale_expr = '(mu."Intensity" * COALESCE(rsm.scale_factor, 1.0))'
                scale_select = f'mu.* REPLACE ({scale_expr} AS "Intensity")'
                scale_orderby = scale_expr

            query = _build_sqlfirst_load_query(
                where_clause,
                scale_select,
                scale_orderby,
                scale_join_clause,
                query_suffix=_MAXLFQ_DIRECTLFQ_QUERY_SUFFIX,
            )
            query_params = [*where_params, min_aa, min_unique]
            try:
                result = feature.parquet_db.execute(query, query_params)
                reader = result.to_arrow_reader(batch_size=1_000_000)
                return pa.Table.from_batches(reader)
            finally:
                if scale_df is not None:
                    feature.parquet_db.unregister("run_scale_map")
        finally:
            feature.parquet_db.close()

    def convert_maxlfq_to_directlfq_format(self, table: pa.Table) -> pd.DataFrame:
        """Convert MaxLFQ SQL-first Arrow output to DirectLFQ wide log2 format."""
        try:
            import polars as pl
        except ImportError as exc:
            raise ImportError(
                "polars is required for MaxLFQ DirectLFQ-streaming conversion. "
                "Install it with `pip install polars` or "
                "`pip install mokume[directlfq]`."
            ) from exc

        pdf = pl.from_arrow(table)
        if isinstance(pdf, pl.Series):
            pdf = pdf.to_frame()

        summed = pdf.group_by(["protein", "ion", "sample_accession"]).agg(
            pl.col("intensity").cast(pl.Float64).sum()
        )
        wide_pl = summed.pivot(
            index=["protein", "ion"],
            on="sample_accession",
            values="intensity",
            aggregate_function=None,
        ).sort(["protein", "ion"])

        sample_cols = sorted(c for c in wide_pl.columns if c not in ("protein", "ion"))
        wide_pl = wide_pl.select(["protein", "ion", *sample_cols])
        wide_pl = wide_pl.with_columns(
            [
                pl.when(pl.col(c) == 0).then(None).otherwise(pl.col(c).log(2)).alias(c)
                for c in sample_cols
            ]
        )

        wide = wide_pl.to_pandas().set_index(["protein", "ion"])
        wide.index.names = ["protein", "ion"]
        return wide

    def load_for_directlfq(self) -> pa.Table:
        """Load and filter peptide rows for DirectLFQ as a streaming Arrow Table.

        Streams the long-format query result through DuckDB's Arrow reader
        instead of materialising it as a pandas DataFrame. On large parquets
        (e.g. PXD030304: 163M peptide rows × 5798 samples) this avoids the
        ~30+ GB pandas object overhead that previously caused OOM kills, and
        cuts wall time roughly 30% on medium datasets (validated on
        PXD017199: 51.5s → 35.4s, 11.8 GB → 7.1 GB peak RSS).

        Protein-accession parsing and the ``min_unique_peptides`` filter are
        deferred to :meth:`convert_to_directlfq_format`, which performs them
        in-place on the polars frame (zero-copy from this Arrow Table).

        Returns
        -------
        pa.Table
            Long-format Arrow Table with columns
            ``[pg_accessions, sequence, sample_accession, intensity]``.
        """
        filter_builder = SQLFilterBuilder(
            remove_contaminants=self.config.filtering.remove_contaminants,
            min_peptide_length=self.config.filtering.min_aa,
            require_unique=True,
        )

        feature = Feature(
            self.config.input.parquet,
            filter_builder=filter_builder,
            duckdb_memory=self.config.runtime.duckdb_memory,
            duckdb_threads=self.config.runtime.duckdb_threads,
        )

        if self.config.input.sdrf:
            feature.enrich_with_sdrf(self.config.input.sdrf)

        where_clause, where_params = filter_builder.build_where_clause()
        query = "".join(
            [
                "SELECT pg_accessions, sequence, sample_accession, intensity ",
                "FROM parquet_db WHERE ",
                where_clause,
            ]
        )

        try:
            result = feature.parquet_db.execute(query, where_params)
            reader = result.to_arrow_reader(batch_size=1_000_000)
            # ``from_batches`` eagerly consumes the reader, materialising every
            # batch into Arrow buffers that own their memory independently of
            # the DuckDB connection, so it is safe to close the connection in
            # the ``finally`` below even though we return the Table.
            return pa.Table.from_batches(reader)
        finally:
            feature.parquet_db.close()

    def convert_to_directlfq_format(self, table: pa.Table) -> pd.DataFrame:
        """Convert long-format Arrow Table into the DirectLFQ wide log2 frame.

        The pipeline is end-to-end polars (zero-copy from Arrow):
        protein extraction, ``min_unique_peptides`` filter, sum aggregation,
        pivot to wide, replace 0 with null, ``log2``. A single
        :py:meth:`polars.DataFrame.to_pandas` at the end yields the
        ``MultiIndex(['protein', 'ion'])``-indexed pandas frame DirectLFQ
        expects.

        Going via polars avoids materialising the long-format pandas
        DataFrame, which dominated peak RSS on large datasets, and lets the
        pivot run in parallel (PXD017199 wall: 51.5s → 18.3s, -64%).
        """
        try:
            import polars as pl
        except ImportError as exc:
            raise ImportError(
                "polars is required for the DirectLFQ format conversion. "
                "Install it with `pip install polars` or "
                "`pip install mokume[directlfq]`."
            ) from exc

        pdf = pl.from_arrow(table)
        if isinstance(pdf, pl.Series):
            pdf = pdf.to_frame()

        # Protein extraction: take pg_accessions[0]; parse UniProt mid-section
        # ('sp|UNIPROT|name' -> 'UNIPROT') when pipe-separated, else as-is.
        # The same rule is implemented in SQL inside _load_for_mokume_sqlfirst
        # and _compute_run_scale_map (this file). Keep all three in sync; do
        # NOT extract a shared Python helper, because calling it via polars
        # ``map_elements`` would force a row-wise Python call and erase the
        # 64% wall-time win that polars's columnar pivot provides.
        pdf = pdf.with_columns(
            pl.col("pg_accessions")
            .list.get(0, null_on_oob=True)
            .fill_null("")
            .alias("_first_acc")
        )
        pdf = pdf.with_columns(
            pl.when(pl.col("_first_acc").str.contains("|", literal=True))
            .then(pl.col("_first_acc").str.split("|").list.get(1, null_on_oob=True))
            .otherwise(pl.col("_first_acc"))
            .alias("protein")
        ).drop(["_first_acc", "pg_accessions"])

        min_pep = self.config.filtering.min_unique_peptides
        if min_pep > 0:
            counts = pdf.group_by("protein").agg(
                pl.col("sequence").n_unique().alias("_n_pep")
            )
            valid = counts.filter(pl.col("_n_pep") >= min_pep).get_column("protein")
            pdf = pdf.filter(pl.col("protein").is_in(valid.implode()))

        summed = pdf.group_by(["protein", "sequence", "sample_accession"]).agg(
            pl.col("intensity").cast(pl.Float64).sum()
        )

        wide_pl = summed.pivot(
            index=["protein", "sequence"],
            on="sample_accession",
            values="intensity",
            aggregate_function=None,
        ).sort(["protein", "sequence"])

        sample_cols = sorted(
            c for c in wide_pl.columns if c not in ("protein", "sequence")
        )
        wide_pl = wide_pl.select(["protein", "sequence", *sample_cols])
        wide_pl = wide_pl.with_columns(
            [
                pl.when(pl.col(c) == 0).then(None).otherwise(pl.col(c).log(2)).alias(c)
                for c in sample_cols
            ]
        )

        wide = wide_pl.to_pandas()
        wide = wide.set_index(["protein", "sequence"])
        wide.index.names = ["protein", "ion"]
        return wide

    def load_for_ratio(self) -> tuple:
        """Load PSM data and detect references for ratio quantification.

        Returns
        -------
        tuple
            (psm_df, ref_samples, sample_to_plex)
        """
        from mokume.quantification.ratio import load_psm_data
        from mokume.normalization.irs import (
            detect_pooled_from_sdrf,
            detect_reference_by_regex,
            detect_plexes_from_sdrf,
        )

        if not self.config.input.sdrf:
            raise ValueError(
                "Ratio quantification requires an SDRF file for reference "
                "sample detection. Use --sdrf to provide one."
            )

        # Detect reference samples (reuse IRS detection logic)
        ref_samples = detect_pooled_from_sdrf(self.config.input.sdrf)

        if ref_samples is None and self.config.irs.reference_samples:
            ref_samples = self.config.irs.reference_samples
            logger.info("Using explicit reference samples: %s", ref_samples)

        if ref_samples is None:
            ref_samples = detect_reference_by_regex(
                self.config.input.sdrf, self.config.irs.reference_regex
            )

        if not ref_samples:
            raise ValueError(
                "Ratio quantification requires reference samples. "
                "None detected from SDRF. Use --irs-reference-samples to specify."
            )

        # Detect plexes
        sample_to_plex = detect_plexes_from_sdrf(self.config.input.sdrf)

        # Load PSM data
        psm_df = load_psm_data(
            parquet_path=self.config.input.parquet,
            sdrf_path=self.config.input.sdrf,
            min_aa=self.config.filtering.min_aa,
            min_unique_peptides=self.config.filtering.min_unique_peptides,
            remove_contaminants=self.config.filtering.remove_contaminants,
        )

        return psm_df, ref_samples, sample_to_plex


class NormalizationStage:
    """Applies run-level and sample-level normalization."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def apply_hierarchical(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply hierarchical sample normalization."""
        logger.info("Applying hierarchical sample normalization...")

        # Convert to wide format for normalization. observed=True avoids the
        # Cartesian-product blowup when SAMPLE_ID is Categorical.
        wide = df.pivot_table(
            index=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            columns=SAMPLE_ID,
            values=NORM_INTENSITY,
            aggfunc="sum",
            observed=True,
        )

        # Log2 transform for normalization
        wide = wide.replace(0, np.nan)
        wide_log2 = np.log2(wide)

        # Load selected proteins if specified
        selected_proteins = None
        if self.config.normalization.proteins_file:
            with open(self.config.normalization.proteins_file, encoding="utf-8") as f:
                selected_proteins = [line.strip() for line in f if line.strip()]
            logger.info(
                "Using %d selected proteins for normalization", len(selected_proteins)
            )

        # Apply hierarchical normalization
        normalizer = HierarchicalSampleNormalizer(
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
            selected_proteins=selected_proteins,
        )

        normalized_log2 = normalizer.fit_transform(wide_log2)

        # Convert back to linear scale
        normalized_wide = 2**normalized_log2

        # Convert back to long format
        normalized_long = normalized_wide.reset_index().melt(
            id_vars=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            var_name=SAMPLE_ID,
            value_name=NORM_INTENSITY,
        )

        # Remove NaN rows
        normalized_long = normalized_long.dropna(subset=[NORM_INTENSITY])

        logger.info(
            "Hierarchical normalization complete: %d rows", len(normalized_long)
        )

        return normalized_long

    def _apply_dataset_normalizer(
        self,
        df: pd.DataFrame,
        normalizer,
        log_space: bool,
        name: str,
    ) -> pd.DataFrame:
        """Pivot peptide-level long DataFrame to wide, fit_transform, melt back."""
        logger.info("Applying %s sample normalization...", name)
        wide = df.pivot_table(
            index=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            columns=SAMPLE_ID,
            values=NORM_INTENSITY,
            aggfunc="sum",
            observed=True,
        )
        if log_space:
            wide = np.log2(wide.replace(0, np.nan))
        normalized_wide = normalizer.fit_transform(wide)
        if log_space:
            normalized_wide = 2**normalized_wide
        normalized_long = (
            normalized_wide.reset_index()
            .melt(
                id_vars=[PROTEIN_NAME, PEPTIDE_CANONICAL],
                var_name=SAMPLE_ID,
                value_name=NORM_INTENSITY,
            )
            .dropna(subset=[NORM_INTENSITY])
        )
        logger.info("%s normalization complete: %d rows", name, len(normalized_long))
        return normalized_long

    def apply_quantile(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply quantile sample normalization."""
        from mokume.normalization.distribution import QuantileNormalizer

        return self._apply_dataset_normalizer(
            df, QuantileNormalizer(), log_space=False, name="Quantile"
        )

    def apply_median_center(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply sample-level median centering in log2 space."""
        from mokume.normalization.distribution import MedianCenterNormalizer

        return self._apply_dataset_normalizer(
            df, MedianCenterNormalizer(), log_space=True, name="MedianCenter"
        )

    def apply_mean_center(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply sample-level mean centering in log2 space."""
        from mokume.normalization.distribution import MeanCenterNormalizer

        return self._apply_dataset_normalizer(
            df, MeanCenterNormalizer(), log_space=True, name="MeanCenter"
        )

    def apply_rlr(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply Robust Linear Regression normalization (NormalyzerDE)."""
        from mokume.normalization.rlr import RlrNormalizer

        return self._apply_dataset_normalizer(
            df,
            RlrNormalizer(log_transform=False),
            log_space=False,
            name="RLR",
        )

    def apply_loess(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply LOESS regression-based normalization (in log2 space)."""
        from mokume.normalization.loess import LOESSNormalizer

        return self._apply_dataset_normalizer(
            df, LOESSNormalizer(), log_space=True, name="LOESS"
        )

    def apply_tmm(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply Trimmed Mean of M-values normalization (on linear intensities)."""
        from mokume.normalization.tmm import TMMNormalizer

        return self._apply_dataset_normalizer(
            df, TMMNormalizer(), log_space=False, name="TMM"
        )

    def apply_irs(self, protein_df: pd.DataFrame, dataset=None) -> pd.DataFrame:
        """Apply Internal Reference Scaling normalization for multi-plex TMT data.

        The optional ``dataset`` parameter is accepted for compatibility with
        the metadata-driven pipeline runner (``pipeline.runner``); it is not
        used by the IRS algorithm itself.

        Detection priority for reference samples:
        1. characteristics[pooled sample] column in SDRF
        2. Explicit sample names (--irs-reference-samples)
        3. Explicit column + values (--irs-sdrf-column + --irs-sdrf-values)
        4. Regex scan across factor/characteristic columns
        """
        from mokume.normalization.irs import (
            IRSNormalizer,
            detect_pooled_from_sdrf,
            detect_reference_by_column,
            detect_reference_by_regex,
            detect_plexes_from_sdrf,
        )

        if not self.config.input.sdrf:
            raise ValueError("IRS normalization requires an SDRF file (--sdrf)")

        # Detect reference samples (priority order)
        ref_samples = None

        # 1. Check for characteristics[pooled sample] column
        ref_samples = detect_pooled_from_sdrf(self.config.input.sdrf)

        # 2. Explicit sample names override
        if ref_samples is None and self.config.irs.reference_samples:
            ref_samples = self.config.irs.reference_samples
            logger.info("Using explicit reference samples: %s", ref_samples)

        # 3. Explicit column + values
        if (
            ref_samples is None
            and self.config.irs.sdrf_column
            and self.config.irs.sdrf_values
        ):
            ref_samples = detect_reference_by_column(
                self.config.input.sdrf,
                self.config.irs.sdrf_column,
                self.config.irs.sdrf_values,
            )

        # 4. Regex fallback
        if ref_samples is None:
            ref_samples = detect_reference_by_regex(
                self.config.input.sdrf, self.config.irs.reference_regex
            )

        if not ref_samples:
            logger.warning("No reference samples detected for IRS, skipping")
            return protein_df

        # Detect plexes
        sample_to_plex = detect_plexes_from_sdrf(self.config.input.sdrf)

        # Apply IRS
        normalizer = IRSNormalizer(
            reference_samples=ref_samples, stat=self.config.irs.stat
        )
        protein_df = normalizer.fit_transform(protein_df, sample_to_plex)

        logger.info("IRS normalization complete: %d proteins", len(protein_df))

        # Optionally remove reference samples from output
        if self.config.irs.remove_reference:
            protein_col = protein_df.columns[0]
            cols_to_keep = [protein_col] + [
                c
                for c in protein_df.columns
                if c == protein_col or c not in ref_samples
            ]
            # Deduplicate while preserving order
            seen = set()
            unique_cols = []
            for c in cols_to_keep:
                if c not in seen:
                    seen.add(c)
                    unique_cols.append(c)
            protein_df = protein_df[unique_cols]
            logger.info("Removed %d reference samples from output", len(ref_samples))

        return protein_df

    def apply_coverage_filter(
        self, protein_df: pd.DataFrame, dataset=None
    ) -> pd.DataFrame:
        """Apply coverage filter using condition mapping from SDRF.

        The optional ``dataset`` parameter is accepted for compatibility with
        the pipeline runner (``pipeline.runner``) and is otherwise unused.
        """
        from mokume.quantification.ratio import apply_coverage_filter
        from mokume.normalization.irs import detect_condition_from_sdrf

        if not self.config.input.sdrf:
            logger.warning("Coverage filter requires SDRF file, skipping")
            return protein_df

        sample_to_condition = detect_condition_from_sdrf(self.config.input.sdrf)

        # Filter out reference/powder conditions for coverage computation
        protein_col = protein_df.columns[0]
        available_samples = [c for c in protein_df.columns if c != protein_col]
        sample_to_condition = {
            s: c
            for s, c in sample_to_condition.items()
            if s in available_samples
            and "powder" not in c.lower()
            and "pool" not in c.lower()
        }

        return apply_coverage_filter(
            protein_df,
            sample_to_condition,
            self.config.quantification.coverage_threshold,
        )


class QuantificationStage:
    """Routes to the appropriate quantification method."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def quantify(self, peptide_df: pd.DataFrame) -> pd.DataFrame:
        """Apply protein quantification method."""
        from mokume.quantification import get_quantification_method, TopNQuantification

        quant_method = self.config.quantification.method.lower()

        if quant_method == "ibaq":
            return self._quantify_ibaq(peptide_df)
        if quant_method in (
            "maxlfq",
            "sum",
            "all",
            "abd",
            "abundance",
            "intensity",
            "reporter",
            "spectral_count",
            "spectralcount",
            "count",
        ):
            # Propagate the global parallelism budget so MaxLFQ -> DirectLFQ
            # does not silently default to cpu_count workers (each forked
            # worker COW-copies the wide pivot table; on PXD030304-scale
            # inputs that path was OOM-killing the 125 GB host).
            # RuntimeConfig.effective_workers() also reasons about
            # ``duckdb_memory`` so the worker count drops automatically when
            # the memory budget cannot absorb cpu_count concurrent fork-COW
            # copies.
            quant_kwargs = {}
            effective = self.config.runtime.effective_workers()
            if effective is not None:
                quant_kwargs["n_jobs"] = effective
                quant_kwargs["num_cores"] = effective
            method = get_quantification_method(quant_method, **quant_kwargs)
            result = method.quantify(
                peptide_df,
                protein_column=PROTEIN_NAME,
                peptide_column=PEPTIDE_CANONICAL,
                intensity_column=NORM_INTENSITY,
                sample_column=SAMPLE_ID,
            )
            return self.to_wide_format(result, quant_method)
        if quant_method.startswith("top"):
            match = re.match(r"top(\d+)", quant_method)
            if match:
                n = int(match.group(1))
            else:
                n = 3
            method = TopNQuantification(n=n)
            # TopN itself is single-process; nothing extra to thread.
            result = method.quantify(
                peptide_df,
                protein_column=PROTEIN_NAME,
                peptide_column=PEPTIDE_CANONICAL,
                intensity_column=NORM_INTENSITY,
                sample_column=SAMPLE_ID,
            )
            return self.to_wide_format(result, quant_method)
        if quant_method == "median":
            return self._quantify_median(peptide_df)
        raise ValueError(f"Unknown quantification method: {quant_method}")

    def to_wide_format(self, long_df: pd.DataFrame, method_name: str) -> pd.DataFrame:
        """Convert long format quantification results to wide format."""
        _ = method_name
        intensity_col = "Intensity"

        if intensity_col not in long_df.columns:
            # Fallback: use the last numeric column
            numeric_cols = long_df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                intensity_col = numeric_cols[-1]

        logger.debug("Using intensity column: %s", intensity_col)

        # Pivot to wide format
        wide_df = long_df.pivot(
            index=PROTEIN_NAME,
            columns=SAMPLE_ID,
            values=intensity_col,
        )

        return wide_df.reset_index()

    def _quantify_ibaq(self, peptide_df: pd.DataFrame) -> pd.DataFrame:
        """Quantify using piBAQ (paralog-aware iBAQ).

        Delegates to :func:`mokume.quantification.ibaq.compute_pibaq`, the
        same core used by the standalone ``peptides2protein`` CLI. The two
        entry points share that core, so they agree on iBAQ values **when
        configured identically** -- the digestion enzyme, ``min_aa`` /
        ``max_aa``, family ``min_shared`` and any YAML family overrides must
        all match. Those knobs are sourced here from
        :class:`~mokume.pipeline.config.QuantificationConfig`
        (``ibaq_enzyme`` / ``ibaq_max_aa`` / ``ibaq_min_shared`` /
        ``ibaq_families_yaml`` / ``ibaq_min_anchors`` /
        ``ibaq_high_anchor_threshold``) and :attr:`FilterConfig.min_aa`.

        TPA is intentionally **not** emitted on this path: the pipeline
        contract returns a single wide protein x sample iBAQ matrix, which
        has no column for a parallel TPA value. Use the ``peptides2protein``
        CLI with ``--tpa`` when a TPA table is needed. Accordingly
        ``mw_map`` is left ``None`` here.
        """
        from mokume.io.fasta import digest_fasta_full
        from mokume.quantification.ibaq import compute_pibaq
        from mokume.quantification.families import (
            discover_families,
            load_families_yaml,
            merge_overrides,
        )

        if not self.config.input.fasta_file:
            raise ValueError(
                "iBAQ quantification requires a FASTA file. Use --fasta to provide one."
            )

        quant_cfg = self.config.quantification
        logger.info(
            "Computing piBAQ for %d proteins using FASTA "
            "(enzyme=%s, max_aa=%d, min_shared=%d)...",
            peptide_df[PROTEIN_NAME].nunique(),
            quant_cfg.ibaq_enzyme,
            quant_cfg.ibaq_max_aa,
            quant_cfg.ibaq_min_shared,
        )

        accession_to_peptides, peptide_to_accessions, _ = digest_fasta_full(
            fasta=self.config.input.fasta_file,
            enzyme=quant_cfg.ibaq_enzyme,
            min_aa=self.config.filtering.min_aa,
            max_aa=quant_cfg.ibaq_max_aa,
            canonicalize_isoforms=True,
            compute_mw=False,
        )
        families = discover_families(
            accession_to_peptides,
            peptide_to_accessions,
            min_shared=quant_cfg.ibaq_min_shared,
        )
        if quant_cfg.ibaq_families_yaml:
            overrides = load_families_yaml(Path(quant_cfg.ibaq_families_yaml))
            families = merge_overrides(families, overrides)

        long_df = compute_pibaq(
            peptide_df,
            accession_to_peptides,
            peptide_to_accessions,
            families,
            mw_map=None,
            min_anchors=quant_cfg.ibaq_min_anchors,
            high_anchor_threshold=quant_cfg.ibaq_high_anchor_threshold,
        )

        if long_df.empty:
            logger.info("iBAQ complete: 0 proteins")
            return pd.DataFrame(columns=[PROTEIN_NAME])

        result_wide = long_df.pivot(
            index=PROTEIN_NAME,
            columns=SAMPLE_ID,
            values="Ibaq",
        )
        logger.info("iBAQ complete: %d proteins", len(result_wide))
        return result_wide.reset_index()

    def _quantify_median(self, peptide_df: pd.DataFrame) -> pd.DataFrame:
        """Quantify using median of peptides."""
        result = (
            peptide_df.groupby([PROTEIN_NAME, SAMPLE_ID], observed=False)[
                NORM_INTENSITY
            ]
            .median()
            .reset_index()
        )

        result_wide = result.pivot(
            index=PROTEIN_NAME, columns=SAMPLE_ID, values=NORM_INTENSITY
        )

        return result_wide.reset_index()


class ImputationStage:
    """Impute missing values on the wide protein-level matrix.

    The matrix is transformed to log2 space before imputation and back to
    linear space afterwards because most imputation methods assume
    log-normal intensities (notably MinProb/MinDet/QRILC).
    """

    _SUPPORTED_METHODS = (
        "none",
        "knn",
        "minprob",
        "mindet",
        "qrilc",
        "missforest",
        "seqknn",
        "impseq",
    )

    def __init__(self, config: PipelineConfig):
        self.config = config

    def impute(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """Apply imputation to the wide protein matrix."""
        if not self.config.imputation.enabled:
            return protein_df

        method = self.config.imputation.method.lower()
        if method == "none":
            return protein_df

        if method not in self._SUPPORTED_METHODS:
            raise ValueError(
                f"Unknown imputation method: {method}. "
                f"Supported: {', '.join(self._SUPPORTED_METHODS)}"
            )

        if PROTEIN_NAME in protein_df.columns:
            protein_col = PROTEIN_NAME
        elif "protein" in protein_df.columns:
            protein_col = "protein"
        else:
            protein_col = protein_df.columns[0]

        sample_cols = [c for c in protein_df.columns if c != protein_col]
        wide = protein_df.set_index(protein_col)[sample_cols]

        wide_log2 = np.log2(wide.replace(0, np.nan))
        n_before = int(wide_log2.isna().sum().sum())
        if n_before == 0:
            logger.info("Imputation (%s): no missing values, skipping", method)
            return protein_df

        imputed_log2 = self._dispatch(method, wide_log2)
        n_after = int(imputed_log2.isna().sum().sum())
        logger.info(
            "Imputation (%s): %d -> %d missing values", method, n_before, n_after
        )

        imputed_linear = 2**imputed_log2
        return imputed_linear.reset_index()

    def _dispatch(self, method: str, wide: pd.DataFrame) -> pd.DataFrame:
        """Dispatch to the matching imputation backend."""
        cfg = self.config.imputation

        if method == "knn":
            from sklearn.impute import KNNImputer

            imputer = KNNImputer(n_neighbors=cfg.n_neighbors, keep_empty_features=True)
            arr = imputer.fit_transform(wide)
            return pd.DataFrame(arr, index=wide.index, columns=wide.columns)

        if method == "minprob":
            from mokume.imputation.censored import impute_minprob

            return impute_minprob(
                wide, quantile=cfg.quantile, shift=cfg.shift, scale=cfg.scale
            )

        if method == "mindet":
            from mokume.imputation.censored import impute_mindet

            return impute_mindet(wide, quantile=cfg.quantile)

        if method == "qrilc":
            from mokume.imputation.qrilc import impute_qrilc

            return impute_qrilc(wide)

        if method == "missforest":
            from mokume.imputation.missforest import impute_missforest

            return impute_missforest(wide)

        if method == "seqknn":
            from mokume.imputation.seqknn import impute_seqknn

            return impute_seqknn(wide, k=cfg.n_neighbors)

        if method == "impseq":
            from mokume.imputation.impseq import impute_impseq

            return impute_impseq(wide)

        raise ValueError(f"Unknown imputation method: {method}")


class PostprocessingStage:
    """Batch correction, DE, plotting, reports."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def _protein_and_sample_columns(self, protein_df: pd.DataFrame) -> tuple[str, list]:
        """Return the protein identifier column and quantitative sample columns."""
        if PROTEIN_NAME in protein_df.columns:
            protein_col = PROTEIN_NAME
        elif "protein" in protein_df.columns:
            protein_col = "protein"
        else:
            protein_col = protein_df.columns[0]

        return protein_col, [c for c in protein_df.columns if c != protein_col]

    def _batch_method(self) -> BatchDetectionMethod:
        """Resolve the configured batch-detection method."""
        try:
            return BatchDetectionMethod.from_str(self.config.batch.method)
        except ValueError:
            logger.warning(
                "Unknown batch method '%s', using sample_prefix",
                self.config.batch.method,
            )
            return BatchDetectionMethod.SAMPLE_PREFIX

    def _detect_batch_indices(self, sample_cols: list) -> Optional[list]:
        """Detect and validate batch assignments for sample columns."""
        batch_indices = detect_batches(
            sample_ids=sample_cols,
            method=self._batch_method(),
            batch_column_values=(
                self._get_batch_column_values(sample_cols)
                if self.config.batch.column
                else None
            ),
        )

        unique_batches = len(set(batch_indices))
        logger.info("Detected %d batches for batch correction", unique_batches)
        if unique_batches < 2:
            logger.warning("Only 1 batch detected, skipping batch correction")
            return None

        min_samples = min(Counter(batch_indices).values())
        if min_samples < 2:
            logger.warning(
                "Some batches have fewer than 2 samples (min=%d), "
                "skipping batch correction",
                min_samples,
            )
            return None
        return batch_indices

    def _batch_covariates(self, sample_cols: list) -> Optional[list]:
        """Extract configured SDRF covariates for batch correction."""
        if not (self.config.input.sdrf and self.config.batch.covariates):
            return None

        covariates = extract_covariates_from_sdrf(
            self.config.input.sdrf,
            sample_cols,
            self.config.batch.covariates,
        )
        if covariates:
            logger.info(
                "Extracted %d covariates to preserve biological signal",
                len(self.config.batch.covariates),
            )
        return covariates

    def _complete_batch_matrix(
        self,
        intensity_matrix: pd.DataFrame,
    ) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:
        """Split complete and incomplete rows before ComBat correction."""
        has_nan = intensity_matrix.isna().any(axis=1)
        complete_matrix = intensity_matrix[~has_nan]
        incomplete_matrix = intensity_matrix[has_nan]

        if len(complete_matrix) == 0:
            logger.warning("No proteins with complete data for ComBat, skipping")
            return None

        if len(incomplete_matrix) > 0:
            logger.info(
                "ComBat: %d complete proteins, %d with missing values "
                "(kept uncorrected)",
                len(complete_matrix),
                len(incomplete_matrix),
            )
        return complete_matrix, incomplete_matrix

    def _rebuild_corrected_df(
        self,
        corrected_matrix: pd.DataFrame,
        incomplete_matrix: pd.DataFrame,
        protein_col: str,
    ) -> pd.DataFrame:
        """Recombine corrected complete rows with uncorrected incomplete rows."""
        if len(incomplete_matrix) > 0:
            corrected_matrix = pd.concat([corrected_matrix, incomplete_matrix])

        corrected_df = corrected_matrix.reset_index()
        return corrected_df.rename(columns={"index": protein_col})

    def apply_batch_correction(
        self, protein_df: pd.DataFrame, dataset=None
    ) -> pd.DataFrame:
        """Apply batch correction after protein quantification.

        Batch = technical variation to REMOVE (from runs/files)
        Covariates = biological signal to PRESERVE (from SDRF characteristics)

        The optional ``dataset`` parameter is accepted for compatibility with
        the pipeline runner (``pipeline.runner``) and is otherwise unused.
        """
        if not is_batch_correction_available():
            raise ImportError(
                "Batch correction requires inmoose package. "
                "Install with: pip install mokume[batch-correction]"
            )

        protein_col, sample_cols = self._protein_and_sample_columns(protein_df)
        if len(sample_cols) < 2:
            logger.warning("Not enough samples for batch correction, skipping")
            return protein_df

        intensity_matrix = protein_df.set_index(protein_col)[sample_cols]
        batch_indices = self._detect_batch_indices(sample_cols)
        if batch_indices is None:
            return protein_df

        matrices = self._complete_batch_matrix(intensity_matrix)
        if matrices is None:
            return protein_df
        complete_matrix, incomplete_matrix = matrices

        logger.info("Applying ComBat batch correction...")
        try:
            corrected_matrix = apply_batch_correction(
                df=complete_matrix,
                batch=batch_indices,
                covs=self._batch_covariates(sample_cols),
                kwargs={
                    "par_prior": self.config.batch.parametric,
                    "mean_only": self.config.batch.mean_only,
                    "ref_batch": self.config.batch.ref_batch,
                },
            )

            corrected_df = self._rebuild_corrected_df(
                corrected_matrix,
                incomplete_matrix,
                protein_col,
            )
            logger.info(
                "Batch correction complete: %d proteins, %d samples",
                len(corrected_df),
                len(sample_cols),
            )
            return corrected_df

        except Exception as e:
            logger.error("Batch correction failed: %s", e)
            logger.warning("Returning uncorrected protein intensities")
            return protein_df

    def _parse_de_contrasts(self, sample_to_condition: dict) -> list[tuple[str, str]]:
        """Parse configured DE contrasts into condition pairs."""
        if not self.config.de.contrasts:
            conditions = sorted(set(sample_to_condition.values()))
            raise ValueError(
                "Differential expression requires explicit contrasts "
                "via --de-contrasts.\n"
                f"Available conditions ({len(conditions)}): {conditions}\n"
                "Format: --de-contrasts 'GroupA vs GroupB,GroupA vs GroupC'"
            )

        contrasts = []
        for contrast_text in self.config.de.contrasts:
            if " vs " in contrast_text:
                parts = contrast_text.split(" vs ", 1)
            elif "-" in contrast_text:
                parts = contrast_text.split("-", 1)
            else:
                logger.warning(
                    "Invalid contrast format '%s', expected 'A vs B' or 'A-B'",
                    contrast_text,
                )
                continue
            contrasts.append((parts[0].strip(), parts[1].strip()))
        return contrasts

    def _resolve_de_method(self) -> str:
        """Resolve the configured DE method, including auto selection."""
        de_method = self.config.de.method
        if de_method != "auto":
            return de_method

        quant = self.config.quantification.method.lower()
        de_method = "deqms" if quant == "directlfq" else "limrots"
        logger.info("Auto-selected DE method: %s (quant=%s)", de_method, quant)
        return de_method

    def _load_de_peptide_counts(self, de_method: str) -> Optional[pd.Series]:
        """Load protein-level unique peptide counts when DEqMS needs them."""
        needs_counts = de_method == "deqms" or (
            de_method == "ensemble"
            and self.config.de.ensemble_methods
            and "deqms" in self.config.de.ensemble_methods
        )
        if not (needs_counts and self.config.input.parquet):
            return None

        try:
            pep_df = pd.read_parquet(
                self.config.input.parquet,
                columns=["anchor_protein", "sequence"],
            )
            return pep_df.groupby("anchor_protein")["sequence"].nunique()
        except (FileNotFoundError, KeyError, ValueError) as exc:
            logger.warning("Could not load peptide counts for DEqMS: %s", exc)
            return None

    def _de_output_path(self, contrast_key: str, n_contrasts: int) -> Optional[str]:
        """Return the output path for a DE contrast, if output is enabled."""
        if not self.config.de.output:
            return None
        if n_contrasts == 1:
            return self.config.de.output

        base, ext = os.path.splitext(self.config.de.output)
        return f"{base}_{contrast_key}{ext}" if ext else f"{base}_{contrast_key}.csv"

    def _save_de_result(
        self,
        result: pd.DataFrame,
        contrast_key: str,
        n_contrasts: int,
    ) -> None:
        """Save a DE result table when configured."""
        output_path = self._de_output_path(contrast_key, n_contrasts)
        if output_path is None:
            return
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_path, index=False)
        logger.info("DE results saved to %s", output_path)

    def _run_ensemble_de(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: dict,
        contrasts: list[tuple[str, str]],
        peptide_counts: Optional[pd.Series],
    ) -> dict:
        """Run ensemble DE for all configured contrasts."""
        from mokume.analysis.ensemble import run_ensemble

        ensemble_methods = self.config.de.ensemble_methods or [
            "limrots",
            "deqms",
            "proda",
        ]
        all_results = {}
        for contrast in contrasts:
            key = f"{contrast[0]}-{contrast[1]}"
            logger.info(
                "Running ensemble DE: %s (methods=%s, min_k=%d)",
                key,
                ensemble_methods,
                self.config.de.ensemble_min_k,
            )
            result = run_ensemble(
                protein_df,
                sample_to_condition,
                contrast,
                methods=ensemble_methods,
                min_k=self.config.de.ensemble_min_k,
                fdr_method=self.config.de.fdr_method,
                fdr_threshold=self.config.de.fdr_threshold,
                log2fc_threshold=self.config.de.log2fc_threshold,
                peptide_counts=peptide_counts,
            )
            all_results[key] = {"df": result, "conditions": contrast}
            self._save_de_result(result, key, len(contrasts))
        return all_results

    def _run_single_de_method(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: dict,
        contrasts: list[tuple[str, str]],
        de_method: str,
        peptide_counts: Optional[pd.Series],
    ) -> dict:
        """Run a non-ensemble DE method for all configured contrasts."""
        from mokume.analysis.differential_expression import DifferentialExpression

        de = DifferentialExpression(
            method=de_method,
            log2fc_threshold=self.config.de.log2fc_threshold,
            fdr_threshold=self.config.de.fdr_threshold,
            fdr_method=self.config.de.fdr_method,
            skip_log2=(self.config.quantification.method.lower() == "ratio"),
            peptide_counts=peptide_counts,
        )

        all_results = {}
        for contrast in contrasts:
            key = f"{contrast[0]}-{contrast[1]}"
            logger.info("Running DE: %s", key)
            result = de.run(protein_df, sample_to_condition, contrast)
            all_results[key] = {"df": result, "conditions": contrast}
            self._save_de_result(result, key, len(contrasts))
        return all_results

    def run_differential_expression(
        self, protein_df: pd.DataFrame, dataset=None
    ) -> Optional[dict]:
        """Run differential expression analysis.

        The optional ``dataset`` parameter is accepted for compatibility with
        the pipeline runner (``pipeline.runner``) and is otherwise unused.
        """
        from mokume.normalization.irs import detect_condition_from_sdrf

        if not self.config.input.sdrf:
            raise ValueError("Differential expression requires an SDRF file (--sdrf)")

        sample_to_condition = detect_condition_from_sdrf(self.config.input.sdrf)
        contrasts = self._parse_de_contrasts(sample_to_condition)
        if not contrasts:
            return None

        de_method = self._resolve_de_method()
        peptide_counts = self._load_de_peptide_counts(de_method)

        if de_method == "ensemble":
            return self._run_ensemble_de(
                protein_df,
                sample_to_condition,
                contrasts,
                peptide_counts,
            )

        return self._run_single_de_method(
            protein_df,
            sample_to_condition,
            contrasts,
            de_method,
            peptide_counts,
        )

    def _plot_sample_conditions(self, protein_df: pd.DataFrame) -> Optional[dict]:
        """Return sample-to-condition mapping for plotting when SDRF is available."""
        from mokume.normalization.irs import detect_condition_from_sdrf

        if not self.config.input.sdrf:
            return None

        sample_to_condition = detect_condition_from_sdrf(self.config.input.sdrf)
        if not self.config.irs.remove_reference:
            return sample_to_condition

        protein_col = protein_df.columns[0]
        available_samples = [c for c in protein_df.columns if c != protein_col]
        return {
            sample: condition
            for sample, condition in sample_to_condition.items()
            if sample in available_samples
        }

    def _generate_volcano_plots(
        self, de_results: Optional[dict], plot_dir: Path
    ) -> None:
        """Generate configured volcano plots."""
        if not (self.config.output.plot_volcano and de_results):
            return

        from mokume.plotting.differential_expression import plot_volcano

        for contrast_name, de_entry in de_results.items():
            output_file = str(plot_dir / f"volcano_{contrast_name}.png")
            plot_volcano(
                de_entry["df"],
                log2fc_threshold=self.config.de.log2fc_threshold,
                fdr_threshold=self.config.de.fdr_threshold,
                highlight_genes=self.config.output.highlight_genes,
                title=f"Volcano Plot: {contrast_name}",
                output_file=output_file,
            )
            logger.info("Volcano plot saved to %s", output_file)

    def _significant_de_proteins(self, de_df: pd.DataFrame) -> pd.DataFrame:
        """Return significant DE rows under configured thresholds."""
        return de_df[
            (de_df["adj_pvalue"] < self.config.de.fdr_threshold)
            & (de_df["log2FC"].abs() > self.config.de.log2fc_threshold)
        ]

    def _generate_heatmap_plots(
        self,
        protein_df: pd.DataFrame,
        de_results: Optional[dict],
        sample_to_condition: Optional[dict],
        plot_dir: Path,
    ) -> None:
        """Generate configured per-contrast DE heatmaps."""
        if not (self.config.output.plot_heatmap and sample_to_condition and de_results):
            return

        from mokume.plotting.differential_expression import plot_heatmap

        protein_col = protein_df.columns[0]
        for contrast_name, de_entry in de_results.items():
            sig = self._significant_de_proteins(de_entry["df"])
            if sig.empty:
                logger.info(
                    "Heatmap skipped for %s: no significant proteins",
                    contrast_name,
                )
                continue

            cond_a, cond_b = de_entry["conditions"]
            sig_proteins = (
                sig.sort_values(by="log2FC", key=lambda x: x.abs(), ascending=False)[
                    "ProteinName"
                ]
                .head(50)
                .tolist()
            )
            contrast_samples = [
                sample
                for sample in protein_df.columns
                if sample != protein_col
                and sample_to_condition.get(sample) in (cond_a, cond_b)
            ]
            contrast_mapping = {
                sample: condition
                for sample, condition in sample_to_condition.items()
                if condition in (cond_a, cond_b)
            }
            output_file = str(plot_dir / f"heatmap_{contrast_name}.png")
            plot_heatmap(
                protein_df[[protein_col] + contrast_samples],
                contrast_mapping,
                proteins=sig_proteins,
                title=(
                    f"DE Heatmap: {contrast_name} "
                    f"(top {len(sig_proteins)}/{len(sig)} sig)"
                ),
                output_file=output_file,
            )
            logger.info("Heatmap saved to %s", output_file)

    def _generate_pca_plot(
        self,
        protein_df: pd.DataFrame,
        sample_to_condition: Optional[dict],
        plot_dir: Path,
    ) -> None:
        """Generate configured PCA-by-condition plot."""
        if not (self.config.output.plot_pca and sample_to_condition):
            return

        from mokume.plotting.differential_expression import plot_pca_conditions

        output_file = str(plot_dir / "pca_conditions.png")
        plot_pca_conditions(
            protein_df,
            sample_to_condition,
            title="PCA by Condition",
            output_file=output_file,
        )
        logger.info("PCA plot saved to %s", output_file)

    def generate_plots(
        self,
        protein_df: pd.DataFrame,
        de_results: Optional[dict] = None,
        dataset=None,
    ):
        """Generate requested plots.

        The optional ``dataset`` parameter is accepted for compatibility with
        the pipeline runner (``pipeline.runner``) and is otherwise unused.
        """
        from mokume.plotting import is_plotting_available

        if not is_plotting_available():
            logger.warning(
                "Plotting dependencies not available. "
                "Install with: pip install mokume[plotting]"
            )
            return

        plot_dir = Path(self.config.output.plot_dir)
        plot_dir.mkdir(parents=True, exist_ok=True)

        sample_to_condition = self._plot_sample_conditions(protein_df)
        self._generate_volcano_plots(de_results, plot_dir)
        self._generate_heatmap_plots(
            protein_df,
            de_results,
            sample_to_condition,
            plot_dir,
        )
        self._generate_pca_plot(protein_df, sample_to_condition, plot_dir)

    def generate_interactive_report(
        self, protein_df: pd.DataFrame, de_results: dict, dataset=None
    ):
        """Generate interactive HTML report for DE results.

        The optional ``dataset`` parameter is accepted for compatibility with
        the pipeline runner (``pipeline.runner``) and is otherwise unused.
        """
        from mokume.reports import is_interactive_available

        if not is_interactive_available():
            logger.warning(
                "Interactive report dependencies (plotly) not available. "
                "Install with: pip install mokume[reports]"
            )
            return

        from mokume.reports.interactive import generate_de_report
        from mokume.normalization.irs import detect_condition_from_sdrf

        if not self.config.input.sdrf:
            logger.warning("Interactive report requires SDRF file, skipping")
            return

        sample_to_condition = detect_condition_from_sdrf(self.config.input.sdrf)

        for contrast_name, de_entry in de_results.items():
            de_df = de_entry["df"] if isinstance(de_entry, dict) else de_entry
            if self.config.output.report_output:
                if len(de_results) == 1:
                    output_html = self.config.output.report_output
                else:
                    base = self.config.output.report_output.rsplit(".", 1)
                    output_html = f"{base[0]}_{contrast_name}.html"
            else:
                plot_dir = self.config.output.plot_dir or "."
                output_html = str(Path(plot_dir) / f"report_{contrast_name}.html")

            generate_de_report(
                de_results=de_df,
                protein_df=protein_df,
                sample_to_condition=sample_to_condition,
                output_html=output_html,
                title=f"DE Report: {contrast_name}",
                highlight_genes=self.config.output.highlight_genes,
                log2fc_threshold=self.config.de.log2fc_threshold,
                fdr_threshold=self.config.de.fdr_threshold,
            )
            logger.info("Interactive report saved to %s", output_html)

    def _get_batch_column_values(self, sample_ids: list) -> Optional[list]:
        """Get batch values from SDRF for explicit batch column."""
        if not self.config.input.sdrf or not self.config.batch.column:
            return None

        try:
            from mokume.core.constants import load_sdrf as _load_sdrf

            sdrf = _load_sdrf(self.config.input.sdrf)

            batch_col = self.config.batch.column.lower()
            if batch_col not in sdrf.columns:
                logger.warning(
                    "Batch column '%s' not in SDRF",
                    self.config.batch.column,
                )
                return None

            # Map sample IDs to batch values
            sample_to_batch = dict(zip(sdrf["source name"], sdrf[batch_col]))
            return [sample_to_batch.get(s, "unknown") for s in sample_ids]

        except Exception as e:
            logger.warning("Failed to extract batch column: %s", e)
            return None
