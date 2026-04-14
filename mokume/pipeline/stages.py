"""
Pipeline stage classes for the quantification pipeline.

Each stage handles a distinct phase of the proteomics quantification workflow:
- LoadingStage: Data loading and filtering
- NormalizationStage: Run-level and sample-level normalization
- QuantificationStage: Protein quantification methods
- PostprocessingStage: Batch correction, DE, plotting, reports
"""

import os

import numpy as np
import pandas as pd
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
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
    INTENSITY,
    PARQUET_COLUMNS,
    AGGREGATION_LEVEL_SAMPLE,
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


class LoadingStage:
    """Loads and filters data from parquet + SDRF."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def load_for_mokume(self) -> pd.DataFrame:
        """Load data and apply normalization for mokume quantification methods."""
        filter_builder = SQLFilterBuilder(
            remove_contaminants=self.config.filtering.remove_contaminants,
            min_peptide_length=self.config.filtering.min_aa,
            require_unique=True,
        )

        feature = Feature(self.config.input.parquet, filter_builder=filter_builder)

        if self.config.input.sdrf:
            feature.enrich_with_sdrf(self.config.input.sdrf)
            technical_repetitions, label, sample_names, choice = analyse_sdrf(
                self.config.input.sdrf
            )
        else:
            (
                technical_repetitions,
                label,
                sample_names,
                choice,
            ) = feature.experimental_inference

        # Get normalization factors if needed
        med_map = {}
        sample_norm = self.config.normalization.sample_method.lower()
        sample_norm_method = PeptideNormalizationMethod.from_str(sample_norm)

        if sample_norm == "globalmedian":
            med_map = feature.get_median_map()
        elif sample_norm == "conditionmedian":
            med_map = feature.get_median_map_to_condition()

        # Process samples
        all_peptides = []

        for samples, df in feature.iter_samples():
            df.dropna(subset=["pg_accessions"], inplace=True)

            for sample in samples:
                dataset_df = df[df["sample_accession"] == sample].copy()
                dataset_df = dataset_df[dataset_df["unique"] == 1]
                dataset_df = dataset_df[PARQUET_COLUMNS]

                dataset_df = reformat_quantms_feature_table_quant_labels(
                    dataset_df, label, choice
                )
                dataset_df = apply_initial_filtering(
                    dataset_df, self.config.filtering.min_aa, AGGREGATION_LEVEL_SAMPLE
                )

                # Filter by min unique peptides
                dataset_df = dataset_df.groupby(PROTEIN_NAME).filter(
                    lambda x: len(set(x[PEPTIDE_CANONICAL]))
                    >= self.config.filtering.min_unique_peptides
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

        # Combine all peptides
        combined_df = pd.concat(all_peptides, ignore_index=True)

        # Apply dataset-level normalization if selected
        if sample_norm_method == PeptideNormalizationMethod.Hierarchical:
            combined_df = NormalizationStage(self.config).apply_hierarchical(combined_df)
        elif sample_norm_method == PeptideNormalizationMethod.TMM:
            combined_df = NormalizationStage(self.config).apply_tmm(combined_df)

        return combined_df

    def load_for_directlfq(self) -> pd.DataFrame:
        """Load and filter data for DirectLFQ processing."""
        filter_builder = SQLFilterBuilder(
            remove_contaminants=self.config.filtering.remove_contaminants,
            min_peptide_length=self.config.filtering.min_aa,
            require_unique=True,
        )

        feature = Feature(self.config.input.parquet, filter_builder=filter_builder)

        if self.config.input.sdrf:
            feature.enrich_with_sdrf(self.config.input.sdrf)

        # Build query with filters
        where_clause, where_params = filter_builder.build_where_clause()
        query = "".join([
            "SELECT pg_accessions, sequence, sample_accession, intensity",
            " FROM parquet_db WHERE ", where_clause,
        ])

        df = feature.parquet_db.execute(query, where_params).df()

        # Parse protein accessions
        # Extract first element from pg_accessions list, then parse UniProt ID
        first_acc = df["pg_accessions"].str[0].fillna("")
        df["protein"] = np.where(
            first_acc.str.contains("|", regex=False),
            first_acc.str.split("|").str[1],
            first_acc,
        )

        # Filter by min unique peptides
        peptide_counts = df.groupby("protein")["sequence"].nunique()
        valid_proteins = peptide_counts[
            peptide_counts >= self.config.filtering.min_unique_peptides
        ].index
        df = df[df["protein"].isin(valid_proteins)]

        return df

    def convert_to_directlfq_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert to DirectLFQ expected format (wide, log2, MultiIndex)."""
        # Pivot to wide format
        wide = df.pivot_table(
            index=["protein", "sequence"],
            columns="sample_accession",
            values="intensity",
            aggfunc="sum",
        )

        # Replace 0 with NaN and log2 transform
        wide = wide.replace(0, np.nan)
        wide = np.log2(wide)

        # Set index names for DirectLFQ
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
            logger.info(f"Using explicit reference samples: {ref_samples}")

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

        # Convert to wide format for normalization
        wide = df.pivot_table(
            index=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            columns=SAMPLE_ID,
            values=NORM_INTENSITY,
            aggfunc="sum",
        )

        # Log2 transform for normalization
        wide = wide.replace(0, np.nan)
        wide_log2 = np.log2(wide)

        # Load selected proteins if specified
        selected_proteins = None
        if self.config.normalization.proteins_file:
            with open(self.config.normalization.proteins_file) as f:
                selected_proteins = [line.strip() for line in f if line.strip()]
            logger.info(f"Using {len(selected_proteins)} selected proteins for normalization")

        # Apply hierarchical normalization
        normalizer = HierarchicalSampleNormalizer(
            num_samples_quadratic=self.config.quantification.directlfq_num_samples_quadratic,
            selected_proteins=selected_proteins,
        )

        normalized_log2 = normalizer.fit_transform(wide_log2)

        # Convert back to linear scale
        normalized_wide = 2 ** normalized_log2

        # Convert back to long format
        normalized_long = normalized_wide.reset_index().melt(
            id_vars=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            var_name=SAMPLE_ID,
            value_name=NORM_INTENSITY,
        )

        # Remove NaN rows
        normalized_long = normalized_long.dropna(subset=[NORM_INTENSITY])

        logger.info(f"Hierarchical normalization complete: {len(normalized_long)} rows")

        return normalized_long

    def apply_tmm(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply TMM (Trimmed Mean of M-values) sample normalization."""
        from mokume.normalization.tmm import TMMNormalizer

        logger.info("Applying TMM sample normalization...")

        # Convert to wide format for normalization
        wide = df.pivot_table(
            index=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            columns=SAMPLE_ID,
            values=NORM_INTENSITY,
            aggfunc="sum",
        )

        # Apply TMM normalization (works on linear-scale intensities)
        normalizer = TMMNormalizer()
        normalized_wide = normalizer.fit_transform(wide)

        # Convert back to long format
        normalized_long = normalized_wide.reset_index().melt(
            id_vars=[PROTEIN_NAME, PEPTIDE_CANONICAL],
            var_name=SAMPLE_ID,
            value_name=NORM_INTENSITY,
        )

        # Remove NaN rows
        normalized_long = normalized_long.dropna(subset=[NORM_INTENSITY])

        logger.info(f"TMM normalization complete: {len(normalized_long)} rows")

        return normalized_long

    def apply_irs(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """Apply Internal Reference Scaling normalization for multi-plex TMT data.

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
            logger.info(f"Using explicit reference samples: {ref_samples}")

        # 3. Explicit column + values
        if ref_samples is None and self.config.irs.sdrf_column and self.config.irs.sdrf_values:
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

        logger.info(f"IRS normalization complete: {len(protein_df)} proteins")

        # Optionally remove reference samples from output
        if self.config.irs.remove_reference:
            protein_col = protein_df.columns[0]
            cols_to_keep = [protein_col] + [
                c for c in protein_df.columns
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
            logger.info(
                f"Removed {len(ref_samples)} reference samples from output"
            )

        return protein_df

    def apply_coverage_filter(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """Apply coverage filter using condition mapping from SDRF."""
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
            s: c for s, c in sample_to_condition.items()
            if s in available_samples
            and "powder" not in c.lower()
            and "pool" not in c.lower()
        }

        return apply_coverage_filter(
            protein_df, sample_to_condition, self.config.quantification.coverage_threshold
        )


class QuantificationStage:
    """Routes to the appropriate quantification method."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def quantify(self, peptide_df: pd.DataFrame) -> pd.DataFrame:
        """Apply protein quantification method."""
        import re
        from mokume.quantification import get_quantification_method, TopNQuantification

        quant_method = self.config.quantification.method.lower()

        if quant_method == "ibaq":
            return self._quantify_ibaq(peptide_df)
        elif quant_method in ("maxlfq", "sum", "all"):
            method = get_quantification_method(quant_method)
            result = method.quantify(
                peptide_df,
                protein_column=PROTEIN_NAME,
                peptide_column=PEPTIDE_CANONICAL,
                intensity_column=NORM_INTENSITY,
                sample_column=SAMPLE_ID,
            )
            return self.to_wide_format(result, quant_method)
        elif quant_method.startswith("top"):
            match = re.match(r"top(\d+)", quant_method)
            if match:
                n = int(match.group(1))
            else:
                n = 3
            method = TopNQuantification(n=n)
            result = method.quantify(
                peptide_df,
                protein_column=PROTEIN_NAME,
                peptide_column=PEPTIDE_CANONICAL,
                intensity_column=NORM_INTENSITY,
                sample_column=SAMPLE_ID,
            )
            return self.to_wide_format(result, quant_method)
        elif quant_method == "median":
            return self._quantify_median(peptide_df)
        else:
            raise ValueError(f"Unknown quantification method: {quant_method}")

    def to_wide_format(self, long_df: pd.DataFrame, method_name: str) -> pd.DataFrame:
        """Convert long format quantification results to wide format."""
        intensity_col = "Intensity"

        if intensity_col not in long_df.columns:
            # Fallback: use the last numeric column
            numeric_cols = long_df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 0:
                intensity_col = numeric_cols[-1]

        logger.debug(f"Using intensity column: {intensity_col}")

        # Pivot to wide format
        wide_df = long_df.pivot(
            index=PROTEIN_NAME,
            columns=SAMPLE_ID,
            values=intensity_col,
        )

        return wide_df.reset_index()

    def _quantify_ibaq(self, peptide_df: pd.DataFrame) -> pd.DataFrame:
        """Quantify using iBAQ method."""
        from mokume.quantification.ibaq import extract_fasta

        if not self.config.input.fasta_file:
            raise ValueError(
                "iBAQ quantification requires a FASTA file. "
                "Use --fasta to provide one."
            )

        proteins = peptide_df[PROTEIN_NAME].unique().tolist()

        logger.info(f"Computing iBAQ for {len(proteins)} proteins using FASTA...")

        unique_peptide_counts, mw_dict, found_proteins = extract_fasta(
            fasta=self.config.input.fasta_file,
            enzyme="Trypsin",
            proteins=proteins,
            min_aa=self.config.filtering.min_aa,
            max_aa=50,
            tpa=False,
        )

        logger.info(f"Found {len(found_proteins)} proteins in FASTA")

        peptide_df = peptide_df[peptide_df[PROTEIN_NAME].isin(found_proteins)]

        protein_intensities = (
            peptide_df.groupby([PROTEIN_NAME, SAMPLE_ID], observed=False)[NORM_INTENSITY]
            .sum()
            .reset_index()
        )

        num_peptides = protein_intensities[PROTEIN_NAME].map(unique_peptide_counts).fillna(1)
        protein_intensities["iBAQ"] = np.where(
            num_peptides > 0,
            protein_intensities[NORM_INTENSITY] / num_peptides,
            0,
        )

        result_wide = protein_intensities.pivot(
            index=PROTEIN_NAME,
            columns=SAMPLE_ID,
            values="iBAQ",
        )

        logger.info(f"iBAQ complete: {len(result_wide)} proteins")

        return result_wide.reset_index()

    def _quantify_median(self, peptide_df: pd.DataFrame) -> pd.DataFrame:
        """Quantify using median of peptides."""
        result = (
            peptide_df.groupby([PROTEIN_NAME, SAMPLE_ID], observed=False)[NORM_INTENSITY]
            .median()
            .reset_index()
        )

        result_wide = result.pivot(
            index=PROTEIN_NAME, columns=SAMPLE_ID, values=NORM_INTENSITY
        )

        return result_wide.reset_index()


class PostprocessingStage:
    """Batch correction, DE, plotting, reports."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def apply_batch_correction(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """Apply batch correction after protein quantification.

        Batch = technical variation to REMOVE (from runs/files)
        Covariates = biological signal to PRESERVE (from SDRF characteristics)
        """
        if not is_batch_correction_available():
            raise ImportError(
                "Batch correction requires inmoose package. "
                "Install with: pip install mokume[batch-correction]"
            )

        # Identify protein column and sample columns
        if PROTEIN_NAME in protein_df.columns:
            protein_col = PROTEIN_NAME
        elif "protein" in protein_df.columns:
            protein_col = "protein"
        else:
            protein_col = protein_df.columns[0]

        sample_cols = [c for c in protein_df.columns if c != protein_col]

        if len(sample_cols) < 2:
            logger.warning("Not enough samples for batch correction, skipping")
            return protein_df

        # Create intensity matrix for batch correction (features x samples)
        intensity_matrix = protein_df.set_index(protein_col)[sample_cols]

        # 1. Detect batches (technical variation to remove)
        try:
            batch_method = BatchDetectionMethod.from_str(self.config.batch.method)
        except ValueError:
            logger.warning(
                f"Unknown batch method '{self.config.batch.method}', "
                "using sample_prefix"
            )
            batch_method = BatchDetectionMethod.SAMPLE_PREFIX

        batch_indices = detect_batches(
            sample_ids=sample_cols,
            method=batch_method,
            batch_column_values=(
                self._get_batch_column_values(sample_cols)
                if self.config.batch.column else None
            ),
        )

        unique_batches = len(set(batch_indices))
        logger.info(f"Detected {unique_batches} batches for batch correction")

        if unique_batches < 2:
            logger.warning("Only 1 batch detected, skipping batch correction")
            return protein_df

        # Check minimum samples per batch
        from collections import Counter
        batch_counts = Counter(batch_indices)
        min_samples = min(batch_counts.values())
        if min_samples < 2:
            logger.warning(
                f"Some batches have fewer than 2 samples (min={min_samples}), "
                "skipping batch correction"
            )
            return protein_df

        # 2. Extract covariates from SDRF (biological signal to preserve)
        covariates = None
        if self.config.input.sdrf and self.config.batch.covariates:
            covariates = extract_covariates_from_sdrf(
                self.config.input.sdrf,
                sample_cols,
                self.config.batch.covariates,
            )
            if covariates:
                logger.info(
                    f"Extracted {len(self.config.batch.covariates)} "
                    f"covariates to preserve biological signal"
                )

        # 3. Apply ComBat batch correction
        has_nan = intensity_matrix.isna().any(axis=1)
        complete_matrix = intensity_matrix[~has_nan]
        incomplete_matrix = intensity_matrix[has_nan]

        if len(complete_matrix) == 0:
            logger.warning("No proteins with complete data for ComBat, skipping")
            return protein_df

        n_incomplete = len(incomplete_matrix)
        if n_incomplete > 0:
            logger.info(
                f"ComBat: {len(complete_matrix)} complete proteins, "
                f"{n_incomplete} with missing values (kept uncorrected)"
            )

        logger.info("Applying ComBat batch correction...")
        try:
            corrected_matrix = apply_batch_correction(
                df=complete_matrix,
                batch=batch_indices,
                covs=covariates,
                kwargs={
                    "par_prior": self.config.batch.parametric,
                    "mean_only": self.config.batch.mean_only,
                    "ref_batch": self.config.batch.ref_batch,
                },
            )

            # Recombine corrected + uncorrected proteins
            if n_incomplete > 0:
                corrected_matrix = pd.concat([corrected_matrix, incomplete_matrix])

            # Reconstruct DataFrame with protein column
            corrected_df = corrected_matrix.reset_index()
            corrected_df = corrected_df.rename(columns={"index": protein_col})

            logger.info(
                f"Batch correction complete: {len(corrected_df)} proteins, "
                f"{len(sample_cols)} samples"
            )
            return corrected_df

        except Exception as e:
            logger.error(f"Batch correction failed: {e}")
            logger.warning("Returning uncorrected protein intensities")
            return protein_df

    def run_differential_expression(
        self, protein_df: pd.DataFrame
    ) -> Optional[dict]:
        """Run differential expression analysis."""
        from mokume.analysis.differential_expression import DifferentialExpression
        from mokume.normalization.irs import detect_condition_from_sdrf

        if not self.config.input.sdrf:
            raise ValueError("Differential expression requires an SDRF file (--sdrf)")

        # Get condition mapping from SDRF
        sample_to_condition = detect_condition_from_sdrf(self.config.input.sdrf)

        # Parse contrasts
        contrasts = []
        if self.config.de.contrasts:
            for c in self.config.de.contrasts:
                # Prefer " vs " delimiter to support hyphenated condition names
                if " vs " in c:
                    parts = c.split(" vs ", 1)
                    contrasts.append((parts[0].strip(), parts[1].strip()))
                elif "-" in c:
                    parts = c.split("-", 1)
                    contrasts.append((parts[0].strip(), parts[1].strip()))
                else:
                    logger.warning(f"Invalid contrast format '{c}', expected 'A vs B' or 'A-B'")
        else:
            # Auto-detect: compare all conditions pairwise
            conditions = sorted(set(sample_to_condition.values()))
            exp_conditions = [
                c for c in conditions
                if "powder" not in c.lower() and "pool" not in c.lower()
            ]
            if len(exp_conditions) == 2:
                contrasts = [(exp_conditions[0], exp_conditions[1])]
            else:
                logger.warning(
                    f"Cannot auto-detect contrasts from {exp_conditions}. "
                    "Use --de-contrasts to specify."
                )
                return None

        if not contrasts:
            return None

        # Auto-select DE method based on quantification if "auto"
        de_method = self.config.de.method
        if de_method == "auto":
            quant = self.config.quantification.method.lower()
            de_method = "deqms" if quant == "directlfq" else "limrots"
            logger.info(f"Auto-selected DE method: {de_method} (quant={quant})")

        # Load peptide counts for DEqMS
        peptide_counts = None
        if de_method == "deqms" and self.config.input.parquet:
            try:
                pep_df = pd.read_parquet(
                    self.config.input.parquet,
                    columns=["anchor_protein", "sequence"],
                )
                peptide_counts = pep_df.groupby("anchor_protein")["sequence"].nunique()
            except (FileNotFoundError, KeyError, ValueError) as exc:
                logger.warning("Could not load peptide counts for DEqMS: %s", exc)

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
            logger.info(f"Running DE: {key}")
            result = de.run(protein_df, sample_to_condition, contrast)
            all_results[key] = result

            # Save DE results
            if self.config.de.output:
                if len(contrasts) == 1:
                    output_path = self.config.de.output
                else:
                    base, ext = os.path.splitext(self.config.de.output)
                    output_path = f"{base}_{key}{ext}" if ext else f"{base}_{key}.csv"
                result.to_csv(output_path, index=False)
                logger.info(f"DE results saved to {output_path}")

        return all_results

    def generate_plots(
        self, protein_df: pd.DataFrame, de_results: Optional[dict] = None
    ):
        """Generate requested plots."""
        from mokume.normalization.irs import detect_condition_from_sdrf
        from mokume.plotting import is_plotting_available

        if not is_plotting_available():
            logger.warning(
                "Plotting dependencies not available. "
                "Install with: pip install mokume[plotting]"
            )
            return

        from mokume.plotting.differential_expression import (
            plot_volcano,
            plot_heatmap,
            plot_pca_conditions,
        )

        plot_dir = Path(self.config.output.plot_dir)
        plot_dir.mkdir(parents=True, exist_ok=True)

        # Get condition mapping
        sample_to_condition = None
        if self.config.input.sdrf:
            sample_to_condition = detect_condition_from_sdrf(self.config.input.sdrf)
            if self.config.irs.remove_reference:
                protein_col = protein_df.columns[0]
                available_samples = [
                    c for c in protein_df.columns if c != protein_col
                ]
                sample_to_condition = {
                    s: c for s, c in sample_to_condition.items()
                    if s in available_samples
                }

        # Volcano plot
        if self.config.output.plot_volcano and de_results:
            for contrast_name, de_df in de_results.items():
                output_file = str(plot_dir / f"volcano_{contrast_name}.png")
                plot_volcano(
                    de_df,
                    log2fc_threshold=self.config.de.log2fc_threshold,
                    fdr_threshold=self.config.de.fdr_threshold,
                    highlight_genes=self.config.output.highlight_genes,
                    title=f"Volcano Plot: {contrast_name}",
                    output_file=output_file,
                )
                logger.info(f"Volcano plot saved to {output_file}")

        # Heatmap
        if self.config.output.plot_heatmap and sample_to_condition:
            output_file = str(plot_dir / "heatmap.png")
            plot_heatmap(
                protein_df,
                sample_to_condition,
                proteins=self.config.output.highlight_genes,
                title="Protein Heatmap",
                output_file=output_file,
            )
            logger.info(f"Heatmap saved to {output_file}")

        # PCA by condition
        if self.config.output.plot_pca and sample_to_condition:
            output_file = str(plot_dir / "pca_conditions.png")
            plot_pca_conditions(
                protein_df,
                sample_to_condition,
                title="PCA by Condition",
                output_file=output_file,
            )
            logger.info(f"PCA plot saved to {output_file}")

    def generate_interactive_report(
        self, protein_df: pd.DataFrame, de_results: dict
    ):
        """Generate interactive HTML report for DE results."""
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

        for contrast_name, de_df in de_results.items():
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
            logger.info(f"Interactive report saved to {output_html}")

    def _get_batch_column_values(self, sample_ids: list) -> Optional[list]:
        """Get batch values from SDRF for explicit batch column."""
        if not self.config.input.sdrf or not self.config.batch.column:
            return None

        try:
            from mokume.core.constants import load_sdrf as _load_sdrf

            sdrf = _load_sdrf(self.config.input.sdrf)

            batch_col = self.config.batch.column.lower()
            if batch_col not in sdrf.columns:
                logger.warning(f"Batch column '{self.config.batch.column}' not in SDRF")
                return None

            # Map sample IDs to batch values
            sample_to_batch = dict(zip(sdrf["source name"], sdrf[batch_col]))
            return [sample_to_batch.get(s, "unknown") for s in sample_ids]

        except Exception as e:
            logger.warning(f"Failed to extract batch column: {e}")
            return None
