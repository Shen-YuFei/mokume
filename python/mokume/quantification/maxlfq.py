"""
MaxLFQ protein quantification method.

This module provides the MaxLFQ algorithm for label-free quantification.

**Implementation Strategy:**
By default, this module uses DirectLFQ (if installed) for maximum accuracy.
If DirectLFQ is not available, it falls back to a built-in implementation
that uses peptide trace alignment (inspired by DirectLFQ).

To install DirectLFQ for best results:
    pip install mokume[directlfq]

The built-in fallback implementation:
- Aligns peptide intensity traces within each protein using median shifts
- Aggregates aligned traces using median
- Scales results to preserve total peptide intensity
- Uses parallelization via joblib for performance

References:
    Cox J, et al. Accurate Proteome-wide Label-free Quantification by Delayed
    Normalization and Maximal Peptide Ratio Extraction, Termed MaxLFQ.
    Mol Cell Proteomics. 2014;13(9):2513-26.

    Ammar C, et al. Accurate label-free quantification by directLFQ to compare
    unlimited numbers of proteomes. Mol Cell Proteomics. 2023.
"""

import importlib
import os
import warnings
from typing import Optional

import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from mokume.quantification.base import ProteinQuantificationMethod
from mokume.core.logger import get_logger
from mokume.core.registry import PluginRegistry
from mokume.core.constants import (
    PROTEIN_NAME,
    PEPTIDE_CANONICAL,
    NORM_INTENSITY,
    SAMPLE_ID,
)

logger = get_logger("mokume.quantification.maxlfq")


def _is_directlfq_available() -> bool:
    """Check if DirectLFQ package is installed and importable."""
    try:
        importlib.import_module("directlfq")
        return True
    except (ImportError, ModuleNotFoundError):
        return False


def _resolve_directlfq_num_cores(threads: int) -> Optional[int]:
    """Translate the joblib-style ``threads`` sentinel into a DirectLFQ core count.

    DirectLFQ's ``num_cores`` expects an actual positive count (or ``None`` to
    fall back to its own default); it does not understand joblib's negative
    "all cores" sentinels. The builtin MaxLFQ path passes ``threads`` straight
    to ``joblib.Parallel(n_jobs=...)``, where ``-1`` means all cores and ``-k``
    means all-but-(k-1); mirror that here so ``threads`` behaves the same on
    both paths instead of silently collapsing ``-1`` to sequential.
    """
    if threads > 0:
        return threads
    if threads < 0:
        return max(1, (os.cpu_count() or 1) + 1 + threads)
    return None


def _select_reference_peptide(
    log_matrix: np.ndarray,
    valid_counts: np.ndarray,
    peptide_ids: Optional[np.ndarray] = None,
) -> int:
    """Pick the reference peptide deterministically, independent of row order.

    The reference anchors the median-shift alignment, so choosing a different one
    changes every aligned trace and therefore the protein quantity. ``np.argmax``
    breaks ties by row position, and ties are the common case (e.g. every peptide
    observed in every sample), so the result depended on the order rows happened to
    arrive in -- which is not guaranteed by a parallel/unordered upstream read. Two
    identical runs could then differ by more than 2 log2 units on a protein.

    With peptide identifiers, sorting each trace before summing and using the
    identifier for final ties makes selection independent of both matrix axes.
    Without identifiers, the value-based fallback remains row-order independent.
    """
    candidates = np.flatnonzero(valid_counts == valid_counts.max())
    if candidates.size == 1:
        return int(candidates[0])

    # Prefer the most intense trace among those with the most measurements.
    ordered_values = np.sort(log_matrix[candidates, :], axis=1)
    totals = np.nansum(ordered_values, axis=1)
    candidates = candidates[np.flatnonzero(totals == totals.max())]
    if candidates.size == 1:
        return int(candidates[0])

    if peptide_ids is not None:
        return int(min(candidates, key=lambda idx: str(peptide_ids[idx])))

    # Final deterministic tiebreak: lexicographically smallest trace. NaNs are
    # mapped to +inf so they sort last and never compare equal to a real value.
    keys = np.where(
        np.isnan(log_matrix[candidates, :]), np.inf, log_matrix[candidates, :]
    )
    order = sorted(range(len(candidates)), key=lambda i: tuple(keys[i]))
    return int(candidates[order[0]])


def _maxlfq_solve_protein(
    peptide_matrix: np.ndarray,
    peptide_ids: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Solve the MaxLFQ optimization problem for a single protein (built-in fallback).

    Uses peptide trace alignment inspired by DirectLFQ for optimal accuracy.
    The algorithm:
    1. Aligns peptide intensity traces using median shifts
    2. Takes median of aligned traces per sample
    3. Scales to preserve total peptide intensity

    Parameters
    ----------
    peptide_matrix : np.ndarray
        Matrix of shape (n_peptides, n_samples) with peptide intensities.
        NaN values indicate missing measurements.
    peptide_ids : np.ndarray, optional
        Stable peptide identifiers aligned with the matrix rows.

    Returns
    -------
    np.ndarray
        Array of protein intensities for each sample.
    """
    n_peptides, n_samples = peptide_matrix.shape

    if n_samples == 0:
        return np.array([])

    if n_samples == 1:
        valid_values = peptide_matrix[~np.isnan(peptide_matrix)]
        if len(valid_values) == 0:
            return np.array([np.nan])
        return np.array([np.median(valid_values)])

    if n_peptides == 1:
        # Single peptide: return its intensities directly
        return peptide_matrix[0, :].copy()

    # Store original sum for scaling
    original_sum = np.nansum(peptide_matrix)
    if original_sum <= 0:
        return np.full(n_samples, np.nan)

    # Log-transform for ratio calculations
    with np.errstate(divide="ignore", invalid="ignore"):
        log_matrix = np.log2(peptide_matrix)

    # Step 1: Align peptide traces
    # Use peptide with most valid values as reference
    valid_counts = np.sum(~np.isnan(log_matrix), axis=1)
    if valid_counts.max() == 0:
        return np.full(n_samples, np.nan)

    ref_peptide_idx = _select_reference_peptide(log_matrix, valid_counts, peptide_ids)
    ref_trace = log_matrix[ref_peptide_idx, :]

    # Align other peptides to reference using median shift
    aligned_matrix = log_matrix.copy()
    for pep_idx in range(n_peptides):
        if pep_idx == ref_peptide_idx:
            continue

        pep_trace = log_matrix[pep_idx, :]

        # Find samples measured in both reference and current peptide
        valid = ~np.isnan(ref_trace) & ~np.isnan(pep_trace)
        if np.sum(valid) > 0:
            # Compute median shift to align this peptide to reference
            shift = np.nanmedian(ref_trace[valid] - pep_trace[valid])
            aligned_matrix[pep_idx, :] = pep_trace + shift

    # Step 2: Take median of aligned traces per sample
    # Suppress warning for samples with no peptides (all-NaN columns)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        log_intensities = np.nanmedian(aligned_matrix, axis=0)

    # Step 3: Scale to preserve total peptide intensity
    intensities = np.power(2, log_intensities)

    # Handle NaN samples using fallback
    for i in range(n_samples):
        if np.isnan(intensities[i]):
            sample_peptides = peptide_matrix[:, i]
            valid_peptides = sample_peptides[~np.isnan(sample_peptides)]
            if len(valid_peptides) > 0:
                intensities[i] = np.median(valid_peptides)

    # Scale to preserve total intensity sum
    current_sum = np.nansum(intensities)
    if current_sum > 0:
        scale_factor = original_sum / current_sum
        intensities = intensities * scale_factor

    return intensities


def _process_protein(
    protein: str,
    protein_data: pd.DataFrame,
    peptide_column: str,
    intensity_column: str,
    sample_column: str,
    samples: np.ndarray,
    min_peptides: int,
    run_column: Optional[str] = None,
) -> list:
    """
    Process a single protein for MaxLFQ quantification (built-in fallback).

    Parameters
    ----------
    protein : str
        Protein identifier.
    protein_data : pd.DataFrame
        Peptide data for this protein.
    peptide_column : str
        Column name for peptides.
    intensity_column : str
        Column name for intensities.
    sample_column : str
        Column name for samples.
    samples : np.ndarray
        Array of all sample names.
    min_peptides : int
        Minimum peptides required for MaxLFQ.

    Returns
    -------
    list
        List of result dictionaries for this protein.
    """
    _ = run_column
    results = []
    # Sorted, not order-of-appearance: ``Series.unique()`` preserves the order rows
    # happen to arrive in, which an unordered/parallel upstream read does not fix.
    # That order feeds the pivot below and therefore the matrix handed to MaxLFQ, so
    # leaving it unsorted makes the whole protein quantification input-order dependent.
    peptides = np.sort(protein_data[peptide_column].unique())

    if len(peptides) < min_peptides:
        # Fall back to median for proteins with few peptides
        for sample in samples:
            sample_data = protein_data[protein_data[sample_column] == sample]
            if len(sample_data) > 0:
                intensity = sample_data[intensity_column].median()
                results.append(
                    {
                        "protein": protein,
                        "sample": sample,
                        "intensity": intensity,
                    }
                )
        return results

    # Create peptide x sample matrix via pivot_table (vectorized, sums duplicates).
    # observed=True avoids materialising every Categorical level of sample_column
    # inside groupby; the reindex below re-introduces the full sample axis with
    # NaN where absent, so the downstream matrix shape is unchanged.
    pivot = protein_data.pivot_table(
        index=peptide_column,
        columns=sample_column,
        values=intensity_column,
        aggfunc="sum",
        observed=True,
    )
    # Reindex to ensure consistent ordering
    pivot = pivot.reindex(index=peptides, columns=samples)
    peptide_matrix = pivot.values

    # Run MaxLFQ algorithm
    intensities = _maxlfq_solve_protein(peptide_matrix, peptides)

    # Store results
    for i, sample in enumerate(samples):
        if not np.isnan(intensities[i]) and intensities[i] > 0:
            results.append(
                {
                    "protein": protein,
                    "sample": sample,
                    "intensity": intensities[i],
                }
            )

    return results


@PluginRegistry.register("quantification", "maxlfq")
class MaxLFQQuantification(ProteinQuantificationMethod):
    """
    MaxLFQ protein quantification with automatic DirectLFQ integration.

    This class provides MaxLFQ-style label-free quantification. By default,
    it uses DirectLFQ (if installed) for maximum accuracy. If DirectLFQ is
    not available, it falls back to a built-in implementation.

    **Recommended:** Install DirectLFQ for best results:
        pip install mokume[directlfq]

    Parameters
    ----------
    min_peptides : int
        Minimum number of peptides required for MaxLFQ calculation.
        Proteins with fewer peptides will use median aggregation.
        Default is 2.
    threads : int
        Number of parallel threads. Use -1 for all available cores,
        1 for single-threaded execution. Default is -1.
    verbose : int
        Verbosity level for parallel processing (0=silent, 10=verbose).
        Default is 0.
    force_builtin : bool
        If True, always use the built-in implementation even if DirectLFQ
        is available. Useful for testing or comparison. Default is False.

    Attributes
    ----------
    using_directlfq : bool
        True if DirectLFQ is being used, False if using built-in fallback.

    Examples
    --------
    >>> from mokume.quantification import MaxLFQQuantification
    >>> maxlfq = MaxLFQQuantification(min_peptides=2, threads=4)
    >>> result = maxlfq.quantify(
    ...     peptide_df,
    ...     protein_column="ProteinName",
    ...     peptide_column="PeptideSequence",
    ...     intensity_column="Intensity",
    ...     sample_column="SampleID"
    ... )
    >>> # Check which implementation was used
    >>> print(f"Used DirectLFQ: {maxlfq.using_directlfq}")

    Notes
    -----
    DirectLFQ typically provides slightly better accuracy than the built-in
    implementation. If you need the most accurate results, install DirectLFQ:

        pip install mokume[directlfq]

    The built-in implementation uses peptide trace alignment (inspired by
    DirectLFQ) and achieves ~0.95 correlation with DIA-NN's MaxLFQ values.

    References
    ----------
    Cox J, et al. Accurate Proteome-wide Label-free Quantification by Delayed
    Normalization and Maximal Peptide Ratio Extraction, Termed MaxLFQ.
    Mol Cell Proteomics. 2014;13(9):2513-26.

    Ammar C, et al. Accurate label-free quantification by directLFQ to compare
    unlimited numbers of proteomes. Mol Cell Proteomics. 2023.
    """

    def __init__(
        self,
        min_peptides: int = 2,
        threads: int = -1,
        verbose: int = 0,
        force_builtin: bool = False,
    ):
        """
        Initialize MaxLFQ quantification.

        Parameters
        ----------
        min_peptides : int
            Minimum number of peptides required for MaxLFQ calculation.
        threads : int
            Number of parallel threads (-1 for all cores, 1 for single-threaded).
        verbose : int
            Verbosity level for parallel processing.
        force_builtin : bool
            If True, use built-in implementation even if DirectLFQ is available.
        """
        self.min_peptides = min_peptides
        self.force_builtin = force_builtin
        self.threads = threads
        self.verbose = verbose

        # Determine which implementation to use
        self._directlfq_available = _is_directlfq_available()
        self.using_directlfq = self._directlfq_available and not force_builtin

        if self.using_directlfq:
            logger.info("MaxLFQ: Using DirectLFQ for quantification")
        else:
            if force_builtin:
                logger.info("MaxLFQ: Using built-in implementation (forced)")
            else:
                logger.info(
                    "MaxLFQ: Using built-in implementation "
                    "(install 'directlfq' for better accuracy: pip install mokume[directlfq])"
                )

    @property
    def name(self) -> str:
        if self.using_directlfq:
            return "MaxLFQ (DirectLFQ)"
        return "MaxLFQ (built-in)"

    def _quantify_with_directlfq(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str,
        peptide_column: str,
        intensity_column: str,
        sample_column: str,
    ) -> pd.DataFrame:
        """Run quantification using DirectLFQ."""
        from mokume.quantification.directlfq import DirectLFQQuantification

        directlfq = DirectLFQQuantification(
            min_nonan=self.min_peptides,
            num_cores=_resolve_directlfq_num_cores(self.threads),
        )

        result_df = directlfq.quantify(
            peptide_df,
            protein_column=protein_column,
            peptide_column=peptide_column,
            intensity_column=intensity_column,
            sample_column=sample_column,
        )

        # DirectLFQ already outputs 'Intensity' column - no rename needed

        return result_df

    def _quantify_builtin(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str,
        peptide_column: str,
        intensity_column: str,
        sample_column: str,
    ) -> pd.DataFrame:
        """Run quantification using built-in implementation."""
        # Get unique samples and proteins
        samples = np.sort(peptide_df[sample_column].unique())
        proteins = peptide_df[protein_column].unique()

        logger.info(
            "Processing %d proteins across %d samples",
            len(proteins),
            len(samples),
        )
        logger.info("Threads: %s", self.threads)

        # Group data by protein for efficient access
        grouped = peptide_df.groupby(protein_column)

        # Process proteins in parallel
        all_results = Parallel(n_jobs=self.threads, verbose=self.verbose)(
            delayed(_process_protein)(
                protein,
                group,
                peptide_column,
                intensity_column,
                sample_column,
                samples,
                self.min_peptides,
            )
            for protein, group in grouped
        )

        # Flatten results
        results = []
        for protein_results in all_results:
            results.extend(protein_results)

        # Create result DataFrame
        result_df = pd.DataFrame(results)

        if len(result_df) > 0:
            result_df = result_df.rename(
                columns={
                    "protein": protein_column,
                    "sample": sample_column,
                    "intensity": "Intensity",
                }
            )

        return result_df

    def quantify(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str = PROTEIN_NAME,
        peptide_column: str = PEPTIDE_CANONICAL,
        intensity_column: str = NORM_INTENSITY,
        sample_column: str = SAMPLE_ID,
        run_column: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Quantify proteins using the MaxLFQ algorithm.

        Uses DirectLFQ if available, otherwise falls back to built-in
        implementation. Check `self.using_directlfq` to see which
        implementation is being used.

        Parameters
        ----------
        peptide_df : pd.DataFrame
            DataFrame containing peptide-level data.
        protein_column : str
            Column name for protein identifiers.
        peptide_column : str
            Column name for peptide sequences.
        intensity_column : str
            Column name for intensity values.
        sample_column : str
            Column name for sample identifiers.
        run_column : str, optional
            Column name for run identifiers. If provided, quantification
            is performed at the run level instead of sample level.
            Note: DirectLFQ delegation does not support run_column yet,
            so the built-in implementation will be used when run_column is provided.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: protein_column, sample_column,
            (run_column if provided), 'Intensity'.
        """
        logger.info("Running MaxLFQ quantification (%s)", self.name)

        # If run_column is provided, use built-in implementation
        # DirectLFQ delegation doesn't support run-level aggregation yet
        if run_column is not None and run_column in peptide_df.columns:
            logger.info("Using built-in implementation for run-level quantification")
            result_df = self._quantify_builtin_with_runs(
                peptide_df,
                protein_column,
                peptide_column,
                intensity_column,
                sample_column,
                run_column,
            )
        elif self.using_directlfq:
            result_df = self._quantify_with_directlfq(
                peptide_df,
                protein_column,
                peptide_column,
                intensity_column,
                sample_column,
            )
        else:
            result_df = self._quantify_builtin(
                peptide_df,
                protein_column,
                peptide_column,
                intensity_column,
                sample_column,
            )

        n_proteins = result_df[protein_column].nunique() if len(result_df) > 0 else 0
        n_samples = result_df[sample_column].nunique() if len(result_df) > 0 else 0
        logger.info("MaxLFQ complete: %d proteins, %d samples", n_proteins, n_samples)

        return result_df

    def _quantify_builtin_with_runs(
        self,
        peptide_df: pd.DataFrame,
        protein_column: str,
        peptide_column: str,
        intensity_column: str,
        sample_column: str,
        run_column: str,
    ) -> pd.DataFrame:
        """
        Run quantification at run level using built-in implementation.

        This processes each (sample, run) combination separately, similar
        to how DIA-NN performs MaxLFQ at the run level.
        """
        # Create a combined grouping column for sample+run
        # Use a separator unlikely to appear in sample/run names
        sep = "|||"
        peptide_df = peptide_df.copy()
        peptide_df["_sample_run"] = (
            peptide_df[sample_column].astype(str)
            + sep
            + peptide_df[run_column].astype(str)
        )

        # Get unique sample-run combinations
        sample_runs = np.sort(peptide_df["_sample_run"].unique())
        proteins = peptide_df[protein_column].unique()

        logger.info(
            "Processing %d proteins across %d sample-run combinations",
            len(proteins),
            len(sample_runs),
        )
        logger.info("Threads: %s", self.threads)

        # Group data by protein for efficient access
        grouped = peptide_df.groupby(protein_column)

        # Process proteins in parallel
        all_results = Parallel(n_jobs=self.threads, verbose=self.verbose)(
            delayed(_process_protein)(
                protein,
                group,
                peptide_column,
                intensity_column,
                "_sample_run",  # Use combined column for grouping
                sample_runs,
                self.min_peptides,
            )
            for protein, group in grouped
        )

        # Flatten results
        results = []
        for protein_results in all_results:
            results.extend(protein_results)

        # Create result DataFrame
        result_df = pd.DataFrame(results)

        if len(result_df) > 0:
            # Split sample_run back into sample and run using the same separator
            # Use regex=False to treat separator as literal string
            result_df[[sample_column, run_column]] = result_df["sample"].str.split(
                sep, n=1, expand=True, regex=False
            )
            result_df = result_df.drop(columns=["sample"])
            result_df = result_df.rename(
                columns={
                    "protein": protein_column,
                    "intensity": "Intensity",
                }
            )

        return result_df
