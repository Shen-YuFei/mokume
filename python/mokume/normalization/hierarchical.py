"""
Hierarchical clustering-based sample normalization.

This module provides a native mokume implementation of DirectLFQ-style
normalization algorithms, enabling the use of hierarchical sample alignment
with any quantification method (iBAQ, Top3, MaxLFQ, etc.).

The implementation is inspired by DirectLFQ's approach but does not depend
on the directlfq package.

References
----------
- Thielert et al. (2024). directLFQ: A protein intensity estimation
  algorithm for DDA, DIA and targeted proteomics data.
  https://github.com/MannLabs/directlfq
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.optimize import least_squares
from typing import Optional, Union
from enum import Enum, auto

from mokume.core.logger import get_logger

logger = get_logger("mokume.normalization.hierarchical")


class DistanceMetric(Enum):
    """
    Distance metrics for computing sample similarity.

    Attributes
    ----------
    MEDIAN : auto
        Absolute difference of medians between overlapping values.
        Most robust to outliers.
    VARIANCE : auto
        Absolute difference of variances between overlapping values.
    OVERLAP : auto
        Inverse of the overlap coefficient between sample ranges.
    """

    MEDIAN = auto()
    VARIANCE = auto()
    OVERLAP = auto()

    @classmethod
    def from_str(cls, name: str) -> "DistanceMetric":
        """Get metric from string name."""
        if name is None:
            return cls.MEDIAN
        name_lower = name.lower()
        for member in cls:
            if member.name.lower() == name_lower:
                return member
        raise ValueError(
            f"Unknown distance metric: {name}. Options: median, variance, overlap"
        )


class HierarchicalSampleNormalizer:
    """
    Hierarchical clustering-based sample normalization.

    Aligns sample distributions using:
    1. Distance matrix computation between all sample pairs
    2. Hierarchical clustering to find optimal alignment order
    3. Iterative shifting of samples to align distributions

    This is mokume's native implementation, inspired by DirectLFQ's approach.
    Use this with any quantification method (iBAQ, Top3, MaxLFQ, etc.).

    Parameters
    ----------
    num_samples_quadratic : int, default=50
        Use quadratic optimization for datasets with fewer samples.
        Above this threshold, use faster linear optimization.
    selected_proteins : list[str], optional
        If provided, compute normalization factors using only these proteins.
        Useful for normalizing on housekeeping genes or stable proteins.
    distance_metric : DistanceMetric or str, default=DistanceMetric.MEDIAN
        Metric for computing sample distances.
        Options: 'median', 'variance', 'overlap'
    min_overlap : int, default=10
        Minimum number of overlapping non-NaN values required between
        two samples to compute their distance. If fewer, distance is set to inf.

    Attributes
    ----------
    normalization_factors_ : dict[str, float]
        After fitting, contains the shift factor for each sample column.
        Positive values mean the sample was shifted up.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> # Create sample intensity data (log2 scale)
    >>> data = pd.DataFrame({
    ...     'sample1': [10.0, 11.0, 12.0, np.nan],
    ...     'sample2': [10.5, 11.5, np.nan, 13.0],
    ...     'sample3': [9.5, 10.5, 11.5, 12.5],
    ... }, index=['pep1', 'pep2', 'pep3', 'pep4'])
    >>>
    >>> normalizer = HierarchicalSampleNormalizer()
    >>> normalized = normalizer.fit_transform(data)
    >>> print(normalizer.normalization_factors_)
    {'sample1': 0.0, 'sample2': 0.25, 'sample3': -0.5}
    """

    def __init__(
        self,
        num_samples_quadratic: int = 50,
        selected_proteins: Optional[list[str]] = None,
        distance_metric: Union[DistanceMetric, str] = DistanceMetric.MEDIAN,
        min_overlap: int = 10,
    ):
        self.num_samples_quadratic = num_samples_quadratic
        self.selected_proteins = selected_proteins
        self.min_overlap = min_overlap

        if isinstance(distance_metric, str):
            self.distance_metric = DistanceMetric.from_str(distance_metric)
        else:
            self.distance_metric = distance_metric

        self.normalization_factors_: Optional[dict[str, float]] = None
        self._leaf_order: Optional[np.ndarray] = None

    def fit(self, intensity_df: pd.DataFrame) -> "HierarchicalSampleNormalizer":
        """
        Compute normalization factors from intensity data.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            DataFrame with features as rows and samples as columns.
            Can have a MultiIndex [protein, peptide/ion] or single index.
            Values should be log2-transformed intensities.
            NaN values are handled gracefully.

        Returns
        -------
        HierarchicalSampleNormalizer
            Returns self for method chaining.
        """
        df = intensity_df.copy()

        # Filter to selected proteins if specified
        if self.selected_proteins is not None:
            if isinstance(df.index, pd.MultiIndex):
                protein_level = df.index.get_level_values(0)
                mask = protein_level.isin(self.selected_proteins)
                df = df[mask]
                logger.info(
                    "Using %d features from %d selected proteins for normalization",
                    len(df),
                    len(self.selected_proteins),
                )
            else:
                logger.warning(
                    "selected_proteins specified but DataFrame doesn't have MultiIndex. "
                    "Ignoring selected_proteins filter."
                )

        if len(df) == 0:
            raise ValueError(
                "No features remaining after filtering. Check selected_proteins."
            )

        n_samples = len(df.columns)
        logger.info("Computing normalization factors for %d samples", n_samples)

        # Handle edge cases
        if n_samples == 1:
            self.normalization_factors_ = {df.columns[0]: 0.0}
            logger.info("Single sample - no normalization needed.")
            return self

        if n_samples == 2:
            # Simple case: shift second sample to match first
            col1, col2 = df.columns[0], df.columns[1]
            mask = ~(df[col1].isna() | df[col2].isna())
            if mask.sum() >= self.min_overlap:
                shift = np.median(df.loc[mask, col1]) - np.median(df.loc[mask, col2])
            else:
                shift = 0.0
            self.normalization_factors_ = {col1: 0.0, col2: shift}
            self._leaf_order = np.array([0, 1])
            logger.info(
                "Normalization factors computed. Range: [%.3f, %.3f]",
                min(self.normalization_factors_.values()),
                max(self.normalization_factors_.values()),
            )
            return self

        # Compute distance matrix
        dist_matrix = self._compute_distance_matrix(df)

        # Handle case where all distances are inf (no overlap)
        if np.all(np.isinf(dist_matrix)):
            logger.warning("No overlapping values between samples. Using zero shifts.")
            self.normalization_factors_ = {col: 0.0 for col in df.columns}
            return self

        # Replace inf with large value for clustering
        max_finite = np.nanmax(dist_matrix[np.isfinite(dist_matrix)])
        dist_matrix_for_clustering = np.where(
            np.isinf(dist_matrix), max_finite * 10, dist_matrix
        )

        # Hierarchical clustering
        condensed = squareform(dist_matrix_for_clustering)
        linkage_matrix = linkage(condensed, method="average")
        self._leaf_order = leaves_list(linkage_matrix)

        # Compute shift factors
        self.normalization_factors_ = self._compute_shifts(df, self._leaf_order)

        logger.info(
            "Normalization factors computed. Range: [%.3f, %.3f]",
            min(self.normalization_factors_.values()),
            max(self.normalization_factors_.values()),
        )

        return self

    def transform(self, intensity_df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply normalization factors to intensity data.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            DataFrame with features as rows and samples as columns.
            Must have the same column names as the data used for fitting.

        Returns
        -------
        pd.DataFrame
            Normalized intensity matrix with same structure as input.

        Raises
        ------
        ValueError
            If fit() has not been called first.
        """
        if self.normalization_factors_ is None:
            raise ValueError("Must call fit() before transform()")

        normalized = intensity_df.copy()

        for col in normalized.columns:
            if col in self.normalization_factors_:
                normalized[col] = normalized[col] + self.normalization_factors_[col]
            else:
                logger.warning(
                    "Column '%s' not found in normalization factors. Leaving unchanged.",
                    col,
                )

        return normalized

    def fit_transform(self, intensity_df: pd.DataFrame) -> pd.DataFrame:
        """
        Fit and transform in one step.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            DataFrame with features as rows and samples as columns.

        Returns
        -------
        pd.DataFrame
            Normalized intensity matrix.
        """
        return self.fit(intensity_df).transform(intensity_df)

    def _compute_distance_matrix(self, df: pd.DataFrame) -> np.ndarray:
        """Compute pairwise distances between all samples."""
        n_samples = len(df.columns)
        dist_matrix = np.zeros((n_samples, n_samples))

        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                dist = self._sample_distance(df.iloc[:, i].values, df.iloc[:, j].values)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

        return dist_matrix

    def _sample_distance(self, s1: np.ndarray, s2: np.ndarray) -> float:
        """
        Compute distance between two samples based on overlapping values.

        Parameters
        ----------
        s1, s2 : np.ndarray
            Sample intensity vectors (may contain NaN).

        Returns
        -------
        float
            Distance between samples. Returns inf if insufficient overlap.
        """
        # Get overlapping (non-NaN in both) values
        mask = ~(np.isnan(s1) | np.isnan(s2))
        n_overlap = np.sum(mask)

        if n_overlap < self.min_overlap:
            return np.inf

        s1_overlap = s1[mask]
        s2_overlap = s2[mask]

        if self.distance_metric == DistanceMetric.MEDIAN:
            return abs(np.median(s1_overlap) - np.median(s2_overlap))

        if self.distance_metric == DistanceMetric.VARIANCE:
            return abs(np.var(s1_overlap) - np.var(s2_overlap))

        # Overlap coefficient based on shared intensity range
        min_max = min(np.max(s1_overlap), np.max(s2_overlap))
        max_min = max(np.min(s1_overlap), np.min(s2_overlap))
        overlap = max(0, min_max - max_min)
        total_range = max(np.max(s1_overlap), np.max(s2_overlap)) - min(
            np.min(s1_overlap), np.min(s2_overlap)
        )
        if total_range == 0:
            return 0.0
        return 1.0 - (overlap / total_range)

    def _compute_shifts(
        self, df: pd.DataFrame, leaf_order: np.ndarray
    ) -> dict[str, float]:
        """Compute shift factors for each sample."""
        n_samples = len(df.columns)

        if n_samples <= self.num_samples_quadratic:
            logger.debug(
                "Using quadratic optimization (%d <= %d)",
                n_samples,
                self.num_samples_quadratic,
            )
            return self._compute_shifts_quadratic(df, leaf_order)

        logger.debug(
            "Using linear optimization (%d > %d)",
            n_samples,
            self.num_samples_quadratic,
        )
        return self._compute_shifts_linear(df, leaf_order)

    def _compute_shifts_linear(
        self, df: pd.DataFrame, leaf_order: np.ndarray
    ) -> dict[str, float]:
        """
        Linear optimization: shift samples along clustering order.

        This is O(n) in the number of samples, suitable for large datasets.
        """
        shifts = {df.columns[leaf_order[0]]: 0.0}
        cumulative_shift = 0.0

        for i in range(1, len(leaf_order)):
            prev_col = df.columns[leaf_order[i - 1]]
            curr_col = df.columns[leaf_order[i]]

            # Compute shift to align current to previous
            mask = ~(df[prev_col].isna() | df[curr_col].isna())
            if mask.sum() >= self.min_overlap:
                shift = np.median(df.loc[mask, prev_col]) - np.median(
                    df.loc[mask, curr_col]
                )
            else:
                shift = 0.0
                logger.debug(
                    "Insufficient overlap between %s and %s. Using zero shift.",
                    prev_col,
                    curr_col,
                )

            cumulative_shift += shift
            shifts[curr_col] = cumulative_shift

        return shifts

    def _compute_shifts_quadratic(
        self, df: pd.DataFrame, leaf_order: np.ndarray
    ) -> dict[str, float]:
        """
        Quadratic optimization: consider all pairwise comparisons.

        Uses least squares to find optimal shifts that minimize the sum of
        squared median differences between all pairs. This is O(n^2) but
        gives better results for small datasets.
        """
        n = len(df.columns)
        cols = [df.columns[i] for i in leaf_order]

        # Precompute pairwise median differences
        median_diffs = {}
        overlap_counts = {}

        for i in range(n):
            for j in range(i + 1, n):
                mask = ~(df[cols[i]].isna() | df[cols[j]].isna())
                overlap_counts[(i, j)] = mask.sum()
                if mask.sum() >= self.min_overlap:
                    median_diffs[(i, j)] = np.median(df.loc[mask, cols[i]]) - np.median(
                        df.loc[mask, cols[j]]
                    )
                else:
                    median_diffs[(i, j)] = 0.0

        def residuals(x):
            """Residuals: differences between (s_i + shift_i) and (s_j + shift_j)."""
            shifts = np.concatenate([[0.0], x])  # First sample fixed at 0
            res = []
            for i in range(n):
                for j in range(i + 1, n):
                    if overlap_counts[(i, j)] >= self.min_overlap:
                        # Expected difference after shifting should be zero
                        expected_diff = (shifts[i] - shifts[j]) - median_diffs[(i, j)]
                        # Weight by overlap count (more overlap = more confidence)
                        weight = np.sqrt(overlap_counts[(i, j)])
                        res.append(expected_diff * weight)
            return np.array(res) if res else np.array([0.0])

        # Initial guess from linear method
        linear_shifts = self._compute_shifts_linear(df, leaf_order)
        x0 = np.array(
            [linear_shifts[c] for c in cols[1:]]
        )  # Exclude first (fixed at 0)

        # Optimize
        try:
            result = least_squares(residuals, x0, method="lm", max_nfev=1000)
            optimal_shifts = np.concatenate([[0.0], result.x])

            if not result.success:
                logger.warning(
                    "Quadratic optimization did not converge: %s. "
                    "Using linear shifts instead.",
                    result.message,
                )
                return linear_shifts

        except Exception as e:
            logger.warning(
                "Quadratic optimization failed: %s. Using linear shifts instead.",
                e,
            )
            return linear_shifts

        return {cols[i]: optimal_shifts[i] for i in range(n)}


class HierarchicalIonAligner:
    """
    Align ions within each protein group using hierarchical approach.

    This aligns peptide/ion intensity traces within a single protein,
    useful for LFQ-style quantification methods that need aligned ions
    before computing protein intensity.

    Parameters
    ----------
    num_samples_quadratic : int, default=10
        Use quadratic optimization for proteins with fewer ions.
        Since proteins typically have few peptides, default is lower.
    min_overlap : int, default=3
        Minimum overlapping samples for ion alignment.

    Examples
    --------
    >>> # Align ions for a single protein
    >>> protein_df = pd.DataFrame({
    ...     'sample1': [10.0, 10.5, 11.0],
    ...     'sample2': [10.2, 10.7, np.nan],
    ...     'sample3': [9.8, 10.3, 10.8],
    ... }, index=['ion1', 'ion2', 'ion3'])
    >>>
    >>> aligner = HierarchicalIonAligner()
    >>> aligned = aligner.align_protein_ions(protein_df)
    """

    def __init__(
        self,
        num_samples_quadratic: int = 10,
        min_overlap: int = 3,
    ):
        self.num_samples_quadratic = num_samples_quadratic
        self.min_overlap = min_overlap

    def align_protein_ions(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """
        Align ion intensities within a single protein group.

        Parameters
        ----------
        protein_df : pd.DataFrame
            DataFrame with ions as rows and samples as columns.
            Values should be log2-transformed intensities.

        Returns
        -------
        pd.DataFrame
            Aligned ion intensities with same structure.
        """
        if len(protein_df) <= 1:
            return protein_df

        # Transpose: treat ions as "samples" to align their traces
        # Each row (ion) becomes a column, each column (sample) becomes a row
        transposed = protein_df.T

        normalizer = HierarchicalSampleNormalizer(
            num_samples_quadratic=self.num_samples_quadratic,
            min_overlap=self.min_overlap,
        )

        try:
            aligned_transposed = normalizer.fit_transform(transposed)
            return aligned_transposed.T
        except Exception as e:
            logger.debug("Ion alignment failed for protein: %s. Returning original.", e)
            return protein_df

    def align_all_proteins(
        self,
        intensity_df: pd.DataFrame,
        protein_column: str = None,
    ) -> pd.DataFrame:
        """
        Align ions for all proteins in a dataset.

        Parameters
        ----------
        intensity_df : pd.DataFrame
            DataFrame with MultiIndex [protein, ion] and samples as columns,
            OR a long-format DataFrame with protein and ion columns.
        protein_column : str, optional
            If long format, the column name for protein IDs.

        Returns
        -------
        pd.DataFrame
            Aligned intensities with same structure as input.
        """
        _ = protein_column
        if not isinstance(intensity_df.index, pd.MultiIndex):
            raise NotImplementedError(
                "Long format alignment not yet implemented. "
                "Please provide wide format with MultiIndex [protein, ion]."
            )

        result_dfs = []
        proteins = intensity_df.index.get_level_values(0).unique()

        for protein in proteins:
            protein_data = intensity_df.loc[protein]
            aligned = self.align_protein_ions(protein_data)
            aligned.index = pd.MultiIndex.from_tuples(
                [(protein, ion) for ion in aligned.index],
                names=intensity_df.index.names,
            )
            result_dfs.append(aligned)

        return pd.concat(result_dfs)
