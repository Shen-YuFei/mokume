"""
TMM (Trimmed Mean of M-values) normalization.

This module provides TMM normalization for proteomics and RNA-seq data.
TMM is a robust normalization method that handles composition bias by
computing scaling factors based on a trimmed weighted mean of log fold-changes.

References
----------
- Robinson, M.D., Oshlack, A. (2010). A scaling normalization method for
  differential expression analysis of RNA-seq data. Genome Biology, 11(3), R25.
- RNAnorm Python implementation: https://github.com/genialis/RNAnorm
"""

import numpy as np
import pandas as pd
from typing import Optional

from mokume.core.logger import get_logger

logger = get_logger("mokume.normalization.tmm")


class TMMNormalizer:
    """
    Trimmed Mean of M-values (TMM) normalization.

    TMM is a robust normalization method that computes scaling factors
    based on a weighted trimmed mean of log fold-changes, filtering out
    genes/proteins with extreme expression or extreme fold-changes.

    This addresses composition bias where a small number of highly
    expressed features can skew simple normalization methods.

    Parameters
    ----------
    m_trim : float, default=0.3
        Proportion of M-values to trim from each tail.
        Higher values = more robust, but uses less data.
    a_trim : float, default=0.05
        Proportion of A-values to trim from each tail.
        Removes proteins with extreme average expression.
    ref_sample : str or None, default=None
        Reference sample for computing M-values.
        If None, automatically selects sample with upper quartile
        closest to mean upper quartile.
    log_transform : bool, default=False
        If True, input data is already log2-transformed.
        If False (default), expects raw intensities and will log-transform.

    Attributes
    ----------
    norm_factors_ : pd.Series
        Normalization factors per sample after fitting.
        Effective library size = raw_library_size * norm_factor.
    ref_sample_ : str
        The reference sample used for normalization.
    geometric_mean_ : float
        Geometric mean of raw factors (used for rescaling).

    References
    ----------
    Robinson, M.D., Oshlack, A. (2010). A scaling normalization method
    for differential expression analysis of RNA-seq data.
    Genome Biology, 11(3), R25.

    Examples
    --------
    >>> import pandas as pd
    >>> from mokume.normalization.tmm import TMMNormalizer
    >>>
    >>> # Protein intensities (proteins x samples)
    >>> data = pd.DataFrame({
    ...     'Sample1': [100, 200, 300, 10000],
    ...     'Sample2': [110, 190, 310, 12000],
    ...     'Sample3': [95, 205, 295, 8000],
    ... }, index=['ProtA', 'ProtB', 'ProtC', 'ProtD'])
    >>>
    >>> tmm = TMMNormalizer()
    >>> normalized = tmm.fit_transform(data)
    """

    def __init__(
        self,
        m_trim: float = 0.3,
        a_trim: float = 0.05,
        ref_sample: Optional[str] = None,
        log_transform: bool = False,
    ):
        self.m_trim = m_trim
        self.a_trim = a_trim
        self.ref_sample = ref_sample
        self.log_transform = log_transform

        # Fitted attributes
        self.norm_factors_: Optional[pd.Series] = None
        self.ref_sample_: Optional[str] = None
        self.geometric_mean_: Optional[float] = None

    def _select_reference_sample(self, X: pd.DataFrame) -> str:
        """
        Select reference sample as the one with upper quartile
        closest to mean upper quartile across samples.

        Parameters
        ----------
        X : pd.DataFrame
            Intensity matrix (proteins x samples).

        Returns
        -------
        str
            Sample name to use as reference.
        """
        # Compute upper quartile for each sample (ignoring zeros and NaN)
        upper_quartiles = {}
        for col in X.columns:
            vals = X[col].replace(0, np.nan).dropna()
            if len(vals) > 0:
                upper_quartiles[col] = np.percentile(vals, 75)
            else:
                upper_quartiles[col] = np.nan

        uq_series = pd.Series(upper_quartiles)
        mean_uq = uq_series.mean()

        # Select sample with UQ closest to mean UQ
        ref_sample = (uq_series - mean_uq).abs().idxmin()
        logger.debug(f"Selected reference sample: {ref_sample} (UQ={upper_quartiles[ref_sample]:.2f})")
        return ref_sample

    def _compute_tmm_factor(
        self,
        X: pd.DataFrame,
        sample: str,
        lib_sizes: pd.Series,
    ) -> float:
        """
        Compute TMM normalization factor for a single sample.

        Parameters
        ----------
        X : pd.DataFrame
            Intensity matrix (proteins x samples).
        sample : str
            Sample to compute factor for.
        lib_sizes : pd.Series
            Library sizes (sum of intensities) per sample.

        Returns
        -------
        float
            TMM factor for this sample (before geometric mean rescaling).
        """
        ref = self.ref_sample_

        if sample == ref:
            return 1.0

        # Get intensities for sample and reference
        y_sample = X[sample].values
        y_ref = X[ref].values

        # Get library sizes
        n_sample = lib_sizes[sample]
        n_ref = lib_sizes[ref]

        # Only use proteins with non-zero values in both samples
        valid = (y_sample > 0) & (y_ref > 0) & np.isfinite(y_sample) & np.isfinite(y_ref)

        if valid.sum() < 10:
            logger.warning(f"Sample {sample}: Only {valid.sum()} valid proteins for TMM. Returning factor=1.")
            return 1.0

        y_s = y_sample[valid]
        y_r = y_ref[valid]

        # Compute M-values: log2 fold-change
        # M = log2(Y_sample / N_sample) - log2(Y_ref / N_ref)
        #   = log2(Y_sample * N_ref / (Y_ref * N_sample))
        M = np.log2((y_s / n_sample) / (y_r / n_ref))

        # Compute A-values: average log expression
        # A = 0.5 * (log2(Y_sample / N_sample) + log2(Y_ref / N_ref))
        A = 0.5 * (np.log2(y_s / n_sample) + np.log2(y_r / n_ref))

        # Compute weights (inverse asymptotic variance)
        # w = (N_sample - Y_sample) / (N_sample * Y_sample) + (N_ref - Y_ref) / (N_ref * Y_ref)
        # This is the delta method approximation on binomial variance
        w = (n_sample - y_s) / (n_sample * y_s) + (n_ref - y_r) / (n_ref * y_r)

        # Avoid division by zero for weights
        w = np.where(w > 0, 1.0 / w, 0)

        # Double trimming
        n_proteins = len(M)
        m_lo = int(np.floor(n_proteins * self.m_trim))
        m_hi = int(np.ceil(n_proteins * (1 - self.m_trim)))
        a_lo = int(np.floor(n_proteins * self.a_trim))
        a_hi = int(np.ceil(n_proteins * (1 - self.a_trim)))

        # Get indices of proteins to keep after M-trimming
        m_order = np.argsort(M)
        keep_m = np.zeros(n_proteins, dtype=bool)
        keep_m[m_order[m_lo:m_hi]] = True

        # Get indices of proteins to keep after A-trimming
        a_order = np.argsort(A)
        keep_a = np.zeros(n_proteins, dtype=bool)
        keep_a[a_order[a_lo:a_hi]] = True

        # Keep proteins that pass both filters
        keep = keep_m & keep_a

        if keep.sum() < 5:
            logger.warning(f"Sample {sample}: Only {keep.sum()} proteins after trimming. Using less stringent trim.")
            keep = keep_m  # Fall back to just M-trimming

        if keep.sum() < 5:
            logger.warning(f"Sample {sample}: Still too few proteins. Returning factor=1.")
            return 1.0

        # Compute weighted mean of M-values for kept proteins
        M_keep = M[keep]
        w_keep = w[keep]

        if w_keep.sum() > 0:
            tmm_factor = np.sum(w_keep * M_keep) / np.sum(w_keep)
        else:
            tmm_factor = np.mean(M_keep)

        # Convert from log2 to factor
        factor = 2 ** tmm_factor

        return factor

    def fit(self, X: pd.DataFrame, y=None) -> "TMMNormalizer":
        """
        Compute TMM normalization factors.

        Parameters
        ----------
        X : pd.DataFrame
            Intensity matrix (proteins x samples).
            Should contain raw intensities (not log-transformed).
        y : None
            Ignored. Present for scikit-learn compatibility.

        Returns
        -------
        self
            Fitted normalizer.
        """
        if X.empty:
            raise ValueError("Input DataFrame is empty")

        # Handle log-transformed input
        if self.log_transform:
            X_raw = 2 ** X  # Convert back to raw for TMM computation
        else:
            X_raw = X

        # Replace negative values and zeros with NaN for library size computation
        X_raw = X_raw.replace([np.inf, -np.inf], np.nan)

        # Select reference sample
        if self.ref_sample is not None:
            if self.ref_sample not in X_raw.columns:
                raise ValueError(f"Reference sample '{self.ref_sample}' not found in data")
            self.ref_sample_ = self.ref_sample
        else:
            self.ref_sample_ = self._select_reference_sample(X_raw)

        # Compute library sizes (sum of non-zero values per sample)
        lib_sizes = X_raw.replace(0, np.nan).sum(skipna=True)

        # Compute TMM factors for each sample
        factors = {}
        for sample in X_raw.columns:
            factors[sample] = self._compute_tmm_factor(X_raw, sample, lib_sizes)

        factors_series = pd.Series(factors)

        # Rescale so geometric mean = 1
        log_factors = np.log(factors_series)
        self.geometric_mean_ = np.exp(log_factors.mean())
        self.norm_factors_ = factors_series / self.geometric_mean_

        logger.info(f"TMM normalization fitted: {len(X_raw.columns)} samples, ref={self.ref_sample_}")
        logger.debug(f"Normalization factors: {self.norm_factors_.to_dict()}")

        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Apply TMM normalization to data.

        Computes: normalized = X / effective_library_size * scale_factor

        For proteomics, this scales each sample by its TMM factor,
        effectively adjusting for composition bias.

        Parameters
        ----------
        X : pd.DataFrame
            Intensity matrix (proteins x samples).

        Returns
        -------
        pd.DataFrame
            TMM-normalized intensities.
        """
        if self.norm_factors_ is None:
            raise ValueError("TMMNormalizer has not been fitted. Call fit() first.")

        # For proteomics, we simply divide by the norm factor
        # (Higher factor = sample had more signal from outliers = divide to reduce)
        # This is equivalent to multiplying by 1/factor

        result = X.copy()
        for sample in X.columns:
            if sample in self.norm_factors_.index:
                result[sample] = X[sample] / self.norm_factors_[sample]
            else:
                logger.warning(f"Sample {sample} not in fitted factors, leaving unchanged")

        return result

    def fit_transform(self, X: pd.DataFrame, y=None) -> pd.DataFrame:
        """Fit and transform in one step."""
        return self.fit(X, y).transform(X)

    def get_norm_factors(self) -> pd.Series:
        """
        Get the normalization factors.

        Returns
        -------
        pd.Series
            Normalization factor per sample.
        """
        if self.norm_factors_ is None:
            raise ValueError("TMMNormalizer has not been fitted. Call fit() first.")
        return self.norm_factors_.copy()

    def get_effective_library_sizes(self, X: pd.DataFrame) -> pd.Series:
        """
        Compute effective library sizes.

        effective = raw_library_size * norm_factor

        Parameters
        ----------
        X : pd.DataFrame
            Intensity matrix.

        Returns
        -------
        pd.Series
            Effective library size per sample.
        """
        if self.norm_factors_ is None:
            raise ValueError("TMMNormalizer has not been fitted. Call fit() first.")

        raw_lib_sizes = X.replace(0, np.nan).sum(skipna=True)
        return raw_lib_sizes * self.norm_factors_


def tmm_normalize(
    df: pd.DataFrame,
    m_trim: float = 0.3,
    a_trim: float = 0.05,
    ref_sample: Optional[str] = None,
    log_transform: bool = False,
) -> pd.DataFrame:
    """
    Convenience function to apply TMM normalization.

    Parameters
    ----------
    df : pd.DataFrame
        Intensity matrix (proteins x samples).
    m_trim : float, default=0.3
        Proportion of M-values to trim.
    a_trim : float, default=0.05
        Proportion of A-values to trim.
    ref_sample : str, optional
        Reference sample. Auto-selected if None.
    log_transform : bool, default=False
        If True, input data is already log2-transformed.

    Returns
    -------
    pd.DataFrame
        TMM-normalized data.

    Examples
    --------
    >>> from mokume.normalization.tmm import tmm_normalize
    >>> normalized = tmm_normalize(raw_intensities)
    """
    return TMMNormalizer(
        m_trim=m_trim,
        a_trim=a_trim,
        ref_sample=ref_sample,
        log_transform=log_transform,
    ).fit_transform(df)
