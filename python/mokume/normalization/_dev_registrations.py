"""
Registration of dev-only normalizers absent from main's registry.

Main already registers feature normalizers
(global/iqr/max/max_min/mean/median/none) and sample normalizers
(conditionmedian/globalmedian/hierarchical/irs/none/tmm). The dev branch
additionally implements RLR, LOESS, and quantile normalization but never
registered them with :class:`~mokume.core.registry.PluginRegistry`.

This module supplies thin adapters that DELEGATE to the dev algorithm
implementations without touching their math:

- ``rlr``      -> :class:`mokume.normalization.rlr.RlrNormalizer`
                  (sample-level; per-protein-median pseudo-reference).
- ``loess``    -> :class:`mokume.normalization.loess.LOESSNormalizer`
                  (feature group per the registration plan; MA-plot LOWESS
                  bias correction across samples).
- ``quantile`` -> :class:`mokume.normalization.distribution.QuantileNormalizer`
                  (feature group per the registration plan; wraps
                  :func:`mokume.normalization.protein.quantile_normalize`).

Each adapter pivots the long feature/protein table to a wide
(proteins x samples) matrix, calls the dev normalizer's ``fit_transform``
/ ``transform``, and melts the result back to long format, mirroring the
existing dataset-level adapters in
:mod:`mokume.normalization.sample_normalizers`.
"""

import numpy as np
import pandas as pd

from mokume.core.constants import (
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PROTEIN_NAME,
    SAMPLE_ID,
)
from mokume.core.logger import get_logger
from mokume.core.registry import PluginRegistry
from mokume.normalization.base import FeatureNormalizer, SampleNormalizer
from mokume.normalization.distribution import QuantileNormalizer
from mokume.normalization.loess import LOESSNormalizer
from mokume.normalization.rlr import RlrNormalizer

logger = get_logger("mokume.normalization._dev_registrations")


def _pivot_long_to_wide(
    df: pd.DataFrame,
    intensity_col: str,
    sample_col: str,
    protein_col: str,
    peptide_col: str,
) -> pd.DataFrame:
    """Pivot a long feature table to a wide (rows x samples) matrix.

    Uses ``[protein_col, peptide_col]`` as the row index when both are
    present, falling back to whichever identifier columns exist so the
    adapters also work on protein-level tables lacking a peptide column.
    """
    index_cols = [c for c in (protein_col, peptide_col) if c in df.columns]
    if not index_cols:
        raise ValueError(
            "Cannot pivot to wide format: neither "
            f"'{protein_col}' nor '{peptide_col}' present in columns "
            f"{list(df.columns)}"
        )
    wide = df.pivot_table(
        index=index_cols,
        columns=sample_col,
        values=intensity_col,
        aggfunc="sum",
    )
    return wide, index_cols


def _melt_wide_to_long(
    wide: pd.DataFrame,
    index_cols: list,
    intensity_col: str,
    sample_col: str,
) -> pd.DataFrame:
    """Melt a wide (rows x samples) matrix back to long format."""
    result = wide.reset_index().melt(
        id_vars=index_cols,
        var_name=sample_col,
        value_name=intensity_col,
    )
    return result.dropna(subset=[intensity_col])


@PluginRegistry.register("normalization.sample", "rlr")
class RlrSampleNormalizerPlugin(SampleNormalizer):
    """Adapter for :class:`~mokume.normalization.rlr.RlrNormalizer`.

    Robust Linear Regression normalization (NormalyzerDE-style). Fits a
    per-sample robust regression against the per-protein median in log2
    space, so it needs the full dataset at once (``is_dataset_level``).

    Operates on the full dataset: long -> wide -> fit_transform -> long.
    """

    def __init__(self, log_transform: bool = False):
        self.log_transform = log_transform

    @property
    def name(self) -> str:
        return "rlr"

    @property
    def is_dataset_level(self) -> bool:
        return True

    def normalize(
        self,
        df,
        intensity_col=NORM_INTENSITY,
        sample_col=SAMPLE_ID,
        condition_col=None,
        **kwargs,
    ):
        protein_col = kwargs.get("protein_col", PROTEIN_NAME)
        peptide_col = kwargs.get("peptide_col", PEPTIDE_CANONICAL)
        log_transform = kwargs.get("log_transform", self.log_transform)

        logger.info("Applying RLR sample normalization...")

        wide, index_cols = _pivot_long_to_wide(
            df, intensity_col, sample_col, protein_col, peptide_col
        )
        wide = wide.replace(0, np.nan)

        normalizer = RlrNormalizer(log_transform=log_transform)
        normalized_wide = normalizer.fit_transform(wide)

        result = _melt_wide_to_long(
            normalized_wide, index_cols, intensity_col, sample_col
        )
        logger.info("RLR normalization complete: %d rows", len(result))
        return result


@PluginRegistry.register("normalization.feature", "loess")
class LoessFeatureNormalizerPlugin(FeatureNormalizer):
    """Adapter for :class:`~mokume.normalization.loess.LOESSNormalizer`.

    LOWESS-based intensity-dependent bias correction on MA plots. The dev
    implementation operates on a wide (proteins x samples) log-scale
    matrix, so this adapter overrides ``normalize`` to pivot long -> wide,
    delegate to the dev normalizer, and melt back. ``transform_series`` is
    implemented only to satisfy the abstract base; the per-run
    orchestration in :meth:`FeatureNormalizer.normalize` is bypassed.
    """

    def __init__(self, frac: float = 0.75, reference: str = "median"):
        self.frac = frac
        self.reference = reference

    @property
    def name(self) -> str:
        return "loess"

    def transform_series(self, series: pd.Series) -> float:
        # LOESS is a cross-sample method; the base-class per-run pathway
        # (which would call this) is not used because normalize() is
        # overridden. Return a harmless scalar matching the base signature.
        return 1.0

    def normalize(
        self,
        df,
        intensity_col=NORM_INTENSITY,
        group_col=None,
        sample_col=SAMPLE_ID,
        **kwargs,
    ):
        protein_col = kwargs.get("protein_col", PROTEIN_NAME)
        peptide_col = kwargs.get("peptide_col", PEPTIDE_CANONICAL)
        frac = kwargs.get("frac", self.frac)
        reference = kwargs.get("reference", self.reference)

        logger.info("Applying LOESS feature normalization...")

        wide, index_cols = _pivot_long_to_wide(
            df, intensity_col, sample_col, protein_col, peptide_col
        )
        wide = wide.replace(0, np.nan)
        wide_log2 = np.log2(wide)

        normalizer = LOESSNormalizer(frac=frac, reference=reference)
        normalized_log2 = normalizer.fit_transform(wide_log2)

        normalized_wide = 2**normalized_log2

        result = _melt_wide_to_long(
            normalized_wide, index_cols, intensity_col, sample_col
        )
        logger.info("LOESS normalization complete: %d rows", len(result))
        return result


@PluginRegistry.register("normalization.feature", "quantile")
class QuantileFeatureNormalizerPlugin(FeatureNormalizer):
    """Adapter for
    :class:`~mokume.normalization.distribution.QuantileNormalizer`.

    Quantile normalization aligns every sample onto a shared reference
    distribution. The dev implementation operates on a wide
    (proteins x samples) matrix, so this adapter overrides ``normalize``
    to pivot long -> wide, delegate to the dev normalizer, and melt back.
    ``transform_series`` is implemented only to satisfy the abstract base.
    """

    @property
    def name(self) -> str:
        return "quantile"

    def transform_series(self, series: pd.Series) -> float:
        # Quantile normalization is a cross-sample method; the base-class
        # per-run pathway is not used because normalize() is overridden.
        # Return a harmless scalar matching the base signature.
        return 1.0

    def normalize(
        self,
        df,
        intensity_col=NORM_INTENSITY,
        group_col=None,
        sample_col=SAMPLE_ID,
        **kwargs,
    ):
        protein_col = kwargs.get("protein_col", PROTEIN_NAME)
        peptide_col = kwargs.get("peptide_col", PEPTIDE_CANONICAL)

        logger.info("Applying quantile feature normalization...")

        wide, index_cols = _pivot_long_to_wide(
            df, intensity_col, sample_col, protein_col, peptide_col
        )
        wide = wide.replace(0, np.nan)

        normalizer = QuantileNormalizer()
        normalized_wide = normalizer.fit_transform(wide)

        result = _melt_wide_to_long(
            normalized_wide, index_cols, intensity_col, sample_col
        )
        logger.info("Quantile normalization complete: %d rows", len(result))
        return result
