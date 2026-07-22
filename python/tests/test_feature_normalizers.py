"""Unit tests for feature-level (within-run) normalizers.

The base ``FeatureNormalizer.normalize()`` sums a per-run SCALAR metric
from ``transform_series()``, averages it across runs, and rescales each
run so its metric equals the sample average. These tests guard that
contract: every concrete ``transform_series`` must return a scalar float
(not a Series), and ``normalize()`` must produce finite intensities.
"""

import numbers

import numpy as np
import pandas as pd
import pytest

import mokume._plugins  # noqa: F401  # side-effect: populate the registry
from mokume.core.constants import (
    NORM_INTENSITY,
    PEPTIDE_CANONICAL,
    PROTEIN_NAME,
    SAMPLE_ID,
    TECHREPLICATE,
)
from mokume.core.registry import PluginRegistry

# Concrete scalar-metric normalizers exercised through normalize().
SCALAR_METRIC_NORMALIZERS = ["mean", "median", "max", "iqr", "max_min", "global"]


@pytest.fixture
def multi_run_df():
    """Tiny long feature table: 2 samples x 2 runs, 3 shared peptides.

    Run "2" carries a different scale/offset from run "1" within each
    sample, giving ``normalize()`` genuine per-run metrics to balance.
    """
    rows = []
    for sample in ("S1", "S2"):
        for run in ("1", "2"):
            for i, pep in enumerate(("PEPA", "PEPB", "PEPC")):
                base = 100.0 * (i + 1)
                scale = 1.0 if run == "1" else 2.0
                offset = 10.0 if sample == "S2" else 0.0
                rows.append(
                    {
                        PROTEIN_NAME: f"P{i}",
                        PEPTIDE_CANONICAL: pep,
                        SAMPLE_ID: sample,
                        TECHREPLICATE: run,
                        NORM_INTENSITY: base * scale + offset,
                    }
                )
    return pd.DataFrame(rows)


@pytest.mark.parametrize("name", SCALAR_METRIC_NORMALIZERS)
def test_transform_series_returns_scalar(name, multi_run_df):
    normalizer = PluginRegistry.get("normalization.feature", name)
    metric = normalizer.transform_series(multi_run_df[NORM_INTENSITY])
    assert isinstance(metric, numbers.Real)
    assert not isinstance(metric, pd.Series)
    assert np.isfinite(float(metric))


@pytest.mark.parametrize("name", SCALAR_METRIC_NORMALIZERS)
def test_normalize_produces_finite_output(name, multi_run_df):
    normalizer = PluginRegistry.get("normalization.feature", name)
    result = normalizer.normalize(multi_run_df.copy())
    values = result[NORM_INTENSITY]
    assert values.notna().all()
    assert np.isfinite(values.to_numpy()).all()
    assert len(result) == len(multi_run_df)


def test_none_transform_series_is_harmless_scalar():
    normalizer = PluginRegistry.get("normalization.feature", "none")
    metric = normalizer.transform_series(pd.Series([1.0, 2.0, 3.0]))
    assert isinstance(metric, numbers.Real)
    assert metric == 1.0


@pytest.mark.parametrize("name", ["max_min", "iqr"])
def test_zero_range_metric_is_guarded(name):
    """A constant run has no spread; the metric must not divide by zero."""
    normalizer = PluginRegistry.get("normalization.feature", name)
    metric = normalizer.transform_series(pd.Series([5.0, 5.0, 5.0]))
    assert metric == 1.0


class TestLegacyEnumReplicateMetric:
    """The FeatureNormalizationMethod enum owes the same scalar contract.

    ``normalize_runs`` divides each run's intensities by that run's metric
    relative to the sample average, so a Series here misaligns on index, turns
    the running total into NaN, and blanks every intensity for the sample.
    """

    @pytest.mark.parametrize("name", ["Mean", "Median", "Max", "Global", "IQR"])
    def test_supported_methods_return_a_scalar(self, name):
        from mokume.model.normalization import FeatureNormalizationMethod

        values = pd.Series([1.0, 2.0, 3.0, 4.0])
        metric = getattr(FeatureNormalizationMethod, name)._replicate_metric(values)
        assert isinstance(metric, numbers.Real)
        assert not isinstance(metric, pd.Series)
        assert np.isfinite(metric)

    def test_rescaling_methods_say_so_instead_of_returning_a_series(self):
        """Max_Min has no single scale factor; it must refuse, not emit NaN."""
        from mokume.model.normalization import FeatureNormalizationMethod

        values = pd.Series([1.0, 2.0, 3.0, 4.0])
        with pytest.raises(NotImplementedError, match="scale factor"):
            FeatureNormalizationMethod.Max_Min._replicate_metric(values)

    def test_multi_run_sample_stays_finite(self):
        """The end-to-end path that used to produce an all-NaN sample."""
        from mokume.core.constants import NORM_INTENSITY, SAMPLE_ID, TECHREPLICATE
        from mokume.model.normalization import FeatureNormalizationMethod

        df = pd.DataFrame(
            {
                SAMPLE_ID: ["S1"] * 6,
                TECHREPLICATE: [1, 1, 1, 2, 2, 2],
                NORM_INTENSITY: [10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
            }
        )
        out = FeatureNormalizationMethod.Median.normalize_runs(df.copy(), 2)
        assert out[NORM_INTENSITY].notna().all()
        assert np.isfinite(out[NORM_INTENSITY]).all()
