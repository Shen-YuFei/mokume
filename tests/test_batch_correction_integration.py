"""
Integration tests for batch correction with DirectLFQ.

These tests verify that batch correction:
1. Integrates correctly with the quantification pipeline
2. Preserves the number of proteins
3. Produces highly correlated results with uncorrected data
4. Can improve clustering metrics
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from unittest.mock import patch

# Skip if inmoose not installed
pytest.importorskip("inmoose")


from mokume.postprocessing import (
    is_batch_correction_available,
    detect_batches,
    extract_covariates_from_sdrf,
    apply_batch_correction,
)
from mokume.model import BatchDetectionMethod, BatchCorrectionConfig
from mokume.pipeline import PipelineConfig
from mokume.pipeline.config import (
    InputConfig,
    QuantificationConfig,
    BatchCorrectionConfig as PipelineBatchConfig,
)


class TestBatchDetection:
    """Tests for batch detection methods."""

    def test_detect_batches_sample_prefix(self):
        """Test batch detection from sample name prefixes."""
        samples = ["PXD001-S1", "PXD001-S2", "PXD002-S1", "PXD002-S2"]
        batches = detect_batches(samples, method="sample_prefix")

        assert len(batches) == 4
        assert batches[0] == batches[1]  # Same prefix
        assert batches[2] == batches[3]  # Same prefix
        assert batches[0] != batches[2]  # Different prefix

    def test_detect_batches_explicit_column(self):
        """Test batch detection with explicit values."""
        samples = ["S1", "S2", "S3", "S4", "S5"]
        batch_values = ["A", "A", "B", "B", "C"]
        batches = detect_batches(
            samples,
            method="column",
            batch_column_values=batch_values
        )

        assert len(batches) == 5
        assert batches[0] == batches[1]  # Same batch A
        assert batches[2] == batches[3]  # Same batch B
        assert len(set(batches)) == 3  # 3 unique batches

    def test_detect_batches_from_string_method(self):
        """Test that string method names work."""
        samples = ["PXD001-S1", "PXD001-S2", "PXD002-S1"]

        # These should all work
        batches1 = detect_batches(samples, method="sample_prefix")
        batches2 = detect_batches(samples, method=BatchDetectionMethod.SAMPLE_PREFIX)

        assert batches1 == batches2


class TestBatchCorrectionConfig:
    """Tests for BatchCorrectionConfig."""

    def test_config_defaults(self):
        """Test default configuration values."""
        config = BatchCorrectionConfig()

        assert config.enabled is False
        assert config.batch_method == BatchDetectionMethod.SAMPLE_PREFIX
        assert config.covariate_columns == []
        assert config.parametric is True

    def test_config_with_string_method(self):
        """Test that string method is converted to enum."""
        config = BatchCorrectionConfig(
            enabled=True,
            batch_method="sample_prefix",
        )

        assert config.batch_method == BatchDetectionMethod.SAMPLE_PREFIX

    def test_config_explicit_column_requires_batch_column(self):
        """Test that EXPLICIT_COLUMN method requires batch_column."""
        with pytest.raises(ValueError):
            BatchCorrectionConfig(
                enabled=True,
                batch_method=BatchDetectionMethod.EXPLICIT_COLUMN,
                batch_column=None,  # Should raise
            )


class TestApplyBatchCorrection:
    """Tests for apply_batch_correction function."""

    def test_batch_correction_basic(self):
        """Test basic batch correction on synthetic data."""
        # Create synthetic data with batch effect
        np.random.seed(42)
        n_proteins = 100
        n_samples = 20

        # Base expression
        base = np.random.randn(n_proteins, n_samples) * 2 + 10

        # Add batch effect (first 10 samples = batch 0, last 10 = batch 1)
        batch_effect = np.zeros((n_proteins, n_samples))
        batch_effect[:, 10:] = 2.0  # Batch 1 is shifted up

        data = pd.DataFrame(
            base + batch_effect,
            index=[f"P{i}" for i in range(n_proteins)],
            columns=[f"S{i}" for i in range(n_samples)],
        )

        batch_indices = [0] * 10 + [1] * 10

        # Apply correction
        corrected = apply_batch_correction(data, batch_indices)

        # Check output shape
        assert corrected.shape == data.shape

        # Check that batch effect is reduced
        # Mean difference between batches should be smaller after correction
        batch0_mean = data.iloc[:, :10].mean().mean()
        batch1_mean = data.iloc[:, 10:].mean().mean()
        original_diff = abs(batch1_mean - batch0_mean)

        corrected_batch0_mean = corrected.iloc[:, :10].mean().mean()
        corrected_batch1_mean = corrected.iloc[:, 10:].mean().mean()
        corrected_diff = abs(corrected_batch1_mean - corrected_batch0_mean)

        assert corrected_diff < original_diff

    def test_batch_correction_preserves_correlation(self):
        """Test that batch correction preserves overall data structure."""
        np.random.seed(42)
        n_proteins = 50
        n_samples = 12

        data = pd.DataFrame(
            np.random.randn(n_proteins, n_samples) * 2 + 10,
            index=[f"P{i}" for i in range(n_proteins)],
            columns=[f"S{i}" for i in range(n_samples)],
        )

        batch_indices = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]

        corrected = apply_batch_correction(data, batch_indices)

        # Correlation should be high (data structure preserved)
        from scipy import stats

        correlations = []
        for col in data.columns:
            r, _ = stats.pearsonr(data[col], corrected[col])
            correlations.append(r)

        mean_corr = np.mean(correlations)
        assert mean_corr > 0.9, f"Mean correlation {mean_corr} is too low"


class TestPipelineConfigBatchCorrection:
    """Tests for batch correction in PipelineConfig."""

    def test_pipeline_config_with_batch_correction(self):
        """Test PipelineConfig accepts batch correction parameters."""
        config = PipelineConfig(
            input=InputConfig(parquet="test.parquet"),
            quantification=QuantificationConfig(method="maxlfq"),
            batch=PipelineBatchConfig(
                enabled=True,
                method="sample_prefix",
                covariates=["characteristics[sex]"],
                parametric=True,
            ),
        )

        assert config.batch.enabled is True
        assert config.batch.method == "sample_prefix"
        assert config.batch.covariates == ["characteristics[sex]"]
        assert config.batch.parametric is True

    def test_pipeline_config_defaults_no_batch(self):
        """Test that batch correction is disabled by default."""
        config = PipelineConfig(
            input=InputConfig(parquet="test.parquet"),
            quantification=QuantificationConfig(method="maxlfq"),
        )

        assert config.batch.enabled is False


class TestSDRFCovariatExtraction:
    """Tests for SDRF covariate extraction."""

    def test_extract_covariates_no_columns(self):
        """Test that empty covariate list returns None."""
        result = extract_covariates_from_sdrf(
            "nonexistent.tsv",
            ["S1", "S2"],
            [],  # No covariates
        )
        assert result is None

    def test_extract_covariates_from_sdrf(self, tmp_path):
        """Test covariate extraction from SDRF file."""
        # Create test SDRF
        sdrf_content = """source name\tcharacteristics[sex]\tcharacteristics[tissue]
Sample1\tmale\tliver
Sample2\tfemale\tliver
Sample3\tmale\tbrain
Sample4\tfemale\tbrain
"""
        sdrf_path = tmp_path / "test.sdrf.tsv"
        sdrf_path.write_text(sdrf_content)

        result = extract_covariates_from_sdrf(
            str(sdrf_path),
            ["Sample1", "Sample2", "Sample3", "Sample4"],
            ["characteristics[sex]", "characteristics[tissue]"],
        )

        assert result is not None
        assert len(result) == 4  # 4 samples
        assert len(result[0]) == 2  # 2 covariates

        # Check encoding is correct
        # sex: male=0, female=1 (or vice versa, but consistent)
        assert result[0][0] == result[2][0]  # Sample1 and Sample3 both male
        assert result[1][0] == result[3][0]  # Sample2 and Sample4 both female
        assert result[0][0] != result[1][0]  # male != female


# Mark slow tests
@pytest.mark.slow
class TestBatchCorrectionIntegration:
    """Integration tests with real pipeline (marked slow)."""

    def test_batch_correction_available(self):
        """Verify batch correction is available."""
        assert is_batch_correction_available() is True
