"""
Tests for preprocessing filter configurations.

This module tests that all example filter configurations can be loaded
and that the filter system works correctly.
"""

import pytest
from pathlib import Path

from mokume.preprocessing.filters import (
    load_filter_config,
    get_filter_pipeline,
    generate_example_config,
)
from mokume.model.filters import PreprocessingFilterConfig


# Path to example filter configurations
EXAMPLE_FILTERS_DIR = Path(__file__).parent / "example" / "filters"


def get_example_config_files():
    """Get all example filter configuration files."""
    if not EXAMPLE_FILTERS_DIR.exists():
        return []
    return list(EXAMPLE_FILTERS_DIR.glob("*.yaml"))


class TestFilterConfigurations:
    """Tests for filter configuration loading and validation."""

    @pytest.mark.parametrize(
        "config_file",
        get_example_config_files(),
        ids=lambda p: p.stem,
    )
    def test_load_example_config(self, config_file):
        """Test that example configurations can be loaded."""
        config = load_filter_config(config_file)

        assert config is not None
        assert isinstance(config, PreprocessingFilterConfig)
        assert config.name is not None
        assert config.enabled is True

    @pytest.mark.parametrize(
        "config_file",
        get_example_config_files(),
        ids=lambda p: p.stem,
    )
    def test_create_pipeline_from_config(self, config_file):
        """Test that filter pipelines can be created from configurations."""
        config = load_filter_config(config_file)
        pipeline = get_filter_pipeline(config)

        assert pipeline is not None
        assert len(pipeline) >= 0  # Some configs may have no active filters

    def test_basic_qc_config(self):
        """Test the basic QC configuration specifically."""
        config_path = EXAMPLE_FILTERS_DIR / "basic_qc.yaml"
        if not config_path.exists():
            pytest.skip("basic_qc.yaml not found")

        config = load_filter_config(config_path)

        assert config.name == "basic_qc"
        assert config.protein.min_unique_peptides == 2
        assert config.protein.remove_contaminants is True
        assert config.peptide.min_peptide_length == 7

    def test_stringent_filtering_config(self):
        """Test the stringent filtering configuration specifically."""
        config_path = EXAMPLE_FILTERS_DIR / "stringent_filtering.yaml"
        if not config_path.exists():
            pytest.skip("stringent_filtering.yaml not found")

        config = load_filter_config(config_path)

        assert config.name == "stringent_filtering"
        assert config.intensity.min_intensity == 1000.0
        assert config.intensity.cv_threshold == 0.3
        assert config.peptide.allowed_charge_states == [2, 3, 4]
        assert "Oxidation" in config.peptide.exclude_modifications

    def test_generate_example_config(self, tmp_path):
        """Test that example config generation works."""
        yaml_path = tmp_path / "test_config.yaml"
        generate_example_config(yaml_path)

        assert yaml_path.exists()

        # Load and verify
        config = load_filter_config(yaml_path)
        assert config.name == "example_config"

    def test_generate_json_config(self, tmp_path):
        """Test that JSON config generation works."""
        json_path = tmp_path / "test_config.json"
        generate_example_config(json_path, format="json")

        assert json_path.exists()

        # Load and verify
        config = load_filter_config(json_path)
        assert config is not None

    def test_config_apply_overrides(self):
        """Test that CLI overrides work correctly."""
        config = PreprocessingFilterConfig(name="test")

        config.apply_overrides(
            {
                "min_intensity": 500.0,
                "cv_threshold": 0.25,
                "charge_states": [2, 3],
                "min_unique_peptides": 3,
                "max_missing_rate": 0.4,
            }
        )

        assert config.intensity.min_intensity == 500.0
        assert config.intensity.cv_threshold == 0.25
        assert config.peptide.allowed_charge_states == [2, 3]
        assert config.protein.min_unique_peptides == 3
        assert config.run_qc.max_missing_rate == 0.4


class TestFilterPipeline:
    """Tests for filter pipeline functionality."""

    def test_empty_pipeline(self):
        """Test that disabled config creates empty pipeline."""
        config = PreprocessingFilterConfig(name="disabled", enabled=False)
        pipeline = get_filter_pipeline(config)

        assert len(pipeline) == 0

    def test_pipeline_with_intensity_filters(self):
        """Test pipeline with intensity filters."""
        config = PreprocessingFilterConfig(name="test")
        config.intensity.min_intensity = 100.0
        config.intensity.cv_threshold = 0.5

        pipeline = get_filter_pipeline(config)

        # Should have at least MinIntensityFilter and CVThresholdFilter
        filter_names = [f.name for f in pipeline.filters]
        assert "MinIntensityFilter" in filter_names
        assert "CVThresholdFilter" in filter_names

    def test_pipeline_with_protein_filters(self):
        """Test pipeline with protein filters."""
        config = PreprocessingFilterConfig(name="test")
        config.protein.min_unique_peptides = 2
        config.protein.remove_contaminants = True

        pipeline = get_filter_pipeline(config)

        filter_names = [f.name for f in pipeline.filters]
        assert "ContaminantFilter" in filter_names
        assert "MinPeptideFilter" in filter_names
