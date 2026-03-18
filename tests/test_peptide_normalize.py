import logging

import pytest
from mokume.normalization.peptide import (
    peptide_normalization,
    SQLFilterBuilder,
    Feature,
)
from pathlib import Path

TESTS_DIR = Path(__file__).parent

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class TestSQLFilterBuilder:
    """Tests for the SQLFilterBuilder class."""

    def test_default_where_clause(self):
        """Test that default filter builder generates expected WHERE clause."""
        builder = SQLFilterBuilder()
        where_clause, params = builder.build_where_clause()

        # Should include intensity > 0
        if "intensity > 0" not in where_clause:
            raise AssertionError("Missing 'intensity > 0' in where_clause")
        # Should include peptide length filter (parameterized)
        if 'LENGTH("sequence") >= ?' not in where_clause:
            raise AssertionError("Missing LENGTH filter in where_clause")
        if 7 not in params:
            raise AssertionError("Missing 7 in params")
        # Should include unique peptide filter
        if '"unique" = 1' not in where_clause:
            raise AssertionError("Missing unique filter in where_clause")
        # Should include contaminant filters (parameterized with ? placeholders)
        if "NOT LIKE ?" not in where_clause:
            raise AssertionError("Missing NOT LIKE placeholder in where_clause")
        if "%CONTAMINANT%" not in params:
            raise AssertionError("Missing %CONTAMINANT% in params")
        if "%DECOY%" not in params:
            raise AssertionError("Missing %DECOY% in params")
        if "%ENTRAP%" not in params:
            raise AssertionError("Missing %ENTRAP% in params")

    def test_custom_contaminant_patterns(self):
        """Test filter builder with custom contaminant patterns."""
        builder = SQLFilterBuilder(
            contaminant_patterns=["CONTAM", "REV_"],
            min_peptide_length=5,
        )
        where_clause, params = builder.build_where_clause()

        if "%CONTAM%" not in params:
            raise AssertionError("Missing %CONTAM% in params")
        if "%REV_%" not in params:
            raise AssertionError("Missing %REV_% in params")
        if "%DECOY%" in params:
            raise AssertionError("%DECOY% should not be in params")
        if 'LENGTH("sequence") >= ?' not in where_clause:
            raise AssertionError("Missing LENGTH filter")
        if 5 not in params:
            raise AssertionError("Missing 5 in params")

    def test_disable_contaminant_filter(self):
        """Test that contaminant filter can be disabled."""
        builder = SQLFilterBuilder(remove_contaminants=False)
        where_clause, params = builder.build_where_clause()

        if "NOT LIKE" in where_clause:
            raise AssertionError("NOT LIKE should not be in where_clause")
        if any("%" in str(p) for p in params):
            raise AssertionError("No pattern params expected when contaminants disabled")
        # Other filters should still be present
        if "intensity > 0" not in where_clause:
            raise AssertionError("Missing 'intensity > 0'")

    def test_min_intensity_threshold(self):
        """Test that min intensity threshold is applied."""
        builder = SQLFilterBuilder(min_intensity=1000.0)
        where_clause, params = builder.build_where_clause()

        if "intensity >= ?" not in where_clause:
            raise AssertionError("Missing intensity threshold placeholder")
        if 1000.0 not in params:
            raise AssertionError("Missing 1000.0 in params")

    def test_disable_unique_requirement(self):
        """Test that unique peptide requirement can be disabled."""
        builder = SQLFilterBuilder(require_unique=False)
        where_clause, _params = builder.build_where_clause()

        if '"unique" = 1' in where_clause:
            raise AssertionError("unique filter should not be present")


class TestFeatureWideFormat:
    """Tests for Feature class with wide format parquet (quantms.io/qpx)."""

    @pytest.fixture
    def feature_path(self):
        """Path to test feature parquet file (wide format)."""
        return str(TESTS_DIR / "example/feature_wide.parquet")

    @pytest.fixture
    def sdrf_path(self):
        """Path to test SDRF file."""
        return str(TESTS_DIR / "example/PXD020192.sdrf.tsv")

    def test_feature_loads_wide_format(self, feature_path):
        """Test that Feature class can load wide format parquet with UNNEST."""
        feature = Feature(feature_path)

        # Should have samples from UNNEST
        samples = feature.get_unique_samples()
        if len(samples) == 0:
            raise AssertionError("Expected at least one sample")

    def test_feature_with_filter_builder(self, feature_path):
        """Test Feature class accepts and stores filter_builder."""
        builder = SQLFilterBuilder(
            remove_contaminants=True,
            min_peptide_length=7,
        )
        feature = Feature(feature_path, filter_builder=builder)

        if feature.filter_builder is None:
            raise AssertionError("filter_builder should not be None")
        if feature.filter_builder.remove_contaminants is not True:
            raise AssertionError("remove_contaminants should be True")
        if feature.filter_builder.min_peptide_length != 7:
            raise AssertionError(f"Expected min_peptide_length=7, got {feature.filter_builder.min_peptide_length}")

    def test_get_median_map(self, feature_path):
        """Test that get_median_map works with wide format."""
        feature = Feature(feature_path)
        med_map = feature.get_median_map()

        # Should return results
        if len(med_map) == 0:
            raise AssertionError("Expected non-empty median map")
        # All values should be positive
        for sample, factor in med_map.items():
            if factor <= 0:
                raise AssertionError(f"Expected positive factor for {sample}, got {factor}")

    def test_get_median_map_with_filter(self, feature_path):
        """Test that get_median_map uses filter_builder when provided."""
        # Without filter
        feature_no_filter = Feature(feature_path)
        med_map_unfiltered = feature_no_filter.get_median_map()

        # With filter (excluding contaminants)
        builder = SQLFilterBuilder(remove_contaminants=True)
        feature_filtered = Feature(feature_path, filter_builder=builder)
        med_map_filtered = feature_filtered.get_median_map()

        # Both should return results
        if len(med_map_unfiltered) == 0:
            raise AssertionError("Expected non-empty unfiltered median map")
        if len(med_map_filtered) == 0:
            raise AssertionError("Expected non-empty filtered median map")

        # The samples should be the same (filtering applies to features, not samples)
        if set(med_map_unfiltered.keys()) != set(med_map_filtered.keys()):
            raise AssertionError("Sample keys mismatch between filtered and unfiltered")

    def test_get_low_frequency_peptides(self, feature_path):
        """Test that get_low_frequency_peptides works with wide format."""
        feature = Feature(feature_path)
        low_freq = feature.get_low_frequency_peptides()

        # Should return a tuple
        if not isinstance(low_freq, tuple):
            raise AssertionError(f"Expected tuple, got {type(low_freq)}")

    def test_get_median_map_to_condition(self, feature_path):
        """Test that get_median_map_to_condition works with wide format."""
        feature = Feature(feature_path)
        med_map = feature.get_median_map_to_condition()

        # Should return dict of dicts
        if not isinstance(med_map, dict):
            raise AssertionError(f"Expected dict, got {type(med_map)}")
        for condition, samples in med_map.items():
            if not isinstance(samples, dict):
                raise AssertionError(f"Expected dict for condition {condition}, got {type(samples)}")

    def test_enrich_with_sdrf(self, feature_path, sdrf_path):
        """Test that enrich_with_sdrf enriches data with SDRF metadata."""
        feature = Feature(feature_path)

        # Before enrichment, condition should equal sample_accession
        conditions_before = feature.get_unique_conditions()
        _ = feature.get_unique_samples()
        # Conditions default to sample_accession
        if len(conditions_before) == 0:
            raise AssertionError("Expected at least one condition before enrichment")

        # Enrich with SDRF
        feature.enrich_with_sdrf(sdrf_path)

        # After enrichment, conditions should be from SDRF
        conditions_after = feature.get_unique_conditions()
        if len(conditions_after) == 0:
            raise AssertionError("Expected at least one condition after enrichment")

    def test_iter_samples(self, feature_path):
        """Test that iter_samples works with wide format."""
        feature = Feature(feature_path)

        count = 0
        for samples, df in feature.iter_samples(sample_num=5):
            if len(samples) > 5:
                raise AssertionError(f"Batch too large: {len(samples)} > 5")
            if len(df) == 0:
                raise AssertionError("Expected non-empty DataFrame in batch")
            count += 1

        if count == 0:
            raise AssertionError("Expected at least one batch")


class TestFeatureNewQPXFormat:
    """Tests for Feature class with new QPX format parquet (charge, run_file_name, {label, intensity})."""

    @pytest.fixture
    def feature_path(self):
        return str(TESTS_DIR / "example/feature_wide_new_qpx.parquet")

    @pytest.fixture
    def sdrf_path(self):
        return str(TESTS_DIR / "example/sdrf_new_qpx.tsv")

    def test_loads_new_qpx_format(self, feature_path):
        """Test that Feature detects and loads new QPX schema."""
        feature = Feature(feature_path)
        if feature._is_new_qpx is not True:
            raise AssertionError("Expected _is_new_qpx to be True")
        if feature._charge_col != "charge":
            raise AssertionError(f"Expected _charge_col='charge', got '{feature._charge_col}'")
        if feature._run_col != "run_file_name":
            raise AssertionError(f"Expected _run_col='run_file_name', got '{feature._run_col}'")

        samples = feature.get_unique_samples()
        if len(samples) == 0:
            raise AssertionError("Expected at least one sample")

    def test_unnested_columns_present(self, feature_path):
        """Test that unnested view has expected column names."""
        feature = Feature(feature_path)
        df = feature.parquet_db.sql("SELECT * FROM parquet_db LIMIT 1").df()
        for col in ["charge", "run_file_name", "sample_accession", "channel", "intensity", "condition"]:
            if col not in df.columns:
                raise AssertionError(f"Missing column: {col}")

    def test_enrich_with_sdrf_maps_sample_accession(self, feature_path, sdrf_path):
        """Test that enrich_with_sdrf correctly maps (run_file_name, label) -> source name."""
        feature = Feature(feature_path)
        feature.enrich_with_sdrf(sdrf_path)

        samples = feature.get_unique_samples()
        # After SDRF enrichment, samples should be SDRF source names
        if "Sample_A_126" not in samples:
            raise AssertionError(f"'Sample_A_126' not found in samples: {samples}")
        if "Sample_A_127N" not in samples:
            raise AssertionError(f"'Sample_A_127N' not found in samples: {samples}")

    def test_enrich_with_sdrf_maps_condition(self, feature_path, sdrf_path):
        """Test that enrich_with_sdrf maps conditions from SDRF factor values."""
        feature = Feature(feature_path)
        feature.enrich_with_sdrf(sdrf_path)

        conditions = feature.get_unique_conditions()
        if "normal" not in conditions:
            raise AssertionError(f"'normal' not found in conditions: {conditions}")
        if "disease" not in conditions:
            raise AssertionError(f"'disease' not found in conditions: {conditions}")

    def test_get_median_map(self, feature_path):
        """Test get_median_map works with new QPX format."""
        feature = Feature(feature_path)
        med_map = feature.get_median_map()
        if len(med_map) == 0:
            raise AssertionError("Expected non-empty median map")
        for sample, factor in med_map.items():
            if factor <= 0:
                raise AssertionError(f"Expected positive factor for {sample}, got {factor}")


class TestPeptideNormalizationWideFormat:
    """Tests for peptide_normalization with wide format parquet."""

    def test_normalization_without_sdrf(self):
        """Test peptide normalization without SDRF (uses defaults)."""
        args = {
            "parquet": str(TESTS_DIR / "example/feature_wide.parquet"),
            "sdrf": None,
            "min_aa": 7,
            "min_unique": 1,
            "remove_ids": None,
            "remove_decoy_contaminants": True,
            "remove_low_frequency_peptides": False,
            "output": str(TESTS_DIR / "example" / "out" / "PXD020192-no-sdrf.csv"),
            "skip_normalization": False,
            "nmethod": "median",
            "pnmethod": "globalMedian",
            "log2": True,
            "save_parquet": False,
        }

        out = Path(args["output"])
        if out.exists():
            out.unlink()

        # Should complete without errors
        peptide_normalization(**args)

        # Output should exist
        if not out.exists():
            raise AssertionError("Output file was not created")

        # Clean up
        if out.exists():
            out.unlink()

    def test_normalization_with_sdrf(self):
        """Test peptide normalization with SDRF enrichment."""
        args = {
            "parquet": str(TESTS_DIR / "example/feature_wide.parquet"),
            "sdrf": str(TESTS_DIR / "example/PXD020192.sdrf.tsv"),
            "min_aa": 7,
            "min_unique": 1,
            "remove_ids": None,
            "remove_decoy_contaminants": True,
            "remove_low_frequency_peptides": False,
            "output": str(TESTS_DIR / "example" / "out" / "PXD020192-with-sdrf.csv"),
            "skip_normalization": False,
            "nmethod": "median",
            "pnmethod": "globalMedian",
            "log2": True,
            "save_parquet": False,
        }

        out = Path(args["output"])
        if out.exists():
            out.unlink()

        # Should complete without errors
        peptide_normalization(**args)

        # Output should exist
        if not out.exists():
            raise AssertionError("Output file was not created")

        # Clean up
        if out.exists():
            out.unlink()

    def test_normalization_with_filter_config(self):
        """Test peptide normalization with PreprocessingFilterConfig."""
        from mokume.model.filters import PreprocessingFilterConfig

        # Create a filter config
        filter_config = PreprocessingFilterConfig(
            name="test_config",
            enabled=True,
        )
        # Set contaminant removal
        filter_config.protein.remove_contaminants = True
        filter_config.protein.remove_decoys = True

        args = {
            "parquet": str(TESTS_DIR / "example/feature_wide.parquet"),
            "sdrf": str(TESTS_DIR / "example/PXD020192.sdrf.tsv"),
            "min_aa": 7,
            "min_unique": 1,
            "remove_ids": None,
            "remove_decoy_contaminants": True,
            "remove_low_frequency_peptides": False,
            "output": str(TESTS_DIR / "example" / "out" / "PXD020192-filtered.csv"),
            "skip_normalization": False,
            "nmethod": "median",
            "pnmethod": "globalMedian",
            "log2": True,
            "save_parquet": False,
            "filter_config": filter_config,
        }

        out = Path(args["output"])
        if out.exists():
            out.unlink()

        # Should complete without errors
        peptide_normalization(**args)

        # Output should exist
        if not out.exists():
            raise AssertionError("Output file was not created")

        # Clean up
        if out.exists():
            out.unlink()
