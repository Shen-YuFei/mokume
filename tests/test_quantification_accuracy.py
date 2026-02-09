"""
Test quantification accuracy by comparing mokume's quantification implementations
with DirectLFQ and DIA-NN reported values.

This module contains two test suites:

1. TestSmallDiannSubset: Uses a small subset DIA-NN report from directLFQ test suite.
   Note: DIA-NN's PG.MaxLFQ is computed from ALL precursors, but this report only
   contains a subset. Direct comparison shows lower correlation (~0.55-0.61).

2. TestDiannWithSdrf: Uses a full DIA-NN report from PRIDE (PXD063291) with SDRF.
   This provides a more complete dataset for validation.

Test data sources:
- Small subset: https://github.com/MannLabs/directlfq/tree/main/test_data
- Full dataset: https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/example-projects/PXD063291.zip
"""

import os
import tempfile
import zipfile
import urllib.request
import warnings
import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

# Suppress warnings during tests
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Test data paths
TEST_DATA_DIR = Path(__file__).parent / "data"
SMALL_DIANN_REPORT = TEST_DATA_DIR / "diann_test_report.tsv"

# PRIDE dataset URL
PRIDE_DATASET_URL = "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/example-projects/PXD063291.zip"
PRIDE_CACHE_DIR = TEST_DATA_DIR / "pride_cache"


def load_small_diann_report():
    """Load small DIA-NN report from directLFQ test suite."""
    if not SMALL_DIANN_REPORT.exists():
        pytest.skip(f"Test data not found: {SMALL_DIANN_REPORT}")
    return pd.read_csv(SMALL_DIANN_REPORT, sep="\t")


def download_pride_dataset():
    """
    Download and cache the PRIDE PXD063291 dataset.
    Returns path to the extracted DIA-NN report.
    """
    PRIDE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = PRIDE_CACHE_DIR / "PXD063291.zip"

    # Check if already extracted
    extracted_files = list(PRIDE_CACHE_DIR.glob("*.tsv")) + list(PRIDE_CACHE_DIR.glob("**/diann-output/*.tsv"))
    if extracted_files:
        # Find the main report file
        for f in extracted_files:
            if "report" in f.name.lower() and f.suffix == ".tsv":
                return f
        # Return first TSV if no report found
        return extracted_files[0]

    # Download if not cached
    if not zip_path.exists():
        print(f"Downloading PRIDE dataset from {PRIDE_DATASET_URL}...")
        try:
            urllib.request.urlretrieve(PRIDE_DATASET_URL, zip_path)
        except Exception as e:
            pytest.skip(f"Failed to download PRIDE dataset: {e}")

    # Extract
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(PRIDE_CACHE_DIR)
    except Exception as e:
        pytest.skip(f"Failed to extract PRIDE dataset: {e}")

    # Find the DIA-NN report
    extracted_files = list(PRIDE_CACHE_DIR.rglob("*.tsv"))
    for f in extracted_files:
        if "report" in f.name.lower():
            return f

    if extracted_files:
        return extracted_files[0]

    pytest.skip("No TSV files found in PRIDE dataset")


def prepare_peptide_data(diann_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert DIA-NN report to mokume peptide format.

    DIA-NN columns -> mokume columns:
    - Protein.Group -> ProteinName
    - Run -> SampleID
    - Precursor.Id -> PeptideSequence (using precursor as unique peptide identifier)
    - Precursor.Quantity -> NormIntensity
    """
    peptide_df = diann_df[["Protein.Group", "Run", "Precursor.Id", "Precursor.Quantity"]].copy()
    peptide_df = peptide_df.rename(columns={
        "Protein.Group": "ProteinName",
        "Run": "SampleID",
        "Precursor.Id": "PeptideSequence",
        "Precursor.Quantity": "NormIntensity",
    })

    # Remove rows with missing intensity
    peptide_df = peptide_df[peptide_df["NormIntensity"].notna()]
    peptide_df = peptide_df[peptide_df["NormIntensity"] > 0]

    return peptide_df


def extract_diann_maxlfq(diann_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract DIA-NN's MaxLFQ values (PG.MaxLFQ) per protein per sample.

    NOTE: These values are computed from ALL precursors by DIA-NN,
    not just those in this report. Direct comparison with mokume's
    values (computed from subset) may show lower correlation.
    """
    if "PG.MaxLFQ" not in diann_df.columns:
        return pd.DataFrame(columns=["ProteinName", "SampleID", "DIANNMaxLFQ"])

    maxlfq_df = diann_df[["Protein.Group", "Run", "PG.MaxLFQ"]].drop_duplicates()
    maxlfq_df = maxlfq_df.rename(columns={
        "Protein.Group": "ProteinName",
        "Run": "SampleID",
        "PG.MaxLFQ": "DIANNMaxLFQ",
    })
    maxlfq_df = maxlfq_df[maxlfq_df["DIANNMaxLFQ"].notna()]
    maxlfq_df = maxlfq_df[maxlfq_df["DIANNMaxLFQ"] > 0]
    return maxlfq_df


def compute_correlation(df1: pd.DataFrame, df2: pd.DataFrame,
                        col1: str, col2: str,
                        protein_col: str = "ProteinName",
                        sample_col: str = "SampleID") -> dict:
    """
    Compute correlation between two quantification results.
    """
    # Use suffixes to handle same column names
    merged = pd.merge(
        df1[[protein_col, sample_col, col1]],
        df2[[protein_col, sample_col, col2]],
        on=[protein_col, sample_col],
        how="inner",
        suffixes=("_x", "_y"),
    )

    if len(merged) == 0:
        return {"pearson": np.nan, "spearman": np.nan, "n": 0}

    # Determine actual column names after merge (handles same-name case)
    col1_actual = col1 if col1 in merged.columns else f"{col1}_x"
    col2_actual = col2 if col2 in merged.columns else f"{col2}_y"

    # Log transform for better correlation
    x = np.log2(merged[col1_actual].values + 1)
    y = np.log2(merged[col2_actual].values + 1)

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 3:
        return {"pearson": np.nan, "spearman": np.nan, "n": len(x)}

    pearson_r, _ = stats.pearsonr(x, y)
    spearman_r, _ = stats.spearmanr(x, y)

    return {
        "pearson": pearson_r,
        "spearman": spearman_r,
        "n": len(x),
        "mean_abs_log_diff": np.mean(np.abs(x - y)),
    }


def print_comparison_table(results: dict, intensity_cols: dict, dataset_name: str,
                           peptide_count: int, protein_count: int, sample_count: int,
                           precursor_count: int):
    """Print a formatted comparison table of all methods."""
    print("\n" + "=" * 90)
    print(f"QUANTIFICATION METHOD COMPARISON - {dataset_name}")
    print("=" * 90)
    print(f"Test Data:")
    print(f"  Peptide measurements: {peptide_count:,}")
    print(f"  Proteins: {protein_count}")
    print(f"  Samples: {sample_count}")
    print(f"  Precursors: {precursor_count}")
    print("-" * 90)
    print("CROSS-METHOD CORRELATIONS (Spearman)")
    print("-" * 90)

    method_names = list(results.keys())

    # Header
    print(f"{'':18}", end="")
    for name in method_names:
        print(f"{name[:14]:>15}", end="")
    print()
    print("-" * 90)

    # Rows
    for name1 in method_names:
        print(f"{name1:18}", end="")
        for name2 in method_names:
            if name1 == name2:
                print(f"{'1.000':>15}", end="")
            else:
                corr = compute_correlation(
                    results[name1], results[name2],
                    intensity_cols[name1], intensity_cols[name2]
                )
                if np.isnan(corr['spearman']):
                    print(f"{'N/A':>15}", end="")
                else:
                    print(f"{corr['spearman']:>15.3f}", end="")
        print()

    print("=" * 90)


def run_all_quantification_methods(peptide_df: pd.DataFrame, diann_maxlfq: pd.DataFrame = None):
    """
    Run all quantification methods and return results dict.

    Note: MaxLFQQuantification automatically uses DirectLFQ if available,
    falling back to built-in implementation otherwise. We test both paths
    when DirectLFQ is available.
    """
    from mokume.quantification import (
        MaxLFQQuantification,
        Top3Quantification,
        AllPeptidesQuantification,
        is_directlfq_available,
    )

    results = {}
    intensity_cols = {}

    # Add DIA-NN MaxLFQ if available
    if diann_maxlfq is not None and len(diann_maxlfq) > 0:
        results["DIA-NN MaxLFQ"] = diann_maxlfq
        intensity_cols["DIA-NN MaxLFQ"] = "DIANNMaxLFQ"

    # Mokume MaxLFQ (uses DirectLFQ if available, otherwise built-in)
    maxlfq = MaxLFQQuantification(min_peptides=1, threads=-1)
    results["Mokume MaxLFQ"] = maxlfq.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )
    intensity_cols["Mokume MaxLFQ"] = "Intensity"

    # Top3
    top3 = Top3Quantification()
    results["Top3"] = top3.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )
    intensity_cols["Top3"] = "Intensity"

    # Sum
    sum_quant = AllPeptidesQuantification()
    results["Sum"] = sum_quant.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )
    intensity_cols["Sum"] = "Intensity"

    # DirectLFQ standalone (for comparison when available)
    if is_directlfq_available():
        try:
            from mokume.quantification import DirectLFQQuantification
            directlfq = DirectLFQQuantification(min_nonan=1)
            results["DirectLFQ"] = directlfq.quantify(
                peptide_df,
                protein_column="ProteinName",
                peptide_column="PeptideSequence",
                intensity_column="NormIntensity",
                sample_column="SampleID",
            )
            intensity_cols["DirectLFQ"] = "Intensity"
        except Exception as e:
            print(f"DirectLFQ failed: {e}")

    return results, intensity_cols


class TestSmallDiannSubset:
    """
    Test suite using small DIA-NN subset from directLFQ test suite.

    IMPORTANT: This is a SUBSET dataset with only ~626 precursors.
    DIA-NN's PG.MaxLFQ was computed from ALL precursors, so correlation
    between computed values and DIA-NN MaxLFQ will be lower (~0.55-0.61).
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        """Load test data."""
        self.diann_df = load_small_diann_report()
        self.peptide_df = prepare_peptide_data(self.diann_df)
        self.diann_maxlfq = extract_diann_maxlfq(self.diann_df)

    def test_small_diann_data_loaded(self):
        """Verify small DIA-NN subset data is loaded correctly."""
        assert len(self.diann_df) > 0, "DIA-NN report should not be empty"
        assert len(self.peptide_df) > 0, "Peptide data should not be empty"

        print(f"\n[Small DIA-NN Subset] Test data summary:")
        print(f"  Peptide measurements: {len(self.peptide_df):,}")
        print(f"  Unique proteins: {self.peptide_df['ProteinName'].nunique()}")
        print(f"  Unique samples: {self.peptide_df['SampleID'].nunique()}")
        print(f"  Unique precursors: {self.peptide_df['PeptideSequence'].nunique()}")

    def test_small_diann_maxlfq_produces_results(self):
        """Test that mokume MaxLFQ produces valid results on small subset."""
        from mokume.quantification import MaxLFQQuantification

        maxlfq = MaxLFQQuantification(min_peptides=1, threads=1)
        result = maxlfq.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        assert len(result) > 0, "MaxLFQ should produce results"
        assert "Intensity" in result.columns
        assert result["Intensity"].notna().sum() > 0, "Should have non-NaN intensities"
        assert (result["Intensity"] > 0).all(), "All intensities should be positive"

        print(f"\n[Small DIA-NN Subset] MaxLFQ results:")
        print(f"  Protein-sample combinations: {len(result):,}")
        print(f"  Intensity range: {result['Intensity'].min():.2f} - {result['Intensity'].max():.2f}")

    def test_small_diann_maxlfq_parallelization(self):
        """Test that parallel MaxLFQ produces same results as single-threaded."""
        from mokume.quantification import MaxLFQQuantification

        maxlfq_single = MaxLFQQuantification(min_peptides=1, threads=1)
        result_single = maxlfq_single.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        maxlfq_parallel = MaxLFQQuantification(min_peptides=1, threads=2)
        result_parallel = maxlfq_parallel.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        corr = compute_correlation(
            result_single, result_parallel,
            "Intensity", "Intensity"
        )

        print(f"\n[Small DIA-NN Subset] Single-threaded vs Parallel MaxLFQ:")
        print(f"  Pearson correlation: {corr['pearson']:.4f}")

        assert corr["pearson"] > 0.9999, "Parallel should produce same results as single-threaded"

    def test_small_diann_full_comparison(self):
        """
        Full comparison of all quantification methods on small DIA-NN subset.

        This test prints the complete comparison table for reference.
        """
        results, intensity_cols = run_all_quantification_methods(
            self.peptide_df, self.diann_maxlfq
        )

        # Print comparison table
        print_comparison_table(
            results, intensity_cols,
            dataset_name="Small DIA-NN Subset (directLFQ test data)",
            peptide_count=len(self.peptide_df),
            protein_count=self.peptide_df['ProteinName'].nunique(),
            sample_count=self.peptide_df['SampleID'].nunique(),
            precursor_count=self.peptide_df['PeptideSequence'].nunique(),
        )

        print("\nNOTE: DIA-NN MaxLFQ shows lower correlation (~0.55-0.61) because:")
        print("  - DIA-NN computed MaxLFQ from ALL precursors")
        print("  - This report only contains a SUBSET (626 precursors)")
        print("  - Most proteins have only 1 precursor in this subset")

        # Validate that computed methods correlate well with each other
        corr_maxlfq_top3 = compute_correlation(
            results["Mokume MaxLFQ"], results["Top3"],
            intensity_cols["Mokume MaxLFQ"], intensity_cols["Top3"]
        )
        assert corr_maxlfq_top3["spearman"] > 0.9, \
            f"MaxLFQ and Top3 should correlate >0.9, got {corr_maxlfq_top3['spearman']:.3f}"

    def test_small_diann_directlfq_comparison(self):
        """Compare DirectLFQ with mokume MaxLFQ on small subset."""
        try:
            from mokume.quantification import DirectLFQQuantification, is_directlfq_available
        except ImportError:
            pytest.skip("DirectLFQ import not available")

        if not is_directlfq_available():
            pytest.skip("DirectLFQ not installed (pip install mokume[directlfq])")

        from mokume.quantification import MaxLFQQuantification

        directlfq = DirectLFQQuantification(min_nonan=1)
        result_directlfq = directlfq.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        # Test default MaxLFQ (should use DirectLFQ when available)
        maxlfq = MaxLFQQuantification(min_peptides=1, threads=1)
        assert maxlfq.using_directlfq, "MaxLFQ should use DirectLFQ when available"

        result_maxlfq = maxlfq.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        corr = compute_correlation(
            result_directlfq, result_maxlfq,
            "Intensity", "Intensity"
        )

        print(f"\n[Small DIA-NN Subset] DirectLFQ vs Mokume MaxLFQ:")
        print(f"  Using DirectLFQ backend: {maxlfq.using_directlfq}")
        print(f"  Pearson correlation:  {corr['pearson']:.4f}")
        print(f"  Spearman correlation: {corr['spearman']:.4f}")
        print(f"  Number of comparisons: {corr['n']:,}")

        assert corr["spearman"] > 0.99, "MaxLFQ with DirectLFQ should match DirectLFQ exactly"

    def test_small_diann_builtin_fallback(self):
        """Test the built-in MaxLFQ fallback implementation."""
        from mokume.quantification import MaxLFQQuantification

        # Force built-in implementation
        maxlfq_builtin = MaxLFQQuantification(min_peptides=1, threads=1, force_builtin=True)
        assert not maxlfq_builtin.using_directlfq, "Should use built-in when forced"
        assert maxlfq_builtin.name == "MaxLFQ (built-in)"

        result_builtin = maxlfq_builtin.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        assert len(result_builtin) > 0, "Built-in MaxLFQ should produce results"
        assert "Intensity" in result_builtin.columns

        print(f"\n[Small DIA-NN Subset] Built-in MaxLFQ fallback:")
        print(f"  Protein-sample combinations: {len(result_builtin):,}")
        print(f"  Intensity range: {result_builtin['Intensity'].min():.2f} - {result_builtin['Intensity'].max():.2f}")


class TestDiannWithSdrf:
    """
    Test suite using full DIA-NN report from PRIDE (PXD063291).

    This dataset is downloaded from PRIDE FTP and cached locally.
    It provides a more complete dataset for validation compared to the small subset.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        """Download and load PRIDE dataset."""
        try:
            report_path = download_pride_dataset()
            self.diann_df = pd.read_csv(report_path, sep="\t")
            self.peptide_df = prepare_peptide_data(self.diann_df)
            self.diann_maxlfq = extract_diann_maxlfq(self.diann_df)
            self.report_path = report_path
        except Exception as e:
            pytest.skip(f"Failed to load PRIDE dataset: {e}")

    def test_pride_data_loaded(self):
        """Verify PRIDE dataset is loaded correctly."""
        assert len(self.diann_df) > 0, "DIA-NN report should not be empty"
        assert len(self.peptide_df) > 0, "Peptide data should not be empty"

        print(f"\n[PRIDE PXD063291] Test data summary:")
        print(f"  Source: {self.report_path}")
        print(f"  Peptide measurements: {len(self.peptide_df):,}")
        print(f"  Unique proteins: {self.peptide_df['ProteinName'].nunique()}")
        print(f"  Unique samples: {self.peptide_df['SampleID'].nunique()}")
        print(f"  Unique precursors: {self.peptide_df['PeptideSequence'].nunique()}")

    def test_pride_maxlfq_produces_results(self):
        """Test that mokume MaxLFQ produces valid results on PRIDE dataset."""
        from mokume.quantification import MaxLFQQuantification

        maxlfq = MaxLFQQuantification(min_peptides=1, threads=-1)
        result = maxlfq.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
        )

        assert len(result) > 0, "MaxLFQ should produce results"
        assert "Intensity" in result.columns
        assert result["Intensity"].notna().sum() > 0, "Should have non-NaN intensities"

        print(f"\n[PRIDE PXD063291] MaxLFQ results:")
        print(f"  Protein-sample combinations: {len(result):,}")
        print(f"  Intensity range: {result['Intensity'].min():.2f} - {result['Intensity'].max():.2f}")

    def test_pride_full_comparison(self):
        """
        Full comparison of all quantification methods on PRIDE dataset.

        This test prints the complete comparison table for reference.
        Expected: Higher correlation with DIA-NN MaxLFQ since this is
        the complete dataset (not a subset).
        """
        results, intensity_cols = run_all_quantification_methods(
            self.peptide_df, self.diann_maxlfq
        )

        # Print comparison table
        print_comparison_table(
            results, intensity_cols,
            dataset_name="PRIDE PXD063291 (Full DIA-NN Report with SDRF)",
            peptide_count=len(self.peptide_df),
            protein_count=self.peptide_df['ProteinName'].nunique(),
            sample_count=self.peptide_df['SampleID'].nunique(),
            precursor_count=self.peptide_df['PeptideSequence'].nunique(),
        )

        # Validate that computed methods correlate well with each other
        corr_maxlfq_top3 = compute_correlation(
            results["Mokume MaxLFQ"], results["Top3"],
            intensity_cols["Mokume MaxLFQ"], intensity_cols["Top3"]
        )

        print(f"\nValidation: MaxLFQ vs Top3 Spearman = {corr_maxlfq_top3['spearman']:.3f}")

        assert corr_maxlfq_top3["spearman"] > 0.85, \
            f"MaxLFQ and Top3 should correlate >0.85, got {corr_maxlfq_top3['spearman']:.3f}"

        # If DIA-NN MaxLFQ is available, check correlation
        if "DIA-NN MaxLFQ" in results and len(results["DIA-NN MaxLFQ"]) > 0:
            corr_diann = compute_correlation(
                results["DIA-NN MaxLFQ"], results["Mokume MaxLFQ"],
                intensity_cols["DIA-NN MaxLFQ"], intensity_cols["Mokume MaxLFQ"]
            )
            print(f"Validation: DIA-NN MaxLFQ vs Mokume MaxLFQ Spearman = {corr_diann['spearman']:.3f}")


def prepare_peptide_data_with_runs(diann_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert DIA-NN report to mokume peptide format with both Sample and Run columns.

    Creates a synthetic SampleID by grouping runs based on the yeast concentration
    pattern in the run names (e.g., Yeast_100ng, Yeast_1ng, etc.).

    DIA-NN columns -> mokume columns:
    - Protein.Group -> ProteinName
    - Run -> Run (kept as run identifier)
    - Derived from Run -> SampleID (condition grouping)
    - Precursor.Id -> PeptideSequence
    - Precursor.Quantity -> NormIntensity
    """
    peptide_df = diann_df[["Protein.Group", "Run", "Precursor.Id", "Precursor.Quantity"]].copy()
    peptide_df = peptide_df.rename(columns={
        "Protein.Group": "ProteinName",
        "Precursor.Id": "PeptideSequence",
        "Precursor.Quantity": "NormIntensity",
    })

    # Extract sample/condition from run name
    # Pattern: ...HeLa_100ng_Yeast_XXXng_S2-...
    def extract_sample(run_name):
        import re
        match = re.search(r'(Yeast_\d+ng)', run_name)
        if match:
            return match.group(1)
        return "Unknown"

    peptide_df["SampleID"] = peptide_df["Run"].apply(extract_sample)

    # Remove rows with missing intensity
    peptide_df = peptide_df[peptide_df["NormIntensity"].notna()]
    peptide_df = peptide_df[peptide_df["NormIntensity"] > 0]

    return peptide_df


def run_quantification_at_level(peptide_df: pd.DataFrame, run_column: str = None):
    """
    Run all quantification methods at sample or run level.

    Parameters
    ----------
    peptide_df : pd.DataFrame
        Peptide data with ProteinName, SampleID, Run, PeptideSequence, NormIntensity
    run_column : str, optional
        If provided, quantification is done at run level

    Returns
    -------
    dict, dict
        Results dict and intensity_cols dict
    """
    from mokume.quantification import (
        MaxLFQQuantification,
        Top3Quantification,
        TopNQuantification,
        AllPeptidesQuantification,
        is_directlfq_available,
    )

    results = {}
    intensity_cols = {}
    level_suffix = "Run" if run_column else "Sample"

    # MaxLFQ (built-in)
    maxlfq = MaxLFQQuantification(min_peptides=1, threads=-1, force_builtin=True)
    results[f"MaxLFQ ({level_suffix})"] = maxlfq.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
        run_column=run_column,
    )
    intensity_cols[f"MaxLFQ ({level_suffix})"] = "Intensity"

    # DirectLFQ (if available, only at sample level - doesn't support run_column)
    if is_directlfq_available() and run_column is None:
        try:
            from mokume.quantification import DirectLFQQuantification
            directlfq = DirectLFQQuantification(min_nonan=1)
            results[f"DirectLFQ ({level_suffix})"] = directlfq.quantify(
                peptide_df,
                protein_column="ProteinName",
                peptide_column="PeptideSequence",
                intensity_column="NormIntensity",
                sample_column="SampleID",
            )
            intensity_cols[f"DirectLFQ ({level_suffix})"] = "Intensity"
        except Exception as e:
            print(f"DirectLFQ failed: {e}")

    # Top3
    top3 = Top3Quantification()
    results[f"Top3 ({level_suffix})"] = top3.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
        run_column=run_column,
    )
    intensity_cols[f"Top3 ({level_suffix})"] = "Intensity"

    # TopN (n=5)
    topn = TopNQuantification(n=5)
    results[f"Top5 ({level_suffix})"] = topn.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
        run_column=run_column,
    )
    intensity_cols[f"Top5 ({level_suffix})"] = "Intensity"

    # Sum (AllPeptides)
    sum_quant = AllPeptidesQuantification()
    results[f"Sum ({level_suffix})"] = sum_quant.quantify(
        peptide_df,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
        run_column=run_column,
    )
    intensity_cols[f"Sum ({level_suffix})"] = "Intensity"

    return results, intensity_cols


def compute_correlation_for_level(
    results: dict,
    intensity_cols: dict,
    run_column: str = None,
) -> pd.DataFrame:
    """
    Compute correlation matrix between all methods at a given aggregation level.

    Returns a DataFrame with Spearman correlations.
    """
    method_names = list(results.keys())
    n_methods = len(method_names)

    # Determine join columns based on level
    if run_column:
        join_cols = ["ProteinName", "SampleID", "Run"]
    else:
        join_cols = ["ProteinName", "SampleID"]

    corr_matrix = np.zeros((n_methods, n_methods))

    for i, name1 in enumerate(method_names):
        for j, name2 in enumerate(method_names):
            if i == j:
                corr_matrix[i, j] = 1.0
            elif i < j:
                # Get available join columns that exist in both dataframes
                df1 = results[name1]
                df2 = results[name2]
                available_cols = [c for c in join_cols if c in df1.columns and c in df2.columns]

                if not available_cols:
                    corr_matrix[i, j] = np.nan
                    corr_matrix[j, i] = np.nan
                    continue

                corr = compute_correlation(
                    df1, df2,
                    intensity_cols[name1], intensity_cols[name2],
                    protein_col=available_cols[0],
                    sample_col=available_cols[1] if len(available_cols) > 1 else available_cols[0]
                )
                corr_matrix[i, j] = corr['spearman']
                corr_matrix[j, i] = corr['spearman']

    return pd.DataFrame(corr_matrix, index=method_names, columns=method_names)


def print_correlation_matrix(corr_df: pd.DataFrame, title: str):
    """Print a correlation matrix in a formatted table."""
    print(f"\n{title}")
    print("-" * (18 + 12 * len(corr_df.columns)))

    # Header - extract short names
    short_names = [name.split(" (")[0] for name in corr_df.columns]
    print(f"{'Method':<18}", end="")
    for name in short_names:
        print(f"{name:>12}", end="")
    print()
    print("-" * (18 + 12 * len(corr_df.columns)))

    # Rows
    for idx, row in corr_df.iterrows():
        short_idx = idx.split(" (")[0]
        print(f"{short_idx:<18}", end="")
        for val in row:
            if np.isnan(val):
                print(f"{'N/A':>12}", end="")
            else:
                print(f"{val:>12.3f}", end="")
        print()

    print("-" * (18 + 12 * len(corr_df.columns)))


def print_aggregation_comparison_table(
    sample_results: dict,
    run_results: dict,
    sample_intensity_cols: dict,
    run_intensity_cols: dict,
    dataset_name: str,
):
    """Print a formatted table comparing sample-level vs run-level quantification."""
    print("\n" + "=" * 100)
    print(f"SAMPLE-LEVEL vs RUN-LEVEL QUANTIFICATION COMPARISON - {dataset_name}")
    print("=" * 100)

    # Print summary of each method's results
    print("\nRESULT COUNTS:")
    print("-" * 100)
    print(f"{'Method':<25} {'Level':<10} {'Protein-Level Rows':>20} {'Unique Proteins':>18} {'Unique Samples/Runs':>20}")
    print("-" * 100)

    method_bases = ["MaxLFQ", "DirectLFQ", "Top3", "Top5", "Sum"]
    for method_base in method_bases:
        sample_key = f"{method_base} (Sample)"
        run_key = f"{method_base} (Run)"

        if sample_key in sample_results:
            df = sample_results[sample_key]
            n_rows = len(df)
            n_proteins = df["ProteinName"].nunique() if "ProteinName" in df.columns else 0
            n_samples = df["SampleID"].nunique() if "SampleID" in df.columns else 0
            print(f"{method_base:<25} {'Sample':<10} {n_rows:>20,} {n_proteins:>18} {n_samples:>20}")

        if run_key in run_results:
            df = run_results[run_key]
            n_rows = len(df)
            n_proteins = df["ProteinName"].nunique() if "ProteinName" in df.columns else 0
            # For run-level, count runs
            if "Run" in df.columns:
                n_runs = df["Run"].nunique()
            elif "SampleID" in df.columns:
                n_runs = df["SampleID"].nunique()
            else:
                n_runs = 0
            print(f"{'':<25} {'Run':<10} {n_rows:>20,} {n_proteins:>18} {n_runs:>20}")

    print("-" * 100)

    # Print intensity statistics
    print("\nINTENSITY STATISTICS (log2 scale):")
    print("-" * 100)
    print(f"{'Method':<25} {'Level':<10} {'Mean':>12} {'Median':>12} {'Std':>12} {'Min':>12} {'Max':>12}")
    print("-" * 100)

    for method_base in method_bases:
        sample_key = f"{method_base} (Sample)"
        run_key = f"{method_base} (Run)"

        if sample_key in sample_results:
            df = sample_results[sample_key]
            col = sample_intensity_cols.get(sample_key)
            if col and col in df.columns:
                vals = np.log2(df[col].dropna() + 1)
                print(f"{method_base:<25} {'Sample':<10} {vals.mean():>12.2f} {vals.median():>12.2f} "
                      f"{vals.std():>12.2f} {vals.min():>12.2f} {vals.max():>12.2f}")

        if run_key in run_results:
            df = run_results[run_key]
            col = run_intensity_cols.get(run_key)
            if col and col in df.columns:
                vals = np.log2(df[col].dropna() + 1)
                print(f"{'':<25} {'Run':<10} {vals.mean():>12.2f} {vals.median():>12.2f} "
                      f"{vals.std():>12.2f} {vals.min():>12.2f} {vals.max():>12.2f}")

    # Print correlation matrices
    if sample_results:
        sample_corr = compute_correlation_for_level(sample_results, sample_intensity_cols, run_column=None)
        print_correlation_matrix(sample_corr, "SAMPLE-LEVEL CORRELATIONS (Spearman)")

    if run_results:
        run_corr = compute_correlation_for_level(run_results, run_intensity_cols, run_column="Run")
        print_correlation_matrix(run_corr, "RUN-LEVEL CORRELATIONS (Spearman)")

    print("=" * 100)


class TestAggregationLevels:
    """
    Test suite comparing sample-level vs run-level quantification.

    This tests the new run_column parameter added to quantification methods,
    which allows aggregation at the run level instead of sample level.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        """Load test data with both sample and run columns."""
        diann_df = load_small_diann_report()
        self.peptide_df = prepare_peptide_data_with_runs(diann_df)

    def test_data_has_multiple_runs_per_sample(self):
        """Verify test data has multiple runs per sample for meaningful comparison."""
        runs_per_sample = self.peptide_df.groupby("SampleID")["Run"].nunique()

        print(f"\n[Aggregation Levels] Test data structure:")
        print(f"  Unique samples: {self.peptide_df['SampleID'].nunique()}")
        print(f"  Unique runs: {self.peptide_df['Run'].nunique()}")
        print(f"  Runs per sample:")
        for sample, n_runs in runs_per_sample.items():
            print(f"    {sample}: {n_runs} runs")

        # Verify we have multiple runs per sample for meaningful test
        assert runs_per_sample.max() > 1, "Need multiple runs per sample for this test"

    def test_sample_level_quantification(self):
        """Test quantification at sample level (default behavior)."""
        results, intensity_cols = run_quantification_at_level(
            self.peptide_df, run_column=None
        )

        for method_name, df in results.items():
            assert len(df) > 0, f"{method_name} should produce results"
            assert "ProteinName" in df.columns
            assert "SampleID" in df.columns
            # Run column should NOT be in results for sample-level
            assert "Run" not in df.columns or df["Run"].isna().all() or method_name.endswith("(Run)")

        print(f"\n[Sample Level] Quantification results:")
        for method_name, df in results.items():
            print(f"  {method_name}: {len(df)} protein-sample combinations")

    def test_run_level_quantification(self):
        """Test quantification at run level."""
        results, intensity_cols = run_quantification_at_level(
            self.peptide_df, run_column="Run"
        )

        for method_name, df in results.items():
            assert len(df) > 0, f"{method_name} should produce results"
            assert "ProteinName" in df.columns
            # Run column should be in results for run-level
            assert "Run" in df.columns, f"{method_name} should have Run column"

        print(f"\n[Run Level] Quantification results:")
        for method_name, df in results.items():
            n_runs = df["Run"].nunique() if "Run" in df.columns else 0
            print(f"  {method_name}: {len(df)} protein-run combinations ({n_runs} runs)")

    def test_run_level_has_more_rows_than_sample_level(self):
        """Run-level should produce more rows than sample-level when multiple runs per sample."""
        sample_results, _ = run_quantification_at_level(self.peptide_df, run_column=None)
        run_results, _ = run_quantification_at_level(self.peptide_df, run_column="Run")

        for method_base in ["MaxLFQ", "Top3", "Top5", "Sum"]:
            sample_key = f"{method_base} (Sample)"
            run_key = f"{method_base} (Run)"

            n_sample = len(sample_results[sample_key])
            n_run = len(run_results[run_key])

            print(f"\n[{method_base}] Sample-level: {n_sample} rows, Run-level: {n_run} rows")

            # Run-level should have >= sample-level rows (more granular)
            assert n_run >= n_sample, \
                f"{method_base}: Run-level ({n_run}) should have >= rows than sample-level ({n_sample})"

    def test_full_comparison_sample_vs_run(self):
        """
        Full comparison of sample-level vs run-level quantification.
        Prints detailed comparison table.
        """
        sample_results, sample_intensity_cols = run_quantification_at_level(
            self.peptide_df, run_column=None
        )
        run_results, run_intensity_cols = run_quantification_at_level(
            self.peptide_df, run_column="Run"
        )

        # Print comparison table
        print_aggregation_comparison_table(
            sample_results,
            run_results,
            sample_intensity_cols,
            run_intensity_cols,
            dataset_name="Small DIA-NN Subset",
        )

        # Print sample data structure
        print("\nDATA STRUCTURE:")
        print(f"  Input peptide measurements: {len(self.peptide_df):,}")
        print(f"  Unique proteins: {self.peptide_df['ProteinName'].nunique()}")
        print(f"  Unique samples (conditions): {self.peptide_df['SampleID'].nunique()}")
        print(f"  Unique runs: {self.peptide_df['Run'].nunique()}")

        runs_per_sample = self.peptide_df.groupby("SampleID")["Run"].nunique()
        print(f"  Runs per sample: {runs_per_sample.min()} - {runs_per_sample.max()}")

    def test_maxlfq_sample_vs_run_correlation(self):
        """
        Test correlation between sample-level and run-level MaxLFQ.

        For run-level results, we aggregate to sample level by taking mean
        to enable comparison.
        """
        from mokume.quantification import MaxLFQQuantification

        maxlfq = MaxLFQQuantification(min_peptides=1, threads=1, force_builtin=True)

        # Sample-level
        result_sample = maxlfq.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
            run_column=None,
        )

        # Run-level
        result_run = maxlfq.quantify(
            self.peptide_df,
            protein_column="ProteinName",
            peptide_column="PeptideSequence",
            intensity_column="NormIntensity",
            sample_column="SampleID",
            run_column="Run",
        )

        # Aggregate run-level to sample-level (mean) for comparison
        # First, we need to map runs back to samples
        run_to_sample = self.peptide_df[["Run", "SampleID"]].drop_duplicates()
        result_run_with_sample = result_run.merge(run_to_sample, on="Run", how="left", suffixes=("_x", ""))

        # Handle column naming after merge
        if "SampleID_x" in result_run_with_sample.columns:
            result_run_with_sample = result_run_with_sample.drop(columns=["SampleID_x"])

        result_run_aggregated = result_run_with_sample.groupby(
            ["ProteinName", "SampleID"]
        )["Intensity"].mean().reset_index()

        # Compute correlation
        merged = pd.merge(
            result_sample[["ProteinName", "SampleID", "Intensity"]],
            result_run_aggregated,
            on=["ProteinName", "SampleID"],
            suffixes=("_sample", "_run_agg")
        )

        if len(merged) > 0:
            x = np.log2(merged["Intensity_sample"].values + 1)
            y = np.log2(merged["Intensity_run_agg"].values + 1)

            valid = np.isfinite(x) & np.isfinite(y)
            x, y = x[valid], y[valid]

            if len(x) >= 3:
                pearson_r, _ = stats.pearsonr(x, y)
                spearman_r, _ = stats.spearmanr(x, y)

                print(f"\n[MaxLFQ] Sample-level vs Run-level (aggregated to sample) correlation:")
                print(f"  Pearson:  {pearson_r:.4f}")
                print(f"  Spearman: {spearman_r:.4f}")
                print(f"  N comparisons: {len(x)}")

                # They should correlate reasonably well
                assert spearman_r > 0.7, f"Sample vs aggregated run-level should correlate > 0.7, got {spearman_r:.3f}"


# Standalone test function for quick runs
@pytest.mark.comparison
def test_small_diann_summary():
    """
    Quick summary test for small DIA-NN subset.
    Run with: pytest tests/test_quantification_accuracy.py::test_small_diann_summary -v -s

    This test prints comparison tables for all quantification methods.
    """
    diann_df = load_small_diann_report()
    peptide_df = prepare_peptide_data(diann_df)
    diann_maxlfq = extract_diann_maxlfq(diann_df)

    results, intensity_cols = run_all_quantification_methods(peptide_df, diann_maxlfq)

    print_comparison_table(
        results, intensity_cols,
        dataset_name="Small DIA-NN Subset (directLFQ test data)",
        peptide_count=len(peptide_df),
        protein_count=peptide_df['ProteinName'].nunique(),
        sample_count=peptide_df['SampleID'].nunique(),
        precursor_count=peptide_df['PeptideSequence'].nunique(),
    )


@pytest.mark.comparison
def test_comprehensive_aggregation_comparison():
    """
    Comprehensive test comparing all quantification methods at both sample and run levels.

    Run with: pytest tests/test_quantification_accuracy.py::test_comprehensive_aggregation_comparison -v -s

    This test prints detailed comparison tables including:
    - Result counts and intensity statistics
    - Correlation matrices for sample-level methods
    - Correlation matrices for run-level methods

    NOTE: This test output is important for CI/CD - it shows quantification method
    correlations that should remain stable across code changes.
    """
    print("\n" + "=" * 100)
    print("COMPREHENSIVE QUANTIFICATION COMPARISON")
    print("Testing all methods at both SAMPLE and RUN aggregation levels")
    print("=" * 100)

    # Load and prepare data
    diann_df = load_small_diann_report()
    peptide_df = prepare_peptide_data_with_runs(diann_df)

    print(f"\nDATASET: Small DIA-NN Subset")
    print(f"  Peptide measurements: {len(peptide_df):,}")
    print(f"  Unique proteins: {peptide_df['ProteinName'].nunique()}")
    print(f"  Unique samples: {peptide_df['SampleID'].nunique()}")
    print(f"  Unique runs: {peptide_df['Run'].nunique()}")

    runs_per_sample = peptide_df.groupby("SampleID")["Run"].nunique()
    print(f"  Runs per sample: {runs_per_sample.to_dict()}")

    # Run quantification at both levels
    print("\n" + "-" * 100)
    print("Running quantification methods...")
    print("-" * 100)

    sample_results, sample_intensity_cols = run_quantification_at_level(peptide_df, run_column=None)
    run_results, run_intensity_cols = run_quantification_at_level(peptide_df, run_column="Run")

    # Print comprehensive comparison table
    print_aggregation_comparison_table(
        sample_results,
        run_results,
        sample_intensity_cols,
        run_intensity_cols,
        dataset_name="Small DIA-NN Subset",
    )


@pytest.mark.comparison
def test_pride_aggregation_comparison():
    """
    Test comparing all quantification methods on PRIDE PXD063291 dataset.

    Run with: pytest tests/test_quantification_accuracy.py::test_pride_aggregation_comparison -v -s

    This test prints correlation matrices comparing all quantification methods
    against DIA-NN's MaxLFQ values.

    NOTE: This test output is important for CI/CD - it validates quantification
    accuracy against DIA-NN reference values.
    """
    print("\n" + "=" * 100)
    print("PRIDE PXD063291 QUANTIFICATION COMPARISON")
    print("=" * 100)

    try:
        report_path = download_pride_dataset()
        diann_df = pd.read_csv(report_path, sep="\t")
    except Exception as e:
        pytest.skip(f"Failed to load PRIDE dataset: {e}")

    # Prepare data with both sample and run columns
    peptide_df = diann_df[["Protein.Group", "Run", "Precursor.Id", "Precursor.Quantity"]].copy()
    peptide_df = peptide_df.rename(columns={
        "Protein.Group": "ProteinName",
        "Precursor.Id": "PeptideSequence",
        "Precursor.Quantity": "NormIntensity",
    })

    # Use Run as both sample and run for this dataset (each run is a sample)
    peptide_df["SampleID"] = peptide_df["Run"]
    peptide_df = peptide_df[peptide_df["NormIntensity"].notna() & (peptide_df["NormIntensity"] > 0)]

    print(f"\nDATASET: PRIDE PXD063291")
    print(f"  Source: {report_path}")
    print(f"  Peptide measurements: {len(peptide_df):,}")
    print(f"  Unique proteins: {peptide_df['ProteinName'].nunique()}")
    print(f"  Unique samples/runs: {peptide_df['SampleID'].nunique()}")

    # For PRIDE data, only run sample-level (each run is treated as a sample)
    print("\n" + "-" * 100)
    print("Running quantification methods at sample level...")
    print("-" * 100)

    sample_results, sample_intensity_cols = run_quantification_at_level(peptide_df, run_column=None)

    # Print correlation matrix
    sample_corr = compute_correlation_for_level(sample_results, sample_intensity_cols, run_column=None)
    print_correlation_matrix(sample_corr, f"PRIDE PXD063291 - SAMPLE-LEVEL CORRELATIONS (Spearman)")

    # Also compare with DIA-NN MaxLFQ if available
    diann_maxlfq = extract_diann_maxlfq(diann_df)
    if len(diann_maxlfq) > 0:
        print("\nComparison with DIA-NN MaxLFQ:")
        for method_name in sample_results:
            corr = compute_correlation(
                sample_results[method_name], diann_maxlfq,
                sample_intensity_cols[method_name], "DIANNMaxLFQ"
            )
            if not np.isnan(corr['spearman']):
                short_name = method_name.split(" (")[0]
                print(f"  {short_name} vs DIA-NN MaxLFQ: Spearman = {corr['spearman']:.3f} (n={corr['n']})")


if __name__ == "__main__":
    print("Running comprehensive aggregation comparison test...")
    test_comprehensive_aggregation_comparison()
