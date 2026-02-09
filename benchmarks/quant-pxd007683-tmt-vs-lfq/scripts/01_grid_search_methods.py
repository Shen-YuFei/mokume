#!/usr/bin/env python3
"""
Grid search over all quantification and normalization method combinations.

For each combination, computes:
- CV within conditions
- Number of proteins quantified
- Missing value rate

Outputs results to results/grid_search_results.csv
"""

import sys
from pathlib import Path
from itertools import product
import warnings

import pandas as pd
import numpy as np

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_DIR, RESULTS_DIR,
    QUANTIFICATION_METHODS, NORMALIZATION_METHODS, IMPUTATION_METHODS,
    SAMPLE_CONDITIONS, MIN_PEPTIDES,
)

# Mokume imports
from mokume.quantification import get_quantification_method, is_directlfq_available
from mokume.normalization.protein import quantile_normalize, median_center, total_intensity_normalize
from mokume.normalization.hierarchical import HierarchicalSampleNormalizer
from mokume.imputation import impute_missing_values
from mokume.postprocessing import pivot_wider

warnings.filterwarnings("ignore")


def load_parquet_data(technology: str) -> pd.DataFrame:
    """Load parquet file for a technology (tmt or lfq)."""
    parquet_path = DATA_DIR / f"PXD007683-{technology.upper()}.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(
            f"Data file not found: {parquet_path}\n"
            "Run 00_download_data.py first."
        )
    return pd.read_parquet(parquet_path)


def extract_sample_id(sample_name: str) -> str:
    """Extract sample ID from full sample name."""
    # Handle various naming conventions
    for pattern in ["Sample-", "sample-", "S"]:
        if pattern in sample_name:
            parts = sample_name.split(pattern)
            if len(parts) > 1:
                num = "".join(c for c in parts[-1] if c.isdigit())
                if num:
                    return f"Sample-{num}"
    return sample_name


def get_condition(sample_id: str) -> str:
    """Get condition for a sample ID."""
    return SAMPLE_CONDITIONS.get(sample_id, "Unknown")


def apply_quantification(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply quantification method to peptide-level data."""
    if method == "directlfq" and not is_directlfq_available():
        return None

    try:
        quant_method = get_quantification_method(method, min_peptides=MIN_PEPTIDES)

        # Identify columns
        protein_col = "ProteinName" if "ProteinName" in df.columns else "protein_accessions"
        peptide_col = "PeptideSequence" if "PeptideSequence" in df.columns else "peptide_canonical"
        intensity_col = "Intensity" if "Intensity" in df.columns else "intensity"
        sample_col = "SampleID" if "SampleID" in df.columns else "sample_accession"

        # Check required columns exist
        for col in [protein_col, peptide_col, intensity_col, sample_col]:
            if col not in df.columns:
                print(f"    Warning: Column {col} not found. Available: {list(df.columns)[:10]}")
                return None

        result = quant_method.quantify(
            peptide_df=df,
            protein_column=protein_col,
            peptide_column=peptide_col,
            intensity_column=intensity_col,
            sample_column=sample_col,
        )
        return result
    except Exception as e:
        print(f"    Error in {method}: {e}")
        return None


def apply_normalization(df_wide: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply normalization method to wide-format protein data."""
    if method == "none":
        return df_wide

    # Work with numeric columns only
    numeric_df = df_wide.select_dtypes(include=[np.number])

    if method == "median":
        # Log transform, median center, back transform
        log_df = np.log10(numeric_df.replace(0, np.nan))
        centered = median_center(log_df, axis=0)
        return 10 ** centered

    elif method == "quantile":
        return quantile_normalize(numeric_df)

    elif method == "total_intensity":
        return total_intensity_normalize(numeric_df)

    elif method == "hierarchical":
        try:
            normalizer = HierarchicalSampleNormalizer()
            # Hierarchical expects log-transformed data
            log_df = np.log10(numeric_df.replace(0, np.nan))
            normalized = normalizer.fit_transform(log_df)
            return 10 ** normalized
        except Exception as e:
            print(f"    Hierarchical normalization failed: {e}")
            return numeric_df

    return df_wide


def apply_imputation(df_wide: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply imputation method to wide-format protein data."""
    if method == "none":
        return df_wide

    numeric_df = df_wide.select_dtypes(include=[np.number])

    if method == "min":
        # Replace NaN with min/2 (MNAR assumption)
        min_val = numeric_df.min().min()
        if pd.isna(min_val) or min_val <= 0:
            min_val = 1.0
        return numeric_df.fillna(min_val / 2)

    elif method in ["knn", "median", "mean"]:
        try:
            sklearn_method = "median" if method == "median" else method
            return impute_missing_values(numeric_df, method=sklearn_method)
        except Exception as e:
            print(f"    Imputation failed: {e}")
            return numeric_df

    return df_wide


def compute_cv_within_conditions(df_wide: pd.DataFrame, sample_conditions: dict) -> float:
    """Compute mean CV within conditions."""
    cvs = []

    # Map columns to conditions
    col_conditions = {}
    for col in df_wide.columns:
        sample_id = extract_sample_id(str(col))
        condition = get_condition(sample_id)
        if condition != "Unknown":
            col_conditions[col] = condition

    # Group columns by condition
    condition_cols = {}
    for col, cond in col_conditions.items():
        if cond not in condition_cols:
            condition_cols[cond] = []
        condition_cols[cond].append(col)

    # Compute CV for each protein in each condition
    for condition, cols in condition_cols.items():
        if len(cols) < 2:
            continue

        subset = df_wide[cols]
        means = subset.mean(axis=1)
        stds = subset.std(axis=1)

        # CV = std / mean, only for proteins with non-zero mean
        valid = means > 0
        cv_values = (stds[valid] / means[valid]).dropna()
        cvs.extend(cv_values.tolist())

    return np.median(cvs) if cvs else np.nan


def run_grid_search(technology: str) -> pd.DataFrame:
    """Run grid search for a single technology."""
    print(f"\n{'='*60}")
    print(f"Grid Search: {technology.upper()}")
    print("=" * 60)

    # Load data
    print("\nLoading data...")
    df = load_parquet_data(technology)
    print(f"  Loaded {len(df):,} rows")

    results = []

    # Get available methods
    quant_methods = QUANTIFICATION_METHODS.copy()
    if not is_directlfq_available():
        quant_methods = [m for m in quant_methods if m != "directlfq"]
        print("  Note: DirectLFQ not available, skipping")

    total_combinations = len(quant_methods) * len(NORMALIZATION_METHODS) * len(IMPUTATION_METHODS)
    print(f"\nTesting {total_combinations} combinations...")

    combo_num = 0
    for quant_method, norm_method, impute_method in product(
        quant_methods, NORMALIZATION_METHODS, IMPUTATION_METHODS
    ):
        combo_num += 1
        combo_name = f"{quant_method}_{norm_method}_{impute_method}"
        print(f"\n  [{combo_num}/{total_combinations}] {combo_name}")

        try:
            # Step 1: Quantification
            protein_df = apply_quantification(df, quant_method)
            if protein_df is None:
                print("    Skipped (quantification failed)")
                continue

            # Step 2: Pivot to wide format
            # Identify intensity column from quantification output
            intensity_col = None
            for col in ["intensity", "Intensity", "NormIntensity", f"Top{quant_method[-1]}Intensity" if quant_method.startswith("top") else ""]:
                if col in protein_df.columns:
                    intensity_col = col
                    break

            if intensity_col is None:
                # Try to find any intensity-like column
                for col in protein_df.columns:
                    if "intensity" in col.lower():
                        intensity_col = col
                        break

            if intensity_col is None:
                print(f"    Skipped (no intensity column found)")
                continue

            protein_col = "ProteinName" if "ProteinName" in protein_df.columns else "protein_accessions"
            sample_col = "SampleID" if "SampleID" in protein_df.columns else "sample_accession"

            df_wide = pivot_wider(
                protein_df,
                row_name=protein_col,
                col_name=sample_col,
                values=intensity_col,
                fillna=False,
            )

            # Step 3: Normalization
            df_norm = apply_normalization(df_wide, norm_method)

            # Step 4: Imputation
            df_final = apply_imputation(df_norm, impute_method)

            # Compute metrics
            n_proteins = len(df_final)
            missing_rate = df_final.isna().sum().sum() / df_final.size if df_final.size > 0 else 1.0
            cv_within = compute_cv_within_conditions(df_final, SAMPLE_CONDITIONS)

            results.append({
                "technology": technology,
                "quant_method": quant_method,
                "norm_method": norm_method,
                "impute_method": impute_method,
                "n_proteins": n_proteins,
                "missing_rate": missing_rate,
                "cv_within_condition": cv_within,
            })

            print(f"    Proteins: {n_proteins}, Missing: {missing_rate:.1%}, CV: {cv_within:.3f}")

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    return pd.DataFrame(results)


def main():
    """Run grid search for both technologies."""
    print("=" * 60)
    print("PXD007683 Benchmark: Grid Search")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    for technology in ["tmt", "lfq"]:
        try:
            results = run_grid_search(technology)
            all_results.append(results)
        except FileNotFoundError as e:
            print(f"\nERROR: {e}")
            continue

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)

        # Save results
        output_path = RESULTS_DIR / "grid_search_results.csv"
        combined.to_csv(output_path, index=False)
        print(f"\n\nResults saved to: {output_path}")

        # Print summary
        print("\n" + "=" * 60)
        print("Summary: Top 10 Methods by CV")
        print("=" * 60)

        for tech in ["tmt", "lfq"]:
            tech_results = combined[combined["technology"] == tech]
            if len(tech_results) > 0:
                print(f"\n{tech.upper()}:")
                top10 = tech_results.nsmallest(10, "cv_within_condition")
                for _, row in top10.iterrows():
                    print(f"  {row['quant_method']:10s} + {row['norm_method']:12s} + {row['impute_method']:8s}: "
                          f"CV={row['cv_within_condition']:.3f}, n={row['n_proteins']}")


if __name__ == "__main__":
    main()
