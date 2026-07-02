#!/usr/bin/env python3
"""
Grid search over all quantification and normalization method combinations.

For each combination, computes:
- CV within conditions
- Number of proteins quantified
- Missing value rate

Outputs results to results/grid_search_results.csv
"""

import warnings
from itertools import product
from pathlib import Path

import pandas as pd
import numpy as np
import mokume

from config import (
    DATA_DIR, RESULTS_DIR,
    QUANTIFICATION_METHODS, NORMALIZATION_METHODS, IMPUTATION_METHODS,
    SAMPLE_CONDITIONS, MIN_PEPTIDES,
)
from benchmark_utils import (
    quantile_normalize, median_center, total_intensity_normalize,
    hierarchical_normalize,
    impute_values, pivot_wider,
)

warnings.filterwarnings("ignore")

# Path for temporary files used by mokume.peptides2protein
_TEMP_DIR = Path(__file__).parent.parent / "_temp"

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
    """Apply quantification method to peptide-level data via mokume API."""
    try:
        _TEMP_DIR.mkdir(parents=True, exist_ok=True)

        # Write peptide data to a temporary file for mokume
        input_path = _TEMP_DIR / f"peptides_{method}.parquet"
        output_path = _TEMP_DIR / f"quantified_{method}.parquet"
        df.to_parquet(input_path, index=False)

        mokume.peptides2protein(
            peptides=str(input_path),
            method=method,
            output=str(output_path),
        )

        if output_path.exists():
            result = pd.read_parquet(output_path)
            # Clean up temp files
            input_path.unlink(missing_ok=True)
            output_path.unlink(missing_ok=True)
            return result
        else:
            return None

    except Exception as e:
        print(f"    Error in {method}: {e}")
        return None

def apply_normalization(df_wide: pd.DataFrame, method: str) -> pd.DataFrame:
    """Apply normalization method to wide-format protein data."""
    if method == "none":
        return df_wide

    numeric_df = df_wide.select_dtypes(include=[np.number])

    if method == "median":
        log_df = np.log10(numeric_df.replace(0, np.nan))
        centered = median_center(log_df, axis=0)
        return 10 ** centered

    elif method == "quantile":
        return quantile_normalize(numeric_df)

    elif method == "total_intensity":
        return total_intensity_normalize(numeric_df)

    elif method == "hierarchical":
        try:
            log_df = np.log10(numeric_df.replace(0, np.nan))
            normalized = hierarchical_normalize(log_df)
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
        min_val = numeric_df.min().min()
        if pd.isna(min_val) or min_val <= 0:
            min_val = 1.0
        return numeric_df.fillna(min_val / 2)

    elif method in ("knn", "median", "mean"):
        try:
            return impute_values(numeric_df, method=method)
        except Exception as e:
            print(f"    Imputation failed: {e}")
            return numeric_df

    return df_wide

def compute_cv_within_conditions(df_wide: pd.DataFrame, sample_conditions: dict) -> float:
    """Compute mean CV within conditions."""
    cvs = []

    col_conditions = {}
    for col in df_wide.columns:
        sample_id = extract_sample_id(str(col))
        condition = get_condition(sample_id)
        if condition != "Unknown":
            col_conditions[col] = condition

    condition_cols = {}
    for col, cond in col_conditions.items():
        if cond not in condition_cols:
            condition_cols[cond] = []
        condition_cols[cond].append(col)

    for condition, cols in condition_cols.items():
        if len(cols) < 2:
            continue

        subset = df_wide[cols]
        means = subset.mean(axis=1)
        stds = subset.std(axis=1)

        valid = means > 0
        cv_values = (stds[valid] / means[valid]).dropna()
        cvs.extend(cv_values.tolist())

    return np.median(cvs) if cvs else np.nan

def run_grid_search(technology: str) -> pd.DataFrame:
    """Run grid search for a single technology."""
    print(f"\n{'='*60}")
    print(f"Grid Search: {technology.upper()}")
    print("=" * 60)

    print("\nLoading data...")
    df = load_parquet_data(technology)
    print(f"  Loaded {len(df):,} rows")

    results = []

    quant_methods = QUANTIFICATION_METHODS.copy()
    # Check directlfq availability
    try:
        import directlfq  # noqa: F401
    except ImportError:
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
            intensity_col = None
            for col in protein_df.columns:
                if "intensity" in col.lower():
                    intensity_col = col
                    break

            if intensity_col is None:
                print("    Skipped (no intensity column found)")
                continue

            protein_col = "ProteinName" if "ProteinName" in protein_df.columns else "protein_accessions"
            sample_col = "SampleID" if "SampleID" in protein_df.columns else "sample_accession"

            df_wide = pivot_wider(
                protein_df,
                row_name=protein_col,
                col_name=sample_col,
                values=intensity_col,
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

        output_path = RESULTS_DIR / "grid_search_results.csv"
        combined.to_csv(output_path, index=False)
        print(f"\n\nResults saved to: {output_path}")

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
