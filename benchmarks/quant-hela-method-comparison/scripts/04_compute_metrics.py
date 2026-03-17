#!/usr/bin/env python3
"""
Phase 4: Compute Evaluation Metrics

Calculate benchmarking metrics for all quantification methods:
- Coefficient of Variation (CV) within experiments
- Cross-experiment correlation
- Expression profile stability (MAD, IQR)
- Rank consistency (Spearman correlation)
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    ALL_DATASETS,
    HELA_DATASETS,
    PROTEIN_QUANT_DIR,
    ANALYSIS_DIR,
    QUANTIFICATION_METHODS,
    PROTEINS_OF_INTEREST,
    DatasetInfo,
)


# Intensity column names for each method
INTENSITY_COLUMNS = {
    "ibaq": "IbaqLog",      # Log-transformed iBAQ (computed by library)
    "ibaq_raw": "Ibaq",     # Raw iBAQ values (not log-transformed)
    "maxlfq": "MaxLFQIntensity",
    "directlfq": "DirectLFQIntensity",
    "top3": "Top3Intensity",
    "topn": "Top10Intensity",
    "top10": "Top10Intensity",
    "top5": "Top5Intensity",
    "sum": "SumIntensity",
}


def load_quantification_result(
    dataset_id: str,
    method: str,
) -> Optional[pd.DataFrame]:
    """Load quantification result for a dataset and method."""
    quant_dir = PROTEIN_QUANT_DIR / dataset_id

    # Handle topN naming
    if method.startswith("top") and method != "top3":
        filename = f"{method}.parquet"
    else:
        filename = f"{method}.parquet"

    path = quant_dir / filename

    if not path.exists():
        return None

    return pd.read_parquet(path)


def compute_cv_per_protein(
    df: pd.DataFrame,
    intensity_col: str,
    protein_col: str = "ProteinName",
    sample_col: str = "SampleID",
) -> pd.DataFrame:
    """
    Compute coefficient of variation per protein across samples.

    CV = std / mean

    Returns DataFrame with columns: ProteinName, CV, Mean, Std, N
    """
    # Pivot to wide format
    wide = df.pivot_table(
        index=protein_col,
        columns=sample_col,
        values=intensity_col,
        aggfunc="first"
    )

    # Compute statistics
    result = pd.DataFrame({
        "ProteinName": wide.index,
        "Mean": wide.mean(axis=1, skipna=True).values,
        "Std": wide.std(axis=1, skipna=True).values,
        "N": wide.notna().sum(axis=1).values,
    })

    # CV = std/mean
    result["CV"] = result["Std"] / result["Mean"]

    # Filter: require at least 3 samples
    result = result[result["N"] >= 3]

    return result


def compute_cv_summary(
    dataset_id: str,
    method: str,
) -> Optional[Dict]:
    """
    Compute CV summary statistics for a dataset and method.

    Returns dict with: mean_cv, median_cv, n_proteins, method, dataset
    """
    df = load_quantification_result(dataset_id, method)
    if df is None:
        return None

    intensity_col = INTENSITY_COLUMNS.get(method)
    if intensity_col is None or intensity_col not in df.columns:
        # Try to find intensity column
        for col in df.columns:
            if "intensity" in col.lower() or "ibaq" in col.lower():
                intensity_col = col
                break

    if intensity_col not in df.columns:
        print(f"  WARNING: No intensity column found for {method}")
        return None

    cv_df = compute_cv_per_protein(df, intensity_col)

    return {
        "dataset": dataset_id,
        "method": method,
        "mean_cv": cv_df["CV"].mean(),
        "median_cv": cv_df["CV"].median(),
        "cv_q25": cv_df["CV"].quantile(0.25),
        "cv_q75": cv_df["CV"].quantile(0.75),
        "n_proteins": len(cv_df),
    }


def compute_cross_experiment_correlation(
    datasets: Dict[str, DatasetInfo],
    method: str,
) -> pd.DataFrame:
    """
    Compute pairwise Pearson correlation between dataset median expressions.

    Returns correlation matrix as DataFrame.
    """
    # Collect median expression per protein per dataset
    medians = {}

    intensity_col = INTENSITY_COLUMNS.get(method)

    for dataset_id, dataset in datasets.items():
        df = load_quantification_result(dataset_id, method)
        if df is None:
            continue

        if intensity_col not in df.columns:
            for col in df.columns:
                if "intensity" in col.lower() or "ibaq" in col.lower():
                    intensity_col = col
                    break

        if intensity_col not in df.columns:
            continue

        # Compute median per protein
        median_expr = df.groupby("ProteinName")[intensity_col].median()
        medians[dataset_id] = median_expr

    if len(medians) < 2:
        return pd.DataFrame()

    # Create expression matrix
    all_proteins = set()
    for m in medians.values():
        all_proteins.update(m.index)

    expr_matrix = pd.DataFrame(index=sorted(all_proteins))
    for dataset_id, median_expr in medians.items():
        expr_matrix[dataset_id] = median_expr

    # Compute correlation matrix
    # Use log-transformed values if not already log
    if "log" not in method.lower() and method != "ibaq":
        expr_matrix = np.log2(expr_matrix + 1)

    corr_matrix = expr_matrix.corr(method="pearson")

    return corr_matrix


def compute_tmt_lfq_correlation(method: str) -> Optional[Dict]:
    """
    Compute correlation between TMT and LFQ quantification (PXD007683).

    Returns correlation statistics and raw values for plotting.
    """
    tmt_df = load_quantification_result("PXD007683-TMT", method)
    lfq_df = load_quantification_result("PXD007683-LFQ", method)

    if tmt_df is None or lfq_df is None:
        return None

    intensity_col = INTENSITY_COLUMNS.get(method)

    # Find intensity columns
    for df_name, df in [("TMT", tmt_df), ("LFQ", lfq_df)]:
        if intensity_col not in df.columns:
            for col in df.columns:
                if "intensity" in col.lower() or "ibaq" in col.lower():
                    intensity_col = col
                    break

    if intensity_col not in tmt_df.columns or intensity_col not in lfq_df.columns:
        return None

    # Get median per protein
    tmt_median = tmt_df.groupby("ProteinName")[intensity_col].median()
    lfq_median = lfq_df.groupby("ProteinName")[intensity_col].median()

    # Merge on common proteins
    common = tmt_median.index.intersection(lfq_median.index)

    if len(common) < 10:
        return None

    # Check if values should be log-transformed
    # iBAQ (IbaqLog) is already log-transformed by the library
    # ibaq_raw, DirectLFQ, Top3, TopN, Sum are raw intensities - need log transformation
    if method == "ibaq":
        # IbaqLog is already log-transformed
        tmt_vals = tmt_median.loc[common].values
        lfq_vals = lfq_median.loc[common].values
    else:
        # All other methods need log transformation
        tmt_vals = np.log2(tmt_median.loc[common].values + 1)
        lfq_vals = np.log2(lfq_median.loc[common].values + 1)

    # Remove NaN
    valid = np.isfinite(tmt_vals) & np.isfinite(lfq_vals)
    proteins = np.array(list(common))[valid]
    tmt_vals = tmt_vals[valid]
    lfq_vals = lfq_vals[valid]

    pearson_r, pearson_p = stats.pearsonr(tmt_vals, lfq_vals)
    spearman_r, spearman_p = stats.spearmanr(tmt_vals, lfq_vals)

    return {
        "method": method,
        "n_proteins": len(tmt_vals),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "proteins": proteins,
        "tmt_values": tmt_vals,
        "lfq_values": lfq_vals,
    }


def compute_rank_consistency(
    datasets: Dict[str, DatasetInfo],
    method: str,
) -> pd.DataFrame:
    """
    Compute rank consistency (Spearman correlation) across datasets.

    Tests if highly/lowly expressed proteins maintain relative ranks.
    """
    # Get median expression per protein per dataset
    medians = {}
    intensity_col = INTENSITY_COLUMNS.get(method)

    for dataset_id, dataset in datasets.items():
        df = load_quantification_result(dataset_id, method)
        if df is None:
            continue

        if intensity_col not in df.columns:
            for col in df.columns:
                if "intensity" in col.lower() or "ibaq" in col.lower():
                    intensity_col = col
                    break

        if intensity_col not in df.columns:
            continue

        median_expr = df.groupby("ProteinName")[intensity_col].median()
        medians[dataset_id] = median_expr.rank()

    if len(medians) < 2:
        return pd.DataFrame()

    # Create rank matrix
    all_proteins = set()
    for m in medians.values():
        all_proteins.update(m.index)

    rank_matrix = pd.DataFrame(index=sorted(all_proteins))
    for dataset_id, ranks in medians.items():
        rank_matrix[dataset_id] = ranks

    # Compute Spearman correlation matrix
    corr_matrix = rank_matrix.corr(method="spearman")

    return corr_matrix


def compute_expression_stability(
    datasets: Dict[str, DatasetInfo],
    method: str,
    proteins: List[str] = None,
) -> pd.DataFrame:
    """
    Compute expression stability metrics (MAD, IQR) for selected proteins.

    Parameters
    ----------
    datasets : dict
        Datasets to analyze
    method : str
        Quantification method
    proteins : list, optional
        Proteins to analyze. Defaults to PROTEINS_OF_INTEREST.

    Returns
    -------
    pd.DataFrame
        Stability metrics per protein: MAD, IQR, range, n_datasets
    """
    if proteins is None:
        proteins = PROTEINS_OF_INTEREST

    intensity_col = INTENSITY_COLUMNS.get(method)

    # Collect expression values per protein
    protein_values = {p: [] for p in proteins}

    for dataset_id, dataset in datasets.items():
        df = load_quantification_result(dataset_id, method)
        if df is None:
            continue

        if intensity_col not in df.columns:
            for col in df.columns:
                if "intensity" in col.lower() or "ibaq" in col.lower():
                    intensity_col = col
                    break

        if intensity_col not in df.columns:
            continue

        # Get median expression per protein
        median_expr = df.groupby("ProteinName")[intensity_col].median()

        for protein in proteins:
            if protein in median_expr.index:
                val = median_expr.loc[protein]
                if pd.notna(val) and val > 0:
                    protein_values[protein].append(np.log2(val + 1))

    # Compute stability metrics
    results = []
    for protein, values in protein_values.items():
        if len(values) >= 2:
            values = np.array(values)
            results.append({
                "ProteinName": protein,
                "MAD": np.median(np.abs(values - np.median(values))),
                "IQR": np.percentile(values, 75) - np.percentile(values, 25),
                "Range": values.max() - values.min(),
                "Mean": values.mean(),
                "Std": values.std(),
                "N_datasets": len(values),
            })

    return pd.DataFrame(results)


def run_all_metrics(
    datasets: dict = None,
    methods: List[str] = None,
) -> dict:
    """
    Compute all benchmarking metrics.

    Returns dict with all results.
    """
    if datasets is None:
        datasets = HELA_DATASETS
    if methods is None:
        methods = QUANTIFICATION_METHODS

    print("=" * 70)
    print("Computing Benchmarking Metrics")
    print("=" * 70)

    results = {
        "cv_summary": [],
        "cross_experiment_corr": {},
        "rank_consistency": {},
        "expression_stability": {},
        "tmt_lfq": {},
    }

    # 1. CV Summary
    print("\n1. Computing CV per method/dataset...")
    for method in methods:
        for dataset_id in datasets:
            cv = compute_cv_summary(dataset_id, method)
            if cv:
                results["cv_summary"].append(cv)
                print(f"   {dataset_id}/{method}: mean_cv={cv['mean_cv']:.3f}")

    # 2. Cross-experiment correlation
    print("\n2. Computing cross-experiment correlations...")
    for method in methods:
        corr_matrix = compute_cross_experiment_correlation(datasets, method)
        if not corr_matrix.empty:
            results["cross_experiment_corr"][method] = corr_matrix
            mean_corr = corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].mean()
            print(f"   {method}: mean_corr={mean_corr:.3f}")

    # 3. Rank consistency
    print("\n3. Computing rank consistency...")
    for method in methods:
        rank_corr = compute_rank_consistency(datasets, method)
        if not rank_corr.empty:
            results["rank_consistency"][method] = rank_corr
            mean_rank = rank_corr.values[np.triu_indices_from(rank_corr.values, k=1)].mean()
            print(f"   {method}: mean_rank_corr={mean_rank:.3f}")

    # 4. Expression stability for proteins of interest
    print("\n4. Computing expression stability for housekeeping proteins...")
    for method in methods:
        stability = compute_expression_stability(datasets, method)
        if not stability.empty:
            results["expression_stability"][method] = stability
            mean_mad = stability["MAD"].mean()
            print(f"   {method}: mean_MAD={mean_mad:.3f}")

    # 5. TMT vs LFQ correlation (PXD007683)
    print("\n5. Computing TMT vs LFQ correlation (PXD007683)...")
    for method in methods:
        tmt_lfq = compute_tmt_lfq_correlation(method)
        if tmt_lfq:
            results["tmt_lfq"][method] = tmt_lfq
            print(f"   {method}: R={tmt_lfq['pearson_r']:.3f}, n={tmt_lfq['n_proteins']}")

    return results


def save_results(results: dict, output_dir: Path = None):
    """Save all results to files."""
    if output_dir is None:
        output_dir = ANALYSIS_DIR

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nSaving results to {output_dir}...")

    # CV Summary
    if results["cv_summary"]:
        cv_df = pd.DataFrame(results["cv_summary"])
        cv_df.to_csv(output_dir / "cv_comparison.csv", index=False)
        print("  Saved: cv_comparison.csv")

    # Cross-experiment correlation
    for method, corr_matrix in results["cross_experiment_corr"].items():
        corr_matrix.to_csv(output_dir / f"cross_experiment_corr_{method}.csv")
    print("  Saved: cross_experiment_corr_*.csv")

    # Rank consistency
    for method, rank_matrix in results["rank_consistency"].items():
        rank_matrix.to_csv(output_dir / f"rank_consistency_{method}.csv")
    print("  Saved: rank_consistency_*.csv")

    # Expression stability
    for method, stability in results["expression_stability"].items():
        stability.to_csv(output_dir / f"expression_stability_{method}.csv", index=False)
    print("  Saved: expression_stability_*.csv")

    # TMT vs LFQ correlation
    if results.get("tmt_lfq"):
        tmt_lfq_summary = []
        for method, data in results["tmt_lfq"].items():
            # Save raw values for plotting
            plot_df = pd.DataFrame({
                "ProteinName": data["proteins"],
                "TMT": data["tmt_values"],
                "LFQ": data["lfq_values"],
            })
            plot_df.to_csv(output_dir / f"tmt_lfq_values_{method}.csv", index=False)

            # Summary stats
            tmt_lfq_summary.append({
                "method": method,
                "n_proteins": data["n_proteins"],
                "pearson_r": data["pearson_r"],
                "spearman_r": data["spearman_r"],
            })

        if tmt_lfq_summary:
            pd.DataFrame(tmt_lfq_summary).to_csv(output_dir / "tmt_lfq_comparison.csv", index=False)
            print("  Saved: tmt_lfq_comparison.csv, tmt_lfq_values_*.csv")

    # Summary metrics
    summary = []
    for method in QUANTIFICATION_METHODS:
        row = {"method": method}

        # Mean CV
        cv_data = [x for x in results["cv_summary"] if x["method"] == method]
        if cv_data:
            row["mean_cv"] = np.mean([x["mean_cv"] for x in cv_data])

        # Mean cross-experiment correlation
        if method in results["cross_experiment_corr"]:
            corr = results["cross_experiment_corr"][method]
            row["mean_cross_corr"] = corr.values[np.triu_indices_from(corr.values, k=1)].mean()

        # Mean rank consistency
        if method in results["rank_consistency"]:
            rank = results["rank_consistency"][method]
            row["mean_rank_corr"] = rank.values[np.triu_indices_from(rank.values, k=1)].mean()

        # Mean expression stability
        if method in results["expression_stability"]:
            row["mean_mad"] = results["expression_stability"][method]["MAD"].mean()

        summary.append(row)

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / "summary_metrics.csv", index=False)
    print("  Saved: summary_metrics.csv")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Compute benchmarking metrics"
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=QUANTIFICATION_METHODS,
        help="Compute metrics for specific method only"
    )
    parser.add_argument(
        "--hela-only",
        action="store_true",
        help="Only use HeLa datasets"
    )

    args = parser.parse_args()

    # Select datasets
    datasets = HELA_DATASETS if args.hela_only else ALL_DATASETS

    # Select methods
    methods = [args.method] if args.method else QUANTIFICATION_METHODS

    # Run
    results = run_all_metrics(datasets=datasets, methods=methods)

    # Save
    save_results(results)

    print("\nDone!")


if __name__ == "__main__":
    main()
