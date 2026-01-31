#!/usr/bin/env python3
"""
Generate sample intensity distribution plots before and after normalization
for PXD007683 TMT and LFQ datasets.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# Paths
RAW_DATA_DIR = Path(__file__).parent / "raw_data"
PLOTS_DIR = Path(__file__).parent / "plots"
PLOTS_DIR.mkdir(exist_ok=True)


def load_tmt_data():
    """Load TMT feature data."""
    df = pd.read_csv(RAW_DATA_DIR / "PXD007683-TMT_feature.csv")
    # Create sample ID from Channel
    df['SampleID'] = 'Channel-' + df['Channel'].astype(str)
    df['LogIntensity'] = np.log10(df['Intensity'].replace(0, np.nan))
    return df


def load_lfq_data():
    """Load LFQ feature data."""
    df = pd.read_csv(RAW_DATA_DIR / "PXD007683-LFQ_feature.csv")
    # Create sample ID from BioReplicate
    df['SampleID'] = 'Sample-' + df['BioReplicate'].astype(str)
    df['LogIntensity'] = np.log10(df['Intensity'].replace(0, np.nan))
    return df


def median_normalize(df, intensity_col='LogIntensity', sample_col='SampleID'):
    """Apply median normalization."""
    df = df.copy()

    # Calculate global median
    global_median = df[intensity_col].median()

    # Calculate per-sample median and normalize
    sample_medians = df.groupby(sample_col)[intensity_col].transform('median')
    df['NormIntensity'] = df[intensity_col] - sample_medians + global_median

    return df


def plot_intensity_distribution(df, intensity_col, sample_col, title, output_path,
                                 xlabel='log10(Intensity)', colors=None):
    """Plot intensity distribution for each sample using KDE."""
    fig, ax = plt.subplots(figsize=(10, 7))

    samples = sorted(df[sample_col].unique())

    if colors is None:
        cmap = plt.cm.get_cmap('tab20', len(samples))
        colors = [cmap(i) for i in range(len(samples))]

    for i, sample in enumerate(samples):
        sample_data = df[df[sample_col] == sample][intensity_col].dropna()

        if len(sample_data) > 10:
            # KDE estimation
            try:
                kde = stats.gaussian_kde(sample_data)
                x_range = np.linspace(sample_data.min() - 0.5, sample_data.max() + 0.5, 200)
                y = kde(x_range)
                ax.plot(x_range, y, label=sample, color=colors[i], linewidth=1.5)
            except:
                pass

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(title='Sample', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")


def plot_before_after(df_raw, df_norm, sample_col, dataset_name, output_path):
    """Plot before and after normalization side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    samples = sorted(df_raw[sample_col].unique())
    cmap = plt.cm.get_cmap('tab20', len(samples))
    colors = [cmap(i) for i in range(len(samples))]

    # Before normalization
    ax = axes[0]
    for i, sample in enumerate(samples):
        sample_data = df_raw[df_raw[sample_col] == sample]['LogIntensity'].dropna()
        if len(sample_data) > 10:
            try:
                kde = stats.gaussian_kde(sample_data)
                x_range = np.linspace(sample_data.min() - 0.5, sample_data.max() + 0.5, 200)
                y = kde(x_range)
                ax.plot(x_range, y, label=sample, color=colors[i], linewidth=1.5)
            except:
                pass

    ax.set_xlabel('log10(Intensity)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f'{dataset_name} - Before Normalization', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # After normalization
    ax = axes[1]
    for i, sample in enumerate(samples):
        sample_data = df_norm[df_norm[sample_col] == sample]['NormIntensity'].dropna()
        if len(sample_data) > 10:
            try:
                kde = stats.gaussian_kde(sample_data)
                x_range = np.linspace(sample_data.min() - 0.5, sample_data.max() + 0.5, 200)
                y = kde(x_range)
                ax.plot(x_range, y, label=sample, color=colors[i], linewidth=1.5)
            except:
                pass

    ax.set_xlabel('log10(Intensity)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f'{dataset_name} - After Median Normalization', fontsize=12, fontweight='bold')
    ax.legend(title='Sample', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {output_path}")


def main():
    print("=" * 60)
    print("Generating Intensity Distribution Plots for PXD007683")
    print("=" * 60)

    # Load TMT data
    print("\nLoading TMT data...")
    tmt_df = load_tmt_data()
    print(f"  Loaded {len(tmt_df):,} features, {tmt_df['SampleID'].nunique()} samples")

    # Normalize TMT
    tmt_norm = median_normalize(tmt_df)

    # Plot TMT
    plot_before_after(
        tmt_df, tmt_norm, 'SampleID',
        'PXD007683-TMT',
        PLOTS_DIR / 'PXD007683_TMT_intensity_distribution.png'
    )

    # Load LFQ data
    print("\nLoading LFQ data...")
    lfq_df = load_lfq_data()
    print(f"  Loaded {len(lfq_df):,} features, {lfq_df['SampleID'].nunique()} samples")

    # Normalize LFQ
    lfq_norm = median_normalize(lfq_df)

    # Plot LFQ
    plot_before_after(
        lfq_df, lfq_norm, 'SampleID',
        'PXD007683-LFQ',
        PLOTS_DIR / 'PXD007683_LFQ_intensity_distribution.png'
    )

    # Combined plot with all 4 panels
    print("\nGenerating combined plot...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    datasets = [
        (tmt_df, tmt_norm, 'SampleID', 'TMT'),
        (lfq_df, lfq_norm, 'SampleID', 'LFQ'),
    ]

    for row, (df_raw, df_norm, sample_col, name) in enumerate(datasets):
        samples = sorted(df_raw[sample_col].unique())
        cmap = plt.cm.get_cmap('tab20', len(samples))
        colors = [cmap(i) for i in range(len(samples))]

        # Before
        ax = axes[row, 0]
        for i, sample in enumerate(samples):
            sample_data = df_raw[df_raw[sample_col] == sample]['LogIntensity'].dropna()
            if len(sample_data) > 10:
                try:
                    kde = stats.gaussian_kde(sample_data)
                    x_range = np.linspace(sample_data.min() - 0.5, sample_data.max() + 0.5, 200)
                    y = kde(x_range)
                    ax.plot(x_range, y, label=sample, color=colors[i], linewidth=1.2)
                except:
                    pass
        ax.set_xlabel('log10(Intensity)', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(f'{name} - Before Normalization', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # After
        ax = axes[row, 1]
        for i, sample in enumerate(samples):
            sample_data = df_norm[df_norm[sample_col] == sample]['NormIntensity'].dropna()
            if len(sample_data) > 10:
                try:
                    kde = stats.gaussian_kde(sample_data)
                    x_range = np.linspace(sample_data.min() - 0.5, sample_data.max() + 0.5, 200)
                    y = kde(x_range)
                    ax.plot(x_range, y, label=sample, color=colors[i], linewidth=1.2)
                except:
                    pass
        ax.set_xlabel('log10(Intensity)', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(f'{name} - After Median Normalization', fontsize=11, fontweight='bold')
        ax.legend(title='Sample', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)
        ax.grid(True, alpha=0.3)

    plt.suptitle('PXD007683: Sample Intensity Distribution', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / 'PXD007683_intensity_distribution_combined.png',
                dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {PLOTS_DIR / 'PXD007683_intensity_distribution_combined.png'}")

    print("\nDone!")


if __name__ == "__main__":
    main()
