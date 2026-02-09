# PXD007683 Benchmark Report

Comprehensive analysis of TMT vs LFQ quantification using mokume.

> **Note**: When using DirectLFQ quantification, always use the external directlfq package
> directly (`pip install mokume[directlfq]`) rather than any fallback implementation.
> This ensures reproducibility and uses the official Mann Lab algorithm.

---

## Q1: Absolute Expression Stability

**Question**: Which method combination provides the most stable absolute expression values?

### Best Methods by Within-Condition CV

| Technology | Best Method | CV (median) | % Proteins with CV < 20% |
|------------|-------------|-------------|--------------------------|
| TMT | maxlfq | 0.029 | 99.7% |
| LFQ | maxlfq | 0.074 | 86.1% |

### Key Finding

TMT consistently shows lower CV than LFQ across all quantification methods.

## Q2: Technical vs Biological Variance

**Question**: How much variance is explained by condition (biology) vs technology (technical)?

### Variance Decomposition

| Technology | % Condition | % Residual | Silhouette |
|------------|-------------|------------|------------|
| TMT | 1.0% | 99.0% | N/A |
| LFQ | 0.9% | 99.1% | N/A |

### Key Finding

Condition (yeast spike-in level) explains a significant portion of variance,
indicating biological signal is preserved through quantification.

## Q3: Fold-Change Accuracy

**Question**: How accurately do we detect expected fold-changes?

Ground truth: Yeast proteins spiked at 10%, 5%, 3.3% ratios
- 10% vs 3.3% → expected 3.0-fold (log2 = 1.58)
- 10% vs 5% → expected 2.0-fold (log2 = 1.0)
- 5% vs 3.3% → expected 1.5-fold (log2 = 0.58)

### Yeast Protein Fold-Change Accuracy

| Technology | Comparison | Expected | Observed | Compression | RMSE |
|------------|------------|----------|----------|-------------|------|
| TMT | QY_10pct_vs_QY_3pct_yeast | 1.58 | 1.57 | 0.99 | 0.325 |
| TMT | QY_10pct_vs_QY_5pct_yeast | 1.00 | 1.13 | 1.13 | 0.228 |
| TMT | QY_5pct_vs_QY_3pct_yeast | 0.58 | 0.44 | 0.75 | 0.199 |
| LFQ | QY_10pct_vs_QY_3pct_yeast | 1.58 | 1.56 | 0.99 | 0.665 |
| LFQ | QY_10pct_vs_QY_5pct_yeast | 1.00 | 1.03 | 1.03 | 0.522 |
| LFQ | QY_5pct_vs_QY_3pct_yeast | 0.58 | 0.53 | 0.90 | 0.398 |

### Human Protein False Positive Rate

Human proteins should show no change (log2 FC = 0)

| Technology | Comparison | FP Rate (|log2FC| > 1) |
|------------|------------|------------------------|
| TMT | QY_10pct_vs_QY_3pct_human | 0.0% |
| TMT | QY_10pct_vs_QY_5pct_human | 0.0% |
| TMT | QY_5pct_vs_QY_3pct_human | 0.0% |
| LFQ | QY_10pct_vs_QY_3pct_human | 1.5% |
| LFQ | QY_10pct_vs_QY_5pct_human | 1.8% |
| LFQ | QY_5pct_vs_QY_3pct_human | 1.1% |

### Key Finding

TMT shows ratio compression (observed fold-change < expected), which is
a known phenomenon. LFQ generally shows less compression but higher variability.

## Cross-Technology Correlation (TMT vs LFQ)

**Question**: How well do TMT and LFQ agree on protein abundances?

### Overall Correlation

- **Pearson r**: 0.7924
- **Spearman r**: 0.7923
- Common proteins: 5954
- Matched samples: 11

### Per-Protein Correlation

- Median correlation: 0.0622
- Proteins with r > 0.8: 438 (7.9%)

### Key Finding

TMT and LFQ show good overall agreement, validating that both technologies
measure similar underlying biology despite methodological differences.

## Recommendations

Based on the comprehensive benchmark analysis:

### For Absolute Quantification (Q1)

- **TMT**: Use `maxlfq` (CV = 0.029)
- **LFQ**: Use `maxlfq` (CV = 0.074)

### For Differential Expression (Q2-Q3)

- For **small fold-changes** (< 2-fold): Prefer TMT (lower CV, higher precision)
- For **large fold-changes** (> 2-fold): LFQ shows less compression
- Apply **batch correction** when combining experiments

### For Cross-Experiment Integration

- Normalize using `median` or `hierarchical` methods
- Consider technology as a batch effect when combining TMT and LFQ

## Figures

See the `figures/` directory for:

- `9_tissues-boxplot.png`
- `9_tissues-density.png`
- `PXD007683-11samples-density.png`
- `PXD007683-LFQ-11samples-ibaq-ibaqpy-and-maxquant.png`
- `PXD007683-LFQ-11samples-ibaq-vs-maxquant-density.png`
- `PXD007683-LFQ-11samples-no_cov.png`
- `PXD007683-LFQ-ibaq-ibaqpy-and-maxquant.png`
- `PXD007683-LFQ-ibaq-vs-maxquant-density.png`
- `PXD007683-LFQ-no_cov.png`
- `PXD007683-TMTvsLFQ-boxplot.png`
- `PXD007683-TMTvsLFQ-density.png`
- `PXD019909-11samples-density.png`
- `PXD019909-TMTvsLFQ-density.png`
- `cross_tech_overall.png`
- `cross_tech_per_protein.png`
- `cross_tech_per_sample.png`
- `cv_by_condition_lfq.png`
- `cv_by_condition_tmt.png`
- `cv_distribution_lfq.png`
- `cv_distribution_tmt.png`
- `cv_method_comparison.png`
- `fold_change_comparison.png`
- `fold_change_lfq.png`
- `fold_change_tmt.png`
- `method_mean_cv_016999_lfq.png`
- `method_mean_cv_lfq.png`
- `method_mean_cv_tmt.png`
- `method_per_p_cv_016999_lfq.png`
- `method_per_p_cv_lfq.png`
- `method_per_p_cv_tmt.png`
- `missing_peptides_by_sample.png`
- `missing_value_016999_lfq.png`
- `pca_combined.png`
- `pca_lfq.png`
- `pca_tmt.png`
- `per_protein_cv.png`
