# Visualization & Reports

mokume provides both static plots (matplotlib/seaborn) and interactive HTML reports (plotly) for quality control and differential expression analysis.

!!! note "Optional dependencies"
    ```bash
    pip install mokume[plotting]   # matplotlib + seaborn
    pip install mokume[reports]    # plotly interactive reports
    pip install mokume[all]        # everything
    ```

## Integrated Plotting (features2proteins)

The easiest way to generate plots is through the unified pipeline:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --de --de-contrasts "NASH-HL" \
    --plot-dir plots/ \
    --plot-volcano \
    --plot-heatmap \
    --plot-pca \
    --highlight-genes "COL10A1,FN1,ALB"
```

This generates:

- `plots/volcano_NASH-HL.png` -- Volcano plot for each contrast
- `plots/heatmap_NASH-HL.png` -- Heatmap of top DE proteins
- `plots/pca_conditions.png` -- PCA colored by experimental condition

## Interactive HTML Reports

Generate a comprehensive interactive QC report:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --de --de-contrasts "NASH-HL" \
    --interactive-report \
    --report-output qc_report.html
```

The interactive report includes:

- **PCA plot** -- Interactive scatter with hover annotations
- **t-SNE plot** -- Non-linear dimensionality reduction
- **Sample correlation heatmap** -- Pairwise Pearson correlations
- **CV distribution** -- Per-condition coefficient of variation
- **Missing value analysis** -- Missing rate per sample
- **Distribution box plots** -- Intensity distributions
- **Silhouette score** -- Condition clustering quality
- **Variance decomposition** -- PVCA-style analysis (condition vs batch vs residual)
- **Volcano plot** -- DE results with searchable gene names
- **Expression heatmap** -- Top DE proteins

### Python API

```python
from mokume.reports.qc_report import generate_qc_report, compute_qc_metrics

# Generate full report
generate_qc_report(
    protein_df=protein_df,
    sample_to_condition=sample_to_condition,
    sample_to_plex=sample_to_plex,      # optional
    de_results=de_results,               # optional
    output_html="qc_report.html",
    title="My Experiment QC Report",
    is_log2=False,
    highlight_genes=["COL10A1", "FN1"],
)

# Or just compute metrics
metrics = compute_qc_metrics(
    protein_df, sample_to_condition,
    sample_to_plex=sample_to_plex,
    de_results=de_results,
)
```

## DE-Specific Report

A focused differential expression report with interactive volcano and expression plots:

```python
from mokume.reports.interactive import generate_de_report

generate_de_report(
    de_results=de_results,
    protein_df=protein_df,
    sample_to_condition=sample_to_condition,
    output_html="de_report.html",
    title="NASH vs HL",
    highlight_genes=["COL10A1", "FN1"],
    log2fc_threshold=0.5,
    fdr_threshold=0.05,
)
```

## Workflow Comparison Report

Compare multiple quantification workflows side by side:

```python
from mokume.reports.workflow_comparison import generate_comparison_report

workflows = {
    "IRS_TMM": {"protein_df": df1, "de_results": de1, "metrics": metrics1},
    "Ratio_ComBat": {"protein_df": df2, "de_results": de2, "metrics": metrics2},
}

generate_comparison_report(
    workflows=workflows,
    output_html="comparison.html",
    title="Method Comparison",
    marker_genes=["COL10A1", "FN1", "ALB"],
)
```

This generates a report with:

- Summary metrics table (silhouette, CV, DE count, pi1)
- PCA grids for each workflow
- DE concordance analysis (Jaccard, log2FC correlation, CAT curves)
- Marker gene results across workflows

## Standalone Plots

### Distribution Plots

```python
from mokume.plotting.distributions import plot_distributions, plot_box_plot

# KDE density plot
fig = plot_distributions(
    dataset=df, field="Intensity", class_field="Condition",
    title="Intensity Distribution", log2=True,
)

# Box plot
fig = plot_box_plot(
    dataset=df, field="Intensity", class_field="SampleID",
    log2=True, rotation=45, violin=False,
)
```

### PCA and t-SNE

```python
from mokume.plotting.pca import compute_pca_with_plot, plot_pca, plot_tsne

# PCA with variance explained plot
df_pca = compute_pca_with_plot(df_wide, n_components=5)

# PCA scatter plot
plot_pca(
    df_pca, output_file="pca.png",
    x_col="PC1", y_col="PC2", hue_col="Condition",
)
```

### Volcano Plot

```python
from mokume.plotting.differential_expression import (
    plot_volcano, plot_heatmap, plot_pca_conditions,
)

# Volcano plot
fig = plot_volcano(
    de_results,
    log2fc_threshold=0.5,
    fdr_threshold=0.05,
    highlight_genes=["COL10A1", "FN1"],
    output_file="volcano.png",
)

# Heatmap of top DE proteins
fig = plot_heatmap(
    protein_df, sample_to_condition,
    top_n=50,
    output_file="heatmap.png",
)

# PCA colored by condition
fig = plot_pca_conditions(
    protein_df, sample_to_condition,
    output_file="pca_conditions.png",
)
```

## t-SNE Visualization Command

A standalone command for t-SNE visualization from protein files:

```bash
mokume tsne-visualization \
    -f protein_folder/ \
    -o proteins.tsv
```
