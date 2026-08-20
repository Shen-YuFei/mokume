# Visualization & Reports

mokume's plotting and reporting live in the Python periphery of the wheel: static
plots (matplotlib/seaborn) and interactive HTML reports (plotly) for quality
control and differential expression analysis. The Rust kernel is pure compute —
it writes the protein matrix and the per-contrast DE result CSVs, and the
periphery commands **read those tables** to render figures, so the cells in the
plots always match the cells in the kernel output (the periphery never
recomputes the numbers).

!!! note "Periphery extras"
    Plotting and reports are **not** part of the Rust CLI. They are wheel-only
    commands; install the extra you need:
    ```bash
    pip install "mokume[plotting]"   # matplotlib + seaborn (DE plots, t-SNE, piBAQ QC)
    pip install "mokume[reports]"    # plotly interactive DE report
    pip install "mokume[analysis]"   # single-matrix QC + workflow-comparison reports
    pip install "mokume[all]"        # everything
    ```

!!! warning "The `--plot-*` / `--interactive-report` flags are gone"
    `features2proteins` no longer accepts `--plot-dir`, `--plot-volcano`,
    `--plot-heatmap`, `--plot-pca`, `--highlight-genes`, `--interactive-report`,
    or `--report-output`. The Rust kernel returns `NotImplemented` for plotting
    and report output. Run the kernel to produce the protein matrix and DE CSVs,
    then call the wheel periphery (`mokume.de_plots` / `mokume.interactive_report`)
    on those files.

## DE Plots from features2proteins Output

First run the kernel to write the protein matrix and one DE result CSV per
contrast:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --de --de-contrasts "NASH vs HL,NASH vs Control" \
    --de-output de_results
```

Then render the plots from those CSVs with the periphery command. `de_plots`
takes an explicit argument list because the per-contrast `--contrast KEY A B CSV`
flag repeats (keyword arguments cannot express a repeated 4-tuple):

```python
import mokume

mokume.de_plots([
    "--protein-matrix", "proteins.csv",
    "--plot-dir", "plots",
    "--sdrf", "experiment.sdrf.tsv",
    "--volcano", "--heatmap", "--pca",
    "--highlight-genes", "COL10A1,FN1,ALB",
    "--contrast", "NASH-HL", "NASH", "HL", "de_results_NASH-HL.csv",
    "--contrast", "NASH-Control", "NASH", "Control", "de_results_NASH-Control.csv",
])
```

This generates:

- `plots/volcano_NASH-HL.png` -- Volcano plot for each contrast
- `plots/heatmap_NASH-HL.png` -- Per-contrast heatmap showing top 50 significant proteins (by |log2FC|) and only the two compared conditions. Skipped if no significant proteins exist for that contrast.
- `plots/pca_conditions.png` -- PCA colored by experimental condition (all samples)

The `--volcano` / `--heatmap` / `--pca` flags select which plots to render;
`--log2fc-threshold` (default 0.5) and `--fdr-threshold` (default 0.05) mirror the
kernel's `--de-log2fc` / `--de-fdr` so the significance cutoffs match. The same
command is runnable as `python -m mokume.commands.de_plots ...`.

## Interactive HTML Report

Generate a comprehensive interactive QC + DE report from the same kernel CSVs.
`interactive_report` also takes an explicit argument list for the repeated
`--contrast` flag:

```python
import mokume

mokume.interactive_report([
    "--protein-matrix", "proteins.csv",
    "--sdrf", "experiment.sdrf.tsv",
    "--report-output", "qc_report.html",
    "--highlight-genes", "COL10A1,FN1",
    "--contrast", "NASH-HL", "NASH", "HL", "de_results/NASH_vs_HL.csv",
])
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

## Single-Matrix QC Report

For a QC report built directly from a protein matrix CSV (with optional DE
results and SDRF grouping), use `mokume.qc_report` (`analysis` extra). It computes
PCA / t-SNE / silhouette / variance decomposition / CV / missing-value and
DE-quality metrics and writes an interactive HTML report:

```python
import mokume

mokume.qc_report(
    protein_matrix="proteins.csv",
    sdrf="experiment.sdrf.tsv",   # optional sample -> condition grouping
    de_results="de_results/NASH_vs_HL.csv",  # optional
    output="qc_report.html",
    title="My Experiment QC Report",
    is_log2=False,
)
```

For volcano gene-highlighting, call
`mokume.reports.qc_report.generate_qc_report` directly.

## Workflow Comparison Report

Compare multiple quantification workflows side by side (`analysis` extra). Pass a
list of workflow dicts, each with a `name` and either a `protein_df` DataFrame or
a `protein_matrix` CSV path (plus optional `de_results`, `sample_to_condition` /
`sdrf`, and `is_log2`):

```python
import mokume

mokume.workflow_comparison(
    workflows=[
        {"name": "IRS_RLR",     "protein_matrix": "irs_rlr.csv",
         "de_results": "irs_rlr_de.csv",     "sdrf": "experiment.sdrf.tsv"},
        {"name": "Ratio_ComBat", "protein_matrix": "ratio_combat.csv",
         "de_results": "ratio_combat_de.csv", "sdrf": "experiment.sdrf.tsv"},
    ],
    output="comparison.html",
    title="Method Comparison",
    marker_genes=["COL10A1", "FN1", "ALB"],
)
```

This generates a report with:

- Summary metrics table (silhouette, CV, DE count, pi1)
- PCA grids for each workflow
- DE concordance analysis (Jaccard, log2FC correlation, CAT curves)
- Marker gene results across workflows

## t-SNE Visualization Command

A standalone periphery command for t-SNE visualization from a folder of protein
files (`plotting` extra). It is **not** a Rust CLI subcommand — it is the wheel
command `mokume.tsne_visualization(**kwargs)` (or
`python -m mokume.commands.visualize`):

```python
import mokume

mokume.tsne_visualization(
    folder="protein_folder/",
    pattern="proteins.tsv",   # default
)
```

## Plotting Library API

Besides the canned commands above, the `mokume.plotting` module exposes the
individual figure functions as a callable library (`plotting` extra). Use these
when you want to build figures directly from a kernel CSV in your own script
instead of running `de_plots`. Each takes a pandas `DataFrame` and either
returns a matplotlib `Figure` or writes to an `output_file` / `file_name` path.

| Function | Purpose |
|----------|---------|
| `plot_volcano(de_results, log2fc_threshold=0.5, fdr_threshold=0.05, highlight_genes=None, output_file=None, ...)` | Volcano plot from a per-contrast DE result frame |
| `plot_heatmap(protein_df, sample_to_condition, proteins=None, top_n=50, output_file=None, ...)` | Clustered heatmap of protein intensities |
| `plot_pca_conditions(protein_df, sample_to_condition, output_file=None, ...)` | PCA scatter coloured by experimental condition |
| `plot_pca(df_pca, output_file, x_col="PC1", y_col="PC2", hue_col="batch", ...)` | PCA scatter from a precomputed PC table |
| `compute_pca_with_plot(df, n_components=5)` | Fit PCA, show the variance-explained scree, return the PC frame |
| `plot_tsne(df, x_col, y_col, hue_col, file_name)` | t-SNE scatter from a frame with embedding columns |
| `plot_distributions(dataset, field, class_field, log2=True, ...)` | Per-class quantile/distribution plot |
| `plot_box_plot(dataset, field, class_field, violin=False, ...)` | Box (or violin) plot of a field grouped by class |
| `is_plotting_available()` | `True` if matplotlib/seaborn are installed |

```python
import pandas as pd
from mokume.plotting import plot_volcano, plot_heatmap

de = pd.read_csv("proteins.de.Treatment_vs_Control.csv")
plot_volcano(de, fdr_threshold=0.05, log2fc_threshold=1.0,
             output_file="volcano.png")

proteins = pd.read_csv("proteins.csv", index_col=0)
sample_to_condition = {"S1": "Treatment", "S2": "Control"}  # ...
plot_heatmap(proteins, sample_to_condition, top_n=50, output_file="heatmap.png")
```

Functions that return a `Figure` (volcano, heatmap, PCA-by-condition,
distributions, box plot) also save to `output_file`/`file_name` when given one;
without it they return the figure for you to compose or save yourself.
