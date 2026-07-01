# tissuemap: Tissue Proteome Atlas

`tissuemap` is a Python periphery command in the mokume wheel. It builds a
per-dataset tissue proteome atlas from QPX parquet outputs. It is intended for
atlas-style tissue exploration and tissue-specificity scoring, not for standard
protein quantification from a single experiment.

!!! note "Periphery command — not a CLI subcommand"
    `tissuemap` is **not** one of the four Rust CLI subcommands
    (`features2proteins`, `features2peptides`, `peptides2protein`,
    `correct-batches`). It lives only in the Python wheel as
    `mokume.tissuemap(**kwargs)` (or `python -m mokume.commands.tissuemap`) and
    requires the `tissuemap` extra:
    ```bash
    pip install "mokume-rs[tissuemap]"
    ```

## When to Use TissueMap

Use `tissuemap` when you want to:

- build a tissue atlas from one dataset or a directory of datasets
- compute AdaTiSS-based tissue-specificity scores
- generate PCA/t-SNE embeddings and tissue-level marker plots
- save AnnData outputs for downstream exploration

If your goal is standard LFQ or TMT protein quantification, start with
[`features2proteins`](../user-guide/features2proteins.md) instead.

## Expected Input Layout

`scan_dir` can point to either:

- a single dataset directory containing `qpx_output/`
- a parent directory containing multiple dataset directories, each with `qpx_output/`

If a dataset is TMT and auto-detection is not sufficient, specify it explicitly
with repeated `tmt_dataset` values.

## Quick Start

```python
import mokume

# Generate a default YAML template
mokume.tissuemap(generate_config="tissuemap.yaml")

# Run a single dataset directory
mokume.tissuemap(
    scan_dir="QPX_data/tissues-mq/PXD016999",
    output_dir="./results",
)

# Run a parent directory containing multiple datasets
mokume.tissuemap(
    scan_dir="QPX_data/tissues-mq",
    tmt_dataset=["PXD016999"],   # a list repeats --tmt-dataset
    output_dir="./results",
    n_jobs=8,
)

# Run with a custom YAML configuration
mokume.tissuemap(
    scan_dir="QPX_data/tissues-mq",
    config="tissuemap.yaml",
    output_dir="./results",
)
```

The wrapper maps keyword arguments to the command's flags (`key=value` →
`--key value` with `_` rewritten to `-`; a list repeats the flag). The same
command is runnable from a shell as `python -m mokume.commands.tissuemap
--scan-dir ... --output-dir ...`.

## What the Pipeline Does

The current TissueMap workflow is organized around per-dataset processing:

1. Discover dataset directories from `scan_dir`
2. Load QPX-derived protein matrices and metadata
3. Apply log2 + median normalization
4. Harmonize tissue labels
5. Filter proteins with excessive missingness or contaminants
6. Apply batch correction
7. Compute AdaTiSS tissue-specificity scores
8. Build PCA + t-SNE embeddings
9. Generate atlas, marker, and tissue-specificity plots
10. Save AnnData and CSV outputs

## Main Options

| Option | Default | Description |
|--------|---------|-------------|
| `scan_dir` | required unless `generate_config` | Dataset directory or parent directory containing datasets |
| `output_dir` | `tissuemap_output` | Output directory for all results |
| `config` | none | YAML configuration file |
| `generate_config` | none | Write a default YAML template and exit |
| `tmt_dataset` | auto | Mark one or more dataset IDs as TMT (list repeats the flag) |
| `n_jobs` | `8` | Threads for dataset processing and embedding |
| `dpi` | `250` | Plot resolution override |
| `imputation_method` | config default | Imputation method used before embedding |
| `embedding_method` | config default | `tsne` or `umap` |

## Configuration Workflow

Use `generate_config` when you want to tune the pipeline before running it
repeatedly:

```python
import mokume

mokume.tissuemap(generate_config="tissuemap.yaml")
```

The generated YAML exposes the main configuration groups:

- `input` — dataset discovery, TMT overrides, minimum tissue sample settings
- `filtering` — NaN threshold and contaminant filtering
- `tissue_specificity` — AdaTiSS thresholds and scoring controls
- `embedding` — PCA/t-SNE parameters
- `plotting` — DPI, PDF export, marker plot controls
- `output` — output directory

Argument values such as `scan_dir`, `output_dir`, `n_jobs`, and `dpi` override
the YAML file.

## Output Files

Each processed dataset gets its own output directory.

| File | Description |
|------|-------------|
| `<ds_id>.corrected.h5ad` | Batch-corrected sample-level AnnData with embeddings and metadata |
| `<ds_id>.ts_scores.h5ad` | Tissue-specificity score matrix as AnnData |
| `protein_ts_scores.csv` | Per-protein tissue-specificity scores and enrichment categories |
| `plots/` | All figures (PNG, plus a PDF copy of each when `save_pdf` is enabled) — see below |

Every figure in `plots/` is written as a PNG and, when `plotting.save_pdf` is
enabled (the default), an accompanying PDF:

| Plot file | Shows |
|-----------|-------|
| `tissue_atlas.png` | t-SNE tissue atlas: samples coloured by organ, with proteomic-group hulls |
| `tissue_dendrogram.png` | Hierarchical dendrogram of tissues (Ward linkage, correlation distance) |
| `slide_atlas_dendrogram.png` | Combined atlas + dendrogram in one 16:9 slide |
| `marker_tsne.png` | t-SNE panels, each coloured by a top tissue marker's expression |
| `marker_heatmap.png` | Heatmap of top tissue markers across tissues |
| `marker_dotplot.png` | Dot plot of marker expression and detection fraction per tissue |
| `ts_distribution.png` | AdaTiSS TS-score distribution + tissue-specific protein counts |
| `specific_per_tissue.png` | Bar chart of tissue-specific protein counts per tissue |
| `pca_scree.png` | PCA scree plot (variance explained per component) |

## Example Outputs

Figures below are real `tissuemap` output on **PXD030304** (a 949-cell-line
proteomic panel), rendered by the plotting code described above.

![Tissue atlas: cell-line proteomes embedded and grouped by tissue of origin](../assets/pxd030304_tissue_atlas.png){ width="100%" }

*`tissue_atlas.png` — the cell-line proteomes embedded by t-SNE and coloured by
organ of origin.*

![t-SNE panels coloured by top tissue-marker expression](../assets/pxd030304_marker_tsne.png){ width="100%" }

*`marker_tsne.png` — the same embedding, each panel coloured by a top
tissue-marker's expression.*

![AdaTiSS tissue-specificity score distribution and tissue-specific protein counts](../assets/pxd030304_ts_distribution.png){ width="100%" }

*`ts_distribution.png` — AdaTiSS tissue-specificity score distribution (with the
GMM-fitted thresholds) and the tissue-specific protein count per tissue.*

## Practical Tips

- Use `features2proteins` for standard quantification workflows; use `tissuemap`
  for atlas-style tissue analysis.
- Start with the generated YAML template if you expect to rerun multiple datasets.
- Pass a list to `tmt_dataset` when one or more datasets should be treated as TMT
  explicitly.
- Keep the default plotting enabled for the first run so you can inspect atlas
  quality and marker behavior.
