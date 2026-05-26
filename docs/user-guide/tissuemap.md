# tissuemap: Tissue Proteome Atlas

The `tissuemap` command builds a per-dataset tissue proteome atlas from QPX parquet outputs. It is intended for atlas-style tissue exploration and tissue-specificity scoring, not for standard protein quantification from a single experiment.

!!! note "Optional dependency"
    ```bash
    pip install mokume[tissuemap]
    ```

## When to Use TissueMap

Use `tissuemap` when you want to:

- build a tissue atlas from one dataset or a directory of datasets
- compute AdaTiSS-based tissue-specificity scores
- generate PCA/t-SNE embeddings and tissue-level marker plots
- save AnnData outputs for downstream exploration

If your goal is standard LFQ or TMT protein quantification, start with [`features2proteins`](features2proteins.md) instead.

## Expected Input Layout

`--scan-dir` can point to either:

- a single dataset directory containing `qpx_output/`
- a parent directory containing multiple dataset directories, each with `qpx_output/`

If a dataset is TMT and auto-detection is not sufficient, specify it explicitly with repeated `--tmt-dataset` flags.

## Quick Start

```bash
# Generate a default YAML template
mokume tissuemap --generate-config tissuemap.yaml

# Run a single dataset directory
mokume tissuemap \
    --scan-dir QPX_data/tissues-mq/PXD016999 \
    --output-dir ./results

# Run a parent directory containing multiple datasets
mokume tissuemap \
    --scan-dir QPX_data/tissues-mq \
    --tmt-dataset PXD016999 \
    --output-dir ./results \
    --n-jobs 8

# Run with a custom YAML configuration
mokume tissuemap \
    --scan-dir QPX_data/tissues-mq \
    --config tissuemap.yaml \
    --output-dir ./results
```

## What the Pipeline Does

The current TissueMap workflow is organized around per-dataset processing:

1. Discover dataset directories from `--scan-dir`
2. Load QPX-derived protein matrices and metadata
3. Apply log2 + median normalization
4. Harmonize tissue labels
5. Filter proteins with excessive missingness or contaminants
6. Apply batch correction
7. Compute AdaTiSS tissue-specificity scores
8. Build PCA + t-SNE embeddings
9. Generate atlas, marker, and tissue-specificity plots
10. Save AnnData and CSV outputs

## Main CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `--scan-dir` | required unless `--generate-config` | Dataset directory or parent directory containing datasets |
| `--output-dir` | `tissuemap_output` | Output directory for all results |
| `--config` | none | YAML configuration file |
| `--generate-config` | none | Write a default YAML template and exit |
| `--tmt-dataset` | auto | Mark one or more dataset IDs as TMT |
| `--n-jobs` | `8` | Threads for dataset processing and embedding |
| `--dpi` | `250` | Plot resolution override |

## Configuration Workflow

Use `--generate-config` when you want to tune the pipeline before running it repeatedly:

```bash
mokume tissuemap --generate-config tissuemap.yaml
```

The generated YAML exposes the main configuration groups:

- `input` — dataset discovery, TMT overrides, minimum tissue sample settings
- `filtering` — NaN threshold and contaminant filtering
- `tissue_specificity` — AdaTiSS thresholds and scoring controls
- `embedding` — PCA/t-SNE parameters
- `plotting` — DPI, PDF export, marker plot controls
- `output` — output directory

CLI values such as `--scan-dir`, `--output-dir`, `--n-jobs`, and `--dpi` override the YAML file.

## Output Files

Each processed dataset gets its own output directory.

| File | Description |
|------|-------------|
| `<ds_id>.corrected.h5ad` | Batch-corrected sample-level AnnData with embeddings and metadata |
| `<ds_id>.ts_scores.h5ad` | Tissue-specificity score matrix as AnnData |
| `protein_ts_scores.csv` | Per-protein tissue-specificity scores and enrichment categories |
| `plots/` | PCA scree plot, atlas/dendrogram, marker plots, and TS distribution plots |

## Python API

```python
from pathlib import Path

from mokume.tissuemap.config import InputConfig, OutputConfig, TissueMapConfig, load_config
from mokume.tissuemap.pipeline import TissueMapPipeline

# Programmatic configuration
config = TissueMapConfig(
    n_jobs=8,
    input=InputConfig(scan_dir=Path("QPX_data/tissues-mq/PXD016999")),
    output=OutputConfig(output_dir=Path("./results")),
)
TissueMapPipeline(config).run()

# Load from YAML and override selected fields
config = load_config(
    Path("tissuemap.yaml"),
    overrides={
        "input.scan_dir": "QPX_data/tissues-mq",
        "output.output_dir": "./results",
        "n_jobs": 8,
    },
)
TissueMapPipeline(config).run()
```

## Practical Tips

- Use `features2proteins` for standard quantification workflows; use `tissuemap` for atlas-style tissue analysis.
- Start with the generated YAML template if you expect to rerun multiple datasets.
- Use repeated `--tmt-dataset` flags when a dataset should be treated as TMT explicitly.
- Keep the default plotting enabled for the first run so you can inspect atlas quality and marker behavior.
