# Configuration

## PipelineConfig

The `PipelineConfig` dataclass controls the `QuantificationPipeline`. It uses nested sub-configurations for each pipeline stage.

```python
from mokume.pipeline import PipelineConfig
from mokume.pipeline.config import (
    InputConfig,
    FilterConfig,
    NormalizationConfig,
    QuantificationConfig,
    IRSConfig,
    BatchCorrectionConfig,
    ImputationConfig,
    DEConfig,
    OutputConfig,
    RuntimeConfig,
)
```

### InputConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `parquet` | `str` | required | Input parquet file path |
| `sdrf` | `str \| None` | `None` | SDRF metadata file |
| `fasta_file` | `str \| None` | `None` | FASTA file (for iBAQ) |

### FilterConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `min_aa` | `int` | `7` | Minimum amino acid length |
| `min_unique_peptides` | `int` | `2` | Minimum unique peptides per protein |
| `remove_contaminants` | `bool` | `True` | Remove contaminants and decoys |

### NormalizationConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `run_method` | `str` | `"median"` | Run-level normalization: median, mean, max, global, max_min, iqr, none |
| `sample_method` | `str` | `"globalMedian"` | Sample-level: globalMedian, conditionMedian, hierarchical, tmm, none |
| `proteins_file` | `str \| None` | `None` | File with protein IDs for normalization |

### QuantificationConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `method` | `str` | `"maxlfq"` | Quantification method |
| `ion_alignment` | `str \| None` | `None` | Ion alignment: none or hierarchical |
| `coverage_threshold` | `float \| None` | `None` | Min non-missing fraction per condition |
| `ratio_fraction_merge` | `str` | `"mean"` | Fraction merge: mean or max |
| `directlfq_num_cores` | `int \| None` | `None` | CPU cores for DirectLFQ |
| `directlfq_min_nonan` | `int` | `1` | Min non-NaN values |
| `directlfq_num_samples_quadratic` | `int` | `50` | Quadratic threshold |

### IRSConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | `bool` | `False` | Enable IRS normalization |
| `reference_samples` | `list \| None` | `None` | Reference sample names |
| `sdrf_column` | `str \| None` | `None` | SDRF column for detection |
| `sdrf_values` | `list \| None` | `None` | Reference indicator values |
| `reference_regex` | `str` | `"pool\|powder\|ref\|reference\|bridge"` | Auto-detection regex |
| `stat` | `str` | `"median"` | Plex reference statistic |
| `remove_reference` | `bool` | `False` | Remove reference samples |

### BatchCorrectionConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | `bool` | `False` | Enable ComBat |
| `method` | `str` | `"sample_prefix"` | Batch detection: sample_prefix, run, column |
| `column` | `str \| None` | `None` | SDRF column (for method="column") |
| `covariates` | `list \| None` | `None` | Covariate columns to preserve |
| `parametric` | `bool` | `True` | Use parametric ComBat |
| `mean_only` | `bool` | `False` | Only correct mean (not variance) |
| `ref_batch` | `int \| None` | `None` | Reference batch index |

### ImputationConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | `bool` | `False` | Enable missing value imputation |
| `method` | `str` | `"none"` | Method: none, knn, minprob, mindet, qrilc, missforest, seqknn, mle, mice, nbavg, gms, bpca, impseq, impseqrob |
| `quantile` | `float` | `0.01` | Quantile for MinProb / MinDet / QRILC low-tail draw |
| `shift` | `float` | `1.6` | MinProb shift in standard deviations |
| `scale` | `float` | `0.3` | MinProb scale factor for sigma |
| `n_neighbors` | `int` | `5` | Number of neighbours for KNN / SeqKNN / NBavg |

The `ImputationStage` runs in log2 space (the protein matrix is transformed before imputation and back to linear afterwards) so that censored-aware methods like MinProb / MinDet / QRILC behave correctly. The CLI exposes this through `--impute`, `--impute-method`, and the per-method `--impute-*` parameters.

### DEConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | `bool` | `False` | Enable DE analysis |
| `contrasts` | `list \| None` | `None` | Contrasts (e.g., ["A-B"]) |
| `method` | `str` | `"auto"` | Method: auto, limrots, limma, deqms, proda, rots, msstats, ensemble |
| `ensemble_methods` | `list \| None` | `None` | Member methods used when `method="ensemble"` (default `["limrots", "deqms", "proda"]`) |
| `ensemble_min_k` | `int` | `2` | Minimum ensemble members that must agree on direction |
| `log2fc_threshold` | `float` | `0.5` | Min absolute log2 fold change |
| `fdr_threshold` | `float` | `0.05` | Max FDR |
| `fdr_method` | `str` | `"bh"` | FDR correction: bh or ihw |
| `output` | `str \| None` | `None` | Output file for DE results |

When `method="ensemble"`, each member method runs on the same contrast and the per-protein verdicts are combined via top-k consensus with Fisher-combined p-values.

### OutputConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `export_peptides` | `str \| None` | `None` | Export peptides to file |
| `export_ions` | `str \| None` | `None` | Export ions to file |
| `plot_dir` | `str \| None` | `None` | Plot output directory |
| `plot_volcano` | `bool` | `False` | Generate volcano plots |
| `plot_heatmap` | `bool` | `False` | Generate heatmaps |
| `plot_pca` | `bool` | `False` | Generate PCA plots |
| `highlight_genes` | `list \| None` | `None` | Genes to highlight in plots |
| `interactive_report` | `bool` | `False` | Generate HTML QC report |
| `report_output` | `str \| None` | `None` | Report output path |

### RuntimeConfig

DuckDB resource limits applied via `SET memory_limit=...` / `SET threads=...`
pragmas during `Feature` initialisation. Also drives `effective_workers()`, the
memory-aware worker budget used by DirectLFQ / MaxLFQ fork pools.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `duckdb_memory` | `str \| None` | `None` (DuckDB autoconfig, ~80% of total RAM) | Memory limit string (e.g. `"80GB"`, `"16384MB"`) |
| `duckdb_threads` | `int \| None` | `None` (all cores) | Thread count for DuckDB (must be >= 1) |

!!! warning "`duckdb_memory` is not a hard process cap"
    The field only sizes DuckDB's internal buffer pool. PyArrow, polars, and
    pandas each have their own independent allocators (mimalloc / jemalloc /
    system malloc), so the surrounding Python process can grow to **2-3x**
    `duckdb_memory` on wide pivots. For production hard caps, use cgroup
    `MemoryMax` (systemd), SLURM `--mem`, or container
    `resources.limits.memory`.

The CLI exposes these as `--duckdb-memory` / `--duckdb-threads`.

### Full Example

```python
config = PipelineConfig(
    input=InputConfig(
        parquet="features.parquet",
        sdrf="experiment.sdrf.tsv",
    ),
    filtering=FilterConfig(
        min_aa=7,
        min_unique_peptides=2,
        remove_contaminants=True,
    ),
    normalization=NormalizationConfig(
        run_method="median",
        sample_method="globalMedian",
    ),
    quantification=QuantificationConfig(
        method="median",
    ),
    irs=IRSConfig(
        enabled=True,
        remove_reference=True,
    ),
    batch=BatchCorrectionConfig(
        enabled=True,
        method="sample_prefix",
        covariates=["characteristics[sex]"],
    ),
    de=DEConfig(
        enabled=True,
        contrasts=["NASH-HL"],
        method="auto",
        fdr_method="ihw",
    ),
    output=OutputConfig(
        plot_dir="plots/",
        plot_volcano=True,
        plot_pca=True,
        interactive_report=True,
        report_output="qc_report.html",
    ),
)
```

---

## TissueMapConfig

The `TissueMapConfig` dataclass controls the `TissueMapPipeline`.

```python
from mokume.tissuemap.config import (
    TissueMapConfig,
    InputConfig as TissueMapInputConfig,
    FilteringConfig,
    TissueSpecificityConfig,
    EmbeddingConfig,
    PlottingConfig,
    OutputConfig as TissueMapOutputConfig,
)
```

### Top-Level Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `n_jobs` | `int` | `8` | Threads used for dataset processing and embedding |
| `input` | `TissueMapInputConfig` | default factory | Dataset discovery and input controls |
| `filtering` | `FilteringConfig` | default factory | Protein filtering controls |
| `tissue_specificity` | `TissueSpecificityConfig` | default factory | AdaTiSS tissue-specificity scoring controls |
| `embedding` | `EmbeddingConfig` | default factory | PCA / t-SNE settings |
| `plotting` | `PlottingConfig` | default factory | Plot output controls |
| `output` | `TissueMapOutputConfig` | default factory | Output directory settings |

### TissueMap InputConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `scan_dir` | `Path` | `Path(".")` | Dataset directory or parent directory containing datasets |
| `tmt_datasets` | `list[str]` | `[]` | Dataset IDs that should be treated as TMT |
| `feature_prefix` | `str \| None` | `None` | Optional custom QPX feature parquet prefix |
| `min_tissue_samples` | `int` | `1` | Minimum number of samples required per tissue label |
| `low_sample_warning_threshold` | `int` | `0` | Warning threshold for low-sample tissues |

### FilteringConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `max_nan_frac` | `float` | `0.95` | Maximum allowed missing fraction per protein |
| `remove_contaminants` | `bool` | `True` | Remove contaminants before downstream analysis |
| `contaminant_pattern` | `str` | `"CONTAM\|ENTRAP\|DECOY"` | Regex used to identify contaminants |

### TissueSpecificityConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `use_pure_mad` | `bool` | `True` | Use MAD-based robust scaling for AdaTiSS |
| `sigma_floor` | `float \| None` | `None` | Manual lower bound for fitted sigma |
| `ts_enriched_threshold` | `float \| None` | `None` | Threshold for tissue-enriched proteins |
| `ts_specific_threshold` | `float \| None` | `None` | Threshold for tissue-specific proteins |
| `ts_housekeeping_threshold` | `float \| None` | `None` | Threshold for housekeeping-like proteins |

### EmbeddingConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `max_nan_frac_for_pca` | `float \| None` | `None` | Optional missingness limit for PCA input proteins |
| `pca_components` | `int` | `50` | Number of PCA components |
| `tsne_perplexity` | `float` | `15.0` | t-SNE perplexity |
| `random_state` | `int` | `42` | Random seed for reproducibility |

### PlottingConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `dpi` | `int` | `250` | Plot resolution |
| `save_pdf` | `bool` | `True` | Save PDF copies of plots |
| `n_marker_top` | `int` | `10` | Number of top markers saved per tissue |

### TissueMap OutputConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | `Path` | `Path("tissuemap_output")` | Output directory for dataset results |

### Example

```python
from pathlib import Path

from mokume.tissuemap.config import InputConfig, OutputConfig, TissueMapConfig

config = TissueMapConfig(
    n_jobs=8,
    input=InputConfig(
        scan_dir=Path("QPX_data/tissues-mq"),
        tmt_datasets=["PXD016999"],
        min_tissue_samples=2,
    ),
    output=OutputConfig(output_dir=Path("./tissuemap_results")),
)
```

---

## Preprocessing Filter Configuration

Filter configurations use YAML format. See [Preprocessing Filters](../concepts/preprocessing.md) for details.

### YAML Structure

```yaml
name: my_filters          # Configuration name
enabled: true             # Master enable/disable

intensity:
  min_intensity: 0.0
  remove_zero_intensity: true
  cv_threshold: null      # null = disabled
  min_replicate_agreement: 1
  quantile_lower: 0.0
  quantile_upper: 1.0

peptide:
  min_peptide_length: 7
  max_peptide_length: 50
  allowed_charge_states: null    # e.g., [2, 3, 4]
  exclude_modifications: []
  max_missed_cleavages: null
  min_search_score: null
  exclude_sequence_patterns: []

protein:
  min_unique_peptides: 2
  remove_contaminants: true
  remove_decoys: true
  contaminant_patterns:
    - CONTAMINANT
    - ENTRAP
    - DECOY
  fdr_threshold: 0.01
  min_coverage: 0.0
  razor_peptide_handling: keep

run_qc:
  min_total_intensity: 0.0
  min_identified_features: 0
  max_missing_rate: 1.0
  min_sample_correlation: null
```

### Pre-configured Templates

Available in `tests/example/filters/`:

| Template | Use Case |
|----------|----------|
| `basic_qc.yaml` | Standard experiments |
| `stringent_filtering.yaml` | Publication-quality |
| `tmt_labeling.yaml` | TMT/iTRAQ |
| `dia_analysis.yaml` | DIA workflows |
| `exploratory_analysis.yaml` | Data exploration |
