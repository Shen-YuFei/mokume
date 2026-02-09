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
    DEConfig,
    OutputConfig,
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
| `run_method` | `str` | `"median"` | Run-level normalization: median, mean, iqr, max, max_min, none |
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

### DEConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | `bool` | `False` | Enable DE analysis |
| `contrasts` | `list \| None` | `None` | Contrasts (e.g., ["A-B"]) |
| `method` | `str` | `"ttest"` | Method: ttest or limma |
| `log2fc_threshold` | `float` | `0.5` | Min absolute log2 fold change |
| `fdr_threshold` | `float` | `0.05` | Max FDR |
| `output` | `str \| None` | `None` | Output file for DE results |

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
        method="ttest",
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
