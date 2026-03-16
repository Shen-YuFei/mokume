# Preprocessing Filters

mokume provides a comprehensive filter system for quality control, configurable via YAML/JSON files or CLI options.

## Filter Categories

### Intensity Filters

| Filter | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| MinIntensityFilter | `min_intensity` | 0.0 | Remove features below threshold |
| CVThresholdFilter | `cv_threshold` | null | Max CV across replicates |
| ReplicateAgreementFilter | `min_replicate_agreement` | 1 | Min replicates with detection |
| QuantileFilter | `quantile_lower/upper` | 0.0/1.0 | Remove intensity outliers |

### Peptide Filters

| Filter | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| PeptideLengthFilter | `min/max_peptide_length` | 7/50 | Peptide length range |
| ChargeStateFilter | `allowed_charge_states` | null | Allowed charges (e.g., [2,3,4]) |
| ModificationFilter | `exclude_modifications` | [] | Remove specific modifications |
| MissedCleavageFilter | `max_missed_cleavages` | null | Max missed cleavages |
| SearchScoreFilter | `min_search_score` | null | Min search engine score |
| SequencePatternFilter | `exclude_sequence_patterns` | [] | Regex patterns to exclude |

### Protein Filters

| Filter | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| ContaminantFilter | `remove_contaminants/decoys` | true | Remove contaminants/decoys |
| MinPeptideFilter | `min_unique_peptides` | 2 | Min unique peptides per protein |
| ProteinFDRFilter | `fdr_threshold` | 0.01 | Protein-level FDR |
| CoverageFilter | `min_coverage` | 0.0 | Min sequence coverage |
| RazorPeptideFilter | `razor_peptide_handling` | "keep" | Handle shared peptides |

### Run/Sample QC Filters

| Filter | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| RunIntensityFilter | `min_total_intensity` | 0.0 | Min total intensity per run |
| MinFeaturesFilter | `min_identified_features` | 0 | Min features per run |
| MissingRateFilter | `max_missing_rate` | 1.0 | Max missing value rate |
| SampleCorrelationFilter | `min_sample_correlation` | null | Min replicate correlation |

## Configuration

### YAML Configuration File

Generate an example configuration:

```bash
mokume features2peptides --generate-filter-config filters.yaml
```

Example `basic_qc.yaml`:

```yaml
name: basic_qc
enabled: true

intensity:
  remove_zero_intensity: true

peptide:
  min_peptide_length: 7
  max_peptide_length: 50

protein:
  min_unique_peptides: 2
  remove_contaminants: true
  remove_decoys: true
  contaminant_patterns:
    - CONTAMINANT
    - ENTRAP
    - DECOY
```

### Using Filter Configurations

```bash
# From config file
mokume features2peptides \
    -p features.parquet -s experiment.sdrf.tsv \
    --filter-config filters.yaml \
    --output peptides.csv

# CLI overrides (take precedence over config file)
mokume features2peptides \
    -p features.parquet -s experiment.sdrf.tsv \
    --filter-config filters.yaml \
    --filter-min-intensity 1000 \
    --filter-cv-threshold 0.3 \
    --output peptides.csv

# CLI-only filtering (no config file)
mokume features2peptides \
    -p features.parquet -s experiment.sdrf.tsv \
    --filter-min-intensity 500 \
    --filter-min-unique-peptides 2 \
    --output peptides.csv
```

### Python API

```python
from mokume.preprocessing.filters import (
    load_filter_config,
    get_filter_pipeline,
    generate_example_config,
)
from mokume.model.filters import PreprocessingFilterConfig

# Generate example
generate_example_config("filters.yaml")

# Load from file
config = load_filter_config("filters.yaml")

# Or create programmatically
config = PreprocessingFilterConfig(name="custom", enabled=True)
config.intensity.min_intensity = 1000.0
config.peptide.allowed_charge_states = [2, 3, 4]
config.protein.min_unique_peptides = 2

# Apply filters
pipeline = get_filter_pipeline(config)
filtered_df, results = pipeline.apply(df)

# Check results
for result in results:
    print(f"{result.filter_name}: removed {result.removed_count} ({result.removal_rate:.1%})")
```

## Pre-configured Templates

mokume includes templates for common scenarios in `tests/example/filters/`:

| Configuration | Use Case | Description |
|---------------|----------|-------------|
| `basic_qc.yaml` | General QC | Minimal filtering for standard experiments |
| `stringent_filtering.yaml` | Publication | High-confidence results with strict thresholds |
| `tmt_labeling.yaml` | TMT/iTRAQ | Optimized for multiplexed labeling |
| `dia_analysis.yaml` | DIA | Optimized for DIA-NN, Spectronaut |
| `exploratory_analysis.yaml` | Exploration | Minimal filtering for data exploration |
