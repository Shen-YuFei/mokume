# Preprocessing Filters

mokume provides a configurable filter system for quality control, driven by YAML/JSON files or CLI options on the `features2peptides` command.

## QPX and SDRF sample identity

When an SDRF is supplied for the current QPX list-of-label/intensity schema, every positive intensity must resolve to exactly one SDRF record before preprocessing starts. Run keys use the filename basename, ignore ASCII case, and remove one trailing `.raw`, `.mzml`, `.d`, or `.wiff` extension. Labels ignore surrounding whitespace and ASCII case, and SDRF controlled-vocabulary values use their `NT=` term. `LFQ`, `label-free`, and `label free sample` are the only equivalent label aliases.

Label-free QPX entries whose intensity labels are run filenames use each intensity label as its owning run. Isobaric entries use the row run together with the reporter label. Duplicate normalized SDRF keys, unmatched runs, missing reporter channels, and ambiguous run-only matches are errors; mokume does not fall back to the first SDRF row or a run name when an SDRF mapping is invalid.

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

!!! note "Group-level filters: what runs and what is a no-op"
    The per-row filters (min-intensity floor, peptide length, charge states,
    excluded modifications, missed cleavages) and the per-`(protein, sample)`
    unique-peptide gate are wired and oracle-locked in the Rust kernel. Among the
    **group-level** filters, CV threshold, quantile outlier removal, and the
    run-QC checks (min-features, min-total-intensity, min-proteins) are
    implemented via a pre-pass that applies them before the normalization median.
    Replicate agreement reproduces Python's degenerate per-sample behaviour (a
    `>= 2` threshold empties the output). `max-missing-rate`, `sample-correlation`,
    `min-search-score`, and `min-coverage` are no-ops on QPX inputs — each warns
    and passes rows through, matching Python's per-sample / column-absent skip.
    Only an unknown `razor-peptide-handling` value returns `NotImplemented`.

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

### Python (wheel)

The filter config is built and applied entirely inside the kernel; the wheel
wrapper passes the same flags as the CLI:

```python
import mokume

# Generate an example config
mokume.features2peptides(generate_filter_config="filters.yaml")

# Run with a config file plus CLI overrides
mokume.features2peptides(
    parquet="features.parquet",
    sdrf="experiment.sdrf.tsv",
    filter_config="filters.yaml",
    filter_min_intensity=1000,
    filter_min_unique_peptides=2,
    output="peptides.csv",
)
```

Charge states and excluded modifications are passed as comma-separated strings
(`filter_charge_states="2,3,4"`), matching the CLI flags.

## Pre-configured Templates

mokume includes templates for common scenarios in `tests/example/filters/`:

| Configuration | Use Case | Description |
|---------------|----------|-------------|
| `basic_qc.yaml` | General QC | Minimal filtering for standard experiments |
| `stringent_filtering.yaml` | Publication | High-confidence results with strict thresholds |
| `tmt_labeling.yaml` | TMT/iTRAQ | Optimized for multiplexed labeling |
| `dia_analysis.yaml` | DIA | Optimized for DIA-NN, Spectronaut |
| `exploratory_analysis.yaml` | Exploration | Minimal filtering for data exploration |
