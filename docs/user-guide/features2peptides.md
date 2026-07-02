# features2peptides: Peptide Normalization

The `features2peptides` command normalizes feature-level mass spectrometry data into peptide intensities. This is the first step of the two-step pipeline, giving you fine-grained control over normalization before protein quantification.

## Basic Usage

=== "CLI"

    ```bash
    mokume features2peptides \
        -p features.parquet \
        -s experiment.sdrf.tsv \
        --run-normalization median \
        --sample-normalization globalMedian \
        --output peptides.csv
    ```

=== "Python (wheel)"

    The wheel wrapper maps keyword arguments to CLI flags (`key=value` → `--key value` with `_` rewritten to `-`; `key=True` → `--key`) and runs the same kernel in-process:

    ```python
    import mokume

    mokume.features2peptides(
        parquet="features.parquet",
        sdrf="experiment.sdrf.tsv",
        run_normalization="median",
        sample_normalization="globalmedian",
        output="peptides.csv",
    )
    ```

=== "Python (explicit argv)"

    ```python
    import mokume

    mokume.run([
        "features2peptides",
        "--parquet", "features.parquet",
        "--sdrf", "experiment.sdrf.tsv",
        "--run-normalization", "median",
        "--sample-normalization", "globalmedian",
        "--output", "peptides.csv",
    ])
    ```

## Processing Steps

The command performs these steps in order:

1. Parse protein identifiers and retain unique peptides
2. Remove entries with empty intensity or condition
3. Filter peptides by minimum amino acid length
4. Remove low-confidence proteins (< min unique peptides)
5. Optionally remove decoys, contaminants, and specified proteins
6. Normalize at feature level between MS runs (`--run-normalization`)
7. Merge peptidoforms across fractions and technical replicates
8. Normalize at sample level (`--sample-normalization`)
9. Remove low-frequency peptides (optional)
10. Assemble peptidoforms to peptides
11. Optional log2 transformation

## Normalization Methods

### Feature-Level (`--run-normalization`)

| Method | Description |
|--------|-------------|
| `median` | Normalize by median across MS runs (default) |
| `mean` | Normalize by mean across MS runs |
| `max` | Normalize by the maximum intensity within each run |
| `global` | Normalize by total intensity within each run |
| `max_min` | Apply min-max scaling |
| `iqr` | Normalize by interquartile range |
| `none` | Skip feature normalization |

### Sample-Level (`--sample-normalization`)

| Method | Description |
|--------|-------------|
| `globalMedian` | Adjust all samples to global median (default) |
| `conditionMedian` | Adjust samples within each condition |
| `none` | Skip sample normalization |

!!! note "Only scalar-per-sample methods change the peptide output"
    In the peptide flow only the factor-based normalizers (`globalmedian` /
    `conditionmedian`) are applied. The dataset-level methods
    (`quantile`, `rlr`, `loess`, `hierarchical`,
    `mediancenter`, `meancenter`) are accepted but are a deterministic **no-op**
    here (same result as `--sample-normalization none`): they need the full
    matrix, which the streaming peptide pass does not hold, and Python's
    per-sample loop also leaves them unchanged. All of these methods **are** implemented and
    oracle-verified in [`features2proteins`](features2proteins.md#normalization-options),
    where the full matrix exists — run dataset-level normalization there.

## Filtering Options

```bash
mokume features2peptides \
    -p features.parquet \
    -s experiment.sdrf.tsv \
    --min_aa 7 \
    --min_unique 2 \
    --remove_decoy_contaminants \
    --remove_low_frequency_peptides \
    --output peptides.csv
```

| Option | Default | Description |
|--------|---------|-------------|
| `--min_aa` | 7 | Minimum amino acid length |
| `--min_unique` | 2 | Minimum unique peptides per protein |
| `--keep-shared-peptides` | off | Keep shared/non-unique peptides and skip the unique-peptide gate |
| `--remove_decoy_contaminants` | off | Remove decoys and contaminants |
| `--remove_low_frequency_peptides` | off | Remove peptides in <20% of samples |
| `--remove_ids` | none | File with protein IDs to exclude |

## TMT / ITRAQ Options

For labeled datasets, `features2peptides` also supports IRS-style scaling and control over aggregation level:

| Option | Default | Description |
|--------|---------|-------------|
| `--irs_channel` | none | Explicit pooled/reference channel label |
| `--irs_autodetect_regex` | none | Regex to detect pooled samples from SDRF |
| `--irs_stat` | `median` | IRS per-run statistic: median or mean |
| `--irs_scope` | `global` | IRS scaling scope: global, by_mixture, or two_stage |
| `--aggregation_level` | `sample` | Aggregate intensities at sample or run level |

!!! note "Channel-based IRS"
    The channel IRS path (`--irs_channel` / `--irs_autodetect_regex`) scales on
    the TMT `mixture` / `channel` columns and is implemented for all three scopes
    (`--irs_scope global` / `by_mixture` / `two_stage`). The reference channel is
    taken from `--irs_channel`, or auto-detected from the SDRF when
    `--irs_autodetect_regex` is given. For cross-plex reference scaling driven by
    the SDRF, see
    [`features2proteins`](features2proteins.md#irs-normalization-multi-plex-tmt).

## Preprocessing Filters

For more advanced filtering, use a YAML/JSON configuration file:

```bash
# Generate example configuration
mokume features2peptides --generate-filter-config filters.yaml

# Use filter configuration
mokume features2peptides \
    -p features.parquet \
    -s experiment.sdrf.tsv \
    --filter-config filters.yaml \
    --output peptides.csv

# CLI overrides (take precedence over config file)
mokume features2peptides \
    -p features.parquet \
    -s experiment.sdrf.tsv \
    --filter-config filters.yaml \
    --filter-min-intensity 1000 \
    --filter-cv-threshold 0.3 \
    --filter-charge-states "2,3,4" \
    --output peptides.csv
```

### CLI Filter Overrides

| Option | Description |
|--------|-------------|
| `--filter-min-intensity` | Minimum intensity threshold |
| `--filter-cv-threshold` | Maximum CV across replicates |
| `--filter-charge-states` | Comma-separated allowed charge states |
| `--filter-max-missed-cleavages` | Maximum missed cleavages |
| `--filter-exclude-modifications` | Comma-separated modifications to exclude |
| `--filter-min-unique-peptides` | Minimum unique peptides per protein |
| `--filter-min-features` | Minimum identified features per run |
| `--filter-max-missing-rate` | Maximum missing value rate (0.0-1.0) |

!!! note "Group-level filters: what runs and what is a no-op"
    The per-row filters (min-intensity floor, peptide length, charge states,
    excluded modifications, missed cleavages) and the per-`(protein, sample)`
    unique-peptide gate are wired in the kernel and oracle-locked vs Python.
    Among the group-level filters, CV threshold (`--filter-cv-threshold`),
    quantile outlier removal, and the run-QC checks `--filter-min-features` /
    min-total-intensity / min-proteins are implemented via a pre-pass. Replicate
    agreement reproduces Python's degenerate per-sample behaviour (a threshold
    `>= 2` empties the output, matching the reference). `--filter-max-missing-rate`,
    sample correlation, min-search-score, and min-coverage are no-ops on the QPX
    streaming model — each warns and passes rows through, exactly as Python skips
    them. Only an unknown `razor-peptide-handling` value returns `NotImplemented`.

See [Preprocessing Filters](../concepts/preprocessing.md) for the full filter reference.

## Output Options

```bash
# Standard CSV output
mokume features2peptides -p data.parquet -o peptides.csv

# Parquet output
mokume features2peptides -p data.parquet -o peptides.csv --save_parquet

# Log2 transform
mokume features2peptides -p data.parquet -o peptides.csv --log2

# Skip normalization entirely
mokume features2peptides -p data.parquet -o peptides.csv --skip_normalization
```

## Python API

The wheel exposes the same command as a thin wrapper; keyword arguments map to CLI
flags (`key=value` → `--key value` with `_` rewritten to `-`; `key=True` → `--key`).

```python
import mokume

mokume.features2peptides(
    parquet="features.parquet",
    sdrf="experiment.sdrf.tsv",
    min_aa=7,
    min_unique=2,
    remove_decoy_contaminants=True,
    remove_low_frequency_peptides=True,
    output="peptides-norm.csv",
    run_normalization="median",
    sample_normalization="globalmedian",
    log2=True,
    save_parquet=True,
)
```

### With Preprocessing Filters

Point the wrapper at a YAML/JSON filter config, or pass the per-filter overrides
directly:

```python
import mokume

mokume.features2peptides(
    parquet="features.parquet",
    sdrf="experiment.sdrf.tsv",
    output="peptides.csv",
    run_normalization="median",
    sample_normalization="globalmedian",
    filter_config="filters.yaml",
    filter_min_intensity=1000,
    filter_charge_states="2,3,4",
    filter_min_unique_peptides=2,
)
```
