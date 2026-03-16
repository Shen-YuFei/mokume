# features2peptides: Peptide Normalization

The `features2peptides` command normalizes feature-level mass spectrometry data into peptide intensities. This is the first step of the two-step pipeline, giving you fine-grained control over normalization before protein quantification.

## Basic Usage

=== "CLI"

    ```bash
    mokume features2peptides \
        -p features.parquet \
        -s experiment.sdrf.tsv \
        --nmethod median \
        --pnmethod globalMedian \
        --output peptides.csv
    ```

=== "Python"

    ```python
    from mokume.normalization.peptide import peptide_normalization

    peptide_normalization(
        parquet="features.parquet",
        sdrf="experiment.sdrf.tsv",
        nmethod="median",
        pnmethod="globalMedian",
        output="peptides.csv",
    )
    ```

## Processing Steps

The command performs these steps in order:

1. Parse protein identifiers and retain unique peptides
2. Remove entries with empty intensity or condition
3. Filter peptides by minimum amino acid length
4. Remove low-confidence proteins (< min unique peptides)
5. Optionally remove decoys, contaminants, and specified proteins
6. Normalize at feature level between MS runs (`--nmethod`)
7. Merge peptidoforms across fractions and technical replicates
8. Normalize at sample level (`--pnmethod`)
9. Remove low-frequency peptides (optional)
10. Assemble peptidoforms to peptides
11. Optional log2 transformation

## Normalization Methods

### Feature-Level (`--nmethod`)

| Method | Description |
|--------|-------------|
| `median` | Normalize by median across MS runs (default) |
| `mean` | Normalize by mean across MS runs |
| `iqr` | Normalize by interquartile range |
| `none` | Skip feature normalization |

### Sample-Level (`--pnmethod`)

| Method | Description |
|--------|-------------|
| `globalMedian` | Adjust all samples to global median (default) |
| `conditionMedian` | Adjust samples within each condition |
| `hierarchical` | DirectLFQ-style hierarchical clustering normalization |
| `none` | Skip sample normalization |

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
| `--remove_decoy_contaminants` | off | Remove decoys and contaminants |
| `--remove_low_frequency_peptides` | off | Remove peptides in <20% of samples |
| `--remove_ids` | none | File with protein IDs to exclude |

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

```python
from mokume.normalization.peptide import peptide_normalization

peptide_normalization(
    parquet="features.parquet",
    sdrf="experiment.sdrf.tsv",
    min_aa=7,
    min_unique=2,
    remove_ids=None,
    remove_decoy_contaminants=True,
    remove_low_frequency_peptides=True,
    output="peptides-norm.csv",
    skip_normalization=False,
    nmethod="median",
    pnmethod="globalMedian",
    log2=True,
    save_parquet=False,
)
```

### With Preprocessing Filters

```python
from mokume.normalization.peptide import peptide_normalization
from mokume.model.filters import PreprocessingFilterConfig

config = PreprocessingFilterConfig(name="custom", enabled=True)
config.intensity.min_intensity = 1000.0
config.peptide.allowed_charge_states = [2, 3, 4]
config.protein.min_unique_peptides = 2

peptide_normalization(
    parquet="features.parquet",
    sdrf="experiment.sdrf.tsv",
    output="peptides.csv",
    nmethod="median",
    pnmethod="globalMedian",
    filter_config=config,
)
```
