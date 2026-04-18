# correct-batches: Batch Correction

The `correct-batches` command applies ComBat batch correction to already-quantified protein data. It reads multiple TSV files from a folder, combines them, and removes batch effects.

!!! tip "Prefer the integrated pipeline"
    For most use cases, batch correction is easier to apply via `features2proteins --batch-correction`. Use this standalone command when you have pre-existing protein quantification files that need correction.

This page documents the standalone CLI command `mokume correct-batches`.

## Basic Usage

=== "CLI"

    ```bash
    mokume correct-batches \
        -f ibaq_folder/ \
        -p "*ibaq.tsv" \
        -o corrected_ibaq.tsv
    ```

=== "Python"

    ```python
    from mokume.postprocessing import (
        apply_batch_correction,
        detect_batches,
        extract_covariates_from_sdrf,
        pivot_wider,
    )

    # Reshape to wide format
    df_wide = pivot_wider(
        df, row_name="ProteinName", col_name="SampleID", values="Ibaq"
    )

    # Detect batches from sample names
    batch_indices = detect_batches(
        sample_ids=df_wide.columns.tolist(),
        method="sample_prefix",
    )

    # Apply ComBat
    df_corrected = apply_batch_correction(
        df=df_wide, batch=batch_indices,
    )
    ```

## CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--folder` | required | Folder containing TSV files |
| `-p/--pattern` | `*ibaq.tsv` | File matching pattern |
| `-o/--output` | required | Output file path |
| `-sid/--sample_id_column` | `SampleID` | Sample ID column name |
| `-pid/--protein_id_column` | `ProteinName` | Protein ID column name |
| `-ibaq/--ibaq_raw_column` | `IBAQ` | Raw intensity column |
| `--ibaq_corrected_column` | `IBAQ_BEC` | Corrected intensity column |
| `--comment` | `#` | Comment character in files |
| `--sep` | `\t` | Field separator |
| `--export_anndata` | off | Export to AnnData h5ad format |

## With Covariates (Python API)

To preserve biological signal during batch correction, specify covariates:

```python
from mokume.postprocessing import (
    apply_batch_correction,
    detect_batches,
    extract_covariates_from_sdrf,
)

batch_indices = detect_batches(
    sample_ids=df_wide.columns.tolist(),
    method="sample_prefix",
)

covariates = extract_covariates_from_sdrf(
    "experiment.sdrf.tsv",
    sample_ids=df_wide.columns.tolist(),
    covariate_columns=["characteristics[sex]", "characteristics[tissue]"],
)

df_corrected = apply_batch_correction(
    df=df_wide,
    batch=batch_indices,
    covs=covariates,
)
```

!!! warning
    Without covariates, batch correction may remove biological signal that correlates with batches. See [Batch Correction concepts](../concepts/batch-correction.md) for details.

## AnnData Export

Export corrected data to AnnData format for downstream analysis with scanpy or other single-cell/proteomics tools:

```bash
mokume correct-batches \
    -f ibaq_folder/ \
    -p "*ibaq.tsv" \
    -o corrected_ibaq.tsv \
    --export_anndata
```

This creates a `.h5ad` file alongside the TSV output.
