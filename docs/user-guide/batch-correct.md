# correct-batches: Batch Correction

The `correct-batches` command applies native Rust ComBat batch correction (oracle-verified vs inmoose) to already-quantified protein data. It reads multiple TSV files from a folder, combines them, and removes batch effects. Because ComBat is native in the kernel, no extra dependency is needed.

!!! tip "Prefer the integrated pipeline"
    For most use cases, batch correction is easier to apply via `features2proteins --batch-correction`. Use this standalone command when you have pre-existing protein quantification files that need correction.

This page documents the `mokume correct-batches` CLI subcommand.

The combined long table must contain one finite piBAQ value for every
protein × sample cell. Structural gaps, blank values, `NaN`, and infinities are
rejected with examples of the affected cells; Mokume never turns them into
zero silently. Impute or otherwise resolve missing values explicitly before
running this command. An explicit numeric zero remains a valid observed value.

## Basic Usage

=== "CLI"

    ```bash
    mokume correct-batches \
        -f pibaq_folder/ \
        -p "*pibaq.tsv" \
        -o corrected_pibaq.tsv
    ```

=== "Python (wheel)"

    The wheel wrapper validates documented keyword arguments, maps them to the
    command's exact CLI flags, and runs the same kernel in-process:

    ```python
    import mokume

    mokume.correct_batches(
        folder="pibaq_folder/",
        pattern="*pibaq.tsv",
        output="corrected_pibaq.tsv",
    )
    ```

=== "Python (explicit argv)"

    ```python
    import mokume

    mokume.run([
        "correct-batches",
        "--folder", "pibaq_folder/",
        "--pattern", "*pibaq.tsv",
        "--output", "corrected_pibaq.tsv",
    ])
    ```

## CLI Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--folder` | required | Folder containing TSV files |
| `-p/--pattern` | `*pibaq.tsv` | File matching pattern |
| `-o/--output` | required | Output file path |
| `--sample_id_column` / `--sid` | `SampleID` | Sample ID column name |
| `--protein_id_column` / `--pid` | `ProteinName` | Protein ID column name |
| `--pibaq_raw_column` / `--pibaq` | `PiBAQ` | Raw intensity column |
| `--pibaq_corrected_column` | `PiBAQBec` | Corrected intensity column |
| `--comment` | `#` | Comment character in files |
| `--sep` | `\t` | Field separator |
| `--export_anndata` | off | Export to AnnData h5ad format |

## With Covariates

To preserve biological signal during batch correction, supply covariates. The
standalone `correct-batches` command runs ComBat on the combined piBAQ folder and
does **not** expose batch-method or covariate options. Covariate-aware correction
is driven from the [`features2proteins`](features2proteins.md#batch-correction)
flow, which extracts the covariates from the SDRF:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --batch-correction \
    --batch-method sample_prefix \
    --batch-covariates "characteristics[sex],characteristics[tissue]"
```

```python
import mokume

mokume.features2proteins(
    parquet="features.parquet",
    output="proteins.csv",
    sdrf="experiment.sdrf.tsv",
    quant_method="maxlfq",
    batch_correction=True,
    batch_method="sample_prefix",
    batch_covariates="characteristics[sex],characteristics[tissue]",
)
```

!!! warning
    Without covariates, batch correction may remove biological signal that correlates with batches. See [Batch Correction concepts](../concepts/batch-correction.md) for details.

## AnnData Export

Export corrected data to AnnData format for downstream analysis with scanpy or other single-cell/proteomics tools:

```bash
mokume correct-batches \
    -f pibaq_folder/ \
    -p "*pibaq.tsv" \
    -o corrected_pibaq.tsv \
    --export_anndata
```

This creates a `.h5ad` file alongside the TSV output.
