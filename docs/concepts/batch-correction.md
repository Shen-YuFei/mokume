# Batch Correction

Batch effects are systematic technical variations introduced during sample processing that can obscure true biological differences.

## What Are Batch Effects?

Common sources include:

- Different processing days or times
- Different instruments or operators
- Different reagent lots
- Multi-site studies with different labs

mokume uses the **ComBat algorithm** to remove batch effects while preserving biological signal. ComBat is a native Rust implementation, oracle-verified against [inmoose](https://github.com/epigenelabs/inmoose) (~1e-6 / 1e-9). It needs no third-party Python dependency.

## Key Concepts

| Term | Definition | Example |
|------|------------|---------|
| **Batch** | Technical variation to **remove** | Samples from Lab A vs Lab B |
| **Covariate** | Biological signal to **preserve** | Tissue type, sex, disease status |

!!! warning "Why covariates matter"
    Without covariates, batch correction may accidentally remove biological signal that correlates with batch assignments.

    For example, if all liver samples were processed on Day 1 and all brain samples on Day 2, naive batch correction would remove the tissue-specific signal. By specifying tissue as a covariate, ComBat preserves this biological variation.

    ```
    Without covariates:  Batch effect removed, but tissue signal also reduced
    With covariates:     Batch effect removed, tissue signal preserved
    ```

## Using Batch Correction

### Integrated Pipeline (Recommended)

Run ComBat as part of `features2proteins` with `--batch-correction`. Batches are
detected from the sample-name prefix (or an explicit SDRF column), and
`--batch-covariates` names the SDRF columns whose biological signal to preserve.

=== "CLI"

    ```bash
    mokume features2proteins \
        -p data.parquet -o proteins.csv -s experiment.sdrf.tsv \
        --quant-method maxlfq \
        --batch-correction \
        --batch-method sample_prefix \
        --batch-covariates "characteristics[sex],characteristics[organism part]"
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.features2proteins(
        parquet="data.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="maxlfq",
        batch_correction=True,
        batch_method="sample_prefix",
        batch_covariates="characteristics[sex],characteristics[organism part]",
    )
    ```

The covariate columns are extracted from the SDRF (column match with a
sample-substring fallback, `pd.factorize` encoding, single-value columns
dropped) and fed to the covariate ComBat design. ComBat runs on the proteins
with no missing cells; the rest are kept uncorrected.

### Standalone piBAQ correction

To correct already-written piBAQ tables (e.g. when merging datasets), use the
dedicated `correct-batches` command, which runs the sample-prefix ComBat flow
over a folder of piBAQ TSVs:

=== "CLI"

    ```bash
    mokume correct-batches \
        --folder ./pibaq_outputs --pattern "*pibaq.tsv" \
        --output corrected.tsv
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.correct_batches(folder="./pibaq_outputs", pattern="*pibaq.tsv",
                           output="corrected.tsv")
    ```

## Batch Detection Methods

| Method | Description | Example |
|--------|-------------|---------|
| `sample_prefix` | Extract from sample name prefix | `PXD001-S1` &rarr; batch `PXD001` |
| `column` | Explicit values from SDRF column (`--batch-column`) | User-specified |
| `run` | Use run/reference file name | Each file is a batch |

!!! warning "`--batch-method run` in the protein-matrix flow"
    `--batch-method run` has no run-level mapping in the `features2proteins`
    protein-matrix flow and **errors at runtime** (the same way Python raises
    `run_info required`). Use `sample_prefix` or `column` here. The PCA + HDBSCAN
    outlier-removal pass is unported (HDBSCAN is not reproducible cross-language).

## When to Use Batch Correction

**Recommended scenarios:**

- Combining datasets from multiple studies (e.g., PXD001 + PXD002)
- Samples processed on different days/instruments
- Multi-site studies with different labs

**Requirements:**

- At least 2 samples per batch
- At least 2 batches

## Best Practices

1. **Always specify covariates** when biological groups correlate with batches
2. **Use SDRF characteristics** to identify biological variables to preserve
3. **Apply at protein level** (after quantification) for best results
4. **Verify results** by checking that biological signal (e.g., condition clustering) is preserved
