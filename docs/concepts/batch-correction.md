# Batch Correction

Batch effects are systematic technical variations introduced during sample processing that can obscure true biological differences.

## What Are Batch Effects?

Common sources include:

- Different processing days or times
- Different instruments or operators
- Different reagent lots
- Multi-site studies with different labs

mokume uses the **ComBat algorithm** (via [inmoose](https://github.com/epigenelabs/inmoose)) to remove batch effects while preserving biological signal.

!!! note "Optional dependency"
    ```bash
    pip install mokume[batch-correction]
    ```

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

```python
from mokume.pipeline import QuantificationPipeline, PipelineConfig
from mokume.pipeline.config import (
    InputConfig, BatchCorrectionConfig, QuantificationConfig,
)

config = PipelineConfig(
    input=InputConfig(parquet="data.parquet", sdrf="experiment.sdrf.tsv"),
    quantification=QuantificationConfig(method="maxlfq"),
    batch=BatchCorrectionConfig(
        enabled=True,
        method="sample_prefix",
        covariates=["characteristics[sex]", "characteristics[organism part]"],
    ),
)

pipeline = QuantificationPipeline(config)
proteins = pipeline.run()  # Returns batch-corrected protein matrix
```

### CLI

```bash
mokume features2proteins \
    -p data.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --batch-correction \
    --batch-method sample_prefix \
    --batch-covariates "characteristics[sex],characteristics[organism part]"
```

### Low-Level API

```python
from mokume.postprocessing import (
    apply_batch_correction,
    detect_batches,
    extract_covariates_from_sdrf,
)

# Detect batches from sample names
batch_indices = detect_batches(
    sample_ids=df_wide.columns.tolist(),
    method="sample_prefix",
)

# Extract covariates (biological signal to preserve)
covariates = extract_covariates_from_sdrf(
    "experiment.sdrf.tsv",
    sample_ids=df_wide.columns.tolist(),
    covariate_columns=["characteristics[sex]"],
)

# Apply ComBat
df_corrected = apply_batch_correction(
    df=df_wide, batch=batch_indices, covs=covariates,
)
```

## Batch Detection Methods

| Method | Description | Example |
|--------|-------------|---------|
| `sample_prefix` | Extract from sample name prefix | `PXD001-S1` &rarr; batch `PXD001` |
| `run` | Use run/reference file name | Each file is a batch |
| `column` | Explicit values from SDRF column | User-specified |

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
