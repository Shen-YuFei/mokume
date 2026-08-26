# Quantification

Aggregate feature-level intensities into a protein x sample matrix. mokume
supports several quantification methods; you pick one with a single option
(`--quant-method` on the CLI, `quant_method=` from the wheel, or
`QuantificationConfig(method=...)` in the package).

The examples below run against the shipped fixture
`python/tests/example/feature_wide.parquet` (500 features, 10 samples) with its
SDRF `python/tests/example/PXD020192.sdrf.tsv`. Every method writes a wide-format
CSV: the first column is `ProteinName`, the remaining columns are one per sample.

## Available methods

| `--quant-method` | Description |
|------------------|-------------|
| `maxlfq` | MaxLFQ pairwise ratio estimation (default) |
| `directlfq` | DirectLFQ intensity traces (native Rust) |
| `top3` | Mean of the 3 most intense peptides (Silva et al. 2006) |
| `top<N>` | Mean of the N most intense peptides — `top5`, `top10`, ... |
| `sum` | Sum of all peptide intensities |
| `median` | Median peptide intensity |

`pibaq`, `ratio`, `abd`, `intensity`, and `peptide-count` are also valid for
feature input. True `spectral-count` pairs `--psm` with its matching feature QPX
via `--parquet` and requires SDRF; piBAQ is covered on the
[Absolute Expression](absolute-expression.md) page. Both count methods use
`none` for run/sample intensity normalization and reject IRS.

## Run a quantification

=== "CLI"

    ```bash
    # MaxLFQ (default) — omit --quant-method to get the same result
    mokume quantify features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_maxlfq.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method maxlfq

    # DirectLFQ (native Rust)
    mokume quantify features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_directlfq.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method directlfq

    # Top3
    mokume quantify features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_top3.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method top3

    # Sum
    mokume quantify features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_sum.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method sum

    # Median
    mokume quantify features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_median.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method median
    ```

=== "Python (wheel)"

    ```python
    import mokume

    # The wheel wrapper maps kwargs to CLI flags (quant_method -> --quant-method)
    # and runs the same Rust kernel in-process.
    for method in ["maxlfq", "directlfq", "top3", "sum", "median"]:
        mokume.features2proteins(
            parquet="python/tests/example/feature_wide.parquet",
            output=f"proteins_{method}.csv",
            sdrf="python/tests/example/PXD020192.sdrf.tsv",
            quant_method=method,
        )
    ```

=== "Python (package)"

    ```python
    from mokume.pipeline.features_to_proteins import QuantificationPipeline
    from mokume.pipeline.config import (
        PipelineConfig,
        InputConfig,
        QuantificationConfig,
    )

    for method in ["maxlfq", "directlfq", "top3", "sum", "median"]:
        config = PipelineConfig(
            input=InputConfig(
                parquet="python/tests/example/feature_wide.parquet",
                sdrf="python/tests/example/PXD020192.sdrf.tsv",
            ),
            quantification=QuantificationConfig(method=method),
        )
        proteins = QuantificationPipeline(config).run()
        print(method, proteins.shape)
        proteins.to_csv(f"proteins_{method}.csv", index=False)
    ```

## What you get

Running the fixture prints an aggregation summary and writes the matrix:

```text
INFO mokume_pipeline: features2proteins aggregation finished
  accepted_features=500 accepted_measurements=500 proteins=117 samples=10
  quant_method=maxlfq output=proteins_maxlfq.csv
```

The exact protein count after filtering depends on the method, because MaxLFQ,
DirectLFQ, and the summary methods differ in how many proteins survive their
minimum-evidence requirements on this small slice:

| Method | Output rows (proteins) | Columns |
|--------|------------------------|---------|
| `maxlfq` | 6 | `ProteinName` + 10 samples |
| `directlfq` | 9 | `ProteinName` + 10 samples |
| `top3` | 11 | `ProteinName` + 10 samples |
| `sum` | 11 | `ProteinName` + 10 samples |
| `median` | 11 | `ProteinName` + 10 samples |

The first two rows of the MaxLFQ output look like this (empty cells are proteins
not observed in that sample):

```text
ProteinName,PXD020192-Sample-10,PXD020192-Sample-11,...
P09382,,22623148408.93,21556287461.08,...
P09417,2266760958.62,,2266760958.62,...
```

!!! tip "Choosing a method"

    See [Quantification Methods](../concepts/quantification.md) for the trade-offs
    between MaxLFQ, DirectLFQ, TopN, and the summary statistics, and when each is
    appropriate for DDA vs DIA and labelled vs label-free data.
