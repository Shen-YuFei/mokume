# Quick Start

This guide shows how to go from raw feature data to protein intensities using mokume.

## Prerequisites

You need:

1. A **parquet file** in quantms.io/qpx format (output from quantms pipeline)
2. Optionally, an **SDRF file** for sample metadata

## One-Step Pipeline (Recommended)

The `features2proteins` command handles everything: loading, filtering, normalization, and quantification.

=== "CLI"

    ```bash
    # MaxLFQ quantification (default)
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        -s experiment.sdrf.tsv

    # With TMT IRS normalization + differential expression
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        -s experiment.sdrf.tsv \
        --quant-method median \
        --irs --irs-remove-reference \
        --de --de-contrasts "NASH-HL" \
        --plot-dir plots/ --plot-volcano --plot-pca

    # DirectLFQ (uses directlfq package for everything)
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        --quant-method directlfq

    # iBAQ (requires FASTA)
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        --quant-method ibaq \
        --fasta proteome.fasta
    ```

=== "Python"

    ```python
    from mokume.pipeline import features_to_proteins

    # Simple MaxLFQ
    proteins = features_to_proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="maxlfq",
    )

    # TMT with IRS + DE
    proteins = features_to_proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="median",
        irs=True,
        irs_remove_reference=True,
        differential_expression=True,
        de_contrasts=["NASH-HL"],
    )
    ```

## Individual Quantification Methods

If you already have normalized peptide data, use the quantification classes directly:

```python
import pandas as pd
from mokume.quantification import (
    TopNQuantification,
    MaxLFQQuantification,
    AllPeptidesQuantification,
    get_quantification_method,
)

peptides = pd.read_csv("peptides.csv")

# Top3 quantification
top3 = TopNQuantification(n=3)
result = top3.quantify(
    peptides,
    protein_column="ProteinName",
    peptide_column="PeptideSequence",
    intensity_column="NormIntensity",
    sample_column="SampleID",
)

# MaxLFQ (auto-uses DirectLFQ if installed)
maxlfq = MaxLFQQuantification(min_peptides=2, threads=4)
result = maxlfq.quantify(peptides)

# Factory function (parses method name)
method = get_quantification_method("top5")  # TopNQuantification(n=5)
result = method.quantify(peptides)
```

## Two-Step Pipeline

For more control, use the peptide normalization step separately:

```bash
# Step 1: Normalize peptides
mokume features2peptides \
    -p features.parquet \
    -s experiment.sdrf.tsv \
    --nmethod median \
    --pnmethod globalMedian \
    --output peptides.csv

# Step 2: Quantify proteins
mokume peptides2protein \
    --method maxlfq \
    -p peptides.csv \
    -o proteins.tsv
```

## What's Next?

- [Quantification Methods](concepts/quantification.md) — understand iBAQ, MaxLFQ, TopN, and more
- [Normalization](concepts/normalization.md) — learn about the normalization pipeline
- [Unified Pipeline](user-guide/features2proteins.md) — full reference for features2proteins
