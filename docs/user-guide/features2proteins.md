# features2proteins: Unified Pipeline

The `features2proteins` command is the recommended way to go from raw feature data to protein quantification. It handles loading, filtering, normalization, quantification, batch correction, IRS, differential expression, and visualization in a single step.

## Basic Usage

=== "CLI"

    ```bash
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        -s experiment.sdrf.tsv \
        --quant-method maxlfq
    ```

=== "Python"

    ```python
    from mokume.pipeline import QuantificationPipeline, PipelineConfig
    from mokume.pipeline.config import (
        InputConfig, QuantificationConfig, NormalizationConfig,
    )

    config = PipelineConfig(
        input=InputConfig(parquet="features.parquet", sdrf="experiment.sdrf.tsv"),
        quantification=QuantificationConfig(method="maxlfq"),
    )
    pipeline = QuantificationPipeline(config)
    proteins = pipeline.run()
    ```

=== "Python (functional)"

    ```python
    from mokume.pipeline import features_to_proteins

    proteins = features_to_proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="maxlfq",
    )
    ```

## Quantification Methods

| Method | CLI Flag | FASTA Required | Description |
|--------|----------|:--------------:|-------------|
| MaxLFQ | `--quant-method maxlfq` | No | Delayed normalization (default) |
| DirectLFQ | `--quant-method directlfq` | No | Hierarchical alignment (requires extra) |
| iBAQ | `--quant-method ibaq` | Yes | Absolute quantification |
| TopN | `--quant-method topn` | No | Average of N most intense peptides |
| Sum | `--quant-method sum` | No | Sum of all peptides |
| Median | `--quant-method median` | No | Median peptide intensity |
| Ratio | `--quant-method ratio` | No | Log2 sample/reference (TMT) |

In practice:

- Use `maxlfq` as the default starting point for standard LFQ workflows.
- Use `directlfq` when you explicitly want the DirectLFQ package to handle normalization and quantification together.
- Use `ibaq` when you need absolute-style quantification and have a FASTA file.
- Use `ratio` for TMT PS-style reference-based analysis.

```bash
# iBAQ (requires FASTA)
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method ibaq --fasta proteome.fasta

# TopN (Top5)
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method topn --topn 5

# DirectLFQ (pip install mokume[directlfq])
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq --directlfq-cores 4
```

## Normalization Options

### Run-Level Normalization

Adjusts for intensity differences between MS runs within each sample.

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --run-normalization median  # median, mean, max, global, max_min, iqr, none
```

### Sample-Level Normalization

Adjusts for systematic differences across samples.

```bash
# Global median (default)
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization globalMedian

# Hierarchical (DirectLFQ-style)
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization hierarchical

# TMM normalization
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization tmm

# With specific normalization proteins
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization hierarchical \
    --normalization-proteins housekeeping.txt
```

- `globalMedian` is the default and a good general-purpose starting point.
- `hierarchical` is useful when you want DirectLFQ-style normalization with a non-DirectLFQ quantification method.
- `tmm` is available for composition-bias-aware sample normalization.

## IRS Normalization (Multi-Plex TMT)

For TMT experiments with shared reference channels across plexes:

```bash
# Auto-detect references from SDRF
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method median \
    --irs --irs-remove-reference

# Explicit reference samples
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method median \
    --irs --irs-reference-samples "p1_11,p2_11"

# Custom regex for reference detection
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --irs --irs-reference-regex "pool|bridge|control"
```

| IRS Option | Default | Description |
|------------|---------|-------------|
| `--irs` | off | Enable IRS normalization |
| `--irs-reference-samples` | auto | Comma-separated reference sample names |
| `--irs-sdrf-column` | auto | SDRF column for reference detection |
| `--irs-sdrf-values` | auto | Values indicating reference samples |
| `--irs-reference-regex` | `pool\|powder\|ref\|reference\|bridge` | Regex for auto-detection |
| `--irs-stat` | `median` | Statistic for plex reference: median or mean |
| `--irs-remove-reference` | off | Remove reference samples from output |

## Ratio Quantification (TMT PS Protocol)

For multi-plex TMT with per-plex reference division:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method ratio \
    --coverage-threshold 0.65 \
    --ratio-fraction-merge mean
```

!!! info
    Ratio quantification handles cross-plex normalization inherently via per-plex reference division. The `--irs` flag is ignored in ratio mode.

## Batch Correction

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --batch-correction \
    --batch-method sample_prefix \
    --batch-covariates "characteristics[sex],characteristics[organism part]"
```

=== "Python"

    ```python
    from mokume.pipeline.config import BatchCorrectionConfig

    config = PipelineConfig(
        input=InputConfig(parquet="data.parquet", sdrf="experiment.sdrf.tsv"),
        quantification=QuantificationConfig(method="maxlfq"),
        batch=BatchCorrectionConfig(
            enabled=True,
            method="sample_prefix",
            covariates=["characteristics[sex]", "characteristics[organism part]"],
        ),
    )
    ```

## Differential Expression

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --de \
    --de-contrasts "NASH-HL,NASH-Control" \
    --de-method limrots \
    --de-log2fc 0.5 \
    --de-fdr 0.05 \
    --de-fdr-method ihw \
    --de-output de_results.csv
```

| DE Option | Default | Description |
|-----------|---------|-------------|
| `--de` | off | Enable differential expression |
| `--de-contrasts` | all pairs | Comma-separated contrasts (e.g., "A-B") |
| `--de-method` | `auto` | Method: auto, limrots, deqms, or proda |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh or ihw |
| `--de-output` | auto | Output file for DE results |

!!! tip
    `--de-method auto` chooses `deqms` for `directlfq` quantification and `limrots` for all others. Use `proda` explicitly when dropout-aware modeling is more appropriate for your matrix. See [Differential Expression concepts](../concepts/differential-expression.md) for a detailed comparison of methods.

## Plots and Reports

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --de --de-contrasts "NASH-HL" \
    --plot-dir plots/ \
    --plot-volcano --plot-heatmap --plot-pca \
    --highlight-genes "COL10A1,FN1,ALB" \
    --interactive-report --report-output qc_report.html
```

## Exporting Intermediate Data

```bash
# Export normalized peptides and ions
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq \
    --export-peptides peptides.csv \
    --export-ions ions.csv
```

## Full Example

A complete TMT multi-plex analysis:

```bash
mokume features2proteins \
    -p features.parquet \
    -o proteins.csv \
    -s experiment.sdrf.tsv \
    --quant-method median \
    --run-normalization median \
    --sample-normalization globalMedian \
    --min-unique 2 \
    --remove-contaminants \
    --irs --irs-remove-reference \
    --batch-correction --batch-method sample_prefix \
    --de --de-contrasts "NASH-HL" --de-method limrots --de-fdr-method ihw \
    --plot-dir plots/ --plot-volcano --plot-pca \
    --interactive-report --report-output qc_report.html
```
