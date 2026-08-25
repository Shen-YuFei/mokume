# features2proteins: Unified Pipeline

The `features2proteins` command is the recommended way to go from raw feature data to protein quantification. It handles loading, filtering, normalization, quantification, batch correction, IRS, imputation, and differential expression in a single step. The compute is single-sourced in the Rust kernel; plotting and interactive reports are a separate wheel periphery step (see [Plots and Reports](#plots-and-reports)).

![The features2proteins pipeline stages: source data, quantify, normalize, impute, batch-correct, differential expression](../assets/pipeline.svg){ width="100%" }

## Basic Usage

=== "CLI"

    ```bash
    mokume features2proteins \
        -p features.parquet \
        -o proteins.csv \
        -s experiment.sdrf.tsv \
        --quant-method maxlfq
    ```

=== "Python (wheel)"

    The wheel wrapper maps keyword arguments to CLI flags (`key=value` → `--key value` with `_` rewritten to `-`; `key=True` → `--key`) and runs the same kernel in-process:

    ```python
    import mokume

    mokume.features2proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="maxlfq",
    )
    ```

=== "Python (explicit argv)"

    ```python
    import mokume

    mokume.run([
        "features2proteins",
        "--parquet", "features.parquet",
        "--output", "proteins.csv",
        "--sdrf", "experiment.sdrf.tsv",
        "--quant-method", "maxlfq",
    ])
    ```

## Input Formats

Provide exactly one feature input:

- `--parquet` accepts a quantms.io/QPX feature parquet file.
- `--msstats` accepts a native MSstats CSV and requires `--sdrf` so runs and
  channels can be mapped to samples.

```bash
mokume features2proteins \
    --msstats msstats.csv \
    --sdrf experiment.sdrf.tsv \
    --output proteins.csv \
    --quant-method maxlfq
```

An MSstats CSV must contain `ProteinName`, `PeptideSequence`, `Intensity`,
`Charge` or `PrecursorCharge`, and `Run` or `Reference`. Isobaric data also
requires `Channel`. `--quant-method ratio` requires PSM-level QPX evidence and
therefore cannot use an MSstats feature table.

## Quantification Methods

| Method | CLI Flag | FASTA Required | Description |
|--------|----------|:--------------:|-------------|
| MaxLFQ | `--quant-method maxlfq` | No | Delayed normalization (default) |
| DirectLFQ | `--quant-method directlfq` | No | Hierarchical alignment (native Rust) |
| piBAQ | `--quant-method pibaq` | Yes | Absolute quantification |
| TopN | `--quant-method top3` / `top5` / any `top<N>` | No | Average of the N most intense peptides |
| Sum | `--quant-method sum` | No | Sum of all peptides |
| Median | `--quant-method median` | No | Median peptide intensity |
| Ratio | `--quant-method ratio` | No | Log2 sample/reference (TMT) |
| TMT Abundance | `--quant-method abd` | No | Median of log2 peptide intensities (TMT) |
| TMT Reporter Intensity | `--quant-method intensity` | No | Sum of raw reporter intensities (TMT) |
| Spectral Count | `--quant-method spectral_count` | No | Count of distinct peptides per (protein, sample) |

In practice:

- Use `maxlfq` as the default starting point for standard LFQ workflows.
- Use `directlfq` when you explicitly want the native Rust DirectLFQ estimator
  to handle normalization and quantification together.
- Use `pibaq` when you need absolute-style quantification and have a FASTA
  file. The pipeline delegates to the same piBAQ algorithm as
  `peptides2protein` -- see
  [Quantification Methods → piBAQ](../concepts/quantification.md#pibaq-paralog-aware-ibaq)
  for family discovery and exact shared-peptide allocation. The wide-format
  pipeline output retains per-protein `PiBAQ` values but does not surface the
  `FamilyId` / `EvidenceLevel` metadata columns.
- Use `ratio` for TMT PS-style reference-based analysis.
- Use `top<N>` for the classic Top3-style summary; `top3` is the method from
  Silva et al. 2006, and any other N works the same way (`top5`, `top10`, ...).

```bash
# piBAQ (requires FASTA)
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method pibaq --fasta proteome.fasta

# TopN — the N lives in the method name (top3, top5, top10, ...)
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method top5

# DirectLFQ (native Rust, no extra dependency)
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq --directlfq-cores 4
```

## Memory & Performance for Large Studies

The Rust kernel reads QPX parquet data in Arrow record batches and accumulates
the compact peptide, protein, and sample structures needed by the selected
method. It does not load the input through DuckDB or build a pandas pivot.

Use `--threads` to size the Rayon thread pool used by parallel Rust sections.
`--threads` and the DirectLFQ-specific `--directlfq-cores` are mutually
exclusive, so the requested worker count can never be shadowed:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq \
    --threads 24 \
    --memory 1GB
```

On Linux, `--memory` sets a soft process RSS budget. It reduces the QPX Arrow
batch size, disables read-ahead, and checks RSS after input batches and major
pipeline phases. If in-memory aggregation state cannot fit, Mokume exits with a
clear budget-exceeded error rather than silently ignoring the option. The guard
can observe a transient overshoot only at the next checkpoint, and smaller
synchronous batches may reduce throughput. Use an external cgroup, scheduler,
or container limit when a hard ceiling is required.

`--duckdb-memory` is intentionally not an alias: the Rust path does not use
DuckDB. `--duckdb-threads` remains an alias for the effective Rayon `--threads`
setting. Runtime pyOpenMS FASTA digestion for piBAQ occurs before the Rust
pipeline starts, so it is not covered by `--memory`.

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

# With specific normalization proteins
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization globalmedian \
    --normalization-proteins housekeeping.txt

# Quantile / MedianCenter / MeanCenter / RLR / LOESS / TMM (dataset-level)
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization quantile

mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization loess

mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization tmm
```

!!! note "Dataset-level normalizers"
    `quantile`, `mediancenter`, `meancenter`, `rlr`, `loess`, and `tmm` are
    dataset-level normalizers applied after peptide aggregation, all native in
    the Rust kernel. `tmm` (Trimmed Mean of M-values,
    `mokume.normalization.tmm.TMMNormalizer`) is robust to composition bias from
    highly abundant proteins.

    MaxLFQ and piBAQ currently accept `quantile` as their dataset-level method;
    requesting RLR, LOESS, hierarchical, centering, or TMM with either method
    is rejected before input loading. For MaxLFQ + quantile, Mokume selects the
    built-in MaxLFQ path so the requested normalization is actually applied.

- `globalMedian` is the default and a good general-purpose starting point.
- `hierarchical` is useful when you want DirectLFQ-style normalization with a non-DirectLFQ quantification method.

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
| `--irs-reference-sample` | none | Repeatable single reference sample name; conflicts with `--irs-reference-samples` |
| `--irs-sdrf-column` | auto | SDRF column for reference detection |
| `--irs-sdrf-values` | auto | Values indicating reference samples |
| `--irs-reference-regex` | `pool\|powder\|ref\|reference\|bridge` | Regex for auto-detection |
| `--irs-stat` | `median` | Statistic for plex reference: median or mean |
| `--irs-remove-reference` | off | Remove reference samples from output |

Every IRS sub-option requires `--irs` and an SDRF. If reference detection finds
no usable sample/plex mapping or no finite scale, the command fails rather than
returning an unscaled matrix.

## Sample Correlation QC

Use `--min-sample-correlation` to remove a sample when its mean Pearson
correlation to the other samples in the same SDRF condition falls below a
threshold:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --min-sample-correlation 0.8
```

The comparison is one-shot on the normalized protein matrix: positive finite
linear intensities are log2-transformed, while the already-log2 `abd` and
`ratio` outputs are used directly. Each pair uses its shared proteins, and all
sample scores are computed before any sample is removed. Conditions with fewer
than two samples and pairs with fewer than three usable proteins are rejected
because the requested correlation cannot be evaluated. Pooled/powder reference
conditions are retained but are not scored.

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
    Ratio quantification handles cross-plex normalization inherently via
    per-plex reference division. Combining it with `--irs` is rejected.

## Batch Correction

ComBat runs as a native Rust kernel (oracle-verified vs inmoose), so no extra dependency is needed.

=== "CLI"

    ```bash
    mokume features2proteins \
        -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
        --quant-method maxlfq \
        --batch-correction \
        --batch-method sample_prefix \
        --batch-covariates "characteristics[sex],characteristics[organism part]"
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.features2proteins(
        parquet="features.parquet",
        output="proteins.csv",
        sdrf="experiment.sdrf.tsv",
        quant_method="maxlfq",
        batch_correction=True,
        batch_method="sample_prefix",
        batch_covariates="characteristics[sex],characteristics[organism part]",
    )
    ```

Only `sample_prefix` and `column` (with `--batch-column`) are exposed in this
protein-matrix flow. ComBat requires at least two batches with two samples each
and at least one protein observed in every sample; otherwise it fails instead
of returning unchanged values. PCA + HDBSCAN outlier removal is not ported.

## Differential Expression

Contrasts must be explicitly specified via `--de-contrasts` (inline) or `--de-contrasts-file` (TSV). Both can be combined.

=== "Inline contrasts"

    ```bash
    mokume features2proteins \
        -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
        --quant-method maxlfq \
        --de \
        --de-contrasts "NASH vs HL,NASH vs Control" \
        --de-method deqms \
        --de-fdr-method ihw \
        --de-output de_results.csv
    ```

=== "Contrasts file"

    ```bash
    mokume features2proteins \
        -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
        --quant-method maxlfq \
        --de \
        --de-contrasts-file contrasts.tsv \
        --de-method deqms \
        --de-fdr-method ihw \
        --de-output de_results.csv
    ```

    Where `contrasts.tsv` is a two-column TSV:

    ```
    group1    group2
    NASH      HL
    NASH      Control
    HL        Control
    ```

| DE Option | Default | Description |
|-----------|---------|-------------|
| `--de` | off | Enable differential expression |
| `--de-contrasts` | — | Comma-separated contrasts (e.g., `"A vs B,A vs C"`) |
| `--de-contrasts-file` | — | TSV file with columns `group1`, `group2` |
| `--de-method` | `auto` | Method: auto, limrots, limma, deqms, proda, rots, ensemble |
| `--de-ensemble-methods` | `limrots,deqms,proda` | Comma-separated DE methods used when `--de-method=ensemble` |
| `--de-ensemble-min-k` | 2 | Minimum ensemble members that must agree on direction (ensemble only) |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change, or `auto` for the data-driven mixture gate |
| `--de-effect-size-gate` | — | Explicit data-driven gate: `mixture` or `null_quantile`; a numeric `--de-log2fc` becomes its fallback |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh, ihw, bky, or storey |
| `--de-output` | — | DE results file; with multiple contrasts each is written as `<stem>_<A-B>.<ext>` |

!!! warning "Contrasts are required"
    If `--de` is enabled but no contrasts are provided
    (neither `--de-contrasts` nor `--de-contrasts-file`),
    the pipeline raises an error listing available conditions.
    Use `" vs "` as the delimiter to support hyphenated
    condition names.
    A DE output path is also required, so completed results cannot be discarded.

!!! tip
    `--de-method auto` chooses `deqms` for `directlfq`
    quantification and `limrots` for all others. All methods
    run in the native Rust kernel — no R or rpy2 required.
    BKY and Storey fall back to BH when their pi0 estimate is not reliable.
    ROTS and LimROTS retain their own permutation FDR, so alternative
    `--de-fdr-method` values are rejected for those methods. Ensemble applies
    the selected correction to the combined result.
    See [Differential Expression
    concepts](../concepts/differential-expression.md) for a
    detailed comparison of methods.

```bash
# Top-k consensus across multiple methods
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --de --de-contrasts "NASH vs HL" \
    --de-method ensemble \
    --de-ensemble-methods "limrots,deqms,proda" \
    --de-ensemble-min-k 2 \
    --de-output ensemble_de.csv
```

## Imputation

Proteomics data is sparse. The pipeline can fill missing values on the
protein matrix after coverage filtering and before batch correction.
Imputation runs in log2 space so that censored-aware methods (MinProb,
MinDet, QRILC) behave correctly.

```bash
# MinProb low-tail draw (Perseus style)
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --impute --impute-method minprob \
    --impute-quantile 0.01 --impute-shift 1.6 --impute-scale 0.3

# KNN imputation
mokume features2proteins \
    -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
    --quant-method maxlfq \
    --impute --impute-method knn --impute-n-neighbors 5
```

| Imputation Option | Default | Description |
|-------------------|---------|-------------|
| `--impute` | off | Enable imputation on the protein matrix |
| `--impute-method` | required with `--impute` | mean, median, constant, zero, most_frequent, knn, minprob, mindet, qrilc, seqknn, impseq, gms, bpca, impseqrob |
| `--impute-quantile` | 0.01 | Quantile for MinProb/MinDet |
| `--impute-shift` | 1.6 | MinProb shift in standard deviations |
| `--impute-scale` | 0.3 | MinProb scale factor for sigma |
| `--impute-n-neighbors` | 5 | Neighbours for KNN/SeqKNN |

Imputation tuning options are method-scoped: quantile applies to MinDet/MinProb,
shift and scale to MinProb, and neighbour count to KNN/SeqKNN. Passing them to
another method is rejected.

!!! warning "`missforest` is Python analysis periphery"
    `missforest` is not offered by the Rust compute command because it wraps
    scikit-learn. Install the Rust-backed wheel's `analysis` dependencies and
    call the Python API:

    ```python
    import mokume

    mokume.impute("proteins.csv", method="missforest", output="imputed.csv")
    ```

## Plots and Reports

Plotting and the interactive HTML report are **not** part of the `features2proteins`
CLI — there are no `--plot-*` / `--interactive-report` flags. The kernel writes the
protein matrix and, with `--de-output`, one DE result CSV per contrast; the Python
periphery then reads those CSVs and renders the figures. Install the `plotting`
and/or `reports` extra:

```python
import mokume

# DE plots (volcano / heatmap / PCA) — explicit argv: the per-contrast
# --contrast KEY A B CSV flag repeats, which keyword arguments cannot express.
mokume.de_plots([
    "--protein-matrix", "proteins.csv",
    "--plot-dir", "plots",
    "--sdrf", "experiment.sdrf.tsv",
    "--volcano", "--heatmap", "--pca",
    "--contrast", "NASH-HL", "NASH", "HL", "de_results.csv",
])

# Interactive HTML report from the same kernel CSVs (reports extra).
mokume.interactive_report([
    "--protein-matrix", "proteins.csv",
    "--sdrf", "experiment.sdrf.tsv",
    "--report-output", "qc_report.html",
    "--contrast", "NASH-HL", "NASH", "HL", "de_results.csv",
])
```

The figures read the kernel's output tables, so the cells in the plots match the
cells in the kernel matrix.

## Exporting Intermediate Data

```bash
# Export normalized peptides from a non-DirectLFQ aggregation
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method sum \
    --export-peptides peptides.csv

# Export the normalized ion table used by DirectLFQ
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq \
    --export-ions ions.csv
```

DirectLFQ and Ratio peptide export is not supported because those calculations
operate on method-specific normalized structures. Conversely, `--export-ions` is DirectLFQ-only. For
non-cell-based aggregation methods, peptide export also rejects dataset-level
sample normalization because the exported peptide values would not represent
the normalized protein matrix.

## Full Example

A complete TMT multi-plex analysis (kernel compute; render plots separately with
the wheel periphery as shown in [Plots and Reports](#plots-and-reports)):

```bash
mokume features2proteins \
    -p features.parquet \
    -o proteins.csv \
    -s experiment.sdrf.tsv \
    --quant-method median \
    --run-normalization median \
    --sample-normalization globalMedian \
    --min-unique 2 \
    --irs --irs-remove-reference \
    --batch-correction --batch-method sample_prefix \
    --de --de-contrasts "NASH vs HL" --de-method deqms --de-fdr-method ihw \
    --de-output de_results.csv
```
