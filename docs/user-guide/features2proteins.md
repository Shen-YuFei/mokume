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

## Quantification Methods

| Method | CLI Flag | FASTA Required | Description |
|--------|----------|:--------------:|-------------|
| MaxLFQ | `--quant-method maxlfq` | No | Delayed normalization (default) |
| DirectLFQ | `--quant-method directlfq` | No | Hierarchical alignment (native Rust) |
| iBAQ | `--quant-method ibaq` | Yes | Absolute quantification |
| TopN | `--quant-method topn` | No | Average of N most intense peptides |
| Sum | `--quant-method sum` | No | Sum of all peptides |
| Median | `--quant-method median` | No | Median peptide intensity |
| Ratio | `--quant-method ratio` | No | Log2 sample/reference (TMT) |
| TMT Abundance | `--quant-method abd` | No | Median of log2 peptide intensities (TMT) |
| TMT Reporter Intensity | `--quant-method intensity` | No | Sum of raw reporter intensities (TMT) |
| Spectral Count | `--quant-method spectral_count` | No | Count of distinct peptides per (protein, sample) |

In practice:

- Use `maxlfq` as the default starting point for standard LFQ workflows.
- Use `directlfq` when you explicitly want the DirectLFQ package to handle normalization and quantification together.
- Use `ibaq` when you need absolute-style quantification and have a FASTA file. The pipeline path delegates to the same piBAQ (paralog-aware iBAQ) algorithm as the `peptides2protein` CLI -- see [Quantification Methods → iBAQ](../concepts/quantification.md#ibaq-pibaq-paralog-aware) for the family discovery and fallback semantics; the wide-format pipeline output retains per-protein `Ibaq` values but does not surface the `FamilyId` / `EvidenceLevel` metadata columns.
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

# DirectLFQ (native Rust, no extra dependency)
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq --directlfq-cores 4
```

## Memory & Performance for Large Studies

When the input parquet has thousands of samples (~5000+), the long-form
features must be pivoted into a wide DirectLFQ matrix. mokume streams the
DuckDB result set through Arrow into polars and pivots there, which keeps the
load step's wall time down (cf. PXD030304: ~32 min on 163M long rows pivoted
into 147,374 × 5,798) and avoids the OOM that pandas pivots used to trigger.

The DuckDB engine itself can be size-capped via `--duckdb-memory` /
`--duckdb-threads`:

```bash
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq \
    --duckdb-memory 40GB \
    --duckdb-threads 16
```

!!! warning "`--duckdb-memory` is *not* a hard process cap"
    The flag only sizes DuckDB's internal buffer pool. PyArrow, polars, and
    pandas allocate independently, so peak Python process RSS can grow to
    **2-3x** the DuckDB cap on wide pivots. For production environments
    that need a strict ceiling, layer one of these on top of mokume:

    - **systemd / cgroup**: `systemd-run --scope -p MemoryMax=80G -- mokume features2proteins ...`
    - **SLURM**: `sbatch --mem=80G ...`
    - **Docker / k8s**: `resources.limits.memory: 80Gi`

    The `directlfq-cores` worker count is automatically reduced when
    `--duckdb-memory` is set, so each forked worker has room for its
    COW-amplified copy of the wide matrix.

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
    --sample-normalization hierarchical \
    --normalization-proteins housekeeping.txt

# Quantile / MedianCenter / MeanCenter / RLR / LOESS (dataset-level)
mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization quantile

mokume features2proteins -p data.parquet -o out.csv \
    --sample-normalization loess
```

!!! note "Dataset-level normalizers"
    `quantile`, `mediancenter`, `meancenter`, `rlr`, and `loess` are
    dataset-level normalizers applied after peptide aggregation, all native in the
    Rust kernel.

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

!!! warning "`--batch-method run` is not available in this flow"
    Only `sample_prefix` and `column` (with `--batch-column`) detection are
    supported in the protein-matrix flow. `--batch-method run` has no run-level
    mapping here and errors at runtime; PCA + HDBSCAN outlier removal is not
    ported.

## Differential Expression

Contrasts must be explicitly specified via `--de-contrasts` (inline) or `--de-contrasts-file` (TSV). Both can be combined.

=== "Inline contrasts"

    ```bash
    mokume features2proteins \
        -p features.parquet -o proteins.csv -s experiment.sdrf.tsv \
        --quant-method maxlfq \
        --de \
        --de-contrasts "NASH vs HL,NASH vs Control" \
        --de-method limrots \
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
| `--de-ensemble-min-k` | 2 | Minimum ensemble members that must agree on direction |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh or ihw |
| `--de-output` | — | DE results file; with multiple contrasts each is written as `<stem>_<A-B>.<ext>` |

!!! warning "Contrasts are required"
    If `--de` is enabled but no contrasts are provided
    (neither `--de-contrasts` nor `--de-contrasts-file`),
    the pipeline raises an error listing available conditions.
    Use `" vs "` as the delimiter to support hyphenated
    condition names.

!!! tip
    `--de-method auto` chooses `deqms` for `directlfq`
    quantification and `limrots` for all others. All methods
    run in the native Rust kernel — no R or rpy2 required.
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
| `--impute-method` | `none` | none, knn, minprob, mindet, qrilc, missforest, seqknn, impseq, gms, bpca, impseqrob |
| `--impute-quantile` | 0.01 | Quantile for MinProb/MinDet/QRILC low-tail draw |
| `--impute-shift` | 1.6 | MinProb shift in standard deviations |
| `--impute-scale` | 0.3 | MinProb scale factor for sigma |
| `--impute-n-neighbors` | 5 | Neighbours for KNN/SeqKNN |

!!! warning "`missforest` is wheel-only"
    Every imputation method above runs in the Rust kernel except `missforest`,
    which wraps scikit-learn's `IterativeImputer`, whose tree-building internals
    are not reproducible cross-language. Passing
    `--impute-method missforest` errors with a pointer to the
    pure-Python fallback. Run it on the protein matrix instead (`analysis` extra):

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
# Export normalized peptides and ions
mokume features2proteins \
    -p features.parquet -o proteins.csv \
    --quant-method directlfq \
    --export-peptides peptides.csv \
    --export-ions ions.csv
```

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
    --remove-contaminants \
    --irs --irs-remove-reference \
    --batch-correction --batch-method sample_prefix \
    --de --de-contrasts "NASH vs HL" --de-method limrots --de-fdr-method ihw \
    --de-output de_results.csv
```
