# CLI Reference

mokume provides six main commands. In practice, most users start with `features2proteins` for quantification workflows or `tissuemap` for tissue atlas workflows.
Use `mokume --help` or `mokume <command> --help` for details.

## features2proteins

The unified pipeline: features to protein quantification in one step.

```bash
mokume features2proteins [OPTIONS]
```

### Required Options

| Option | Description |
|--------|-------------|
| `-p/--parquet` | Input parquet file (quantms.io/qpx format) |
| `-o/--output` | Output protein intensities CSV |

### Input & Filtering

| Option | Default | Description |
|--------|---------|-------------|
| `-s/--sdrf` | none | SDRF file for sample metadata |
| `--min-aa` | 7 | Minimum amino acid length |
| `--min-unique` | 2 | Minimum unique peptides per protein |
| `--remove-contaminants` | on | Remove contaminants and decoys |
| `--keep-contaminants` | off | Keep contaminants and decoys |

### Quantification

| Option | Default | Description |
|--------|---------|-------------|
| `--quant-method` | `maxlfq` | Method: maxlfq, directlfq, ibaq, topn, sum, median, ratio, abd (TMT abundance), intensity (TMT reporter), spectral_count |
| `--fasta` | none | FASTA file (required for iBAQ) |
| `--topn` | 3 | N for TopN quantification |
| `--ion-alignment` | none | Ion alignment: none or hierarchical |
| `--directlfq-cores` | auto | CPU cores for DirectLFQ |
| `--directlfq-min-nonan` | 1 | Min non-NaN values for DirectLFQ |

### Normalization

| Option | Default | Description |
|--------|---------|-------------|
| `--run-normalization` | `median` | Run-level: median, mean, max, global, max_min, iqr, none |
| `--sample-normalization` | `globalMedian` | Sample-level: globalMedian, conditionMedian, hierarchical, tmm, quantile, mediancenter, meancenter, rlr, mbqn, loess, none |
| `--normalization-proteins` | none | File with protein IDs for normalization |

VSN is available as a standalone Python API utility (`mokume.normalization.vsn_normalize`) but not exposed via `--sample-normalization` because its glog2 output is incompatible with the pipeline's downstream linear-scale assumptions.

### IRS (Multi-Plex TMT)

| Option | Default | Description |
|--------|---------|-------------|
| `--irs` | off | Enable IRS normalization |
| `--irs-reference-samples` | auto | Comma-separated reference sample names |
| `--irs-sdrf-column` | auto | SDRF column for reference detection |
| `--irs-sdrf-values` | auto | Values indicating reference samples |
| `--irs-reference-regex` | `pool\|powder\|ref\|reference\|bridge` | Regex for reference auto-detection |
| `--irs-stat` | `median` | Plex reference statistic: median or mean |
| `--irs-remove-reference` | off | Remove reference samples from output |

### Coverage Filter

| Option | Default | Description |
|--------|---------|-------------|
| `--coverage-threshold` | none | Min fraction of non-missing values per condition |

### Ratio Quantification

| Option | Default | Description |
|--------|---------|-------------|
| `--ratio-fraction-merge` | `mean` | Fraction merge strategy: mean or max |

### Batch Correction

| Option | Default | Description |
|--------|---------|-------------|
| `--batch-correction` | off | Enable ComBat batch correction |
| `--batch-method` | `sample_prefix` | Detection method: sample_prefix, run, column |
| `--batch-column` | none | SDRF column used when `--batch-method=column` |
| `--batch-covariates` | none | Comma-separated SDRF columns to preserve |
| `--batch-parametric` / `--batch-nonparametric` | parametric | ComBat estimation mode |
| `--batch-mean-only` | off | Only adjust batch means |
| `--batch-ref` | none | Reference batch ID |

### Imputation

| Option | Default | Description |
|--------|---------|-------------|
| `--impute` | off | Enable missing-value imputation on the protein matrix |
| `--impute-method` | `none` | Method: none, knn, minprob, mindet, qrilc, missforest, seqknn, mle, mice, nbavg, gms, bpca, impseq, impseqrob |
| `--impute-quantile` | 0.01 | Quantile for MinProb/MinDet/QRILC low-tail draw |
| `--impute-shift` | 1.6 | MinProb shift in standard deviations |
| `--impute-scale` | 0.3 | MinProb scale factor for the imputation distribution sigma |
| `--impute-n-neighbors` | 5 | Number of neighbours for KNN/SeqKNN/NBavg |

Imputation runs in log2 space (the matrix is transformed before imputation and back to linear afterwards) so that censored-aware methods like MinProb/MinDet/QRILC behave correctly. The imputation step is applied after coverage filtering and before batch correction in the pipeline.

### Differential Expression

| Option | Default | Description |
|--------|---------|-------------|
| `--de` | off | Enable differential expression analysis |
| `--de-contrasts` | — | Comma-separated contrasts (e.g., `"A vs B,A vs C"`) |
| `--de-contrasts-file` | — | TSV file with columns `group1`, `group2` |
| `--de-method` | `auto` | Method: auto, limrots, limma, deqms, proda, rots, msstats, ensemble |
| `--de-ensemble-methods` | `limrots,deqms,proda` | Comma-separated DE methods used when `--de-method=ensemble` |
| `--de-ensemble-min-k` | 2 | Minimum ensemble members that must agree on direction |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh or ihw |
| `--de-output` | auto | Output file for DE results |

Contrasts must be explicitly provided via `--de-contrasts` and/or `--de-contrasts-file`. Both can be combined.

`--de-method auto` selects `deqms` for `directlfq` quantification and `limrots` for other quantification methods. All methods are pure-Python reimplementations — no R or rpy2 required.

`--de-method ensemble` runs each member method on the same contrast and combines the per-protein verdicts via top-k consensus: a protein is called UP/DOWN only when at least `--de-ensemble-min-k` members agree on direction and the Fisher-combined p-value passes the FDR threshold.

### Plots & Reports

| Option | Default | Description |
|--------|---------|-------------|
| `--plot-dir` | none | Output directory for plots |
| `--plot-volcano` | off | Generate volcano plots |
| `--plot-heatmap` | off | Generate per-contrast DE heatmaps (top 50 by \|log2FC\|) |
| `--plot-pca` | off | Generate PCA plots |
| `--highlight-genes` | none | Comma-separated gene names to highlight |
| `--interactive-report` | off | Generate interactive HTML QC report |
| `--report-output` | auto | Output path for HTML report |

### Export

| Option | Default | Description |
|--------|---------|-------------|
| `--export-peptides` | none | Export normalized peptides to file |
| `--export-ions` | none | Export normalized ions (DirectLFQ only) |

### DuckDB Resource Limits

| Option | Default | Description |
|--------|---------|-------------|
| `--duckdb-memory` | DuckDB autoconfig (~80% of total RAM) | DuckDB memory limit (e.g. `80GB`, `16384MB`). See note below. |
| `--duckdb-threads` | all cores | Number of threads DuckDB may use |

!!! warning "`--duckdb-memory` is not a hard process cap"
    The flag only sizes DuckDB's internal buffer pool. PyArrow, polars, and pandas each
    have their own independent allocators, so the surrounding Python process can grow to
    **2-3x** the DuckDB limit on wide pivots (e.g. PXD030304 at 5798 samples peaks
    ~94 GB of process RSS with `--duckdb-memory 40GB`). For production environments
    that require a hard ceiling, layer one of these on top of mokume:

    - **systemd / cgroup**: `systemd-run --scope -p MemoryMax=80G -- mokume features2proteins ...`
    - **SLURM**: `sbatch --mem=80G ...`
    - **Docker / k8s**: `resources.limits.memory: 80Gi`

---

## features2peptides

Feature-level to peptide-level normalization.

```bash
mokume features2peptides [OPTIONS]
```

### Core Options

| Option | Default | Description |
|--------|---------|-------------|
| `-p/--parquet` | required | Input parquet file |
| `-s/--sdrf` | none | SDRF file for metadata |
| `-o/--output` | required | Output peptide intensity file |
| `--min_aa` | 7 | Minimum amino acid length |
| `--min_unique` | 2 | Minimum unique peptides per protein |
| `--remove_ids` | none | File with protein IDs to exclude |
| `--remove_decoy_contaminants` | off | Remove decoys and contaminants |
| `--remove_low_frequency_peptides` | off | Remove peptides in <20% of samples |

### Normalization

| Option | Default | Description |
|--------|---------|-------------|
| `--run-normalization` | `median` | Feature normalization: median, mean, max, global, max_min, iqr, none |
| `--sample-normalization` | `globalMedian` | Sample normalization: globalMedian, conditionMedian, hierarchical, tmm, none |
| `--skip_normalization` | off | Skip all normalization |
| `--log2` | off | Log2 transform output |
| `--save_parquet` | off | Save output as parquet |

### TMT / ITRAQ

| Option | Default | Description |
|--------|---------|-------------|
| `--irs_channel` | none | Explicit pooled/reference channel label |
| `--irs_autodetect_regex` | none | Regex to detect pooled samples from SDRF |
| `--irs_stat` | `median` | IRS per-run statistic |
| `--irs_scope` | `global` | IRS scaling scope: global, by_mixture, or two_stage |
| `--aggregation_level` | `sample` | Aggregate at sample or run level |

### Filter Configuration

| Option | Default | Description |
|--------|---------|-------------|
| `--filter-config` | none | YAML/JSON filter configuration file |
| `--generate-filter-config` | none | Generate example config and exit |
| `--filter-min-intensity` | none | Min intensity threshold (override) |
| `--filter-cv-threshold` | none | Max CV across replicates (override) |
| `--filter-charge-states` | none | Comma-separated charge states (override) |
| `--filter-max-missed-cleavages` | none | Max missed cleavages (override) |
| `--filter-exclude-modifications` | none | Comma-separated modifications (override) |
| `--filter-min-unique-peptides` | none | Min unique peptides (override) |
| `--filter-min-features` | none | Min features per run (override) |
| `--filter-max-missing-rate` | none | Max missing rate (override) |

---

## peptides2protein

Protein quantification from normalized peptide data.

```bash
mokume peptides2protein [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p/--peptides` | required | Input peptide intensity file |
| `-f/--fasta` | none | FASTA file (required for iBAQ) |
| `--method` | `ibaq` | Method: ibaq, top3, topn, maxlfq, sum, directlfq |
| `-e/--enzyme` | `Trypsin` | Enzyme for in-silico digestion |
| `-n/--normalize` | off | Normalize quantification values |
| `--min_aa` | 7 | Min amino acid length |
| `--max_aa` | 30 | Max amino acid length |
| `-t/--tpa` | off | Calculate TPA (iBAQ only) |
| `-r/--ruler` | off | ProteomicRuler (iBAQ only) |
| `-i/--ploidy` | 2 | Ploidy number |
| `-m/--organism` | `human` | Organism for histone data |
| `-c/--cpc` | 200 | Cellular protein concentration (g/L) |
| `--topn_n` | 3 | N for TopN quantification |
| `--threads` | -1 | Threads for MaxLFQ (-1 = all cores) |
| `--min_nonan` | 1 | Min non-NaN for DirectLFQ |
| `-o/--output` | none | Output file path |
| `--verbose` | off | Print distribution info |
| `--qc_report` | QCprofile.pdf | QC report PDF path |

Use `--method topn --topn_n 5` or `--method topn --topn_n 10` for Top5 or Top10-style quantification. `-o/--output` is effectively required for `ibaq`; for the other methods, omitting it prints the result table to stdout.

---

## correct-batches

Standalone batch correction for pre-quantified data.

```bash
mokume correct-batches [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--folder` | required | Folder with TSV files |
| `-p/--pattern` | `*ibaq.tsv` | File matching pattern |
| `-o/--output` | required | Output file path |
| `-sid/--sample_id_column` | `SampleID` | Sample ID column |
| `-pid/--protein_id_column` | `ProteinName` | Protein ID column |
| `--ibaq_raw_column` | `IBAQ` | Raw intensity column |
| `--ibaq_corrected_column` | `IBAQ_BEC` | Corrected intensity column |
| `--comment` | `#` | Comment character |
| `--sep` | `\t` | Field separator |
| `--export_anndata` | off | Export to AnnData h5ad |

---

## tissuemap

Per-dataset tissue proteome atlas analysis from QPX outputs.

```bash
mokume tissuemap [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--scan-dir` | required unless `--generate-config` | Dataset directory or parent directory containing datasets |
| `--output-dir` | `tissuemap_output` | Output directory for results |
| `--config` | none | YAML configuration file |
| `--generate-config` | none | Generate a default YAML template and exit |
| `--tmt-dataset` | auto | Mark one or more dataset IDs as TMT |
| `--n-jobs` | `8` | Threads for dataset processing and embedding |
| `--dpi` | `250` | Plot resolution override |

!!! note
    Install the optional dependencies first with `pip install mokume[tissuemap]`.

---

## tsne-visualization

t-SNE dimensionality reduction visualization.

```bash
mokume tsne-visualization [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--folder` | required | Folder with protein files |
| `-o/--pattern` | `proteins.tsv` | File matching pattern |
