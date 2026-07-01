# CLI Reference

The `mokume` binary is the Rust compute kernel's command-line entry point. It exposes exactly four compute subcommands — `features2proteins`, `features2peptides`, `peptides2protein`, and `correct-batches` — and shells out to nothing. Most users start with `features2proteins` for quantification workflows.
Use `mokume --help` or `mokume <command> --help` for details.

The same `mokume` command is available two ways: as a standalone Rust binary (built with `cargo`, no Python) and from the `pip install mokume-rs` wheel, which runs the identical kernel in-process through the `mokume._mokume` extension. The flag surface below is single-sourced in Rust and is identical for both.

!!! note "Plotting, tissue maps, and reports live in the Python wheel"
    The visualization periphery — t-SNE, tissue-proteome maps, DE plots, interactive HTML reports, iBAQ QC — is **not** part of this CLI. It ships in the `pip install mokume-rs` wheel as `mokume.tsne_visualization(...)`, `mokume.tissuemap(...)`, `mokume.de_plots([...])`, `mokume.interactive_report([...])`, and `mokume.peptides2protein_qc(...)`. See the [Python API](python-api.md). These commands read the tables the kernel wrote and never recompute the numbers.

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
| `--quant-method` | `maxlfq` | Method: maxlfq, directlfq, ibaq, topn, top3, sum, median, ratio, abd (TMT abundance), intensity (TMT reporter), spectral_count |
| `--fasta` | none | FASTA file (required for iBAQ) |
| `--topn` | 3 | N for TopN quantification |
| `--ion-alignment` | none | Ion alignment: none or hierarchical |
| `--directlfq-cores` | auto | CPU cores for DirectLFQ |
| `--directlfq-min-nonan` | 1 | Min non-NaN values for DirectLFQ |

### Normalization

| Option | Default | Description |
|--------|---------|-------------|
| `--run-normalization` | `median` | Run-level: median, mean, max, global, max_min, iqr, none |
| `--sample-normalization` | `globalMedian` | Sample-level: globalMedian, conditionMedian, hierarchical, quantile, mediancenter, meancenter, rlr, loess, none |
| `--normalization-proteins` | none | File with protein IDs for normalization |

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

ComBat is a native Rust implementation, oracle-verified against inmoose (parametric ~1e-6, covariate / non-parametric / `ref_batch` / `mean_only` paths included). It runs on the proteins with no missing cells; the rest are kept uncorrected.

!!! warning "`--batch-method run` is unsupported in the protein-matrix flow"
    The protein-matrix flow has no run-level mapping, so `--batch-method run` errors at runtime (the same gate Python's `_detect_batch_indices` enforces). Use `sample_prefix` (default) or `column` (with `--batch-column`). The PCA + HDBSCAN outlier-removal pass is not ported.

### Imputation

| Option | Default | Description |
|--------|---------|-------------|
| `--impute` | off | Enable missing-value imputation on the protein matrix |
| `--impute-method` | `none` | Method: none, knn, minprob, mindet, qrilc, missforest, seqknn, impseq, gms, bpca, impseqrob |
| `--impute-quantile` | 0.01 | Quantile for MinProb/MinDet/QRILC low-tail draw |
| `--impute-shift` | 1.6 | MinProb shift in standard deviations |
| `--impute-scale` | 0.3 | MinProb scale factor for the imputation distribution sigma |
| `--impute-n-neighbors` | 5 | Number of neighbours for KNN/SeqKNN |

Imputation runs in log2 space (the matrix is transformed before imputation and back to linear afterwards) so that censored-aware methods like MinProb/MinDet/QRILC behave correctly. The imputation step is applied after coverage filtering and before batch correction in the pipeline.

!!! note "`missforest` is wheel-only"
    Every method above runs in the native Rust kernel except `missforest`, which wraps scikit-learn's `IterativeImputer` (`RandomForestRegressor`) and cannot be reproduced bit-for-bit cross-language. The kernel accepts it but returns a clear error pointing to the wheel's `mokume.impute(matrix, method="missforest")` (the `analysis` extra), which runs the pure-Python imputer.

### Differential Expression

| Option | Default | Description |
|--------|---------|-------------|
| `--de` | off | Enable differential expression analysis |
| `--de-contrasts` | — | Comma-separated contrasts (e.g., `"A vs B,A vs C"`) |
| `--de-contrasts-file` | — | TSV file with columns `group1`, `group2` |
| `--de-method` | `auto` | Method: auto, limrots, limma, deqms, proda, rots, ensemble |
| `--de-ensemble-methods` | `limrots,deqms,proda` | Comma-separated DE methods used when `--de-method=ensemble` |
| `--de-ensemble-min-k` | 2 | Minimum ensemble members that must agree on direction |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh or ihw |
| `--de-output` | none | DE results file; with multiple contrasts each is written as `<stem>_<A-B>.<ext>` |

Contrasts must be explicitly provided via `--de-contrasts` and/or `--de-contrasts-file`. Both can be combined.

`--de-method auto` selects `deqms` for `directlfq` quantification and `limrots` for other quantification methods. All methods run in the native Rust kernel — no R or rpy2 required. Deterministic methods (limma / deqms) are cell-exact against the Python reference; RNG/optimizer-driven methods (rots / limrots / proda) match log2 fold change cell-exactly and p-values at rank level.

`--de-method ensemble` runs each member method on the same contrast and combines the per-protein verdicts via top-k consensus: a protein is called UP/DOWN only when at least `--de-ensemble-min-k` members agree on direction and the Fisher-combined p-value passes the FDR threshold.

### Plots & Reports

`features2proteins` is pure compute and writes no figures. The kernel emits the protein-matrix CSV and (with `--de-output`) one DE result CSV per contrast; render plots and HTML reports from those CSVs with the wheel periphery:

```python
import mokume

# DE volcano / heatmap / PCA from the kernel CSVs (plotting extra):
mokume.de_plots(["--protein-matrix", "proteins.csv", "--plot-dir", "plots",
                 "--volcano", "--heatmap", "--pca",
                 "--contrast", "c1", "A", "B", "de.csv"])

# interactive HTML report (reports extra):
mokume.interactive_report(["--protein-matrix", "proteins.csv", "--report-output", "report.html"])
```

Pass `--help` to either (`python -m mokume.commands.de_plots --help`) for the full flag set.

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
| `--keep-shared-peptides` | off | Keep shared/non-unique peptides and skip the unique-peptide gate |
| `--remove_ids` | none | File with protein IDs to exclude |
| `--remove_decoy_contaminants` | off | Remove decoys and contaminants |
| `--remove_low_frequency_peptides` | off | Remove peptides in <20% of samples |

### Normalization

| Option | Default | Description |
|--------|---------|-------------|
| `--run-normalization` | `median` | Feature normalization: median, mean, max, global, max_min, iqr, none |
| `--sample-normalization` | `globalMedian` | Sample normalization: globalMedian, conditionMedian, hierarchical, none |
| `--skip_normalization` | off | Skip all normalization |
| `--log2` | off | Log2 transform output |
| `--save_parquet` | off | Save output as parquet |

!!! note "Only factor-based sample normalization changes the peptide output"
    In `features2peptides`, only the scalar-per-sample methods (`globalMedian` / `conditionMedian`) are applied; the cross-sample distribution-alignment methods (`quantile`, `rlr`, `loess`, `hierarchical`, `mediancenter`, `meancenter`) are accepted but are a deterministic no-op here (same as `none`), because the streaming peptide loop cannot see the full matrix and Python's per-sample loop also leaves them unchanged. All are implemented and oracle-verified in `features2proteins`, where the full matrix exists.

### TMT / ITRAQ

| Option | Default | Description |
|--------|---------|-------------|
| `--irs_channel` | none | Explicit pooled/reference channel label |
| `--irs_autodetect_regex` | none | Regex to detect pooled samples from SDRF |
| `--irs_stat` | `median` | IRS per-run statistic |
| `--irs_scope` | `global` | IRS scaling scope: global, by_mixture, or two_stage |
| `--aggregation_level` | `sample` | Aggregate at sample or run level |

!!! note "Channel-based IRS in `features2peptides`"
    The `--irs_channel` / `--irs_autodetect_regex` channel path scales on the TMT `mixture` / `channel` columns and is implemented for all three `--irs_scope` values (`global` / `by_mixture` / `two_stage`); the reference channel is taken from `--irs_channel` or auto-detected from the SDRF via `--irs_autodetect_regex`. SDRF-driven multi-plex IRS is also available in `features2proteins` (`--irs` with `--irs-reference-samples` / `--irs-sdrf-column`).

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

!!! note "Group-level filters: what runs and what is a no-op"
    The per-row filters (min-intensity floor, peptide length, charge states, excluded modifications, missed cleavages) and the per-`(protein, sample)` unique-peptide gate run natively and are oracle-locked vs Python. Among the **group-level** filters, CV threshold, quantile outlier removal, and the run-QC checks `--filter-min-features` / min-total-intensity / min-proteins are implemented via a pre-pass; replicate agreement reproduces Python's degenerate per-sample behaviour (a `>= 2` threshold empties the output); and `--filter-max-missing-rate`, sample-correlation, min-search-score, and min-coverage are no-ops that warn and pass rows through (matching Python's per-sample / column-absent skip). Only an unknown `razor-peptide-handling` value returns `NotImplemented`.

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
| `--families` | none | YAML file with explicit family overrides (iBAQ only; see [user guide](../user-guide/peptides2protein.md#family-discovery-tuning)) |
| `--min-shared` | 2 | Minimum shared peptides for auto-family discovery (iBAQ only) |
| `--min-anchors` | 1 | Minimum anchor peptides per family member (iBAQ only) |
| `--high-anchor-threshold` | 3 | Anchors required for an `EvidenceLevel` of `high` (iBAQ only) |
| `-o/--output` | required | Output file path |
| `--verbose` | off | Print a pointer to the wheel QC report command |
| `--qc_report` | QCprofile.pdf | QC report path echoed by `--verbose` |

Use `--method topn --topn_n 5` or `--method topn --topn_n 10` for Top5 or Top10-style quantification. `--output` is required for every method.

!!! note "iBAQ for unported enzymes lives in the wheel"
    The native iBAQ path digests proteins for the ported pyOpenMS enzymes (Trypsin[/P], Lys-C[/P], Arg-C[/P], Chymotrypsin[/P], Glu-C, Asp-N, Lys-N, PepsinA, ...), oracle-locked against pyOpenMS. For any other enzyme pyOpenMS knows (CNBr, V8-DE, unspecific cleavage, ...) the kernel has no cleavage rule and fails with an error pointing to the wheel: `mokume.peptides2protein_ibaq(peptides=..., fasta=..., enzyme="CNBr", output=...)` (the `ibaq` extra), which computes the whole iBAQ table in pure Python.

!!! note "`--verbose` no longer draws a QC PDF"
    QC plotting moved to the wheel. On the iBAQ path `--verbose` prints a one-line pointer to `mokume.peptides2protein_qc(protein_table=..., qc_report=...)` (the `plotting` extra), which draws the same density / box plots from the kernel's output table; it writes no PDF itself.

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
| `--sample_id_column` / `--sid` | `SampleID` | Sample ID column |
| `--protein_id_column` / `--pid` | `ProteinName` | Protein ID column |
| `--ibaq_raw_column` / `--ibaq` | `Ibaq` | Raw intensity column |
| `--ibaq_corrected_column` | `IbaqBec` | Corrected intensity column |
| `--comment` | `#` | Comment character |
| `--sep` | `\t` | Field separator |
| `--export_anndata` | off | Export to AnnData h5ad |

ComBat here is the native Rust implementation, oracle-verified against inmoose. The `.h5ad` written by `--export_anndata` is Rust-native and verified to round-trip through `anndata.read_h5ad`. This command does not expose batch-method or covariate options; for those, use `features2proteins --batch-correction`.

---

## Periphery (tissue maps, t-SNE)

Tissue-proteome maps and t-SNE visualization are **not** CLI subcommands — they ship in the `pip install mokume-rs` wheel and read the kernel's output tables:

```python
import mokume

mokume.tissuemap(scan_dir="./data", output_dir="./out")     # tissuemap extra
mokume.tsne_visualization(folder="./proteins", pattern="proteins.tsv")  # plotting extra
```

Each is also runnable as `python -m mokume.commands.tissuemap --help` / `python -m mokume.commands.visualize --help`. See the [Python API](python-api.md) for the full periphery surface and the matching extras.
