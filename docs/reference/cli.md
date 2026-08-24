# CLI Reference

The `mokume` console command is installed by the `mokume` wheel. It exposes exactly four compute subcommands — `features2proteins`, `features2peptides`, `peptides2protein`, and `correct-batches` — and runs the compiled Rust kernel in-process. Most users start with `features2proteins` for quantification workflows.
Use `mokume --help` or `mokume <command> --help` for details.

The flag surface below is single-sourced in Rust and is shared with the wheel's thin Python API.

!!! note "Plotting, tissue maps, and reports live in the Python wheel"
    The visualization periphery — t-SNE, tissue-proteome maps, DE plots,
    interactive HTML reports, piBAQ QC — is **not** part of this CLI. It ships in
    the `pip install mokume` wheel as `mokume.tsne_visualization(...)`,
    `mokume.tissuemap(...)`, `mokume.de_plots([...])`,
    `mokume.interactive_report([...])`, and
    `mokume.peptides2protein_qc(...)`. Plotting and reporting consume kernel
    tables; TissueMap performs its documented downstream analysis from QPX data.
    See the [Python API](python-api.md).

## features2proteins

The unified pipeline: features to protein quantification in one step.

```bash
mokume features2proteins [OPTIONS]
```

### Input & Output

| Option | Default | Description |
|--------|---------|-------------|
| `-p/--parquet` | none | Input quantms.io/QPX feature parquet; mutually exclusive with `--msstats` |
| `--msstats` | none | Input MSstats CSV; mutually exclusive with `--parquet` and requires `--sdrf` |
| `-o/--output` | required | Output protein intensities CSV |
| `--output-format` | `python-compatible` | Protein identifier header: `ProteinName` for `python-compatible`, `protein` for `rust-native` |

Provide exactly one of `--parquet` or `--msstats`. MSstats input requires
`ProteinName`, `PeptideSequence`, `Intensity`, `Charge` or `PrecursorCharge`,
and `Run` or `Reference`; isobaric data also requires `Channel`. Ratio
quantification requires PSM-level QPX evidence and does not accept MSstats
feature tables.

### Metadata & Filtering

| Option | Default | Description |
|--------|---------|-------------|
| `-s/--sdrf` | none | SDRF file for sample metadata |
| `--min-aa` | 7 | Minimum amino acid length |
| `--min-unique` | 2 | Minimum unique peptides per protein |
| `--keep-contaminants` | off | Keep contaminants and decoys instead of removing them |

### Quantification

| Option | Default | Description |
|--------|---------|-------------|
| `--quant-method` | `maxlfq` | Method: maxlfq, directlfq, pibaq, `top<N>` (top3, top5, top10, ...), sum, median, ratio, abd (TMT abundance), intensity (TMT reporter), spectral_count |
| `--fasta` | none | FASTA file (required for piBAQ) |
| `--directlfq-cores` | auto | Fallback Rayon thread count for DirectLFQ when `--threads` is omitted |
| `--directlfq-min-nonan` | 1 | Min non-NaN values for DirectLFQ |
| `--directlfq-num-samples-quadratic` | 50 | Maximum samples in DirectLFQ's quadratic global-alignment subset |
| `--pibaq-enzyme` | `Trypsin` | Protease name from the installed pyOpenMS catalog |
| `--pibaq-max-aa` | 50 | Maximum theoretical peptide length |
| `--pibaq-min-shared` | 2 | Minimum shared peptides for automatic family discovery |
| `--pibaq-families` | none | YAML file with explicit family overrides |
| `--pibaq-min-anchors` | 1 | Minimum unique-peptide anchors required before proportional family allocation |
| `--pibaq-high-anchor-threshold` | 3 | Anchor threshold used to classify high-evidence family members |

!!! note "Write the N in the method name: `--quant-method top<N>`"
    TopN quantification takes its N from the method name — `--quant-method top3`,
    `top5`, `top10`, and so on for any N ≥ 1.

    The `--topn` option has been removed: replace `--quant-method topn --topn 5`
    with `--quant-method top5`. Bare `--quant-method topn` still works and still
    means N = 3, the same as `--quant-method top3`. A `top` name with no arabic
    numeral (`topa`) is rejected rather than silently treated as Top3.

### Normalization

| Option | Default | Description |
|--------|---------|-------------|
| `--run-normalization` | method-dependent | Run-level: median, mean, max, global, max_min, iqr, none |
| `--sample-normalization` | method-dependent | Sample-level: globalMedian, conditionMedian, hierarchical, quantile, mediancenter, meancenter, rlr, loess, tmm, none |
| `--normalization-proteins` | none | File with protein IDs for normalization |

DirectLFQ and Ratio manage normalization internally and therefore default both
normalization layers to `none`; passing a non-`none` value is rejected. Other
methods default to `median` / `globalMedian`.

### IRS (Multi-Plex TMT)

| Option | Default | Description |
|--------|---------|-------------|
| `--irs` | off | Enable IRS normalization |
| `--irs-reference-samples` | auto | Comma-separated reference sample names |
| `--irs-reference-sample` | none | Repeatable single reference sample name; conflicts with `--irs-reference-samples` |
| `--irs-sdrf-column` | auto | SDRF column for reference detection |
| `--irs-sdrf-values` | auto | Values indicating reference samples |
| `--irs-reference-regex` | `pool\|powder\|ref\|reference\|bridge` | Regex for reference auto-detection |
| `--irs-stat` | `median` | Plex reference statistic: median or mean |
| `--irs-remove-reference` | off | Remove reference samples from output |

IRS options require `--irs` and an SDRF. Reference detection must find usable
reference samples and plex assignments; otherwise the command fails.

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
| `--batch-method` | `sample_prefix` | Detection method: sample_prefix or column |
| `--batch-column` | none | SDRF column used when `--batch-method=column` |
| `--batch-covariates` | none | Comma-separated SDRF columns to preserve |
| `--batch-nonparametric` | off | Use non-parametric ComBat instead of the parametric default |
| `--batch-mean-only` | off | Only adjust batch means |
| `--batch-ref` | none | Reference batch ID |

ComBat is a native Rust implementation, oracle-verified against inmoose
(parametric ~1e-6, covariate / non-parametric / `ref_batch` / `mean_only`
paths included). It runs on proteins with no missing cells; the rest are kept
uncorrected. Invalid batch layouts and an empty complete-protein subset are
errors rather than successful uncorrected outputs.

### Imputation

| Option | Default | Description |
|--------|---------|-------------|
| `--impute` | off | Enable missing-value imputation on the protein matrix |
| `--impute-method` | required with `--impute` | Method: mean, median, constant, zero, most_frequent, knn, minprob, mindet, qrilc, seqknn, impseq, gms, bpca, impseqrob |
| `--impute-quantile` | 0.01 | Quantile for MinProb/MinDet |
| `--impute-shift` | 1.6 | MinProb shift in standard deviations |
| `--impute-scale` | 0.3 | MinProb scale factor for the imputation distribution sigma |
| `--impute-n-neighbors` | 5 | Number of neighbours for KNN/SeqKNN |

Imputation runs in log2 space (the matrix is transformed before imputation and back to linear afterwards) so that censored-aware methods like MinProb/MinDet/QRILC behave correctly. The imputation step is applied after coverage filtering and before batch correction in the pipeline.

`missforest` is intentionally absent from this Rust compute command because it
wraps scikit-learn. Install `mokume[analysis]` and call the wheel's Python
`mokume.impute` API when that method is required.

### Differential Expression

| Option | Default | Description |
|--------|---------|-------------|
| `--de` | off | Enable differential expression analysis |
| `--de-contrasts` | — | Comma-separated contrasts (e.g., `"A vs B,A vs C"`) |
| `--de-contrasts-file` | — | TSV file with columns `group1`, `group2` |
| `--de-method` | `auto` | Method: auto, limrots, limma, deqms, proda, rots, ensemble |
| `--de-ensemble-methods` | `limrots,deqms,proda` | Comma-separated DE methods used when `--de-method=ensemble` |
| `--de-ensemble-min-k` | 2 | Minimum ensemble members that must agree on direction |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change, or `auto` for the data-driven mixture gate |
| `--de-effect-size-gate` | none | Explicit data-driven gate: mixture or null_quantile; a numeric log2FC value is its fallback |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh, ihw, bky, or storey |
| `--de-output` | required with `--de` | DE results file; with multiple contrasts each is written as `<stem>_<A-B>.<ext>` |

Contrasts and an output path must be explicitly provided when `--de` is set;
this prevents a completed DE calculation from being discarded. Inline and file
contrasts can be combined.

`--de-method auto` selects `deqms` for `directlfq` quantification and `limrots`
for other quantification methods. All methods run in the native Rust kernel —
no R or rpy2 required. ROTS and LimROTS retain permutation FDR and therefore
reject alternative `--de-fdr-method` values. Deterministic methods (limma /
deqms) are cell-exact against frozen Python-generated compatibility output;
RNG/optimizer-driven methods (rots / limrots / proda) match log2 fold change
cell-exactly and p-values at rank level.

`--de-method ensemble` runs each member method on the same contrast and combines the per-protein verdicts via top-k consensus: a protein is called UP/DOWN only when at least `--de-ensemble-min-k` members agree on direction and the Fisher-combined p-value passes the FDR threshold. Eligible non-ROTS members use the requested correction; ROTS and LimROTS retain their native permutation FDR. The Fisher-combined p-values use BH by default, with BKY or Storey applied when requested and reliable; IHW remains a member-level correction because the combined rows have no IHW covariate.

### Plots & Reports

`features2proteins` is pure compute and writes no figures. The kernel emits the protein-matrix CSV and (with `--de-output`) one DE result CSV per contrast; render plots and HTML reports from those CSVs with the wheel periphery:

```python
import mokume

# DE volcano / heatmap / PCA from the kernel CSVs (plotting extra):
mokume.de_plots(["--protein-matrix", "proteins.csv", "--plot-dir", "plots",
                 "--sdrf", "experiment.sdrf.tsv",
                 "--volcano", "--heatmap", "--pca",
                 "--contrast", "c1", "A", "B", "de.csv"])

# interactive HTML report (reports extra):
mokume.interactive_report([
    "--protein-matrix", "proteins.csv",
    "--sdrf", "experiment.sdrf.tsv",
    "--report-output", "report.html",
    "--contrast", "c1", "A", "B", "de.csv",
])
```

Pass `--help` to either (`python -m mokume.commands.de_plots --help`) for the full flag set.

### Export

| Option | Default | Description |
|--------|---------|-------------|
| `--export-peptides` | none | Export normalized peptides (not supported with DirectLFQ) |
| `--export-ions` | none | Export normalized ions (DirectLFQ only) |

`--export-peptides` also rejects dataset-level sample normalization for
non-cell-based aggregation methods, because those peptide values would not
represent the normalized protein matrix.

### Runtime Resource Controls

| Option | Default | Description |
|--------|---------|-------------|
| `--threads` | Rayon default | Size the Rayon thread pool used by parallel Rust sections; alias: `--duckdb-threads` |

The Rust stream does not expose a memory-limit option because it cannot enforce
one. Use systemd/cgroup, SLURM, or container resource limits for a hard ceiling.

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
| `--sample-normalization` | `globalMedian` | Sample normalization: globalMedian, conditionMedian, none |
| `--skip_normalization` | off | Skip all normalization |
| `--log2` | off | Log2 transform output |
| `--save_parquet` | off | Save output as parquet |

Dataset-level sample normalizers are not accepted by `features2peptides`; run
them through `features2proteins`, where the full matrix exists.

### TMT / ITRAQ

| Option | Default | Description |
|--------|---------|-------------|
| `--irs_channel` | none | Explicit pooled/reference channel label |
| `--irs_autodetect_regex` | none | Regex to detect pooled samples from SDRF |
| `--irs_stat` | `median` | IRS per-run statistic |
| `--irs_scope` | `global` | IRS scaling scope: global, by_mixture, or two_stage |
| `--aggregation_level` | `sample` | Aggregate at sample or run level |

!!! note "Channel-based IRS in `features2peptides`"
    Choose exactly one of `--irs_channel` and `--irs_autodetect_regex`.
    Autodetection requires an SDRF and must match a reference channel; a
    requested channel must produce scaling factors. `--skip_normalization`
    conflicts with IRS. SDRF-driven multi-plex IRS is also available in
    `features2proteins` (`--irs` with `--irs-reference-samples` /
    `--irs-sdrf-column`).

### Filter Configuration

| Option | Default | Description |
|--------|---------|-------------|
| `--filter-config` | none | YAML/JSON filter configuration file |
| `--generate-filter-config` | none | Generate example config and exit |
| `--filter-min-intensity` | none | Min intensity threshold (override) |
| `--filter-cv-threshold` | none | Max CV across replicates (override) |
| `--filter-charge-states` | none | Comma-separated charge states (override) |
| `--filter-max-missed-cleavages` | none | Max missed cleavages (override) |
| `--filter-peptide-fdr` | none | Max QPX peptide q-value (override) |
| `--filter-exclude-modifications` | none | Comma-separated modifications (override) |
| `--filter-min-unique-peptides` | none | Min unique peptides (override) |
| `--filter-protein-fdr` | none | Max QPX protein-group q-value (override) |
| `--filter-min-features` | none | Min features per run (override) |
| `--filter-max-missing-rate` | none | Max missing feature fraction per technical run (override) |

Only implemented filters appear in the generated example. FDR filtering is
disabled by default and requires populated `peptide_qvalue` or
`pg_global_qvalue` input. Missing rate uses the complete distinct
`(protein, peptide)` universe among the surviving technical runs in each sample.
Active unsupported filter-config values are rejected instead of being logged and
ignored.

---

## peptides2protein

Protein quantification from normalized peptide data.

```bash
mokume peptides2protein [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-p/--peptides` | required | Input peptide intensity file |
| `-f/--fasta` | none | FASTA file (required for piBAQ) |
| `--method` | `pibaq` | Method: pibaq, `top<N>` (top3, top5, top10, ...), maxlfq, sum, directlfq |
| `-e/--enzyme` | `Trypsin` | Enzyme for in-silico digestion |
| `-n/--normalize` | off | Normalize quantification values |
| `--min_aa` | 7 | Min amino acid length |
| `--max_aa` | 30 | Max amino acid length |
| `-t/--tpa` | off | Calculate TPA (piBAQ only) |
| `-r/--ruler` | off | ProteomicRuler (piBAQ only) |
| `-i/--ploidy` | 2 with `--ruler` | Positive ploidy number (ruler only) |
| `-m/--organism` | `human` with `--ruler` | Organism for histone data (ruler only) |
| `-c/--cpc` | 200 with `--ruler` | Positive cellular protein concentration in g/L (ruler only) |
| `--threads` | -1 | Rayon threads for MaxLFQ and DirectLFQ; positive values set the pool, `-1` uses all available CPUs, and `-2` leaves one free; zero is rejected |
| `--min_nonan` | 1 | Min non-NaN for DirectLFQ |
| `--families` | none | YAML file with explicit family overrides (piBAQ only; see [user guide](../user-guide/peptides2protein.md#family-discovery-tuning)) |
| `--min-shared` | 2 | Minimum shared peptides for auto-family discovery (piBAQ only) |
| `--min-anchors` | 1 | Anchor threshold; if no member reaches it, shared signal is split equally and evidence is `family_only` (piBAQ only) |
| `--high-anchor-threshold` | 3 | Anchors every member must reach for `EvidenceLevel=high` (piBAQ only) |
| `-o/--output` | required | Output file path |
| `--verbose` | off | Generate the piBAQ QC PDF after native quantification |
| `--qc_report` | QCprofile.pdf | QC report path used by `--verbose` |

Use `--method top5` or `--method top10` for Top5 or Top10-style quantification; `top3` is the named method from [Silva et al. 2006](https://doi.org/10.1074/mcp.M500230-MCP200). `--output` is required for every method.

!!! warning "`--topn_n` has been removed"
    N now comes from the method name only. Replace `--method topn --topn_n 5`
    with `--method top5`. Bare `--method topn` still works and still means N = 3.
    `--topn_n` existed in mokume 0.1.0, so scripts written against that release
    and passing it need updating.

!!! note "piBAQ uses the installed pyOpenMS catalog"
    Both piBAQ commands query the installed pyOpenMS `ProteaseDB` at runtime and support its complete catalog. Python digests the FASTA and passes the full protein-to-theoretical-peptide map into Rust; there is no separate unported-enzyme branch or `pibaq` extra. At `debug` or `info` log level, the run log records the pyOpenMS version, canonical enzyme, catalog SHA-256, peptide-length bounds, and missed-cleavage count.

!!! note "`--verbose` plots the native result table"
    Rust writes the piBAQ table first, then the wheel renders the density and box
    plots from those exact values. Install `mokume[plotting]` to use this option.

---

## correct-batches

Standalone batch correction for pre-quantified data.

```bash
mokume correct-batches [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--folder` | required | Folder with TSV files |
| `-p/--pattern` | `*pibaq.tsv` | File matching pattern |
| `-o/--output` | required | Output file path |
| `--sample_id_column` / `--sid` | `SampleID` | Sample ID column |
| `--protein_id_column` / `--pid` | `ProteinName` | Protein ID column |
| `--pibaq_raw_column` / `--pibaq` | `PiBAQ` | Raw intensity column |
| `--pibaq_corrected_column` | `PiBAQBec` | Corrected intensity column |
| `--comment` | `#` | Comment character |
| `--sep` | `\t` | Field separator |
| `--export_anndata` | off | Export to AnnData h5ad |

ComBat here is the native Rust implementation, oracle-verified against inmoose. The `.h5ad` written by `--export_anndata` is Rust-native and verified to round-trip through `anndata.read_h5ad`. This command does not expose batch-method or covariate options; for those, use `features2proteins --batch-correction`.

---

## Periphery (tissue maps, t-SNE)

Tissue-proteome maps and t-SNE visualization are **not** CLI subcommands; they
ship in the `pip install mokume` wheel. t-SNE visualization reads a protein
matrix, while TissueMap derives a downstream atlas from QPX data:

```python
import mokume

mokume.tissuemap(scan_dir="./data", output_dir="./out")     # tissuemap extra
mokume.tsne_visualization(folder="./proteins", pattern="proteins.tsv")  # plotting extra
```

Each is also runnable as `python -m mokume.commands.tissuemap --help` / `python -m mokume.commands.visualize --help`. See the [Python API](python-api.md) for the full periphery surface and the matching extras.
