# CLI Reference

mokume provides four main commands. Use `mokume --help` or `mokume <command> --help` for details.

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
| `--quant-method` | `maxlfq` | Method: maxlfq, directlfq, ibaq, topn, sum, median, ratio |
| `--fasta` | none | FASTA file (required for iBAQ) |
| `--topn` | 3 | N for TopN quantification |
| `--ion-alignment` | none | Ion alignment: none or hierarchical |
| `--directlfq-cores` | auto | CPU cores for DirectLFQ |
| `--directlfq-min-nonan` | 1 | Min non-NaN values for DirectLFQ |

### Normalization

| Option | Default | Description |
|--------|---------|-------------|
| `--run-normalization` | `median` | Run-level: median, mean, iqr, max, max_min, none |
| `--sample-normalization` | `globalMedian` | Sample-level: globalMedian, conditionMedian, hierarchical, none |
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

### Differential Expression

| Option | Default | Description |
|--------|---------|-------------|
| `--de` | off | Enable differential expression analysis |
| `--de-contrasts` | all pairs | Comma-separated contrasts (e.g., "A-B,A-C") |
| `--de-method` | `auto` | Method: auto, limrots, deqms, or proda |
| `--de-log2fc` | 0.5 | Minimum absolute log2 fold change |
| `--de-fdr` | 0.05 | Maximum FDR threshold |
| `--de-fdr-method` | `bh` | FDR correction: bh or ihw |
| `--de-output` | auto | Output file for DE results |

### Plots & Reports

| Option | Default | Description |
|--------|---------|-------------|
| `--plot-dir` | none | Output directory for plots |
| `--plot-volcano` | off | Generate volcano plots |
| `--plot-heatmap` | off | Generate heatmaps |
| `--plot-pca` | off | Generate PCA plots |
| `--highlight-genes` | none | Comma-separated gene names to highlight |
| `--interactive-report` | off | Generate interactive HTML QC report |
| `--report-output` | auto | Output path for HTML report |

### Export

| Option | Default | Description |
|--------|---------|-------------|
| `--export-peptides` | none | Export normalized peptides to file |
| `--export-ions` | none | Export normalized ions (DirectLFQ only) |

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
| `--run-normalization` | `median` | Feature normalization: mean, median, iqr, none |
| `--sample-normalization` | `globalMedian` | Sample normalization: globalMedian, conditionMedian, hierarchical, none |
| `--skip_normalization` | off | Skip all normalization |
| `--log2` | off | Log2 transform output |
| `--save_parquet` | off | Save output as parquet |

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
| `--method` | `ibaq` | Method: ibaq, top3, top5, top10, topn, maxlfq, sum, directlfq |
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
| `-o/--output` | required | Output file path |
| `--verbose` | off | Print distribution info |
| `--qc_report` | QCprofile.pdf | QC report PDF path |

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
| `-ibaq/--ibaq_raw_column` | `IBAQ` | Raw intensity column |
| `--ibaq_corrected_column` | `IBAQ_BEC` | Corrected intensity column |
| `--comment` | `#` | Comment character |
| `--sep` | `\t` | Field separator |
| `--export_anndata` | off | Export to AnnData h5ad |

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
