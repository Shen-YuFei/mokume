# Computed Values

Reference for all output columns produced by mokume's quantification methods.

## Quantification Columns

| Column | Method | Description |
|--------|--------|-------------|
| `PiBAQ` | piBAQ | Per-protein piBAQ after shared-peptide signal is allocated proportionally or, without anchor signal, equally. |
| `PiBAQNorm` | piBAQ | `PiBAQ / sum(PiBAQ)` per sample |
| `PiBAQLog` | piBAQ | `10 + log10(PiBAQNorm)` |
| `PiBAQPpb` | piBAQ | `PiBAQNorm * 100,000,000` (parts per billion) |
| `PiBAQBec` | piBAQ + ComBat | Batch effect corrected piBAQ |
| `FamilyId` | piBAQ | Canonical accession identifying the protein family used by the piBAQ algorithm |
| `FamilySize` | piBAQ | Number of canonical members in the family (1 = singleton, isolated protein) |
| `EvidenceLevel` | piBAQ | `high` (all members meet the high-anchor threshold), `medium` (at least one member meets the minimum), or `family_only` (none does) |
| `TopNIntensity` | TopN | Average of top N peptides (e.g., Top3Intensity, Top5Intensity) |
| `MaxLFQIntensity` | MaxLFQ | MaxLFQ algorithm result |
| `DirectLFQIntensity` | DirectLFQ | DirectLFQ intensity traces |
| `SumIntensity` | Sum | Sum of all peptide intensities |
| `Intensity` | Unified pipeline | Standard output column (all methods) |

## Derived Values (piBAQ)

| Column | Formula | Description |
|--------|---------|-------------|
| `TPA` | `NormIntensity / MolecularWeight` | Total Protein Approach |
| `CopyNumber` | ProteomicRuler calculation | Protein copies per cell |
| `Concentration[nM]` | ProteomicRuler calculation | Protein concentration |

## Metadata Columns

| Column | Description |
|--------|-------------|
| `ProteinName` | UniProt accession (e.g., P02452) |
| `SampleID` | Sample identifier |
| `Condition` | Experimental condition |
| `BioReplicate` | Biological replicate |
| `PeptideSequence` | Amino acid sequence |
| `NormIntensity` | Normalized peptide intensity |

## Differential Expression Columns

The kernel writes one DE result CSV per contrast. The leading columns are shared across methods; later columns are method-specific (e.g. `t_stat`/`AveExpr`/`B` for limma, `sca_t`/`peptide_count` for deqms, `d_stat` for rots).

| Column | Description |
|--------|-------------|
| `ProteinName` | UniProt accession (`Protein` in the ensemble output) |
| `log2FC` | Log2 fold change between conditions |
| `pvalue` | Raw p-value |
| `adj_pvalue` | FDR-adjusted p-value |
| `log_pvalue` | Natural logarithm of the raw p-value (DEqMS only; remains finite when `pvalue` underflows to zero) |
| `significance` | Whether the protein passes the `--de-log2fc` / `--de-fdr` thresholds |

## Pipeline Output Format

All quantification methods in the unified pipeline (`features2proteins`) produce a **wide-format** output: rows are proteins, columns are samples.

```text
ProteinName,sample1,sample2,sample3,...
P02452,1234.5,5678.9,2345.6,...
P12345,9876.5,4321.0,8765.4,...
```

For ratio quantification, values are in **log2 space** (log2 sample/reference).
