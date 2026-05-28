# Computed Values

Reference for all output columns produced by mokume's quantification methods.

## Quantification Columns

| Column | Method | Description |
|--------|--------|-------------|
| `Ibaq` | iBAQ | Per-protein iBAQ (proportional branch) **or** family-level iBAQ (fallback). Interpretation depends on `EvidenceLevel`. |
| `IbaqNorm` | iBAQ | `ibaq / sum(ibaq)` per sample |
| `IbaqLog` | iBAQ | `10 + log10(IbaqNorm)` |
| `IbaqPpb` | iBAQ | `IbaqNorm * 100,000,000` (parts per billion) |
| `IbaqBec` | iBAQ + ComBat | Batch effect corrected iBAQ |
| `FamilyId` | iBAQ | Canonical accession identifying the protein family used by the piBAQ algorithm |
| `FamilySize` | iBAQ | Number of canonical members in the family (1 = singleton, isolated protein) |
| `EvidenceLevel` | iBAQ | `high` (≥3 unique anchors), `medium` (1-2 anchors), or `family_only` (zero anchors -> aggregated) |
| `TopNIntensity` | TopN | Average of top N peptides (e.g., Top3Intensity, Top5Intensity) |
| `MaxLFQIntensity` | MaxLFQ | MaxLFQ algorithm result |
| `DirectLFQIntensity` | DirectLFQ | DirectLFQ intensity traces |
| `SumIntensity` | Sum | Sum of all peptide intensities |
| `Intensity` | Unified pipeline | Standard output column (all methods) |

## Derived Values (iBAQ)

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

| Column | Description |
|--------|-------------|
| `log2FoldChange` | Log2 fold change between conditions |
| `pValue` | Raw p-value |
| `adjustedPValue` | FDR-adjusted p-value |
| `significant` | Whether the protein passes thresholds |

## Pipeline Output Format

All quantification methods in the unified pipeline (`features2proteins`) produce a **wide-format** output: rows are proteins, columns are samples.

```
ProteinName,sample1,sample2,sample3,...
P02452,1234.5,5678.9,2345.6,...
P12345,9876.5,4321.0,8765.4,...
```

For ratio quantification, values are in **log2 space** (log2 sample/reference).
