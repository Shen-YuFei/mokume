# Absolute Expression

Absolute quantification estimates how much of each protein is present, not just
relative fold changes. Mokume implements piBAQ, a paralog-aware extension of
the original iBAQ scaling, plus the Total Protein Approach (TPA) and the
ProteomicRuler for copy numbers and concentrations.

Two entry points:

- **`features2proteins --quant-method pibaq --fasta ...`** — piBAQ inside the
  unified pipeline (loads features, filters, normalizes, then computes piBAQ).
- **`peptides2protein --method pibaq ...`** — the standalone step that takes an
  already-normalized peptide table and adds all the absolute columns (TPA,
  ProteomicRuler, ...).

Both need a **FASTA** file. piBAQ retains iBAQ's theoretical-peptide scaling
while allocating shared-peptide intensity explicitly and using an owned-peptide
denominator symmetric with that allocation. Both Rust-backed entry points digest
through every protease registered in the installed pyOpenMS runtime; Python passes
the complete theoretical-peptide map to the Rust aggregation kernel.

## piBAQ inside the pipeline

=== "CLI"

    ```bash
    mokume features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_pibaq.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method pibaq \
        --fasta python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.features2proteins(
        parquet="python/tests/example/feature_wide.parquet",
        output="proteins_pibaq.csv",
        sdrf="python/tests/example/PXD020192.sdrf.tsv",
        quant_method="pibaq",
        fasta="python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta",
    )
    ```

=== "Python (package)"

    ```python
    from mokume.pipeline.features_to_proteins import QuantificationPipeline
    from mokume.pipeline.config import (
        PipelineConfig,
        InputConfig,
        QuantificationConfig,
    )

    config = PipelineConfig(
        input=InputConfig(
            parquet="python/tests/example/feature_wide.parquet",
            sdrf="python/tests/example/PXD020192.sdrf.tsv",
            fasta_file="python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta",
        ),
        quantification=QuantificationConfig(method="pibaq"),
    )
    proteins = QuantificationPipeline(config).run()
    ```

!!! note "Empty output on the tiny fixture is expected"

    `feature_wide.parquet` is a 500-feature slice, and piBAQ requires enough unique
    anchor peptides per protein. On this slice the anchor filter removes every
    protein, so the run completes cleanly but the matrix has zero rows. Use the
    richer `peptides2protein` example below to see fully populated piBAQ output.

## Standalone piBAQ + TPA + ProteomicRuler

`peptides2protein` takes a normalized peptide table
(`python/tests/example/PXD017834-peptides.csv`) and the FASTA, and produces a
long-format table with every absolute-expression column. `--tpa` adds the Total
Protein Approach; `--ruler` runs the ProteomicRuler; `--organism human` selects
the histone reference used by the ruler.

=== "CLI"

    ```bash
    mokume peptides2protein \
        --method pibaq \
        --tpa --ruler --organism human \
        -f python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta \
        -p python/tests/example/PXD017834-peptides.csv \
        -o proteins_pibaq_absolute.tsv
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.peptides2protein(
        method="pibaq",
        tpa=True,
        ruler=True,
        organism="human",
        fasta="python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta",
        peptides="python/tests/example/PXD017834-peptides.csv",
        output="proteins_pibaq_absolute.tsv",
    )
    ```

This writes a 14-column table (1670 rows on the fixture). The first row looks
like:

```text
ProteinName  SampleID            Condition     NormIntensity  PiBAQ    FamilyId  EvidenceLevel  FamilySize  MolecularWeight  TPA       CopyNumber   Moles[nmol]   Weight[ng]  Concentration[nM]
A4D1B5       PXD017834-Sample-1  Blood Plasma  50.233559      1.357664 A4D1B5    high           1           97738.851166     0.000514  7.049979e+06 1.170677e-08  0.001144    0.266295
```

## Normalizing piBAQ in Python

The package exposes the piBAQ normalization step directly. `normalize_pibaq` takes a
DataFrame with `ProteinName`, `SampleID`, `Condition`, and `PiBAQ` columns and adds
the PRIDE/ProteomicsDB-normalized columns.

=== "Python (package)"

    ```python
    import pandas as pd
    from mokume.quantification import normalize_pibaq

    df = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3"],
            "SampleID": ["s1", "s1", "s1"],
            "Condition": ["ctrl", "ctrl", "ctrl"],
            "PiBAQ": [10.0, 30.0, 60.0],
        }
    )
    result = normalize_pibaq(df)
    print(result[["ProteinName", "PiBAQ", "PiBAQNorm", "PiBAQLog", "PiBAQPpb"]])
    ```

    Output — `PiBAQNorm` is the per-sample fraction, `PiBAQLog` is
    `10 + log10(PiBAQNorm)`, and `PiBAQPpb` is `PiBAQNorm * 1e8`:

    ```text
      ProteinName  PiBAQ  PiBAQNorm  PiBAQLog     PiBAQPpb
    0          P1  10.0       0.1   9.0000  10000000.0
    1          P2  30.0       0.3   9.4771  30000000.0
    2          P3  60.0       0.6   9.7782  60000000.0
    ```

    The full end-to-end function `mokume.quantification.peptides_to_protein`
    (FASTA + peptide table -> absolute table with an optional QC-report PDF) is
    the machinery behind the `peptides2protein` command; drive it through the CLI
    or wheel wrapper shown above rather than calling it directly.

## Computed columns

`peptides2protein --method pibaq --tpa --ruler` produces these columns. The
`EvidenceLevel` records the strength of member-resolving anchor evidence.
`family_only` means that no family member reaches the minimum anchor threshold;
shared signal is split equally rather than duplicated across the family.

| Column | Meaning |
|--------|---------|
| `PiBAQ` | Per-protein piBAQ after exact shared-peptide allocation |
| `PiBAQNorm` | `PiBAQ / sum(PiBAQ)` within each sample |
| `PiBAQLog` | `10 + log10(PiBAQNorm)` (ProteomicsDB convention) |
| `PiBAQPpb` | `PiBAQNorm * 1e8` (parts per billion, PRIDE convention) |
| `FamilyId` | Canonical accession of the piBAQ protein family |
| `FamilySize` | Number of members in the family (1 = singleton) |
| `EvidenceLevel` | `high` (all members meet the high threshold), `medium` (some member meets the minimum), or `family_only` (none does) |
| `MolecularWeight` | Protein MW used by TPA and the ruler |
| `TPA` | `NormIntensity / MolecularWeight` (Total Protein Approach) |
| `CopyNumber` | ProteomicRuler protein copies per cell |
| `Moles[nmol]` / `Weight[ng]` / `Concentration[nM]` | ProteomicRuler amounts |

The full reference, including the piBAQ-plus-ComBat `PiBAQBec` column, lives in
[Computed Values](../reference/computed-values.md).
