# Absolute Expression

Absolute quantification estimates how much of each protein is present, not just
relative fold changes. mokume implements iBAQ (Intensity-Based Absolute
Quantification) with a paralog-aware family model (piBAQ), plus the Total Protein
Approach (TPA) and the ProteomicRuler for copy numbers and concentrations.

Two entry points:

- **`features2proteins --quant-method ibaq --fasta ...`** — iBAQ inside the
  unified pipeline (loads features, filters, normalizes, then computes iBAQ).
- **`peptides2protein --method ibaq ...`** — the standalone step that takes an
  already-normalized peptide table and adds all the absolute columns (TPA,
  ProteomicRuler, ...).

Both need a **FASTA** file — iBAQ divides observed intensity by the number of
theoretically observable tryptic peptides, which comes from digesting the FASTA.

## iBAQ inside the pipeline

=== "CLI"

    ```bash
    mokume features2proteins \
        -p python/tests/example/feature_wide.parquet \
        -o proteins_ibaq.csv \
        -s python/tests/example/PXD020192.sdrf.tsv \
        --quant-method ibaq \
        --fasta python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.features2proteins(
        parquet="python/tests/example/feature_wide.parquet",
        output="proteins_ibaq.csv",
        sdrf="python/tests/example/PXD020192.sdrf.tsv",
        quant_method="ibaq",
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
        quantification=QuantificationConfig(method="ibaq"),
    )
    proteins = QuantificationPipeline(config).run()
    ```

!!! note "Empty output on the tiny fixture is expected"

    `feature_wide.parquet` is a 500-feature slice, and iBAQ requires enough unique
    anchor peptides per protein. On this slice the anchor filter removes every
    protein, so the run completes cleanly but the matrix has zero rows. Use the
    richer `peptides2protein` example below to see fully populated iBAQ output.

## Standalone iBAQ + TPA + ProteomicRuler

`peptides2protein` takes a normalized peptide table
(`python/tests/example/PXD017834-peptides.csv`) and the FASTA, and produces a
long-format table with every absolute-expression column. `--tpa` adds the Total
Protein Approach; `--ruler` runs the ProteomicRuler; `--organism human` selects
the histone reference used by the ruler.

=== "CLI"

    ```bash
    mokume peptides2protein \
        --method ibaq \
        --tpa --ruler --organism human \
        -f python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta \
        -p python/tests/example/PXD017834-peptides.csv \
        -o proteins_ibaq_absolute.tsv
    ```

=== "Python (wheel)"

    ```python
    import mokume

    mokume.peptides2protein(
        method="ibaq",
        tpa=True,
        ruler=True,
        organism="human",
        fasta="python/tests/example/Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta",
        peptides="python/tests/example/PXD017834-peptides.csv",
        output="proteins_ibaq_absolute.tsv",
    )
    ```

This writes a 14-column table (1688 rows on the fixture). The first row looks
like:

```text
ProteinName  SampleID            Condition     NormIntensity  Ibaq      FamilyId    EvidenceLevel  FamilySize  MolecularWeight  TPA       CopyNumber    Moles[nmol]  Weight[ng]  Concentration[nM]
A0A075B6I0   PXD017834-Sample-1  Blood Plasma  34.218126      8.554531  A0A075B6I0  medium         1           12806.123706     0.002672  1.063371e+10  0.000018     0.226127    116730.137994
```

## Normalizing iBAQ in Python

The package exposes the iBAQ normalization step directly. `normalize_ibaq` takes a
DataFrame with `ProteinName`, `SampleID`, `Condition`, and `Ibaq` columns and adds
the PRIDE/ProteomicsDB-normalized columns.

=== "Python (package)"

    ```python
    import pandas as pd
    from mokume.quantification import normalize_ibaq

    df = pd.DataFrame(
        {
            "ProteinName": ["P1", "P2", "P3"],
            "SampleID": ["s1", "s1", "s1"],
            "Condition": ["ctrl", "ctrl", "ctrl"],
            "Ibaq": [10.0, 30.0, 60.0],
        }
    )
    result = normalize_ibaq(df)
    print(result[["ProteinName", "Ibaq", "IbaqNorm", "IbaqLog", "IbaqPpb"]])
    ```

    Output — `IbaqNorm` is the per-sample fraction, `IbaqLog` is
    `10 + log10(IbaqNorm)`, and `IbaqPpb` is `IbaqNorm * 1e8`:

    ```text
      ProteinName  Ibaq  IbaqNorm  IbaqLog     IbaqPpb
    0          P1  10.0       0.1   9.0000  10000000.0
    1          P2  30.0       0.3   9.4771  30000000.0
    2          P3  60.0       0.6   9.7782  60000000.0
    ```

    The full end-to-end function `mokume.quantification.peptides_to_protein`
    (FASTA + peptide table -> absolute table with an optional QC-report PDF) is
    the machinery behind the `peptides2protein` command; drive it through the CLI
    or wheel wrapper shown above rather than calling it directly.

## Computed columns

`peptides2protein --method ibaq --tpa --ruler` produces these columns. The
`EvidenceLevel` label tells you how to read `Ibaq`: `high`/`medium` are
per-protein proportional iBAQ, while `family_only` means the value was rolled up
across a paralog family.

| Column | Meaning |
|--------|---------|
| `Ibaq` | Per-protein iBAQ (or family-level, per `EvidenceLevel`) |
| `IbaqNorm` | `Ibaq / sum(Ibaq)` within each sample |
| `IbaqLog` | `10 + log10(IbaqNorm)` (ProteomicsDB convention) |
| `IbaqPpb` | `IbaqNorm * 1e8` (parts per billion, PRIDE convention) |
| `FamilyId` | Canonical accession of the piBAQ protein family |
| `FamilySize` | Number of members in the family (1 = singleton) |
| `EvidenceLevel` | `high` (>=3 anchors), `medium` (1-2), or `family_only` (0) |
| `MolecularWeight` | Protein MW used by TPA and the ruler |
| `TPA` | `NormIntensity / MolecularWeight` (Total Protein Approach) |
| `CopyNumber` | ProteomicRuler protein copies per cell |
| `Moles[nmol]` / `Weight[ng]` / `Concentration[nM]` | ProteomicRuler amounts |

The full reference, including the iBAQ-plus-ComBat `IbaqBec` column, lives in
[Computed Values](../reference/computed-values.md).
