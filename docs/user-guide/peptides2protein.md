# peptides2protein: Protein Quantification

The `peptides2protein` command quantifies proteins from normalized peptide data. It supports multiple quantification methods and is the second step of the two-step pipeline.

## Basic Usage

=== "CLI"

    ```bash
    # iBAQ (default, requires FASTA)
    mokume peptides2protein --method ibaq \
        -f proteome.fasta \
        -p peptides.csv \
        -o proteins-ibaq.tsv

    # TopN (no FASTA needed)
    mokume peptides2protein --method top3 \
        -p peptides.csv \
        -o proteins-top3.tsv

    # MaxLFQ with parallelization
    mokume peptides2protein --method maxlfq \
        --threads 4 \
        -p peptides.csv \
        -o proteins-maxlfq.tsv
    ```

=== "Python"

    ```python
    from mokume.quantification import (
        TopNQuantification,
        MaxLFQQuantification,
        AllPeptidesQuantification,
        peptides_to_protein,
    )
    import pandas as pd

    peptides = pd.read_csv("peptides.csv")

    # TopN
    top3 = TopNQuantification(n=3)
    result = top3.quantify(
        peptides,
        protein_column="ProteinName",
        peptide_column="PeptideSequence",
        intensity_column="NormIntensity",
        sample_column="SampleID",
    )
    ```

## Methods

### iBAQ

Intensity-Based Absolute Quantification. Divides summed peptide intensities by the number of theoretically observable peptides. **Requires a FASTA file**.

```bash
mokume peptides2protein --method ibaq \
    -f proteome.fasta \
    -p peptides.csv \
    -e Trypsin \
    --normalize \
    --output proteins-ibaq.tsv
```

#### Full iBAQ with TPA and ProteomicRuler

```bash
mokume peptides2protein \
    -f proteome.fasta \
    -p peptides.csv \
    -e Trypsin \
    --normalize \
    --tpa \
    --ruler \
    --ploidy 2 \
    --cpc 200 \
    --organism human \
    --output proteins-ibaq.tsv \
    --verbose \
    --qc_report QC.pdf
```

```python
from mokume.quantification import peptides_to_protein

peptides_to_protein(
    fasta="proteome.fasta",
    peptides="peptides.csv",
    enzyme="Trypsin",
    normalize=True,
    tpa=True,
    ruler=True,
    ploidy=2,
    cpc=200,
    organism="human",
    output="proteins-ibaq.tsv",
    min_aa=7,
    max_aa=30,
    verbose=True,
    qc_report="QC.pdf",
)
```

### TopN

Averages the N most intense peptides per protein per sample.

```bash
# Top3 (most common)
mokume peptides2protein --method top3 -p peptides.csv -o out.tsv

# Top5
mokume peptides2protein --method topn --topn_n 5 -p peptides.csv -o out.tsv

# Top10
mokume peptides2protein --method topn --topn_n 10 -p peptides.csv -o out.tsv
```

```python
from mokume.quantification import TopNQuantification

top3 = TopNQuantification(n=3)
top5 = TopNQuantification(n=5)
top10 = TopNQuantification(n=10)
```

### MaxLFQ

Delayed normalization with pairwise peptide ratios. Automatically uses DirectLFQ backend when installed, falling back to the built-in parallelized implementation.

```bash
mokume peptides2protein --method maxlfq \
    --threads 4 \
    -p peptides.csv \
    -o proteins-maxlfq.tsv
```

```python
from mokume.quantification import MaxLFQQuantification

maxlfq = MaxLFQQuantification(min_peptides=2, threads=4)
result = maxlfq.quantify(peptides)

# Check which backend is active
print(maxlfq.using_directlfq)  # True/False
print(maxlfq.name)             # "MaxLFQ (DirectLFQ)" or "MaxLFQ (built-in)"

# Force built-in implementation
maxlfq_builtin = MaxLFQQuantification(min_peptides=2, force_builtin=True)
```

### DirectLFQ

Uses hierarchical normalization with variance-guided pairwise alignment. Requires `pip install mokume[directlfq]`.

```bash
mokume peptides2protein --method directlfq \
    -p peptides.csv \
    -o proteins-directlfq.tsv
```

```python
from mokume.quantification import is_directlfq_available

if is_directlfq_available():
    from mokume.quantification import DirectLFQQuantification
    directlfq = DirectLFQQuantification(min_nonan=2)
    result = directlfq.quantify(peptides)
```

### Sum

Sums all peptide intensities per protein per sample.

```bash
mokume peptides2protein --method sum \
    -p peptides.csv \
    -o proteins-sum.tsv
```

## Factory Function

The `get_quantification_method` factory automatically parses method names:

```python
from mokume.quantification import get_quantification_method

method = get_quantification_method("top3")    # TopNQuantification(n=3)
method = get_quantification_method("top5")    # TopNQuantification(n=5)
method = get_quantification_method("maxlfq", min_peptides=2, threads=-1)

# Check available methods
from mokume.quantification import list_quantification_methods
print(list_quantification_methods())
# {'top3': True, 'topn': True, 'maxlfq': True, 'directlfq': False, 'sum': True}
```

## CLI Options Reference

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--fasta` | none | FASTA file (required for iBAQ) |
| `-p/--peptides` | required | Input peptide intensity file |
| `--method` | `ibaq` | Quantification method: ibaq, top3, topn, maxlfq, sum, directlfq |
| `-e/--enzyme` | `Trypsin` | Enzyme for in-silico digestion |
| `-n/--normalize` | off | Normalize quantification values |
| `--min_aa` | 7 | Minimum amino acid length |
| `--max_aa` | 30 | Maximum amino acid length |
| `-t/--tpa` | off | Calculate TPA (iBAQ only) |
| `-r/--ruler` | off | Use ProteomicRuler (iBAQ only) |
| `-i/--ploidy` | 2 | Ploidy number |
| `-m/--organism` | `human` | Organism for histone data |
| `-c/--cpc` | 200 | Cellular protein concentration (g/L) |
| `--topn_n` | 3 | N for TopN quantification |
| `--threads` | -1 | Threads for MaxLFQ (-1 = all cores) |
| `--min_nonan` | 1 | Min non-NaN values (DirectLFQ) |
| `-o/--output` | none | Output file path |
| `--verbose` | off | Print distribution info |
| `--qc_report` | QCprofile.pdf | Path for QC report PDF |

`-o/--output` is effectively required for `--method ibaq`; for the other methods, omitting it prints the result table to stdout.
