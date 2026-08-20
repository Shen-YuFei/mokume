# peptides2protein: Protein Quantification

The `peptides2protein` command quantifies proteins from normalized peptide data. It supports multiple quantification methods and is the second step of the two-step pipeline.

## Basic Usage

=== "CLI"

    ```bash
    # piBAQ (default, requires FASTA)
    mokume peptides2protein --method pibaq \
        -f proteome.fasta \
        -p peptides.csv \
        -o proteins-pibaq.tsv

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

=== "Python (wheel)"

    The wheel wrapper maps keyword arguments to CLI flags (`key=value` → `--key value` with `_` rewritten to `-`; `key=True` → `--key`) and runs the same kernel in-process:

    ```python
    import mokume

    # piBAQ (requires FASTA)
    mokume.peptides2protein(
        method="pibaq",
        fasta="proteome.fasta",
        peptides="peptides.csv",
        output="proteins-pibaq.tsv",
    )

    # TopN (no FASTA needed)
    mokume.peptides2protein(method="top3", peptides="peptides.csv",
                            output="proteins-top3.tsv")
    ```

## Methods

### piBAQ

piBAQ is Mokume's paralog-aware extension of iBAQ. For each sample,
shared-peptide signal is allocated in proportion to mapped members'
proteotypic-peptide intensities; if all mapped members lack that signal, it is
split equally. Each shared intensity is allocated exactly once. See
[Quantification Methods](../concepts/quantification.md#pibaq-paralog-aware-ibaq)
for the underlying algorithm. **Requires a FASTA file**.

```bash
mokume peptides2protein --method pibaq \
    -f proteome.fasta \
    -p peptides.csv \
    -e Trypsin \
    --normalize \
    --output proteins-pibaq.tsv
```

The output adds three metadata columns -- `FamilyId`, `FamilySize`, `EvidenceLevel` -- so users can audit family support. `family_only` means that no family member reaches the minimum anchor threshold; it does not duplicate one family-level piBAQ across all members. Every reported piBAQ remains a per-member estimate with a per-member theoretical-peptide denominator.

#### Family Discovery Tuning

Families are discovered automatically by collapsing UniProt isoform suffixes (`-2`, `-3`, ...) onto the canonical entry and then grouping proteins on the shared-peptide graph. Two CLI flags tune the auto-discovery; both are optional and rarely need adjustment.

```bash
# Lower the shared-peptide threshold for very tightly homologous families
mokume peptides2protein --method pibaq -f proteome.fasta -p peptides.csv \
    --min-shared 1 -o out.tsv

# Pin specific families with an audit-friendly YAML override
mokume peptides2protein --method pibaq -f proteome.fasta -p peptides.csv \
    --families families.yaml -o out.tsv
```

The YAML schema is:

```yaml
families:
  - name: ACT
    members: [P60709, P63261, P68133]   # canonical accessions only
  - name: HIST_H2A
    members: [P0C0S5, Q96QV6, P04908]
```

#### Full piBAQ with TPA and ProteomicRuler

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
    --output proteins-pibaq.tsv \
    --verbose \
    --qc_report QC.pdf
```

```python
import mokume

mokume.peptides2protein(
    fasta="proteome.fasta",
    peptides="peptides.csv",
    enzyme="Trypsin",
    normalize=True,
    tpa=True,
    ruler=True,
    ploidy=2,
    cpc=200,
    organism="human",
    output="proteins-pibaq.tsv",
    min_aa=7,
    max_aa=30,
)
```

!!! note "piBAQ enzyme coverage and the QC report"
    piBAQ reads the complete protease catalog from the installed pyOpenMS runtime,
    including context-dependent rules, unspecific cleavage, and no cleavage. Python
    supplies the theoretical-peptide map; Rust performs shared-peptide allocation,
    family handling, denominators, TPA, normalization, and output. The `--qc_report`
    PDF is plotting periphery: `--verbose` prints a one-line pointer to
    `mokume.peptides2protein_qc` (`plotting` extra).

### TopN

Averages the N most intense peptides per protein per sample.

```bash
# Top3 (most common) -- the named method from Silva et al. 2006
mokume peptides2protein --method top3 -p peptides.csv -o out.tsv

# Top5
mokume peptides2protein --method top5 -p peptides.csv -o out.tsv

# Top10
mokume peptides2protein --method top10 -p peptides.csv -o out.tsv
```

N is spelled in the method name, so any N works the same way and there is no
companion option to keep in sync.

```python
import mokume

mokume.peptides2protein(method="top3", peptides="peptides.csv", output="out.tsv")
mokume.peptides2protein(method="top5", peptides="peptides.csv", output="out.tsv")
```

### MaxLFQ

Delayed normalization with pairwise peptide ratios. `maxlfq` rolls the peptide matrix up with the native Rust DirectLFQ estimator (delegating to it with `min_nonan = 2`).

```bash
mokume peptides2protein --method maxlfq \
    --threads 4 \
    -p peptides.csv \
    -o proteins-maxlfq.tsv
```

```python
import mokume

mokume.peptides2protein(method="maxlfq", threads=4, peptides="peptides.csv",
                        output="proteins-maxlfq.tsv")
```

### DirectLFQ

Uses hierarchical normalization with variance-guided pairwise alignment, native in the Rust kernel (no extra dependency).

```bash
mokume peptides2protein --method directlfq \
    --min_nonan 2 \
    -p peptides.csv \
    -o proteins-directlfq.tsv
```

```python
import mokume

mokume.peptides2protein(method="directlfq", min_nonan=2, peptides="peptides.csv",
                        output="proteins-directlfq.tsv")
```

### Sum

Sums all peptide intensities per protein per sample.

```bash
mokume peptides2protein --method sum \
    -p peptides.csv \
    -o proteins-sum.tsv
```

```python
import mokume

mokume.peptides2protein(method="sum", peptides="peptides.csv",
                        output="proteins-sum.tsv")
```

## CLI Options Reference

| Option | Default | Description |
|--------|---------|-------------|
| `-f/--fasta` | none | FASTA file (required for piBAQ) |
| `-p/--peptides` | required | Input peptide intensity file |
| `--method` | `pibaq` | Quantification method: pibaq, `top<N>` (for example top3 or top5), maxlfq, sum, directlfq |
| `-e/--enzyme` | `Trypsin` | Enzyme for in-silico digestion |
| `-n/--normalize` | off | Normalize quantification values |
| `--min_aa` | 7 | Minimum amino acid length |
| `--max_aa` | 30 | Maximum amino acid length |
| `-t/--tpa` | off | Calculate TPA (piBAQ only) |
| `-r/--ruler` | off | Use ProteomicRuler (piBAQ only) |
| `-i/--ploidy` | 2 | Ploidy number |
| `-m/--organism` | `human` | Organism for histone data |
| `-c/--cpc` | 200 | Cellular protein concentration (g/L) |
| `--threads` | -1 | Threads for MaxLFQ (-1 = all cores) |
| `--min_nonan` | 1 | Min non-NaN values (DirectLFQ) |
| `--families` | none | Optional YAML file with explicit family overrides (piBAQ only) |
| `--min-shared` | 2 | Minimum shared peptides for auto-family discovery (piBAQ only) |
| `-o/--output` | none | Output file path |
| `--verbose` | off | Print distribution info |
| `--qc_report` | QCprofile.pdf | Path for QC report PDF |

`-o/--output` is effectively required for `--method pibaq`; for the other methods, omitting it prints the result table to stdout.
