# Python API

`pip install mokume-rs` gives you a thin Python wheel over the Rust compute kernel, in the PyO3/maturin layout used by projects such as polars and pydantic-core (Python imports a compiled Rust extension). There is **no** rich class-based API any more: the compute numbers are single-sourced in Rust and exposed in-process through the `mokume._mokume` extension, and the Python periphery only reads the tables the kernel writes.

The package has two layers:

- **Compute wrappers** — `mokume.features2proteins(...)`, `mokume.features2peptides(...)`, `mokume.peptides2protein(...)`, `mokume.correct_batches(...)`, plus `mokume.run([...])` and `mokume.version()`. These run the same clap parsing + dispatch the standalone `mokume` binary uses, **in-process, no subprocess**.
- **Periphery** — plotting, tissue maps, DE plots, interactive reports, iBAQ QC, and the pure-Python method fallback (`missforest`). These live in `mokume.commands.*` / `mokume.reports.*` and are reached through the ergonomic wrappers below. Each needs an [install extra](#install-extras).

```python
import mokume

mokume.version()   # the kernel version string
```

---

## Compute wrappers

Each compute wrapper maps keyword arguments to CLI flags and runs the kernel in-process. The flags are exactly those documented in the [CLI Reference](cli.md) — the wrappers add no surface of their own.

```python
import mokume

# feature parquet -> protein matrix (+ optional DE)
mokume.features2proteins(parquet="features.parquet", output="proteins.csv")

# feature parquet -> peptide-level output
mokume.features2peptides(parquet="features.parquet", output="peptides.csv")

# peptide-level input -> protein quantities
mokume.peptides2protein(method="ibaq", peptides="peptides.parquet",
                        fasta="proteome.fasta", output="proteins.tsv")

# ComBat batch-effect correction on iBAQ output
mokume.correct_batches(folder="ibaq_dir", output="corrected.tsv")
```

### kwargs → flags rule

Each wrapper translates `**kwargs` into a CLI argument list:

| keyword form | becomes | example |
|--------------|---------|---------|
| `key=value` | `--key value` (`_` → `-`) | `quant_method="ibaq"` → `--quant-method ibaq` |
| `key=True` | `--key` (a bare flag) | `batch_correction=True` → `--batch-correction` |
| `key=[a, b]` | the flag repeated | `de=[...]` style list → flag once per item |
| `key=None` / `key=False` | skipped | omitted entirely |

```python
mokume.features2proteins(
    parquet="features.parquet",
    output="proteins.csv",
    sdrf="experiment.sdrf.tsv",
    quant_method="maxlfq",
    batch_correction=True,
    batch_covariates="characteristics[sex]",
    de=True,
    de_contrasts="NASH vs HL",
    duckdb_threads=24,
)
```

### `mokume.run([...])` — full control

When you need flags a keyword cannot express (e.g. a repeated `--contrast KEY A B CSV`), pass the argument vector verbatim. `run` accepts the subcommand name as the first element:

```python
mokume.run(["features2proteins", "--parquet", "x.parquet", "--output", "y.csv"])
mokume.run(["correct-batches", "--folder", "ibaq_dir", "--output", "corrected.tsv"])
```

`mokume.run` and the four wrappers raise on a dispatch failure and surface clap's usage errors; they never tear down the hosting interpreter.

---

## Periphery

The periphery reads the tables the kernel wrote — it never recomputes the numbers, so the cells in the plots match the cells in the kernel output. Each command lives in `mokume.commands.<name>` with a `main(argv)` entry point (runnable as `python -m mokume.commands.<name>`) and most have an ergonomic wrapper on the top-level package.

### Visualization

```python
import mokume

# t-SNE over a folder of protein files (plotting extra)
mokume.tsne_visualization(folder="./proteins", pattern="proteins.tsv")

# per-dataset tissue proteome analysis (tissuemap extra)
mokume.tissuemap(scan_dir="./data", output_dir="./out")

# iBAQ QC report from a protein table (plotting extra)
mokume.peptides2protein_qc(protein_table="proteins.tsv", qc_report="QC.pdf")
```

`de_plots` and `interactive_report` take an explicit argv (the per-contrast `--contrast KEY A B CSV` flag repeats, which keyword arguments cannot express):

```python
# DE volcano / heatmap / PCA from kernel-written CSVs (plotting extra)
mokume.de_plots(["--protein-matrix", "proteins.csv", "--plot-dir", "plots",
                 "--volcano", "--contrast", "c1", "A", "B", "de.csv"])

# interactive HTML report from kernel CSVs (reports extra)
mokume.interactive_report(["--protein-matrix", "proteins.csv", "--report-output", "report.html"])
```

Run `python -m mokume.commands.de_plots --help` / `python -m mokume.commands.interactive_report --help` for the flags.

### QC and workflow-comparison reports

```python
# single-matrix QC report: PCA / t-SNE / silhouette / CV / missing-value / DE-quality
path = mokume.qc_report(
    protein_matrix="proteins.csv",
    sdrf="experiment.sdrf.tsv",
    output="qc.html",
    de_results="de.csv",     # optional
)

# compare several quantification workflows in one HTML report
path = mokume.workflow_comparison(
    workflows=[
        {"name": "maxlfq", "protein_matrix": "maxlfq.csv", "sdrf": "x.sdrf.tsv"},
        {"name": "ibaq",   "protein_matrix": "ibaq.csv",   "sdrf": "x.sdrf.tsv"},
    ],
    output="comparison.html",
)
```

Both need the `analysis` extra. For volcano gene-highlighting, call `mokume.reports.qc_report.generate_qc_report` directly.

### Pure-Python method fallbacks

A method not reproducible bit-for-bit in the Rust kernel: the kernel's `features2proteins` errors point here (needs the `analysis` extra):

```python
# missforest — wraps scikit-learn's IterativeImputer
mokume.impute("proteins.csv", method="missforest", output="imputed.csv")
```

`mokume.impute` also reaches every other supported method (`knn`, `minprob`, `qrilc`, ...); it accepts a wide protein-matrix CSV path or a DataFrame and returns the imputed DataFrame, writing `output` if given.

### iBAQ for unported enzymes

The native iBAQ path digests proteins for the ported pyOpenMS enzymes (Trypsin[/P], Lys-C[/P], Arg-C[/P], Chymotrypsin[/P], Glu-C, Asp-N, Lys-N, PepsinA, ...). For any other enzyme pyOpenMS knows (CNBr, V8-DE, unspecific cleavage, ...) the kernel has no cleavage rule and points you here — the whole iBAQ table is then computed in pure Python (the `ibaq` extra):

```python
mokume.peptides2protein_ibaq(peptides="peptides.parquet", fasta="proteome.fasta",
                             enzyme="CNBr", output="proteins.tsv")
```

---

## Install extras { #install-extras }

The compute path (the `mokume._mokume` extension) needs **no** third-party Python dependencies. Install only the extra for the periphery command you run:

```bash
pip install mokume-rs                 # compute kernel + Python API
pip install "mokume-rs[plotting]"     # + t-SNE / DE plots / iBAQ QC report
pip install "mokume-rs[tissuemap]"    # + per-dataset tissue proteome analysis
pip install "mokume-rs[reports]"      # + interactive HTML DE report
pip install "mokume-rs[ibaq]"         # + pure-Python iBAQ for unported enzymes
pip install "mokume-rs[analysis]"     # + QC / comparison reports + missforest
pip install "mokume-rs[all]"          # everything
```

| Wrapper | Extra | Third-party libraries |
|---------|-------|------------------------|
| `mokume.tsne_visualization` | `plotting` | numpy, pandas, scipy, scikit-learn, matplotlib, seaborn |
| `mokume.peptides2protein_qc` | `plotting` | numpy, pandas, matplotlib, seaborn |
| `mokume.de_plots` | `plotting` | numpy, pandas, matplotlib, seaborn, scikit-learn |
| `mokume.interactive_report` | `reports` | numpy, pandas, plotly |
| `mokume.tissuemap` | `tissuemap` | scanpy, anndata, umap-learn, combat, matplotlib, seaborn, pyarrow |
| `mokume.peptides2protein_ibaq` | `ibaq` | pyopenms, pyarrow, PyYAML, numpy, pandas, scipy |
| `mokume.qc_report` / `mokume.workflow_comparison` | `analysis` | numpy, pandas, scipy, scikit-learn |
| `mokume.impute` | `analysis` | numpy, pandas, scipy, scikit-learn |

The exact dependency lists are declared in `pyproject.toml`'s `[project.optional-dependencies]`. The retired `directlfq` and `batch-correction` extras are gone: DirectLFQ and ComBat are now native Rust and need no extra.

!!! note "Agentic workflows are a separate package"
    The agentic / LLM-driven workflow layer is **not** part of this wheel; it lives in the separate `mokume_py` package.
