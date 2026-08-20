# Python API

`pip install mokume` gives you a thin Python wheel over the Rust compute
kernel, in the PyO3/maturin layout used by projects such as polars and
pydantic-core (Python imports a compiled Rust extension). This wheel does not
expose the separately installed pure-Python package's rich class-based API. Its
kernel-supported compute runs in-process through the `mokume._mokume` extension;
plotting and reporting read the kernel's tables, TissueMap derives downstream
atlas outputs from QPX data, and explicitly documented Python-only fallbacks
compute operations the kernel does not provide.

The package has two layers:

- **Compute wrappers** — full commands such as `mokume.features2proteins(...)`
  plus matrix-level `normalize_matrix`, `impute_matrix`, and
  `differential_expression` calls. They run **in-process, with no subprocess**.
- **Periphery and fallbacks** — plotting, tissue maps, DE plots, interactive
  reports, piBAQ QC, and explicit pure-Python fallbacks such as `missforest`.
  These live in `mokume.commands.*` / `mokume.reports.*` and are reached through
  the ergonomic wrappers below. Each needs an
  [install extra](#install-extras).

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
mokume.peptides2protein(method="pibaq", peptides="peptides.parquet",
                        fasta="proteome.fasta", output="proteins.tsv")

# ComBat batch-effect correction on piBAQ output
mokume.correct_batches(folder="pibaq_dir", output="corrected.tsv")
```

### kwargs → flags rule

Each wrapper translates `**kwargs` into a CLI argument list:

| keyword form | becomes | example |
|--------------|---------|---------|
| `key=value` | `--key value` (`_` → `-`) | `quant_method="pibaq"` → `--quant-method pibaq` |
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
mokume.run(["correct-batches", "--folder", "pibaq_dir", "--output", "corrected.tsv"])
```

`mokume.run` and the four wrappers raise on a dispatch failure and surface clap's usage errors; they never tear down the hosting interpreter.

---

## Matrix-level compute

These calls reuse the same Rust implementations as the full pipeline without
rerunning QPX loading or protein aggregation. Matrices are row-major linear
intensities (`values[protein][sample]`); use `None` or a non-finite float for a
missing cell. Every call accepts an explicit thread count.

```python
import mokume

values = [
    [100.0, 120.0, None, 240.0],
    [400.0, 420.0, 800.0, 820.0],
]

normalized = mokume.normalize_matrix(
    values, "median", ["A1", "A2", "B1", "B2"], threads=24
)
imputed = mokume.impute_matrix(normalized, "mindet", threads=24)
de_rows = mokume.differential_expression(
    ["P1", "P2"],
    imputed,
    2,
    2,
    "limma",
    condition_a="A",
    condition_b="B",
    threads=24,
)
```

`normalize_matrix` supports `none`, `median`, `mean`, `quantile`, `rlr`,
`loess`, `hierarchical`, and `tmm`. `impute_matrix` supports every native Rust
imputer; `missforest` remains the explicit Python fallback. The DE call supports
`limma`, `deqms`, `rots`, `limrots`, `proda`, and `ensemble`, returning a list
of dictionaries with the same result-column names as the corresponding command
table. Because matrix-level DE has no quantification context, it requires a
concrete method and rejects `auto`.

---

## Periphery

The plotting and reporting periphery reads the tables the kernel wrote without
recomputing them, so the cells in its plots match the kernel output. TissueMap
derives downstream atlas outputs from QPX data, while explicit fallbacks compute
operations unsupported by the Rust kernel. Each command lives in
`mokume.commands.<name>` with a `main(argv)` entry point (runnable as
`python -m mokume.commands.<name>`) and most have an ergonomic wrapper on the
top-level package.

### Visualization

```python
import mokume

# t-SNE over a folder of protein files (plotting extra)
mokume.tsne_visualization(folder="./proteins", pattern="proteins.tsv")

# per-dataset tissue proteome analysis (tissuemap extra)
mokume.tissuemap(scan_dir="./data", output_dir="./out")

# piBAQ QC report from a protein table (plotting extra)
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
        {"name": "pibaq",   "protein_matrix": "pibaq.csv",   "sdrf": "x.sdrf.tsv"},
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

### piBAQ for unported enzymes

The native piBAQ path digests proteins for the ported pyOpenMS enzymes (Trypsin[/P], Lys-C[/P], Arg-C[/P], Chymotrypsin[/P], Glu-C, Asp-N, Lys-N, PepsinA, ...). For any other enzyme pyOpenMS knows (CNBr, V8-DE, unspecific cleavage, ...) the kernel has no cleavage rule and points you here — the whole piBAQ table is then computed in pure Python (the `pibaq` extra):

```python
mokume.peptides2protein_pibaq(peptides="peptides.parquet", fasta="proteome.fasta",
                             enzyme="CNBr", output="proteins.tsv")
```

---

## Install extras { #install-extras }

The compute path (the `mokume._mokume` extension) needs **no** third-party Python dependencies. Install only the extra for the periphery command you run:

```bash
pip install mokume                 # compute kernel + Python API
pip install "mokume[plotting]"     # + t-SNE / DE plots / piBAQ QC report
pip install "mokume[tissuemap]"    # + per-dataset tissue proteome analysis
pip install "mokume[reports]"      # + interactive HTML DE report
pip install "mokume[pibaq]"         # + pure-Python piBAQ for unported enzymes
pip install "mokume[analysis]"     # + QC / comparison reports + missforest
pip install "mokume[agentic]"     # + local MCP service for the Mokume Plugin
pip install "mokume[all]"          # everything
```

| Wrapper | Extra | Third-party libraries |
|---------|-------|------------------------|
| `mokume.tsne_visualization` | `plotting` | numpy, pandas, scipy, scikit-learn, matplotlib, seaborn |
| `mokume.peptides2protein_qc` | `plotting` | numpy, pandas, matplotlib, seaborn |
| `mokume.de_plots` | `plotting` | numpy, pandas, matplotlib, seaborn, scikit-learn |
| `mokume.interactive_report` | `reports` | numpy, pandas, plotly |
| `mokume.tissuemap` | `tissuemap` | scanpy, anndata, umap-learn, combat, matplotlib, seaborn, pyarrow |
| `mokume.peptides2protein_pibaq` | `pibaq` | pyopenms, pyarrow, PyYAML, numpy, pandas, scipy |
| `mokume.qc_report` / `mokume.workflow_comparison` | `analysis` | numpy, pandas, scipy, scikit-learn |
| `mokume.impute` | `analysis` | numpy, pandas, scipy, scikit-learn |
| Mokume Plugin MCP service | `agentic` | mcp, numpy, pandas, scipy, scikit-learn, statsmodels, PyYAML |

The exact dependency lists are declared in `pyproject.toml`'s `[project.optional-dependencies]`. The retired `directlfq` and `batch-correction` extras are gone: DirectLFQ and ComBat are now native Rust and need no extra.

!!! note "Agentic reasoning belongs to the host"
    The wheel provides deterministic MCP tools, not a model client. Install and
    enable the [Mokume Plugin](../user-guide/agentic-plugin.md); the host starts
    the local service and keeps ownership of its model credentials.
