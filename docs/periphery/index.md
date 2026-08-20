# Python Periphery

The periphery is the Python half of the toolkit: plotting, tissue maps,
interactive reports, the piBAQ QC report, the pure-Python piBAQ path for unported
enzymes, and the pure-Python method fallbacks. It lives in
`rust/python/mokume/commands/` (plus `rust/python/mokume/reports/`,
`rust/python/mokume/normalization/`, `rust/python/mokume/imputation/` for the
analysis fallbacks) and is reached **only** through the
`pip install mokume` wheel. The Rust compute extension does not import this
periphery.

## Reporting reads kernel output; fallbacks are explicit

Plotting, QC, and reporting commands consume tables the **kernel** wrote without
re-running kernel-supported computation. TissueMap instead derives its
documented downstream normalization, batch correction, tissue-specificity, and
atlas outputs from QPX data.

Two explicit fallbacks compute operations the Rust kernel does not provide:
`mokume.peptides2protein_pibaq` handles enzymes outside the Rust-ported set, and
`mokume.impute(method='missforest')` runs the scikit-learn estimator. They are
documented under [Rust Wheel](../rust-wheel.md) and
[Analysis Fallbacks](analysis-fallbacks.md), rather than being presented as
kernel-produced results.

## Extras matrix

The compute kernel (`mokume._mokume`) needs **no** third-party Python
dependencies. The extras pull in only the libraries a given periphery command
needs. Install just the extra for the command you run.

| Command | Extra | Third-party libraries |
| --- | --- | --- |
| `mokume.tsne_visualization` | `plotting` | numpy, pandas, scipy, scikit-learn, matplotlib, seaborn |
| `mokume.peptides2protein_qc` | `plotting` | numpy, pandas, matplotlib, seaborn |
| `mokume.de_plots` | `plotting` | numpy, pandas, matplotlib, seaborn, scikit-learn |
| `mokume.interactive_report` | `reports` | numpy, pandas, plotly |
| `mokume.tissuemap` | `tissuemap` | scanpy, anndata, umap-learn, combat, matplotlib, seaborn, pyarrow |
| `mokume.peptides2protein_pibaq` | `pibaq` | pyopenms, pyarrow, PyYAML, numpy, pandas, scipy |
| `mokume.qc_report` / `mokume.workflow_comparison` | `analysis` | numpy, pandas, scipy, scikit-learn |
| `mokume.impute` | `analysis` | numpy, pandas, scipy, scikit-learn |

```bash
pip install mokume                 # compute kernel + Python API
pip install "mokume[plotting]"     # + t-SNE / DE plots / piBAQ QC report
pip install "mokume[reports]"      # + interactive HTML DE report
pip install "mokume[tissuemap]"    # + per-dataset tissue proteome analysis
pip install "mokume[pibaq]"         # + pure-Python piBAQ for unported enzymes
pip install "mokume[analysis]"     # + QC/comparison reports + missforest
pip install "mokume[all]"          # everything
```

The exact dependency lists are declared in `pyproject.toml`'s
`[project.optional-dependencies]`.

!!! note "Removed extras"
    The old `directlfq` and `batch-correction` extras are **gone**. DirectLFQ and
    ComBat are now native Rust in the kernel and need no extra and no third-party
    Python dependency.

## The import-then-call pattern

Each periphery command lives in `mokume.commands.<name>` with a `main(argv)`
entry point, is runnable as `python -m mokume.commands.<name>`, and most have an
ergonomic wrapper on the top-level `mokume` package. Import the package and call
the wrapper:

```python
import mokume

# kwarg wrappers (the kwargs -> flags rule applies):
mokume.tsne_visualization(folder="./proteins", pattern="proteins.tsv")
mokume.tissuemap(scan_dir="./data", output_dir="./out")
mokume.peptides2protein_qc(protein_table="proteins.tsv", qc_report="QC.pdf")
mokume.peptides2protein_pibaq(peptides="peptides.parquet", fasta="proteome.fasta",
                             enzyme="CNBr", output="proteins.tsv")

# de_plots / interactive_report take an explicit argv (the per-contrast
# --contrast KEY A B CSV flag repeats, which keyword arguments cannot express):
mokume.de_plots(["--protein-matrix", "proteins.csv", "--plot-dir", "plots",
                 "--volcano", "--contrast", "c1", "A", "B", "de.csv"])

# QC / comparison reports + the pure-Python method fallbacks (analysis extra):
mokume.qc_report(protein_matrix="proteins.csv", sdrf="x.sdrf.tsv", output="qc.html")
mokume.impute("proteins.csv", method="missforest", output="imputed.csv")
```

The kwarg wrappers follow the [kwargs &rarr; flags rule](../rust-wheel.md)
(`key=value` &rarr; `--key value` with `_` rewritten to `-`; `key=True` &rarr;
`--key`; a list repeats the flag; `None` / `False` skipped). `de_plots` and
`interactive_report` take a literal argument list because their per-contrast
flag repeats in a shape keyword arguments cannot express.

## What lives where

- **Tissue Proteome Atlas** &rarr; [`mokume.tissuemap`](tissuemap.md) — PCA /
  t-SNE / UMAP / AdaTiSS tissue-specificity scoring, markers, atlas plots.
- **Visualization & Reports** &rarr;
  [t-SNE / DE plots / interactive report](visualization.md).
- **Analysis Fallbacks** &rarr; [`mokume.impute`,
  `mokume.qc_report`, `mokume.workflow_comparison`](analysis-fallbacks.md) —
  methods the Rust kernel does not reproduce, plus the QC / comparison reports.
