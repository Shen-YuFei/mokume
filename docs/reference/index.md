# Reference

Reference documentation for mokume's CLI commands, the wheel's Python API, output columns, and the periphery commands (TissueMap, t-SNE, DE plots, interactive reports).

## [CLI Reference](cli.md)

The four compute subcommands of the `mokume` binary and their options: `features2proteins`, `features2peptides`, `peptides2protein`, and `correct-batches`.

## [Python API](python-api.md)

The thin `pip install mokume-rs` wheel: the in-process compute wrappers (`mokume.features2proteins(...)`, `mokume.run([...])`, ...) and the Python periphery functions (t-SNE, tissue maps, DE plots, interactive reports, QC, the missforest fallback).

## [Computed Values](computed-values.md)

Column names and formulas for all output values (iBAQ, TopN, MaxLFQ, TPA, CopyNumber, etc.).
