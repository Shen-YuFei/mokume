# Rust Refactor Plan

> **Status — historical planning record.** The authoritative current state is the
> "Current Status" section below. The per-phase sections and the
> "Current Conclusion" / "Key Differences" sections further down predate the
> completed DE / ComBat / batch-correction / channel-IRS build-out: they still
> describe as "not yet supported" capabilities that have since shipped and been
> parity-verified. Read them as history, not as the current gap list. This file is
> excluded from the published docs site (`exclude_docs: internal/*`).

This document records the port of the upstream Python `mokume` implementation to
the Rust `mokume` compute kernel (now shipped in a PyO3/maturin wheel) and any
remaining gaps. The goal is better system performance while
keeping behavior verifiable. Agentic recommendation is an installable host
plugin whose local MCP service calls the Rust-backed matrix APIs.

## Current Status (2026-08-19 — authoritative; per-phase sections below may be stale)

A full method-level audit plus the DE/ComBat build-out and a real-data parity pass
align the computation core: every quant/normalization/imputation/DE/batch method is
either real-data-verified against Python or a documented cross-language-infeasible
gap. The per-phase sections further down predate this and under-report progress;
trust this section and `features2proteins-parity.md`.

Aligned now:

- Quant: `sum`, `median`, `topn`/`top3`, `intensity`, `abd`, `spectral_count`,
  `directlfq` and `maxlfq` (both delegate to the Mann-Labs-aligned DirectLFQ path
  by default, matching Python; real-data cell-exact on PXD003539 — 0/0 protein
  sets, max_rel ~6e-15, per-sample Spearman 1.0 across all 120 samples — after the
  maxlfq per-cell min-unique gate fix, the directlfq first-accession naming fix, and
  the ion-row-order / NaN-linear-shift alignment), and the piBAQ
  family-allocation core. `maxlfq --force-builtin` selects the built-in MaxLFQ.
- Sample normalization: `none`, `globalmedian`, `conditionmedian`, `quantile`,
  `mediancenter`, `meancenter`, `rlr`, `loess`, `hierarchical`.
  `mediancenter`/`meancenter`/`hierarchical` are now
  real-data cell-exact on PXD003539 after a dataset-level centering / pivot
  column-order fix.
- Imputation: `none`, `mindet`, `minprob`, `mean`, `median`, `constant`, `zero`,
  `most_frequent`, `knn`, `seqknn`, `qrilc`, `impseq`, `gms`, `bpca`, `impseqrob`
  (`missforest` unimplemented).
- Differential expression: the full catalog is wired and reachable from the
  `features2proteins --de --de-contrasts "A vs B" --de-method <m>` dispatcher
  (conditions resolve from the SDRF `source name` and `comment[data file]` stem,
  so run-level matrices work; real-data validated on PXD004701): `limma`,
  `deqms` (log2FC cell-exact, count-aware sca_t ~5e-3), `rots`,
  `limrots`, `proda` (deterministic kernels cell-exact; the bootstrap/optimizer
  paths are faithful-not-bit-exact per the accepted tiered-tolerance decision —
  Python's own ROTS output is seed-unstable), and `ensemble` (deterministic
  combine cell-exact). Contrasts come from inline `--de-contrasts "A vs B"` and/or
  the TSV `--de-contrasts-file` (group1/group2 columns, appended after the inline
  entries); `--de-method auto` resolves to `deqms` for directlfq, otherwise
  `limrots` (Python's `_resolve_de_method`). Batch correction: parametric **plus covariate
  (covar_mod) and non-parametric (par_prior=false)** ComBat, with `ref_batch`
  and `mean_only`, is wired into both the standalone `correct-batches` command and
  the `features2proteins` pipeline (`--batch-correction`, sample-prefix detection
  with the covariate-free path). Batch *detection* now ports the full Python set —
  `sample_prefix`, `run_name`, `fraction`, `techreplicate`, and explicit-column
  (`mokume-stats` `batch::detection`, matching `model/batch_correction.py` and the
  `stages.py` resolver), plus SDRF covariate extraction and a deterministic
  PCA+HDBSCAN outlier-removal pass (`batch::{pca,hdbscan,outlier}`, PCA bit-aligned
  to sklearn `svd_flip(u_based_decision=False)`, HDBSCAN EOM matching sklearn on
  continuous data; the outlier-removal pass mirrors Python's library function,
  which `combiner.py` imports but the main pipeline does not call).
  `--export-anndata` (h5ad) is implemented in `correct-batches`: a Rust-native
  `.h5ad` writer (`crates/mokume-command/src/h5ad.rs`, via the statically vendored
  `hdf5-metno`) whose layout/dtypes/values match Python `anndata.write_h5ad`.
- `mokume-stats` is now a Cargo dependency of `mokume-pipeline`, so the limma /
  DE catalog / ComBat libraries are reachable (no longer dead code).

Recently closed (2026-06-21): `ratio` plex/reference now uses an SDRF-driven
`sample_to_plex` map (per-plex reference, matching Python); `maxlfq` delegates to
directlfq with a per-(protein,sample) min-unique gate (real-data exact protein
sets); `directlfq` proteins are named by first accession (matching Python
`pg_accessions[0]`, real-data 0/0 protein-set diff); `IHW` multiple-testing is
ported (covariate binning + weight optimization + weighted BH, cell-exact vs
Python); `features2peptides` `--remove_ids` / `--remove_low_frequency_peptides` /
peptide-level normalization are implemented.

Recently closed (2026-06-27 — `features2peptides` C-stage parity, six commits on top
of the IRS-global work): dataset-level sample normalization (`quantile`/`rlr`/
`loess`/`hierarchical`/`mediancenter`/`meancenter`) is accepted and emits the
same deterministic no-op matrix Python produces (the registered replicate fns return
the frame unchanged and there is no post-loop dataset pass in the `features2peptides`
command — the dataset pivot lives only in `features2proteins`); the CV-threshold and
replicate-agreement intensity group filters are ported (joining the already-ported
quantile / Run-QC / per-row / razor filters, so the whole `--filter-config` pipeline
is live, with the single-sample `dataset_df` semantics mirrored — replicate-agreement
is a degenerate whole-sample wipe); channel IRS gained the `by_mixture` / `two_stage`
scopes and the SDRF source-name autodetect path on top of the `global` scope
(`collect_irs_scale` + scope functions + `resolve_irs_autodetect_channel`,
function-level oracle-locked to Python `get_irs_scaling_factors`); the keep-shared
median pre-pass now includes shared peptides (`require_unique = not keep_shared`,
fixing a 1.28e-2 → 2.31e-7 real-data bias); the no-SDRF `Condition` fallback now
mirrors Python's `run_file_name` (not `"Empty"`); and the directlfq / maxlfq high-ion
residual is resolved to cell-exact (see caveats). All six are golden-tested and
PXD003539 real-data verified (0/0 keyed diff, NormIntensity max_rel <= 2.5e-7;
directlfq/maxlfq max_rel ~6e-15, per-sample Spearman 1.0).

Recently closed (2026-06-28 — CDEF-external surface, four commits): the deterministic
pieces were ported to Rust — advanced batch *detection* (run/fraction/techreplicate/
explicit-column + SDRF covariates + PCA+HDBSCAN outlier removal, golden-tested vs
sklearn) and the native `.h5ad`/AnnData writer for `correct-batches --export-anndata`
(`hdf5-metno` static build, X cell-exact vs Python `anndata`, requires a C compiler +
cmake at build time). The visualization/QC surface ships in the `mokume` Python wheel (PyO3/maturin):
the periphery lives in `python/mokume/commands/` and reads the kernel's
TSV/parquet output to render figures, so the numbers stay single-sourced in Rust.
The former `tsne_visualization` / `tissuemap` subcommands and the
`features2proteins --plot-*` / `--interactive-report` flags moved to wheel APIs
(`mokume.tsne_visualization(...)`, `mokume.de_plots(...)`, etc.), so they are no
longer part of the Rust compute command surface. Meanwhile,
`peptides2protein --verbose` writes the numeric TSV while the wheel renders the QC
PDF (`mokume.peptides2protein_qc`). The `tsne_visualization` command carries a
small sklearn-version shim (`n_iter`→`max_iter`, removed upstream in sklearn 1.5+)
so it runs on current sklearn. Agentic reasoning is now host-owned; the wheel
provides a deterministic MCP service rather than another model or compute core.

Recently closed (2026-06-28b — gap-audit follow-ups, three commits): a 5-agent rs<->py
gap audit found the remaining user-facing gaps were peripheral output. They were first
delivered via vendored Python delegation and then moved into the `mokume` wheel by the
PyO3/maturin migration (see the lead paragraph above); the Rust compute command surface no longer carries them.
`features2proteins` DE plotting (`--plot-volcano`/`--plot-heatmap`/`--plot-pca`/
`--highlight-genes`) and the interactive HTML report (`--interactive-report`) became the
wheel's `mokume.de_plots` / `mokume.interactive_report` (Rust writes the protein matrix +
per-contrast DE CSVs, the wheel renders the figures/HTML; verified end-to-end producing real
PNGs + a plotly HTML, multi-contrast file topology byte-matching Python). `peptides2protein
--method pibaq` with an enzyme outside the natively-ported set became the wheel's pure-Python
`mokume.peptides2protein_pibaq` (CNBr / unspecific-cleavage / V8-DE byte-identical to Python),
and the default `Trypsin` path stayed Rust-native (zero regression). That interim
split has since been superseded: both piBAQ commands now digest through the complete
installed pyOpenMS runtime catalog and send the resulting theoretical-peptide map to
the Rust aggregation kernel; there is no unported-enzyme branch or `pibaq` extra.
Before deleting the native cleavage rules, a migration gate confirmed that all 21
legacy enzyme names were present in the 33-entry pyOpenMS catalog and compared their
peptide sets across 582 sequence cases and six length-bound combinations (73,332
digests, zero missed cleavages). Every result matched exactly. The duplicate rules
and the one-time verifier were then removed instead of becoming a second maintained
digestion implementation. A post-migration whole-catalog gate then compared the
Rust-backed and full-Python piBAQ tables for all 33 proteases: protein/sample keys,
condition and family metadata matched exactly, while `NormIntensity` and `PiBAQ`
matched within `1e-12` absolute and relative tolerance.
The audit also re-confirmed `deqms` (see caveats). Intentionally deferred as
library-only (Python exposes no CLI for them either): the iterative PCA+HDBSCAN
outlier-removal pass beyond the implemented kernel, and the multi-study `Combiner`
orchestration.

Remaining gaps (irreducible or niche):

1. Computation methods the Rust kernel does not reproduce natively but that now
   ship in the `mokume` wheel (`pip install mokume[analysis]`): `missforest`
   (wraps scikit-learn `IterativeImputer` + `RandomForestRegressor`, whose
   tree internals + RNG cannot be reproduced cross-language) via
   `mokume.impute(matrix, method=...)`. The kernel's `features2proteins` returns a
   stable `NotImplemented` whose message points at these wheel functions. The
   single-matrix QC report and the multi-workflow comparison report were never a
   porting target and likewise ship in the wheel (`mokume.qc_report` /
   `mokume.workflow_comparison`, `analysis` extra). `postprocessing/reshape.py`'s three library-only
   helpers and the group-level run-QC filters remain unported (see the parity doc).

Honest caveats: `deqms` log2FC and the quantitative matrix are cell-exact vs Python;
the per-protein peptide *counts* (and thus the `sca_t` moderation) diverge, but the
2026-06-28 audit showed Rust is the *more correct* side here — Python's
`_load_de_peptide_counts` groups the *un-parsed* `anchor_protein` (`stages.py:1795`)
while the matrix names are parsed (`aggregation.py:207`), so under the standard
`sp|ACC|NAME` format the count reindex never matches and silently falls back to 1
(degrading DEqMS to plain eBayes); Rust uses the parsed group name on both sides and
keeps real count-aware moderation. `proda` dropout-heavy proteins can diverge end-to-end (EM fixed-point path sensitivity);
RNG-driven DE methods are reproducible run-to-run via an in-crate splitmix64 PRNG but
are not bit-exact against numpy. The earlier `maxlfq`/`directlfq` high-ion residual
tail is resolved — it was not algorithmic chaos but two deterministic divergences (ion
row order keyed on the integer peptide id instead of Python's `sort(["protein",
"sequence"])` lexicographic order, and a skipped NaN linear-shift distance mask). After
aligning both, PXD003539 directlfq and maxlfq are cell-exact py==rust (0/0 protein and
sample sets, max_rel ~6e-15, per-sample Spearman 1.0 on all 120 samples).

Real-data parity (Wave 1-C, `/tmp/mokume-parity/`, not committed): `scripts/run_parity_matrix.sh`
ran the Phase-1 deterministic matrix across all four cell-line datasets
(PXD003539/004701/041421/030304) — 24/24 combos cell-exact (rel max 0.0 for
sum/median/topn/runmedian/conditionmedian; ~5e-8 float-level for globalmedian;
protein/sample sets exact). `maxlfq` and `directlfq` real-data parity on PXD003539:
cell-exact (0/0 protein and sample sets, max_rel ~6e-15, per-sample Spearman 1.0 on
all 120 samples) after the ion-row-order / NaN-linear-shift alignment. DE methods real-data
verified on PXD004701 (limma/deqms cell-exact, 100% significance-call
agreement; rots/limrots rank-level). `ratio` has no real-data parity because the
cell-line datasets are label-free (no plex data) — it is covered by a synthetic
plexed golden test cell-exact vs Python `RatioQuantification`.

## Repository and Commit Status

The Rust workspace is now a local git repository on branch `dev` (not pushed).
Commit history so far:

- `a9ae3b4` — baseline import of the workspace.
- `f74011b` — Wave 0 modularization.
- `0b875ee` — canonical-peptide collapse fix.
- `de5ed65` — quantile-normalization fix.

## Current Conclusion

The Rust implementation is no longer an empty project. The current `mokume`
workspace already provides:

- The `mokume` CLI entry point and the `features2proteins` main command.
- Streaming QPX parquet reading, with support for expanding nested
  `intensities`.
- SDRF reading and `run_file_name + label` sample matching.
- Integer-ID registries for proteins, peptides, samples, runs, and ions.
- `features2proteins` protein-matrix CSV output.
- Quantification methods: `sum`, `median`, `topn`, `maxlfq`, `directlfq`,
  `pibaq`, `ratio`, `abd`, `intensity`, `spectral_count`.
- Run normalization: `none`, `mean`, `median`, `max`, `global`, `max_min`,
  `iqr`.
- Sample normalization: `none`, `globalmedian`, `conditionmedian`, `quantile`,
  `mediancenter`, `meancenter`.
- IRS, coverage filter, and a subset of missing-value imputation.
- Synthetic golden tests covering several `features2proteins` main paths.

The main gaps in the current Rust implementation are:

- Other top-level non-agent commands are still placeholders: `peptides2protein`,
  `correct-batches`, `tsne_visualization`, `tissuemap`; `features2peptides` only
  implements `--generate-filter-config`, while its main flow is still a
  placeholder.
- The `features2proteins` parameter surface is larger than the implemented
  surface; most unimplemented capabilities return a stable `NotImplemented`,
  while a few parameters are currently only parsed/validated, or are skipped
  under specific quantification methods.
- The top-level `--log-file` already supports creating the parent directory and
  writing tracing logs.
- Apart from `PXD003539 sum + none + none`, which already has Python/Rust output
  and a comparison report, real-data results have not yet undergone systematic
  comparison and acceptance.
- Statistics, batch correction, differential expression, plotting, and report
  capabilities have not yet landed in Rust.
- More cross-language shared test fixtures are still missing, as are persisted,
  automated real-data comparison reports.

## Wave 0 — Modularization (DONE)

Wave 0 split single-file crates into modules with zero behavior change:

- `mokume-quant` was split into `maxlfq.rs` / `directlfq.rs`.
- `mokume-normalization` was split into `run` / `sample` / `irs` / `coverage` /
  `math.rs`.
- A new crate `mokume-imputation` was added (per-method files under `methods/`
  plus `support.rs`).
- Shared numeric primitives were moved into `mokume-core::stats`
  (`quantile_linear` / `median` / `median_finite` / `mean_positive` /
  `mean_finite` / `finite_sd` / `median_sorted`).
- `mokume-stats` was reshaped into a differential-expression and
  batch-correction crate. It now holds a `de/` module (the limma moderated
  t-test for a two-group contrast, with its special functions, Student-t
  survival, and BH correction, verified cell-for-cell against a mokume oracle)
  and a `batch/` module (parametric ComBat, verified against the inmoose
  `pycombat_norm` oracle). Both are library functions; wiring them into the
  `correct-batches` command and a `--de-output` path is still pending.
- CI (`.github/workflows/rust.yml`) was updated to add `--jobs 24` and a
  Chinese-character scan over `crates/**/*.rs`.
- `rust-toolchain.toml` was pinned to 1.96.0 (matching CI).
- A new `scripts/run_parity_matrix.sh` runs the 4-dataset by 6-method-combo
  Python-vs-Rust parity matrix.

Gate: 53 tests pass, clippy clean, fmt clean, zero behavior change.

Deliberate deviation: the full split of the 3189-line
`mokume-pipeline/src/lib.rs` was intentionally deferred to Wave 5. Its
validation gates and orchestration share private helpers, so an aggressive split
is high-risk and low immediate benefit, and no upcoming wave is blocked by it.
It will be done when Wave 5 wires in batch correction and differential
expression. (See also the status note in Wave 0 above and the Wave 5 section.)

## Wave 1 — Semantic-Divergence Fixes (DONE)

Every fix below was verified against the Python reference before changing
anything:

- `intensity` / TMTReporter quant: VERIFIED NOT a bug. Python's loader also
  applies a per-ion MAX collapse before quantify, so no change was made.
- `spectral_count` plus `median` / `abd` / `topn`: these were folding
  per-(peptidoform, charge) ions; they now collapse per-(peptidoform, charge)
  MAX values into the canonical peptide by SUM before the rollup, matching
  Python's `get_peptidoform_normalize_intensities` →
  `sum_peptidoform_intensities`. `sum` was already correct (addition is
  associative). This is locked by a golden test.
- `quantile` sample normalization: previously this assigned average RANK
  NUMBERS. The Python reference
  (`mokume/normalization/protein.py:quantile_normalize`) was ITSELF broken (it
  returned ranks and raised KeyError on within-column ties). Both Rust and
  Python were rewritten to standard quantile normalization (rank-fraction
  interpolation into a cross-sample mean reference distribution, limma
  `normalizeQuantiles` style) and verified cell-for-cell identical on a shared
  oracle. `directlfq` has since been rewritten to port the Mann-Labs `directlfq`
  package algorithm (global agglomerative sample normalization plus per-protein
  ion alignment); on PXD003539 it matches the Python output cell-for-cell at the
  median with per-sample Spearman 0.9999, the residual being floating-point
  tie-breaking in the agglomerative clustering. The earlier self-developed
  approximation has been removed.

## Agentic scope

Agentic recommendation lives in `plugins/mokume/` plus
`rust/python/mokume/agentic/`: the plugin host supplies reasoning, the MCP layer
binds policy and evidence, and the Rust matrix APIs perform normalization,
imputation, and DEA. `mokume mcp serve` is an integration entry point rather
than a fifth compute command.

## Python Capability Surface

The non-agent capabilities of the Python implementation include:

- Top-level CLI: `features2proteins`, `features2peptides`, `peptides2protein`,
  `correct-batches`, `tissuemap`, `tsne-visualization`.
- QPX, SDRF, FASTA, CSV, TSV, parquet, and AnnData reading/writing.
- feature to peptide: filtering, run/sample aggregation, IRS, log2, CSV/parquet
  output.
- peptide to protein: piBAQ, Top3, TopN, MaxLFQ, DirectLFQ, Sum, Median,
  Ratio, TMT abundance, TMT reporter intensity, Spectral count.
- piBAQ details: FASTA digestion, enzyme, min/max peptide length, TPA,
  ProteomicRuler, organism, CPC, ploidy, family YAML, shared-peptide family
  rollup, and anchor rules.
- Run normalization: `none`, `mean`, `median`, `max`, `global`, `max_min`,
  `iqr`.
- Sample normalization: `none`, `globalmedian`, `conditionmedian`,
  `hierarchical`, `quantile`, `mediancenter`, `meancenter`, `rlr`,
  `loess`. Rust currently covers `none`, `globalmedian`,
  `conditionmedian`, `quantile`, `mediancenter`, and `meancenter`.
- IRS/TMT: reference-channel auto detection, SDRF column/value selection, regex,
  `median`/`mean`, `global`/`by_mixture`/`two_stage` scope, reference-sample
  removal, ratio PS protocol.
- Missing-value imputation: `none`, `knn`, `minprob`, `mindet`, `qrilc`,
  `missforest`, `seqknn`, `impseq`, `gms`, `bpca`, `impseqrob`.
- Differential expression: `limrots`, `limma`, `deqms`, `proda`, `rots`,
  `ensemble`, plus `BH`, `IHW`, `BKY`, and `Storey` FDR and both data-driven
  effect-size gates.
- Batch correction: standalone `correct-batches` and
  `features2proteins --batch-correction`.
- Output: protein CSV, peptide CSV/parquet, ion CSV, DE CSV, plot PNG/PDF,
  interactive HTML, QC report, AnnData h5ad, TissueMap CSV/h5ad/PNG/PDF.
- TissueMap: dataset scan, TMT labeling, YAML config, filtering, PCA, t-SNE/UMAP,
  tissue-specificity scoring, marker output, and atlas plots.

## Rust Implemented Surface

The current Rust implementation is concentrated in `features2proteins`:

- The CLI already exposes the Python non-agent command surface, but only
  `features2proteins` runs a main flow.
- The QPX reader supports several candidate column names, including
  `sequence`/`Sequence`, `peptidoform`/`modified_sequence`,
  `charge`/`precursor_charge`, `run_file_name`/`reference_file_name`,
  `intensities`/`primary_intensities`, and
  `pg_accessions`/`protein_accessions`/`proteins`.
- The SDRF parser builds `run` and `run + label` indexes.
- `features2proteins` supports streaming QPX reading, filtering, aggregation,
  IRS, coverage filter, partial imputation, and CSV output.
- piBAQ queries the installed pyOpenMS `ProteaseDB` at runtime, digests every
  FASTA protein with `ProteaseDigestion`, and passes the complete theoretical
  peptide map into Rust for family discovery, shared-peptide assignment,
  denominators, and expression-matrix output.
- DirectLFQ ports the Mann-Labs `directlfq` algorithm and matches the Python
  output to a median relative error of 0 with per-sample Spearman 0.9999 on
  PXD003539 (the residual is agglomerative-clustering tie-breaking, a tolerance
  tier). MaxLFQ still has a Rust-native implementation not yet proven to
  reproduce every default detail of the Python package.
- The current performance design centers on parquet batch reading, integer IDs,
  HashMap aggregation, Rayon protein-level parallelism, and a background
  reader thread that overlaps Arrow decode/flatten with the serial aggregation
  consumer on the main `features2proteins` path.

## Key Differences

### CLI Differences

- Rust keeps agentic recommendation out of the clap compute surface; the Python
  wheel dispatcher exposes only the plugin-internal `mcp serve` entry point.
- Rust exposes the non-agent top-level commands;
  `features2peptides --generate-filter-config` is implemented, while the
  `features2peptides` main flow and the other non-`features2proteins` top-level
  commands still return `NotImplemented`.
- Rust `--log-file` now supports writing to a file and no longer returns
  `NotImplemented`.
- Many `features2proteins` parameters in Rust are already parsed, but their
  corresponding functionality is not yet implemented.

### `features2proteins` Differences

- Python supports fully exporting peptide/ion intermediate results; Rust already
  supports non-DirectLFQ `--export-peptides`, and supports a Rust-native ion
  trace wide-table export for DirectLFQ. The DirectLFQ ion table has not yet
  been proven to fully reproduce the normalized output of the Python `directlfq`
  package.
- Python supports `--normalization-proteins`; Rust now validates and accepts
  this file, and uses the protein list to restrict the normalization-factor
  computation inputs in `globalmedian`, `conditionmedian`, `mediancenter`, and
  `meancenter`. `quantile` does not use this list, and the final output matrix
  is not filtered by this list either.
- Python supports `max_min` run normalization; Rust implements a working
  within-run range-alignment version, covered by a synthetic golden test.
- Python supports more sample-normalization methods; Rust supports `none`,
  `globalmedian`, `conditionmedian`, `quantile`, `mediancenter`, and
  `meancenter`. DirectLFQ and ratio do not currently take the ordinary
  run/sample normalization path, so under these methods some normalization
  parameters are skipped rather than returning `NotImplemented` via the ordinary
  path.
- Python supports a full missing-value imputation catalog; Rust supports only a
  basic subset.
- Python supports batch correction, DE, plotting, and an interactive report;
  Rust does not yet support these.
- Python supports more DirectLFQ parameters; Rust has wired in a native-solver
  subset of `directlfq-min-nonan` and `directlfq-num-samples-quadratic`, and
  `directlfq-cores` is wired in as the fallback thread-count parameter for the
  Rayon thread pool.

### Data and Algorithm Differences

- Python's QPX compatibility tests cover both new and legacy formats; Rust
  already has legacy flatten unit tests, and has added parser assertions for the
  new-format struct `pg_accessions`, `anchor_protein`, and label-only intensity.
- Python's SQL-first tests cover filtering, charge aggregation, piBAQ shared
  rows, and normalization consistency with the legacy pandas path; Rust
  currently lacks a same-source shared fixture.
- Python's piBAQ tests cover denominator, isoform, unknown protein, TPA,
  ProteomicRuler, and family allocation; Rust covers only part of these.
- Python's benchmark has rich trend metrics; Rust still lacks stable real-data
  comparison and performance reporting.

## Refactor Principles

- First complete the verifiable non-agent main flows, then extend to advanced
  analysis and reporting.
- For each Rust capability added, first find Python behavioral evidence: source
  code, tests, documentation, or a fixed fixture.
- Prefer porting the small synthetic data from Python tests into Rust golden
  tests.
- Write real-data comparison results to `/tmp/mokume-parity/<PXD>/`, not into
  the repository.
- When a parameter is exposed but not implemented, prefer a clear, stable
  `NotImplemented` over silently ignoring it; for parameters that are only
  parsed/validated or are skipped under specific quantification methods, state
  this explicitly in the documentation.
- Multi-threaded commands explicitly use 24 threads.
- Design preference: retain Rust streaming reading, integer IDs, batch
  processing, and parallel computation, rather than copying the Python
  pandas/DuckDB shape.

## Phase Plan

### Phase 0 — Divergence Baseline and Acceptance Fixtures

Goal: give every subsequent change a comparable baseline.

Tasks:

- Add a cross-language shared fixture directory, using the small QPX/SDRF/FASTA
  data from the Python tests to generate fixed inputs and expected outputs.
  [Partially done: `crates/mokume-golden-tests/fixtures/python/` and
  `features2proteins_sum_python_compatible.csv` exist, but more Python expected
  outputs are still missing.]
- Port the new QPX and legacy QPX structures from Python
  `test_qpx_format_compat.py` into Rust parser tests. [Done:
  `crates/mokume-io/src/qpx.rs` covers legacy `sample_accession`, new-format
  struct `pg_accessions`, `anchor_protein`, and label-only intensities.]
- Port the Python SQL-first synthetic data into Rust golden tests, covering
  filtering, charge aggregation, shared peptide, and basic normalization.
  [Partially done: golden tests cover filtering, charge aggregation,
  canonical-peptide counting, shared-peptide/piBAQ family allocation,
  `globalmedian`, `conditionmedian`, run `median`, and run `max_min`; a
  same-source SQL-first shared fixture is still missing.]
- Expand the usage notes for `scripts/compare_protein_matrices.py` so each
  real-data report includes protein/sample counts, set differences, shared
  cells, absolute error, relative error, per-sample Spearman, and total-intensity
  difference. [Done: the script notes, command examples, `error_metrics`
  grouping, per-sample Spearman, and total-intensity difference have been added.]
- Add negative tests for unimplemented commands and unimplemented parameters,
  confirming they stably return `NotImplemented`. [Done: covers the non-agent
  top-level placeholder commands and a batch of unimplemented `features2proteins`
  parameters; the help option surface of the non-`features2proteins` placeholder
  commands is also covered by tests, but their main flow is still unimplemented;
  the unimplemented-parameter negative tests are not exhaustive.]

Note: the full split of the 3189-line `mokume-pipeline/src/lib.rs` is
intentionally deferred to Wave 5 (see the Wave 0 deliberate-deviation note and
the Phase 5 section). This is a planned deviation, not an oversight.

Acceptance:

```bash
cd /home/shenyufei/Git-repository/Bigbio/mokume
cargo fmt --all --check
cargo test --workspace --all-targets --jobs 24
cargo clippy --workspace --all-targets --jobs 24 -- -D warnings
```

Status: [Historical record: on 2026-06-04 the workspace passed the three Rust
routine checks above; this documentation pass did not re-run these commands, and
the current workspace state should still be reconfirmed by re-running them as
needed.]

### Phase 1 — `features2proteins` Main-Path Acceptance

Goal: first make the most common feature-to-protein path acceptable on real
data.

Tasks:

- Build a fixed comparison matrix for `PXD003539`, `PXD004701`, `PXD041421`, and
  `PXD030304`. [Partially done: `docs/features2proteins-parity.md` records the
  dataset roles, QPX paths, and SDRF paths; `PXD003539 sum + none + none` has
  generated Python/Rust output and `sum_none_report.json` in
  `/tmp/mokume-parity/PXD003539/`, while the remaining real outputs and JSON
  reports have not yet been generated.]
- Compare the following combinations first:
  - `sum + none + none` [Partially done: `PXD003539` has a report; all of its
    6500 proteins, 120 samples, and 410217 shared cells are cell-for-cell
    identical, with absolute error and relative error both zero; the other PXDs
    are not yet done.]
  - `topn + none + none`
  - `median + none + none`
  - `sum + run median + sample none`
  - `sum + run none + globalmedian`
  - `sum + run none + conditionmedian`
- Generate Python output, Rust output, and a JSON comparison report for each
  combination.
- Record the acceptance threshold for each combination. For `sum + none + none`,
  the goal should be near cell-for-cell identity; for normalization and LFQ
  methods, allow algorithmic differences but the source of the difference must
  be explained.
- Fix QPX/SDRF sample-matching, filtering, aggregation, output-ordering, and
  column-name inconsistencies. [Partially done: the Rust parser, SDRF lookup,
  filtering, stable CSV ordering, and the `ProteinName`/`protein` output columns
  are covered by tests; `PXD003539 sum + none + none` is cell-for-cell
  identical, while full real-data parity has not yet been accepted.]
- Document clearly that `--memory` currently only validates and does not limit
  memory, to avoid misleading users. [Done:
  `docs/features2proteins-parity.md` and the performance concerns in this plan
  record this behavior.]

Acceptance:

```bash
mkdir -p /tmp/mokume-parity/PXD003539/python /tmp/mokume-parity/PXD003539/rust

cd /home/shenyufei/Git-repository/Bigbio/mokume/python
python -m mokume.mokume_cli features2proteins \
  --parquet /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/qpx/PXD003539.feature.parquet \
  --sdrf /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/mokume/sdrf/PXD003539.sdrf.tsv \
  --output /tmp/mokume-parity/PXD003539/python/sum_none.csv \
  --quant-method sum \
  --run-normalization none \
  --sample-normalization none \
  --duckdb-threads 24

cd /home/shenyufei/Git-repository/Bigbio/mokume
mokume features2proteins \
  --parquet /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/qpx/PXD003539.feature.parquet \
  --sdrf /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/mokume/sdrf/PXD003539.sdrf.tsv \
  --output /tmp/mokume-parity/PXD003539/rust/sum_none.csv \
  --quant-method sum \
  --run-normalization none \
  --sample-normalization none \
  --threads 24

python scripts/compare_protein_matrices.py \
  --python /tmp/mokume-parity/PXD003539/python/sum_none.csv \
  --rust /tmp/mokume-parity/PXD003539/rust/sum_none.csv \
  --output /tmp/mokume-parity/PXD003539/sum_none_report.json
```

Note: the full real-data Phase 1 parity matrix (Wave 1-C) has now been run across all
four cell-line datasets — 24/24 deterministic combos cell-exact (see the Current Status
section at the top for details).

### Phase 2 — Complete the Non-Statistical Advanced Capabilities of `features2proteins`

Goal: make the Rust main command cover the capabilities of the Python main flow
that do not depend on plotting and differential expression.

This is Wave 2 — NOT yet done. Remaining work: missing normalization methods
(`rlr`, `loess`, `hierarchical`, peptide-level) and
imputation methods (`qrilc`, `impseq`, `missforest`).

Tasks:

- Continue refining `--export-peptides` and `--export-ions`: non-DirectLFQ
  peptide export and DirectLFQ Rust-native ion trace export have landed; the
  DirectLFQ ion table still needs verification against the normalized output of
  the Python `directlfq` package. [Partially done: synthetic golden tests cover
  the non-DirectLFQ peptide CSV and the DirectLFQ ion wide-table export.]
- Refine `--normalization-proteins` so the implemented dataset-level normalizers
  actually use this protein list. [Subset done: `globalmedian`,
  `conditionmedian`, `mediancenter`, and `meancenter` use the list to compute
  factors, and an empty file or no matching feature errors out clearly;
  `quantile` does not use the list, and the full hierarchical and other
  unimplemented normalizers remain to be done.]
- Implement `max_min` run normalization. [Done: a within-run range-alignment
  version is implemented, covered by a synthetic golden test.]
- Implement sample normalization: `quantile`, `rlr`, `loess`,
  `hierarchical`. [Partially done: `quantile` is wired into the
  `features2proteins` main flow and covered by a synthetic golden test;
  `quantile + --export-peptides` is still unimplemented. `rlr`,
  `loess`, and `hierarchical` remain to be done; the related basic items
  `mediancenter` and `meancenter` are implemented and covered by unit tests.]
- Extend imputation: `qrilc`, `impseq`. For
  `missforest`, decide between a Rust-native algorithm and an external
  crate after a small benchmark and stability assessment. [Partially done: basic
  `mindet`, `minprob`, `mean`, `median`, `constant`, `zero`, `knn`, and
  `seqknn` are implemented; of these, `mindet`, `minprob`, `mean`, `median`,
  `constant`, `zero`, and `knn` have synthetic golden tests, while
  `seqknn` still lacks an explicit golden oracle.]
- Extend DirectLFQ parameters: `directlfq-num-samples-quadratic`. [Done: this
  parameter is wired into the Rust-native DirectLFQ solver subset and has
  config-parsing/acceptance tests; `directlfq-cores` is wired in as the fallback
  thread-count parameter for the Rayon thread pool.]
- Extend piBAQ: non-trypsin enzyme, `pibaq-high-anchor-threshold`, TPA,
  ProteomicRuler, organism, CPC, ploidy, and unknown-protein behavior.
  [Done: runtime FASTA digestion for the complete installed pyOpenMS protease
  catalog; family discovery, family YAML, `pibaq-min-shared`,
  `pibaq-min-anchors`, TPA, ProteomicRuler, organism, CPC, ploidy, and
  shared-peptide assignment are implemented/tested.]

Acceptance:

- The small data from the corresponding Python tests has same-source assertions
  in Rust golden tests.
- Real data covers at least `PXD003539` and one larger dataset.
- Implemented parameters no longer return `NotImplemented`, and still-unimplemented
  parameters keep their negative tests.

### Phase 3 — Implement Standalone `features2peptides`

Goal: cover the Python feature-to-peptide command.

This is Wave 3 — NOT yet done.

Tasks:

- Split a peptide-level intermediate structure out of the Rust pipeline, to avoid
  `features2proteins` and `features2peptides` re-implementing filtering and
  normalization.
- Implement `--keep-shared-peptides`, `--remove_ids`,
  `--remove_decoy_contaminants`, and `--remove_low_frequency_peptides`.
- Implement `--skip_normalization`, run/sample normalization, and `--log2`.
- Implement CSV and parquet output.
- Implement `--filter-config` and `--generate-filter-config`, supporting
  YAML/JSON. [Partially done: `features2peptides --generate-filter-config`
  supports YAML/JSON generation; `--filter-config` and the `features2peptides`
  main flow are still unimplemented.]
- Implement the filter-override parameters: intensity, CV, charge, missed
  cleavages, modifications, unique peptides, features, and missing rate.
- Implement the IRS parameters and `aggregation_level` for features2peptides.

Acceptance:

- Port the core assertions of Python `test_features2peptides_cli.py`,
  `test_peptide_normalize.py`, and `test_filter_configs.py`.
- For the same synthetic QPX/SDRF, the Python and Rust peptide outputs have
  identical or explainably-different row counts, samples, peptides, and values.

### Phase 4 — Implement Standalone `peptides2protein`

Goal: cover the Python peptide-to-protein command and reuse the quantification
algorithms from Phases 1/2.

This is Wave 4 — NOT yet done. Together with the differential-expression
methods (`limma`, `rots`, `deqms`, `limrots`, `proda`, `ensemble`,
plus BH/IHW/BKY/Storey and the data-driven effect-size gates) tracked in Phase 6.

Tasks:

- Support peptide CSV/TSV/parquet input.
- Implement `pibaq`, `top3`, `topn`, `maxlfq`, `sum`, and `directlfq`.
- Output parquet for a `.parquet` suffix, otherwise TSV.
- Complete the piBAQ FASTA, enzyme, TPA, ProteomicRuler, organism, CPC, ploidy,
  family YAML, and anchor parameters.
- Support `--normalize`, `--verbose`, and `--qc_report`. If QC report is not
  implemented first, it must keep an explicit `NotImplemented` and documentation
  note.

Acceptance:

- Port the hand-computable assertions of Python `test_pibaq.py`,
  `test_ibaq_denominator_uniqueness.py`, `test_accession_normalization.py`, and
  `test_quantification_tmt.py`.
- A small peptide table produces identical output in Python and Rust.

### Phase 5 — Batch Correction and Post-Processing

Goal: cover the standalone `correct-batches` and the batch correction in the
main flow.

This is Wave 5 — NOT yet done. Wave 5 covers ComBat batch correction,
`correct-batches`, and pipeline post-processing wiring, and it is also where the
deferred `mokume-pipeline/src/lib.rs` module split is done (see Wave 0 and
Phase 0).

Tasks:

- Implement matrix reshape: long/wide conversion, low-protein-count sample
  filtering, missing-value filtering, and expression metrics.
- Implement batch detection: `sample_prefix`, `run`, and `column`.
- Implement covariate extraction and validation.
- Evaluate a ComBat Rust implementation approach. If there is no stable crate,
  first implement basic parametric ComBat, then add nonparametric, mean-only,
  and reference-batch variants.
- Implement standalone `correct-batches`.
- Wire the same post-processing module into
  `features2proteins --batch-correction`.
- Settle the AnnData/h5ad output approach. If Rust-native h5ad maturity is
  insufficient, first output CSV and keep h5ad as `NotImplemented` until it is
  verified usable.
- Perform the deferred split of `mokume-pipeline/src/lib.rs` while wiring in
  batch correction and differential expression.

Acceptance:

- Port the small assertions of Python `test_batch_correction.py` and
  `test_batch_correction_integration.py`.
- Do manual trend validation on `batch-quartet-multilab`, observing batch
  variance, PCA, and correlation metrics.

### Phase 6 — Differential Expression, Plotting, and Reporting

Goal: cover the statistical analysis and user-visible output of the main flow.

Tasks:

- Implement contrast parsing: CLI string and TSV file.
- Implement `limma`, `rots`, `deqms`, `proda`, `limrots`, and
  `ensemble`.
- Implement `auto` method selection: `directlfq` defaults to `deqms`, others
  default to `limrots`.
- Implement BH, IHW, BKY, and Storey FDR, including the adaptive-pi0 reliability
  fallback used by Python.
- Implement fixed and data-driven (`mixture`, `null_quantile`) effect-size gates.
- Implement DE CSV output, with fields including `ProteinName`, `log2FC`,
  `pvalue`, `adj_pvalue`, and `significance`.
- Implement volcano, heatmap, and PCA output. If using a Rust plotting library,
  first fix the image size, colors, and label rules.
- Implement an interactive HTML report. If Plotly compatibility is insufficient,
  first output a static report and keep an unimplemented note for the HTML
  report.

Acceptance:

- Port the differential-expression small assertions of Python
  `test_de_and_imputation.py`.
- Do trend validation on a benchmark with known spike-ins: DE recall, false
  positives, fold-change RMSE, and FDR.

### Phase 7 — TissueMap

Goal: cover the non-agent tissue-atlas flow of the Python `tissuemap` command.

Tasks:

- Implement dataset scan, PXD subdirectory discovery, and TMT-dataset
  specification and auto-detection.
- Implement YAML config reading, generation, and CLI override.
- Implement QPX dataset loading, tissue-field harmonization, contaminant
  filtering, and NaN filtering.
- Implement PCA input imputation, PCA, and t-SNE/UMAP.
- Implement the AdaTiSS tissue-specificity score and GMM threshold.
- Implement marker output, atlas plot, PCA scree, marker heatmap, dotplot, and
  TS distribution.
- Implement h5ad or an alternative output strategy, documented accordingly.

Acceptance:

- Use a small Python TissueMap config to generate reproducible input/output.
- Run the Rust flow on at least one local tissue dataset, with the output file
  list aligned to Python.

### Phase 8 — Performance Acceptance and Release Preparation

Goal: prove that the Rust refactor genuinely delivers better system performance,
and prepare to replace the Python non-agent CLI.

Tasks:

- Add a Rust benchmark harness recording runtime, peak memory, input row count,
  and output matrix size.
- Compare Python and Rust on `PXD003539`, `PXD004701`, `PXD041421`, and
  `PXD030304`.
- Evaluate HashMap memory usage on large datasets, and if needed introduce
  chunked aggregation, external-sort intermediate results, or sparse-matrix
  write-out.
- Evaluate the cost of the normalization double-scan over parquet, and if needed
  implement a single-scan approximation or an intermediate-statistics cache.
- Add a performance report template recording the command, dataset, thread count,
  memory, runtime, and comparison metrics.
- Confirm that all implemented/unimplemented statuses are consistent across the
  README, the parity docs, and the CLI help. [Currently known: the README still
  needs to sync the `quantile` sample normalization, and the implemented status
  of `features2peptides --generate-filter-config`.]

Acceptance:

- Routine validation passes:

```bash
cd /home/shenyufei/Git-repository/Bigbio/mokume
cargo fmt --all --check
cargo test --workspace --all-targets --jobs 24
cargo clippy --workspace --all-targets --jobs 24 -- -D warnings
```

- The real-data parity reports are complete, and differences are either accepted
  or clearly explained.
- Rust has a stable speed or memory advantage on the main flows.

## Test Migration Checklist

Prioritize porting the non-agent parts of these Python tests:

- QPX compatibility: `tests/test_qpx_format_compat.py` [Core parser migration
  done: the legacy and new struct QPX formats have Rust unit tests.]
- SQL-first loading: `tests/test_loading_sqlfirst.py` [Partially done: filtering,
  charge aggregation, shared peptide, and basic normalization are covered in Rust
  golden tests; a same-source SQL-first fixture is still missing.]
- DirectLFQ loading: `tests/test_loading_directlfq.py` [Partially done: the
  DirectLFQ synthetic solver and ion export have Rust golden tests; a comparison
  against the Python `directlfq` package output is still missing.]
- CLI parameters: `tests/test_features2proteins_cli.py`,
  `tests/test_features2peptides_cli.py` [Partially done: the wheel command help and
  `features2proteins` parameter parsing have tests;
  `features2peptides --generate-filter-config` has YAML/JSON generation and
  missing-parameter error tests; the help option surface of `peptides2protein`,
  `correct-batches`, `tsne_visualization`, and `tissuemap` is covered by tests,
  but the main flow of these commands is still `NotImplemented`.]
- piBAQ: `tests/test_pibaq.py`, `tests/test_ibaq_denominator_uniqueness.py`,
  `tests/test_pibaq_tpa_numerator_proteotypic.py` [Partially done: theoretical
  peptide denominators, family shared-peptide allocation, TPA, and basic piBAQ
  synthetic oracles are covered; the remaining Python edge cases are not all
  ported as Rust golden tests.]
- accession: `tests/test_accession_normalization.py` [Partially done: QPX protein
  accession parsing and group output are covered; FASTA accession normalization is
  still missing.]
- TMT/ratio: `tests/test_quantification_tmt.py` [Partially done: the Ratio
  synthetic oracle is covered; IRS SDRF column/value detection is implemented with
  parameter validation, reference-sample removal has golden-test coverage, but SDRF
  column/value detection still lacks explicit golden coverage; full TMT behavior is
  still missing.]
- Normalization: `tests/test_hierarchical_normalization.py`,
  `tests/test_normalization_rlr.py`,
  `tests/test_normalization_distribution.py` [Partially done: run `median`, run
  `max_min`, `globalmedian`, `conditionmedian`, and `quantile` have golden-test
  coverage; `mediancenter` and `meancenter` are implemented and covered by
  `mokume-normalization` unit tests, but golden-test coverage is not claimed for
  them at present; hierarchical, RLR, and others are still missing.]
- Missing-value imputation: `tests/test_imputation_python.py`,
  `tests/test_imputation_statistical.py` [Partially done: `mindet`, `minprob`,
  `mean`, `median`, `constant`, `zero`, and `knn` are covered; `seqknn` is
  implemented but still lacks an explicit golden oracle; `qrilc` and the other
  advanced methods are still missing.]
- Differential expression: `tests/test_de_and_imputation.py`
- Batch correction: `tests/test_batch_correction.py`,
  `tests/test_batch_correction_integration.py`
- Filter config: `tests/test_filter_configs.py`

Plugin policy, evaluation, and MCP registration are tested separately under
`rust/tests/test_agentic_plugin.py`.

## Real-Data Comparison Matrix

The first round uses the data already in the Rust parity docs:

- `PXD003539`: first-pass baseline.
- `PXD004701`: normalization and TopN baseline.
- `PXD041421`: normalization and TopN baseline.
- `PXD030304`: stress baseline.

Current real-data results:

- `PXD003539`'s `sum + none + none` has generated Python/Rust output and
  `/tmp/mokume-parity/PXD003539/sum_none_report.json`; all 6500 proteins, 120
  samples, and 410217 shared cells are cell-for-cell identical, with absolute
  error and relative error both zero.

For each dataset, first run:

- `--quant-method sum --run-normalization none --sample-normalization none`
- `--quant-method topn --topn 3 --run-normalization none --sample-normalization none`
- `--quant-method median --run-normalization none --sample-normalization none`
- `--quant-method sum --run-normalization median --sample-normalization none`
- `--quant-method sum --run-normalization none --sample-normalization globalmedian`
- `--quant-method sum --run-normalization none --sample-normalization conditionmedian`

Then extend to:

- piBAQ + FASTA.
- Ratio + IRS.
- MaxLFQ.
- DirectLFQ.
- coverage filter.
- imputation.
- batch correction.
- DE.

## Performance Concerns

- `--threads` sizes the Rayon thread pool. The `features2proteins` main path now
  overlaps the parallelizable work (Arrow batch decode and flattening, run on a
  background reader thread with Rayon-parallel per-window flattening) with the
  serial consumer (string-id registration and intensity accumulation). Batches
  reach the consumer in reader order and rows keep their order, so the protein
  matrix stays cell-for-cell identical to a serial pass. The serial consumer is
  the remaining wall-clock floor (~40% of CPU work); parallelizing it would
  require sharded accumulation and is deferred.
- Measured on 2026-08-23 with PXD003539 / PXD004701 (`sum`, 24 threads, one
  warm-up, median of three measured runs): Rust/Python wall times are
  8.95/7.17 seconds and 17.84/13.99 seconds, while peak memory is 0.86/4.25 GiB
  and 1.29/8.02 GiB. Protein sets, sample sets, and all 390,540/544,008 matrix
  cells are exact. The Rust path is currently about 1.25-1.28x slower and
  4.9-6.2x lower-memory for this workload.
- Parquet reading is streaming by batch, but many aggregations are still held in
  HashMap/HashSet.
- Enabling normalization first scans parquet once, then runs the main
  aggregation, increasing IO cost.
- FASTA is currently read fully into memory; a large FASTA carries a fixed memory
  pressure.
- `--memory`/`--duckdb-memory` currently only do string validation, not
  process-memory limiting.
- Output ordering guarantees stability, but sorting a large matrix has a cost.
- A real performance report must record both runtime and peak memory, not just
  wall time.

## Risks

- The parameter surface is larger than the implemented surface; users may
  mistakenly assume some features are already usable.
- Synthetic tests cannot prove all variants of real QPX/SDRF.
- Some Python algorithms depend on libraries such as pandas, DuckDB, scipy,
  statsmodels, sklearn, directlfq, combat, and scanpy; Rust needs to confirm,
  item by item, whether to implement natively, find a crate, or keep
  unimplemented for now.
- The algorithm details of DirectLFQ and MaxLFQ may not be cell-for-cell
  reproducible; an acceptable threshold needs to be determined.
- piBAQ family allocation and shared-peptide rules are prone to small numerical
  differences.
- h5ad, the interactive report, and plotting output may require new dependency
  evaluation in Rust.
- The Rust workspace is on a local `dev` branch and not pushed; before pushing,
  the Bigbio repository's review and validation requirements still need to be
  satisfied.

## Completion Criteria

When this goal is complete, the following must hold:

- The Rust non-agent CLI covers the Python non-agent user-visible features, with
  unimplemented exceptions clearly recorded and confirmed acceptable to keep.
- Agentic recommendation is not a `mokume` compute command; it is a plugin over
  the local MCP integration surface.
- The Python/Rust shared tests, Rust golden tests, and real-data parity reports
  can all prove that the main-path behavior is consistent or the difference is
  acceptable.
- On the main datasets, Rust's runtime or memory usage is better than Python's,
  or at least more stable in large-data scenarios.
- The implemented statuses across the README, parity docs, CLI help, tests, and
  this plan are consistent.
