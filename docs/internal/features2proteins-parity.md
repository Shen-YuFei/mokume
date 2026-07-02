# features2proteins Python/Rust Parity Notes

This document tracks the first parity target for the Rust `features2proteins`
implementation. Python remains the behavioral reference, but Rust behavior must
also be checked against the algorithm intent for each parameter.

## Parameter Alignment

| Parameter | Stage | Python Behavior | Rust Behavior | Status |
| --- | --- | --- | --- | --- |
| `--parquet` | Input | Required QPX feature parquet input. | Required QPX feature parquet input. | Implemented |
| `--sdrf` | Input | Optional SDRF sample metadata. Required for ratio and metadata-aware stages. | Optional SDRF sample metadata. Required for ratio and metadata-aware stages. | Implemented subset |
| `--output` | Output | Required protein matrix output path. | Required protein matrix output path. | Implemented |
| `--quant-method` | Quantification | Selects directlfq, ibaq, maxlfq, topn, sum, median, ratio, abd, intensity, or spectral_count. | Selects the same method set for `features2proteins`. Missing protein x sample cells follow Python's per-method convention (`pivot_table(fill_value=0)` over the observed samples): the additive methods (sum / intensity / spectral_count / directlfq) write `0` for a missing cell in an observed sample, the average/ratio methods (median / topn / abd / maxlfq / ratio) leave it empty; a sample with no observations stays empty for every method. PXD003539 sum is cell-exact (0 NaN-pattern mismatches). | Implemented subset |
| `--topn` | Quantification | Sets N for TopN when `--quant-method topn`. | Sets N for TopN when `--quant-method topn`. | Implemented |
| `--min-aa` | Filtering | Removes peptide sequences shorter than the threshold. | Removes peptide sequences shorter than the threshold. | Implemented |
| `--min-unique` | Filtering | Requires a minimum number of unique peptides per protein/sample cell. | Requires a minimum number of unique peptides per protein/sample cell for non-iBAQ methods. | Implemented subset |
| `--remove-contaminants` / `--keep-contaminants` | Filtering | Controls contaminant and decoy removal. | Controls contaminant and decoy removal. | Implemented subset |
| `--run-normalization` | Normalization | Supports none, mean, median, max, global, max_min, and iqr. | Supports none, mean, median, max, global, max_min, and iqr. | Implemented subset |
| `--sample-normalization` | Normalization | Supports none, globalmedian, conditionmedian, hierarchical, quantile, mediancenter, meancenter, rlr, and loess. | Supports none, globalmedian, conditionmedian, quantile, mediancenter, meancenter, rlr, loess, and hierarchical. hierarchical has non-identity real-path golden oracles and is cell-exact; loess is ~2e-3 vs statsmodels lowess; mediancenter/meancenter/hierarchical are real-data cell-exact on PXD003539. | Implemented subset |
| `--threads` / `--duckdb-threads` | Runtime | Caps Python-side DuckDB or method-specific workers depending on the path. | Configures the Rayon global thread pool for parallel Rust sections. | Implemented subset |
| `--memory` / `--duckdb-memory` | Runtime | Caps DuckDB memory in Python. | Parses and validates the memory string; QPX reading is not DuckDB-based. | Implemented subset |
| `--export-peptides` | Output | Writes normalized peptide-level intermediates. | Writes Python-shaped peptide intermediates for non-DirectLFQ methods; DirectLFQ peptide export still returns `NotImplemented`. | Implemented subset |
| `--export-ions` | Output | Writes normalized ion-level intermediates for DirectLFQ. | Writes a Python-shaped Rust-native DirectLFQ ion trace matrix; non-DirectLFQ methods return `NotImplemented`. | Implemented subset |
| `--normalization-proteins` | Normalization | Restricts normalization proteins to the provided list for selected dataset-level normalizers without filtering the final matrix. | Implemented sample normalizers use the file to restrict factor inputs and still apply the factors to the full output matrix. Empty files or lists with no matching features fail clearly. | Implemented subset |
| `--coverage-threshold` | Postprocessing | Drops proteins below per-condition non-missing coverage. | Applies a per-condition coverage filter when SDRF condition metadata is available. | Implemented subset |
| `--impute` and basic imputation options | Postprocessing | Supports a broad imputation catalog. | Supports none, mindet, minprob, mean, median, constant, zero, most_frequent, knn, seqknn, qrilc, impseq, gms, bpca, and impseqrob. `missforest` is accepted by the CLI but returns `NotImplemented` — a documented gap: it wraps scikit-learn's `IterativeImputer` driven by `RandomForestRegressor`, whose output is the artifact of sklearn's exact tree-building internals plus its RNG. Cross-language alignment would require reimplementing the estimator bit-for-bit (the model differs structurally, not just in RNG draws, unlike the bootstrap DE methods), so no meaningful tolerance tier is reachable; a Rust ML crate (linfa) would not align either. | Implemented subset |
| batch correction options | Postprocessing | Runs ComBat when optional dependencies are available. | ComBat (parametric, covariate (covar_mod), non-parametric (par_prior=false), with ref_batch/mean_only, oracle-verified vs inmoose ~1e-6 / 1e-9) is wired into both the standalone `correct-batches` command and the `features2proteins` pipeline (`--batch-correction`). In `features2proteins`, batches come from `sample_prefix` or explicit `column` detection (`--batch-method column` + `--batch-column`, mapped via SDRF `source name -> column`, missing samples `"unknown"`), and `--batch-covariates` are extracted from the SDRF (`extract_covariates_from_sdrf`: column match, sample-substring fallback, `pd.factorize` encoding, single-value columns dropped — oracle-locked) and fed to the covariate ComBat design. ComBat runs on the proteins with no missing cells (the rest are kept uncorrected, mirroring `_complete_batch_matrix`). `--batch-method run` has no run-level mapping in the protein-matrix flow and errors at runtime, the same as Python's `_detect_batch_indices` (`run_info required`); fraction / techreplicate detection are not exposed by the `features2proteins` CLI (Python `click.Choice` is `sample_prefix`/`run`/`column`). PCA+HDBSCAN outlier removal is not ported; `--export-anndata` (h5ad) is implemented only in `correct-batches` (Rust-native `.h5ad` matching Python `anndata.write_h5ad`, verified via `anndata.read_h5ad`), which does not expose batch-method/covariate options. | Implemented subset |
| differential expression options | Postprocessing | Runs DE methods and writes optional DE output. | `features2proteins --de --de-contrasts "A vs B" --de-method <m>` runs the full catalog via a Rust dispatcher: limma, deqms, rots, limrots, proda, ensemble (with BH). Conditions resolve from the SDRF by both `source name` and `comment[data file]` stem (so run-level matrices work). Deterministic kernels are cell-exact on real data (PXD004701: limma/deqms log2FC ~5e-15, 100% significance-call agreement); RNG/optimizer-driven methods (rots/limrots/proda) are faithful-not-bit-exact (log2FC cell-exact, p-value rank-level). Both BH and IHW FDR are supported (`--de-fdr-method`); IHW (covariate-binned weighted BH) is cell-exact vs Python. Contrasts come from inline `--de-contrasts "A vs B"` and/or the TSV `--de-contrasts-file` (group1/group2 columns, appended after the inline entries, mirroring features2proteins.py:768); `--de-method auto` resolves to `deqms` for directlfq, otherwise `limrots` (Python's `_resolve_de_method`). | Implemented (tiered) |
| plotting and report options | Output | Writes plots and optional HTML reports. | Returns `NotImplemented`. | Not implemented |

## Real-Data Parity Matrix

Use the local cell-line QPX data for Rust/Python comparisons.

| Dataset | Role | Feature Parquet | SDRF |
| --- | --- | --- | --- |
| `PXD003539` | First-pass baseline | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/qpx/PXD003539.feature.parquet` | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell_lines_mokume/PXD003539/sdrf/PXD003539.sdrf.tsv` |
| `PXD004701` | Normalization and TopN baseline | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD004701/qpx/PXD004701.feature.parquet` | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell_lines_mokume/PXD004701/sdrf/PXD004701.sdrf.tsv` |
| `PXD041421` | Normalization and TopN baseline | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD041421/qpx/PXD041421.feature.parquet` | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell_lines_mokume/PXD041421/sdrf/PXD041421.sdrf.tsv` |
| `PXD030304` | Stress baseline | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD030304/qpx/PXD030304.feature.parquet` | `/home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell_lines_mokume/PXD030304/sdrf/PXD030304.sdrf.tsv` |

## Required Comparison Metrics

Every real-data parity report must include:

- protein count
- sample count
- protein set differences
- sample set differences
- shared matrix cell count
- median absolute error
- 95th percentile absolute error
- median relative error
- Spearman correlation per sample
- total intensity difference per sample
- largest protein/sample cell differences

## First-Pass Command Shape

Python baseline:

```bash
python -m mokume.mokume_cli features2proteins \
  --parquet /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/qpx/PXD003539.feature.parquet \
  --sdrf /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell_lines_mokume/PXD003539/sdrf/PXD003539.sdrf.tsv \
  --output /tmp/mokume-parity/PXD003539/python/sum_none.csv \
  --quant-method sum \
  --run-normalization none \
  --sample-normalization none \
  --duckdb-threads 24
```

Rust baseline:

```bash
cargo run -p mokume-cli -- features2proteins \
  --parquet /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell-lines/PXD003539/qpx/PXD003539.feature.parquet \
  --sdrf /home/shenyufei/Git-repository/Bigbio/Bigbio_data/cell_lines_mokume/PXD003539/sdrf/PXD003539.sdrf.tsv \
  --output /tmp/mokume-parity/PXD003539/rust/sum_none.csv \
  --quant-method sum \
  --run-normalization none \
  --sample-normalization none \
  --threads 24
```

Comparison:

```bash
python scripts/compare_protein_matrices.py \
  --python /tmp/mokume-parity/PXD003539/python/sum_none.csv \
  --rust /tmp/mokume-parity/PXD003539/rust/sum_none.csv \
  --output /tmp/mokume-parity/PXD003539/sum_none_report.json
```

## Verification Basis and Known Gaps

As of the 2026-06-20 alignment audit, 26 computation-core methods are parity-checked
and considered aligned: the quant methods reachable through `features2proteins`
(sum, median, topn/top3, intensity, abd, spectral_count, directlfq, maxlfq — which
delegates to the directlfq path by default, matching Python — and the piBAQ
family-allocation core), the run-level and sample-level normalizers listed above, and
14 imputation methods.

Real-data parity (2026-06-21, Wave 1-C): `scripts/run_parity_matrix.sh` was run across
all four cell-line datasets (PXD003539/004701/041421/030304) for the Phase-1
deterministic combos — **24/24 combos cell-exact** (rel max 0.0 for
sum/median/topn/runmedian/conditionmedian; ~5e-8 float-level for globalmedian; protein
and sample sets exact py==rust). `maxlfq` and `directlfq` were separately real-data
parity-checked on PXD003539: **protein sets exact (py==rust, 0/0)** after the maxlfq
per-(protein,sample) min-unique gate fix and the directlfq first-accession naming fix.
The DE family was real-data verified on PXD004701 (limma/deqms cell-exact, 100%
significance-call agreement; rots/limrots rank-level).

directlfq/maxlfq are now **cell-exact on PXD003539** within the f32 tolerance tier
(rel max ~1e-14, 0 cells over 2.5e-7, every per-sample Spearman >= 0.999999; maxlfq
adds a uniform ~2e-8 rescale offset, still below the gate). This closed the prior
"agglomerative-clustering tail", which was not float chaos but two deterministic
mismatches against `directlfq` 0.3.3:

- **Ion row order.** DirectLFQ sorts ion rows by `(protein, sequence)` lexically
  (`stages.py:887` `.sort(["protein", "sequence"])`; directlfq
  `sort_input_df_by_protein_and_quant_id`) and `get_normfacts`' agglomerative merge
  breaks variance-`argmin` ties by row position. The Rust path keyed rows on the raw
  `PeptideId` (registration order), feeding `get_normfacts` a different row order and
  rerouting the merge tree on ~21% of proteins. `DirectLfqIon` now carries the bare
  sequence's global alphabetical rank (`ion_seq_rank`); `IonMatrix::build` sorts on it.
- **NaN distance-to-reference in the per-protein linear shifter.**
  `SampleShifterLinear._shift_to_reference_sample` adds the distance to the whole row
  unconditionally (`normalization.py:412`); when the quadratic-subset reference shares
  no finite sample with a linear ion the distance is `NaN`, so `row += NaN` masks every
  cell. The Rust path skipped a non-finite shift, leaving the lone finite value, which
  resurrected a sample Python drops (a 0 -> value flip that rescaled the whole protein).
  `normalize_ion_profiles` now NaNs the row on a non-finite distance. Affected 19
  proteins / 27 cells on PXD003539.

Both stages were verified in isolation: fed the byte-identical post-sample-norm log2
matrices DirectLFQ hands each protein, the Rust per-protein estimator and the global
sample-shift computation match Python to machine epsilon (~1e-15) on PXD003539 and
PXD004701. On PXD004701 the directlfq protein matrix still shows a per-sample residual
(rel max ~5%); this is a **pre-existing upstream feature-intensity aggregation
difference** (per-protein totals differ identically at baseline da05206 and after this
fix, and the sample-shift / per-protein stages are machine-exact on its captured input),
not the agglomerative clustering — tracked separately. `ratio` has no real-data parity
(cell-line data is label-free / no plex); it is covered by a synthetic plexed golden
test cell-exact vs Python `RatioQuantification`.

Since the 2026-06-20 audit: the DE family (limma/deqms/rots/limrots/proda/
ensemble) plus `IHW` and covariate/non-parametric ComBat were ported and wired;
`most_frequent` imputation and the iBAQ extras (TPA/ProteomicRuler/normalize_ibaq)
landed; `peptides2protein`, `correct-batches`, and `features2peptides` (filters +
peptide normalization) are wired; `ratio` switched to an SDRF-driven `sample_to_plex`
map; `maxlfq` and `directlfq` protein-set divergences were fixed and real-data verified.
Deterministic methods are cell-exact; RNG/optimizer DE methods are
faithful-not-bit-exact (tiered tolerance).

Remaining gaps (irreducible or niche):

- Python-delegated computation methods (functionally available via Python mokume;
  the Rust kernel does not reproduce them, and its error points the user to
  Python): `missforest` (wraps scikit-learn `IterativeImputer` +
  `RandomForestRegressor`; tree + RNG internals not reproducible cross-language —
  a Rust ML crate would not align either).
- `features2peptides` runs a real main flow (default `--min_aa`/`--min_unique`/
  contaminant filtering, factor-based run+sample normalization, `--keep-shared-peptides`,
  `--remove_ids`, `--remove_low_frequency_peptides`, `--log2`, `--save_parquet`,
  `--aggregation_level run` all golden-tested vs Python). Three capabilities are
  still `NotImplemented` — channel-based IRS is **not** the only gap:
  - preprocessing filter pipeline (`--filter-config` YAML/JSON + `--filter-*`): the
    per-row filters (min-intensity floor with the `remove_zero_intensity` 1e-10 floor,
    peptide length, charge states, excluded modifications, trypsin missed cleavages)
    and the per-`(protein, sample)` unique-peptide gate are wired. They apply during
    ingest *before* the unique gate, matching Python's pipeline chain order, and each
    primitive is oracle-locked vs Python (`mokume-pipeline::filters`). The group-level
    filters (CV threshold, replicate agreement, quantile, run-QC min-features /
    max-missing-rate / total-intensity / sample-correlation) need a per-sample chain
    the streaming collect-then-export model does not express, so they fail fast with a
    per-filter `NotImplemented` stage rather than apply an approximate result; the
    rarer per-row filters (search-score, sequence patterns, coverage, razor handling,
    custom contaminant patterns) are rejected the same way. Peptide/protein FDR
    thresholds are a no-op without a q-value column (matching Python's apply-time check).
    A known edge: a custom `--filter-min-intensity > 0` combined with normalization may
    shift the normalization median map (Python's SQL applies it; the Rust factor pass
    uses `intensity > 0`); the per-row filtering itself is exact, verified with
    `--skip_normalization`.
  - dataset-level sample normalization inside `features2peptides`
    (`quantile`/`rlr`/`loess`/`hierarchical`/`mediancenter`/`meancenter`)
    is rejected; only factor-based methods run. **Intentional divergence (not a port gap):**
    these are unbuilt no-op placeholders in Python's peptide flow — each registered
    `*_normalization(dataset_df, sample, med_map)` just `return dataset_df`
    (`model/normalization.py:404-455`) and the per-sample streaming loop has no
    post-loop dataset pass (`normalization/peptide.py:370,412-415`), so Python silently
    applies nothing. The streaming model only fits scalar-per-sample methods
    (global/condition median); cross-sample distribution alignment requires the full
    matrix, which Python defers to the protein stage. Rust returns `NotImplemented`
    rather than reproduce a silent no-op, so a requested normalization never looks
    applied when it is not. (These methods *are* implemented and oracle-verified in
    `features2proteins`, where the full matrix exists.)
  - channel-based IRS (`--irs_channel`/`--irs_autodetect_regex`) needs the
    `mixture`/`channel` columns (unparsed), is TMT-only, and is unverifiable on the
    available label-free datasets. SDRF-driven multi-plex IRS *is* implemented in
    `features2proteins`; only the `features2peptides` channel path is missing.
