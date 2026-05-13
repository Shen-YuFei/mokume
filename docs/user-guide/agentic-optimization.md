# Agentic Optimization

`mokume agentic optimize` runs an **LLM-assisted (or rule-based) search**
over differential-expression workflows: normalization × imputation ×
DE method × FDR × log2FC threshold, including **ensemble (top-k
consensus) DE** strategies.

It takes an already-quantified protein matrix and an SDRF, profiles the
data, proposes a handful of candidate configurations, runs each through
the same engine as `features2proteins` does, scores them, reflects, and
iterates until convergence or the experiment budget is exhausted.

## When to Use

- You have a quantified protein matrix and want a defensible
  "best-effort" DE workflow without manually grid-searching.
- You want to **hedge against single-method bias** via ensemble DE.
- You have a benchmark / ground truth and want to know which workflow
  maximises sensitivity at controlled FDR.

If you just want one fixed pipeline, prefer
[`features2proteins`](features2proteins.md).

## Workflow

```text
┌────────────────────┐
│  Data Profile      │  ← n_samples, n_proteins, missingness,
│  (profiler.py)     │     peptide_counts, data_type (LFQ/TMT/DIA)
└────────┬───────────┘
         │
         ▼
┌────────────────────┐
│  Proposer          │  ← LLM (tool-call, strict mode) OR
│  - LLM             │    rule-based fallback (heuristics.yaml +
│  - Rule-based      │    ensemble_strategies)
└────────┬───────────┘
         │  3–18 CandidateConfigs
         ▼
┌────────────────────┐
│  Runner            │  norm → imputation (log2 space) →
│  (runner.py)       │  DifferentialExpression OR run_ensemble
└────────┬───────────┘
         │  one DE table per config
         ▼
┌────────────────────┐
│  Evaluator         │  ground-truth metrics (AUC/TP/FP) OR
│  (evaluator.py)    │  unsupervised score (DE count, CV, missing)
└────────┬───────────┘
         │  scored results
         ▼
┌────────────────────┐
│  Reflector (LLM)   │  → convergence? next_configs? OR
│                    │    rule-based: stop after N rounds
└────────┬───────────┘
         │
         ▼   (loop until converged or budget hit)
   Best config + audit trail
```

## CLI Reference

```bash
mokume agentic optimize \
  --protein-matrix proteins.tsv \
  --sdrf sdrf.tsv \
  --contrasts "treated vs control" \
  [--ground-truth gt.txt] \
  [--expected-fc expected.yaml] \
  [--max-rounds 5] \
  [--max-experiments 30] \
  [--llm-provider deepseek|custom] \
  [--llm-base-url URL] \
  [--llm-api-key KEY] \
  [--llm-model MODEL] \
  [--no-llm] \
  [--output-dir ./optimization]
```

| Option | Description |
|---|---|
|`--protein-matrix`| Protein intensity matrix (TSV/CSV, first column = protein ID). |
|`--sdrf`| SDRF file used to map samples to conditions. |
|`--contrasts`| Comma-separated contrasts, e.g. `"A vs B, A vs C"`. |
|`--ground-truth`| Optional file with true-positive protein IDs (one per line). Enables Mode A (AUC + TP/FP) scoring. |
|`--expected-fc`| Optional YAML with expected log2FC per contrast (for orientation/calibration). |
|`--max-rounds`| Maximum number of propose → run → reflect rounds (default: 5). |
|`--max-experiments`| Hard cap on total configurations evaluated across all rounds (default: 30). |
|`--llm-provider`| `deepseek` (default) or `custom` (any OpenAI-compatible endpoint). |
|`--llm-base-url`| OpenAI-compatible API base URL. Required when `--llm-provider custom`. |
|`--llm-api-key`| API key. Reads `DEEPSEEK_API_KEY` / `OPENAI_API_KEY` from environment if omitted. |
|`--llm-model`| Model name (auto-selected from provider preset). |
|`--no-llm`| Disable the LLM and run the rule-based engine only. |
|`--output-dir`| Where to write the audit trail, per-round configs, and the final report (default: `./optimization`). |

`--llm-api-key` and `--llm-base-url` are persisted to `.env` on first
successful use so subsequent runs do not need them on the command line.

## Method Catalogue Exposed to the Agent

The proposer (LLM and rule-based) and the runner share the same
catalogue. New methods added to the runner are automatically usable from
the rule engine; for the LLM, they are also enumerated in the strict
tool schema (`mokume/agentic/llm_client.py`).

| Axis | Choices |
|---|---|
|`de_method`| `limma`, `rots`, `deqms`, `proda`, `msstats`, `ensemble` |
|`fdr_method`| `bh`, `ihw` |
|`normalization`| `none`, `median`, `quantile`, `mean`, `rlr`, `vsn`, `loess`, `mbqn` |
|`imputation`| `none`, `minprob`, `mindet`, `knn`, `missforest`, `seqknn`, `qrilc`, `mle`, `mice`, `nbavg`, `gms`, `bpca`, `impseq`, `impseqrob` |
|`ensemble`| `none`, `limma,deqms,proda`, `limma,rots,deqms`, `limma,rots,deqms,proda` |
|`ensemble_k`| Integer 1–5 (top-k consensus across the ensemble methods) |
|`log2fc_threshold`| Typically 0.5 or 1.0 |

Notes:

- All DE, normalization, and imputation methods (including `bpca`,
  `impseq`, `impseqrob`, `vsn`, `proda`, `deqms`, `msstats`) are
  **pure-Python reimplementations**. No `rpy2`, R runtime, or
  Bioconductor packages are required.
- Imputation always runs in **log2 space**, and the runner restores the
  original scale afterwards. See
  [concepts/imputation.md](../concepts/imputation.md).
- Ensemble configs use the same FDR method, normalization, and
  imputation as their non-ensemble peers; only the DE step changes.

## Ensemble DE in the Agent Loop

When the proposer picks `de_method=ensemble`, the runner calls
[`run_ensemble`](../reference/python-api.md) with the methods in
`ensemble` and the `ensemble_k` threshold. The top-k consensus is then
forwarded to the evaluator like any other DE result.

Two default ensemble presets ship in
`mokume/agentic/knowledge/heuristics.yaml`:

```yaml
ensemble_strategies:
  top_k_2of3:
    methods: ["limma", "deqms", "proda"]
    min_k: 2
  top_k_3of4:
    methods: ["limma", "rots", "deqms", "proda"]
    min_k: 3
```

Add more strategies there; the rule-based engine will emit one
`CandidateConfig` per strategy, using the data-type-recommended FDR,
normalization, and imputation as the rest of the pipeline.

## Examples

### LLM-driven (DeepSeek)

```bash
mokume agentic optimize \
  --protein-matrix proteins.tsv \
  --sdrf sdrf.tsv \
  --contrasts "treated vs control" \
  --llm-api-key sk-xxx \
  --output-dir runs/deepseek
```

### Custom OpenAI-compatible provider

```bash
mokume agentic optimize \
  --protein-matrix proteins.tsv \
  --sdrf sdrf.tsv \
  --contrasts "treated vs control" \
  --llm-provider custom \
  --llm-base-url https://my-llm.example.com/v1 \
  --llm-model my-model \
  --llm-api-key sk-xxx
```

### Rule-based only (no LLM)

Useful for reproducible CI runs or when an API key is unavailable. The
rule engine reads `heuristics.yaml` and emits a diverse set of
candidates, **including ensemble strategies**:

```bash
mokume agentic optimize \
  --protein-matrix proteins.tsv \
  --sdrf sdrf.tsv \
  --contrasts "treated vs control" \
  --no-llm \
  --output-dir runs/rule_based
```

### With ground truth (Mode A scoring)

```bash
mokume agentic optimize \
  --protein-matrix proteins.tsv \
  --sdrf sdrf.tsv \
  --contrasts "treated vs control" \
  --ground-truth true_positives.txt \
  --expected-fc expected_fc.yaml \
  --output-dir runs/benchmark
```

## Outputs

Under `--output-dir`:

- `audit_trail.json` — every propose / run / reflect step
- `round_<N>/configs.json` — candidates proposed that round
- `round_<N>/results.csv` — per-config DE counts, CV, score, etc.
- `best_config.json` — winning configuration and its DE table
- `report.md` — LLM-written summary (if LLM is enabled)

## Python API

```python
from mokume.agentic.config import AgenticConfig
from mokume.agentic.optimizer import optimize

config = AgenticConfig(
    use_llm=False,
    max_rounds=3,
    max_experiments=12,
    contrasts=[("treated", "control")],
    output_dir="./runs",
)
results = optimize(
    protein_df=protein_df,
    sample_to_condition=s2c,
    config=config,
)
for contrast, state in results.items():
    print(contrast, state.best_config.name, state.best_score)
```

## See Also

- [features2proteins user guide](features2proteins.md) — fixed-pipeline counterpart
- [differential-expression concepts](../concepts/differential-expression.md) — ensemble + MSstats details
- [imputation concepts](../concepts/imputation.md) — method choices
- [normalization concepts](../concepts/normalization.md) — MBQN / LOESS / VSN
