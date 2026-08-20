# MCP input contract

Call `mokume.inspect_dataset` with the two absolute input paths and one optional
metadata object:

```json
{
  "protein_matrix": "/absolute/proteins.tsv",
  "sdrf": "/absolute/project.sdrf.tsv",
  "metadata": {
    "data_type": "LFQ",
    "quantification": "directlfq",
    "upstream_engine": "quantms",
    "factor_column": null
  }
}
```

The metadata fields are optional; do not invent values that were not supplied or
supported by the inputs.

Call `mokume.evaluate_recommendation` with a two-item contrast, the generated
recommendation block, and an execution-options object:

```json
{
  "protein_matrix": "/absolute/proteins.tsv",
  "sdrf": "/absolute/project.sdrf.tsv",
  "contrast": ["control", "treated"],
  "recommendation": {},
  "options": {
    "output_dir": "/absolute/results",
    "ground_truth": null,
    "expected_direction": null,
    "data_type": "LFQ",
    "quantification": "directlfq",
    "upstream_engine": "quantms",
    "factor_column": null,
    "fdr_threshold": 0.05,
    "input_scale": "auto",
    "threads": 24
  }
}
```

The options object accepts only the fields shown above. Build `contrast` from two
keys in the inspection result's `profile.samples_per_condition`; do not replace
those canonical labels with longer raw SDRF values. `ground_truth`, when present,
must be an absolute path to a one-protein-per-line file and requires
`expected_direction` to be `UP` or `DOWN`. When `ground_truth` is null,
`expected_direction` must also be null and carries no biological meaning.
`options.output_dir` must be an absolute path that does not already exist. Repeat
the same declared `data_type`, `quantification`, `upstream_engine`, and
`factor_column` values used for `inspect_dataset` so evaluation rebinds the same
policy context.

For a Score A result, calculate each candidate's tested-universe size as
`TP + FP + FN + TN` and report it with the ranking. Unequal totals mean that the
candidates were scored over different measurable universes; preserve the returned
ranking, but do not describe it as general superiority without this limitation.

## Recommendation block

Pass `mokume.evaluate_recommendation` an object with exactly this shape:

```json
{
  "configs": [
    {
      "name": "candidate_name",
      "de_method": "limma",
      "fdr_method": "bh",
      "normalization": "none",
      "imputation": "none",
      "ensemble": "none",
      "ensemble_k": 2,
      "log2fc_threshold": 0.5,
      "reasoning": "Why this candidate should be tested.",
      "expected_outcome": "What the evaluation is expected to clarify."
    }
  ],
  "evidence_refs": ["allowed-evidence-id"],
  "confidence": "low",
  "limitations": ["Required limitation copied exactly from inspection."],
  "abstain_reason": null
}
```

`de_method` supports `limrots`, `limma`, `deqms`, `proda`, `rots`, and
`ensemble`. `fdr_method` supports `bh`, `ihw`, `bky`, and `storey`.
Normalization supports `none`, `median`, `quantile`, `mean`, `rlr`, and
`loess`. Imputation supports `none`, `minprob`, `mindet`, `knn`, `missforest`,
`seqknn`, `qrilc`, `impseq`, `impseqrob`, `bpca`, and `gms`.

For `de_method="ensemble"`, choose one returned ensemble preset and set
`ensemble_k` no higher than its member count. For every other DE method,
`ensemble` must be `none`. `log2fc_threshold` is either a number from 0 to 10 or
`auto`.

When abstaining, pass an empty `configs` list, no evidence references, low
confidence, and a non-empty `abstain_reason`. Do not add fields to either the
block or its config objects.
