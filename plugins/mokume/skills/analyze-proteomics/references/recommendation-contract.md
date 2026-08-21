# MCP input contract

Call `mokume.inspect_dataset` with the absolute matrix and SDRF paths, an exact
two-condition contrast, an explicit input scale, and optional peptide-count and
acquisition inputs:

```json
{
  "protein_matrix": "/absolute/proteins.tsv",
  "sdrf": "/absolute/project.sdrf.tsv",
  "contrast": ["control", "treated"],
  "options": {
    "input_scale": "linear",
    "peptide_counts": "/absolute/peptide_counts.tsv",
    "data_type": "LFQ",
    "quantification": "directlfq",
    "upstream_engine": "quantms",
    "factor_column": null
  }
}
```

The protein matrix may be comma- or tab-delimited and must satisfy all of these
requirements:

- The first column contains non-empty, unique protein identifiers.
- At least two later columns contain samples; all sample cells are numeric or
  missing (`NaN`). Positive and negative infinity are rejected.
- Column names are non-empty and unique, and every sample column maps to the SDRF.
- A linear matrix contains at least one positive finite intensity; a log2 matrix
  contains at least one finite intensity.

`options.input_scale` is required and must declare whether those intensities are
`linear`
or `log2`; `auto` is not supported. In a linear matrix, non-positive cells are
canonicalized to missing before profiling and evaluation. In a log2 matrix,
finite zero and negative cells remain observed. `contrast` is required, contains
exactly two distinct SDRF factor labels, and scopes both profiling and evaluation
to those matrix columns; unrelated conditions cannot affect diagnostics or
preprocessing. The remaining options fields are optional.
`options.peptide_counts` is optional only for count-independent candidates; do not invent
values that were not supplied or supported by the inputs. Without a declared
`data_type`, Mokume only infers `LFQ`, `DIA`, or `TMT` from explicit sample-name
markers. Generic names such as `S1` or `sample-01` produce `unknown`, force an
abstention, and require a supported declaration before reinspection. Known engine
aliases are normalized to the catalog spelling; for example, `DIANN`, `DIA NN`,
and `DIA-NN` all become `DIA-NN`.

The peptide-count sidecar may be comma- or tab-delimited and must have exactly
these columns:

```text
protein\tpeptide_count
P12345\t7
P67890\t3
```

Protein identifiers must be unique and match the protein-matrix identifier
values. Counts are positive integers representing unique peptides per protein.
The sidecar may include proteins outside the matrix, but at least one identifier
must overlap. A sidecar is mandatory for `deqms` and every ensemble containing
`deqms`; without it, deterministic policy omits those candidates and evaluation
rejects a host-added candidate. Individual matrix proteins absent from a supplied
sidecar use a count of one.

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
    "input_scale": "linear",
    "peptide_counts": "/absolute/peptide_counts.tsv",
    "threads": 24
  }
}
```

The options object accepts only the fields shown above. Repeat the same two
contrast labels accepted by inspection, using the canonical keys returned in
`profile.samples_per_condition`; do not replace those labels with longer raw SDRF
values. `ground_truth`, when present, must be an absolute path to a
one-protein-per-line file and requires
`expected_direction` to be `UP` or `DOWN`. When `ground_truth` is null,
`expected_direction` must also be null and carries no biological meaning.
`options.output_dir` must be an absolute path that does not already exist.
`options.input_scale` is required and must be `linear` or `log2`. When inspection
used a peptide-count sidecar, repeat it in `options.peptide_counts`, together with
the declared scale, `data_type`, `quantification`, `upstream_engine`, and
`factor_column` values used for `inspect_dataset` so evaluation rebinds the same
policy context. `options.peptide_counts` is required if any candidate uses
`deqms` directly or through an ensemble. Mokume writes the round to a sibling
staging directory and publishes `output_dir` only after every artifact succeeds;
on failure, the target remains absent and may be retried.

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
`ensemble` must be `none`. DEqMS and every ensemble preset containing DEqMS
require `options.peptide_counts`. `log2fc_threshold` is either a number from 0 to
10 or `auto`.

When abstaining, pass an empty `configs` list, no evidence references, low
confidence, and a non-empty `abstain_reason`. Do not add fields to either the
block or its config objects.
