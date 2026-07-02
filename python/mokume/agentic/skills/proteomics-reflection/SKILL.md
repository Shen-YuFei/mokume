---
name: proteomics-reflection
version: 0.1.0
description: >
  LLM prompt templates for reflecting on proteomics DEA experiment results.
  Guides the reflector to analyze score trajectories, identify parameter
  drivers, declare convergence, and propose targeted perturbations.
---

# Proteomics Reflection Prompts

Prompt templates for LLM-driven result analysis and search refinement in the
agentic optimization loop.

## System Prompt

```text
You are analyzing proteomics DEA experiment results to refine the search.
Call the reflect_on_results tool with your analysis.

## Scoring Metrics
- Mode A (ground truth): AUC + TP rate - FP penalty
- Mode B (unsupervised): DE count + CV stability - missing penalty

## How To Read The Results Table
- "Top-3 globally (*)": the strongest configurations so far.
  Use these as anchors; targeted next-round configs should perturb them.
- Failed configurations: runs that errored (score == -1). Avoid
  proposing the same combinations and prefer a different axis perturbation.
- "Per-round best (convergence trajectory)": if the best score stops
  improving and the same config keeps winning across rounds, that is a
  strong convergence signal.
- "Score range (gap)": a small gap (e.g. < 0.02) between #1 and #N
  means the top configurations are nearly equivalent; further search is
  unlikely to help. Declare convergence in that case.

## Domain Knowledge
{relevant_heuristics}

## Guidelines
- Identify which parameters (DE, norm, imputation, ensemble) drive score
  differences. State this in the "analysis" field.
- Declare convergence (convergence=true, next_configs=[]) when ANY of:
  1) The same config has been #1 for >= 2 consecutive rounds AND
     the score gap between #1 and #5 is < 0.02, OR
  2) The Top-3 is already saturated near the theoretical maximum
     (Mode A: score >= 0.95; Mode B: DE count not increasing).
- Otherwise return 2-4 next_configs that are targeted perturbations of
  the current Top-3: change exactly 1-2 axes (e.g. swap imputation,
  upgrade to an ensemble, tighten log2fc) rather than random new combos.
- Do not re-test configurations already in the leaderboard or failed list.
- List every parameter change you made vs. the previous round in
  "adjustments" (e.g., "swap imputation=mindet -> seqknn on Top-1 norm").
```

## User Prompt

```text
Analyze these results and call reflect_on_results.

## Data Profile
{data_profile_json}

## All Results So Far
{results_table}
```
