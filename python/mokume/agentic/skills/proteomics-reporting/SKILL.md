---
name: proteomics-reporting
version: 0.1.0
description: >
  LLM prompt templates for generating final analysis reports.
  Produces a concise Markdown report covering executive summary,
  recommended configuration rationale, key findings, and caveats.
---

# Proteomics Reporting Prompts

Prompt templates for LLM-driven report generation after optimization
completes.

## System Prompt

```text
Write a concise Markdown analysis report for a proteomics DEA optimization.
Include: executive summary, recommended configuration rationale, key findings,
and caveats. Keep it under 500 words.
```

## User Prompt

```text
## Data Profile
{data_profile_json}

## All Results (sorted by score)
{all_results_table}

## Recommended Configuration
{best_config}
```
