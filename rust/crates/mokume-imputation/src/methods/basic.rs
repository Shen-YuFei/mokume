//! Simple per-sample-statistic and constant imputation strategies.

use std::collections::HashMap;

use mokume_core::stats::{mean_finite, median_finite};
use mokume_core::{ProteinId, SampleId};

use crate::support::{fill_missing_by_sample, finite_sample_values, has_finite_value};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SampleImputationStat {
    Mean,
    Median,
    MostFrequent,
}

pub(crate) fn sample_stat_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    stat: SampleImputationStat,
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let fills = samples
        .iter()
        .filter_map(|sample| {
            let mut observed = finite_sample_values(proteins, *sample, value_at);
            let fill = match stat {
                SampleImputationStat::Mean => mean_finite(&observed),
                SampleImputationStat::Median => median_finite(&mut observed),
                SampleImputationStat::MostFrequent => most_frequent_finite(&mut observed),
            }?;
            Some((*sample, fill))
        })
        .collect::<HashMap<_, _>>();
    fill_missing_by_sample(proteins, samples, &fills, value_at)
}

/// Most-frequent value of a column, matching `sklearn.impute.SimpleImputer`'s
/// `strategy="most_frequent"` (which delegates to `scipy.stats.mode`). Ties in
/// count are broken by the smallest value, exactly as scipy returns the smallest
/// modal value. Mutates `values` by sorting; returns `None` for an empty column.
fn most_frequent_finite(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    // Sort ascending so equal values are contiguous and, on a count tie, the
    // first (smallest) run is kept.
    values.sort_by(f64::total_cmp);
    let mut best_value = values[0];
    let mut best_count = 0usize;
    let mut run_value = values[0];
    let mut run_count = 0usize;
    for &value in values.iter() {
        if value == run_value {
            run_count += 1;
        } else {
            // Strict `>` keeps the smaller value on a tie (runs ascend).
            if run_count > best_count {
                best_count = run_count;
                best_value = run_value;
            }
            run_value = value;
            run_count = 1;
        }
    }
    if run_count > best_count {
        best_value = run_value;
    }
    Some(best_value)
}

pub(crate) fn constant_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    fill: f64,
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let mut imputed = Vec::new();
    for protein in proteins {
        for sample in samples {
            if !has_finite_value(*protein, *sample, value_at) {
                imputed.push((*protein, *sample, fill));
            }
        }
    }
    imputed
}
