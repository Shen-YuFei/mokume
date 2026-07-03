use std::collections::{HashMap, HashSet};

use mokume_core::{ProteinId, SampleId};

pub fn coverage_filtered_proteins<F>(
    proteins: impl IntoIterator<Item = ProteinId>,
    samples_by_condition: &HashMap<String, Vec<SampleId>>,
    threshold: f64,
    mut value_at: F,
) -> HashSet<ProteinId>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let proteins = proteins.into_iter().collect::<Vec<_>>();
    if samples_by_condition.is_empty() {
        return proteins.into_iter().collect();
    }

    proteins
        .into_iter()
        .filter(|protein| {
            samples_by_condition.values().all(|samples| {
                let valid = samples
                    .iter()
                    .filter(|sample| {
                        value_at(*protein, **sample)
                            .is_some_and(|value| value.is_finite() && value != 0.0)
                    })
                    .count();
                valid as f64 / samples.len() as f64 >= threshold
            })
        })
        .collect()
}
