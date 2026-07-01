use std::collections::HashMap;

use mokume_core::stats::median_finite;
use mokume_core::{ProteinId, SampleId};

pub fn irs_scaling_factors<F>(
    proteins: impl IntoIterator<Item = ProteinId>,
    plexes: &[String],
    refs_by_plex: &HashMap<String, Vec<SampleId>>,
    stat: &str,
    mut value_at: F,
) -> HashMap<ProteinId, HashMap<String, f64>>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let mut factors = HashMap::<ProteinId, HashMap<String, f64>>::new();
    for protein in proteins {
        let mut refs = Vec::<(String, f64)>::new();
        for plex in plexes {
            let Some(samples) = refs_by_plex.get(plex) else {
                continue;
            };
            if let Some(reference) = irs_reference_value(protein, samples, stat, &mut value_at) {
                refs.push((plex.clone(), reference));
            }
        }
        let positive_refs = refs
            .iter()
            .filter_map(|(_, value)| ((*value > 0.0) && value.is_finite()).then_some((*value).ln()))
            .collect::<Vec<_>>();
        if positive_refs.is_empty() {
            continue;
        }
        let global_ref = (positive_refs.iter().sum::<f64>() / positive_refs.len() as f64).exp();
        for (plex, reference) in refs {
            if reference > 0.0 && reference.is_finite() {
                factors
                    .entry(protein)
                    .or_default()
                    .insert(plex, global_ref / reference);
            }
        }
    }
    factors
}

fn irs_reference_value<F>(
    protein: ProteinId,
    samples: &[SampleId],
    stat: &str,
    value_at: &mut F,
) -> Option<f64>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let mut values = samples
        .iter()
        .filter_map(|sample| value_at(protein, *sample))
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    if values.is_empty() {
        return None;
    }
    if stat.eq_ignore_ascii_case("mean") {
        Some(values.iter().sum::<f64>() / values.len() as f64)
    } else {
        median_finite(&mut values)
    }
}
