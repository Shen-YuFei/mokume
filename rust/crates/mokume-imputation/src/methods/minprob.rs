use std::collections::HashMap;

use mokume_core::stats::{finite_sd, quantile_linear};
use mokume_core::{ProteinId, SampleId};

use crate::support::{fill_missing_by_sample, finite_sample_values};

pub(crate) fn minprob_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    quantile: f64,
    shift: f64,
    scale: f64,
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let fills = samples
        .iter()
        .filter_map(|sample| {
            let mut observed = finite_sample_values(proteins, *sample, value_at);
            let q_low = quantile_linear(&mut observed, quantile)?;
            let sd = finite_sd(&observed).map(|sd| sd * scale).unwrap_or(0.1);
            let sd = if sd < 1e-10 { 0.1 } else { sd };
            let fill = q_low - shift * sd;
            fill.is_finite().then_some((*sample, fill))
        })
        .collect::<HashMap<_, _>>();
    fill_missing_by_sample(proteins, samples, &fills, value_at)
}
