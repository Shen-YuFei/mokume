use mokume_core::stats::quantile_linear;
use mokume_core::{ProteinId, SampleId};

use crate::support::{fill_missing_by_sample, finite_sample_values};

pub(crate) fn mindet_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    quantile: f64,
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let fills = samples
        .iter()
        .filter_map(|sample| {
            let mut observed = finite_sample_values(proteins, *sample, value_at);
            quantile_linear(&mut observed, quantile).map(|fill| (*sample, fill))
        })
        .collect::<std::collections::HashMap<_, _>>();
    fill_missing_by_sample(proteins, samples, &fills, value_at)
}
