use mokume_core::stats::mean_finite;
use mokume_core::{ProteinId, SampleId};

pub(crate) fn knn_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    n_neighbors: usize,
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let n_neighbors = n_neighbors.max(1);
    let matrix = proteins
        .iter()
        .map(|protein| {
            samples
                .iter()
                .map(|sample| value_at(*protein, *sample).filter(|value| value.is_finite()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let column_means = samples
        .iter()
        .enumerate()
        .map(|(sample_index, _)| {
            mean_finite(
                &matrix
                    .iter()
                    .filter_map(|row| row[sample_index])
                    .collect::<Vec<_>>(),
            )
            .unwrap_or(0.0)
        })
        .collect::<Vec<_>>();

    let mut imputed = Vec::new();
    for (receiver_index, protein) in proteins.iter().enumerate() {
        for (sample_index, sample) in samples.iter().enumerate() {
            if matrix[receiver_index][sample_index].is_some() {
                continue;
            }
            let donors = knn_donors(&matrix, receiver_index, sample_index, n_neighbors);
            let fill = if donors.is_empty() {
                column_means[sample_index]
            } else {
                donors.iter().map(|(_, value)| value).sum::<f64>() / donors.len() as f64
            };
            if fill.is_finite() {
                imputed.push((*protein, *sample, fill));
            }
        }
    }
    imputed
}

fn knn_donors(
    matrix: &[Vec<Option<f64>>],
    receiver_index: usize,
    sample_index: usize,
    n_neighbors: usize,
) -> Vec<(f64, f64)> {
    let mut donors = matrix
        .iter()
        .enumerate()
        .filter_map(|(donor_index, donor)| {
            let value = donor[sample_index]?;
            let distance = nan_euclidean_distance(&matrix[receiver_index], donor)?;
            Some((distance, donor_index, value))
        })
        .collect::<Vec<_>>();
    donors.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
    });
    donors
        .into_iter()
        .take(n_neighbors)
        .map(|(distance, _, value)| (distance, value))
        .collect()
}

fn nan_euclidean_distance(left: &[Option<f64>], right: &[Option<f64>]) -> Option<f64> {
    let mut observed = 0usize;
    let mut squared = 0.0;
    for (left, right) in left.iter().zip(right) {
        let (Some(left), Some(right)) = (left, right) else {
            continue;
        };
        observed += 1;
        let delta = left - right;
        squared += delta * delta;
    }
    (observed > 0).then(|| (squared * left.len() as f64 / observed as f64).sqrt())
}
