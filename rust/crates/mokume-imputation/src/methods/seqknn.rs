use mokume_core::{MokumeError, ProteinId, Result, SampleId};

pub(crate) fn seqknn_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    n_neighbors: usize,
    value_at: &mut F,
) -> Result<Vec<(ProteinId, SampleId, f64)>>
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
    if matrix.is_empty() || samples.is_empty() {
        return Ok(Vec::new());
    }

    let missing_counts = matrix
        .iter()
        .map(|row| row.iter().filter(|value| value.is_none()).count())
        .collect::<Vec<_>>();
    if missing_counts.iter().all(|count| *count == 0) {
        return Ok(Vec::new());
    }

    let original_missing = matrix
        .iter()
        .map(|row| row.iter().map(Option::is_none).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut order = (0..proteins.len()).collect::<Vec<_>>();
    order.sort_by(|left, right| {
        missing_counts[*left]
            .cmp(&missing_counts[*right])
            .then_with(|| left.cmp(right))
    });

    let mut sorted_matrix = order
        .iter()
        .map(|index| matrix[*index].clone())
        .collect::<Vec<_>>();
    let mut complete_rows = sorted_matrix
        .iter()
        .take_while(|row| row.iter().all(Option::is_some))
        .cloned()
        .collect::<Vec<_>>();
    if complete_rows.is_empty() {
        return Err(MokumeError::InvalidInput {
            message: "SeqKNN requires an initial set of complete protein rows".to_owned(),
        });
    }

    for row in sorted_matrix.iter_mut().skip(complete_rows.len()) {
        *row = seqknn_imputed_row(row, &complete_rows, n_neighbors);
        complete_rows.push(row.clone());
    }

    let mut imputed = Vec::new();
    for (sorted_index, original_index) in order.into_iter().enumerate() {
        for (sample_index, sample) in samples.iter().enumerate() {
            if !original_missing[original_index][sample_index] {
                continue;
            }
            let Some(value) = sorted_matrix[sorted_index][sample_index] else {
                continue;
            };
            if value.is_finite() {
                imputed.push((proteins[original_index], *sample, value));
            }
        }
    }
    Ok(imputed)
}

fn seqknn_imputed_row(
    target: &[Option<f64>],
    complete_rows: &[Vec<Option<f64>>],
    n_neighbors: usize,
) -> Vec<Option<f64>> {
    let missing_indices = target
        .iter()
        .enumerate()
        .filter_map(|(index, value)| value.is_none().then_some(index))
        .collect::<Vec<_>>();
    let observed_indices = target
        .iter()
        .enumerate()
        .filter_map(|(index, value)| value.is_some().then_some(index))
        .collect::<Vec<_>>();
    if missing_indices.is_empty() || complete_rows.is_empty() {
        return target.to_vec();
    }

    let mut neighbors = complete_rows
        .iter()
        .enumerate()
        .filter_map(|(index, row)| {
            let distance = seqknn_distance(target, row, &observed_indices)?;
            Some((distance, index, row))
        })
        .collect::<Vec<_>>();
    neighbors.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
    });

    let neighbors = neighbors.into_iter().take(n_neighbors).collect::<Vec<_>>();
    let mut result = target.to_vec();
    for sample_index in missing_indices {
        let mut weighted_sum = 0.0;
        let mut weight_sum = 0.0;
        for (distance, _, row) in &neighbors {
            let Some(value) = row[sample_index] else {
                continue;
            };
            // Original SeqKnn 1.0.1 nnmiss uses squared Euclidean distance.
            let weight = 1.0 / (distance + 1e-15);
            weighted_sum += value * weight;
            weight_sum += weight;
        }
        if weight_sum > 0.0 {
            result[sample_index] = Some(weighted_sum / weight_sum);
        }
    }
    result
}

fn seqknn_distance(
    target: &[Option<f64>],
    row: &[Option<f64>],
    observed_indices: &[usize],
) -> Option<f64> {
    let mut squared = 0.0;
    for index in observed_indices {
        let target_value = target[*index]?;
        let row_value = row[*index]?;
        let delta = row_value - target_value;
        squared += delta * delta;
    }
    Some(squared)
}
