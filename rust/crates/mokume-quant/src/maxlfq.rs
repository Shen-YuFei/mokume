use std::collections::HashMap;

use mokume_core::{PeptideId, SampleId};

use crate::{median, sample_order, PeptideMeasurement};

pub fn max_lfq(measurements: &[PeptideMeasurement]) -> Vec<(SampleId, f64)> {
    let samples = sample_order(measurements);
    solve_max_lfq_builtin(measurements, &samples)
}

pub fn max_lfq_with_samples(
    measurements: &[PeptideMeasurement],
    samples: &[SampleId],
) -> Vec<(SampleId, f64)> {
    solve_max_lfq_builtin(measurements, samples)
}

fn solve_max_lfq_builtin(
    measurements: &[PeptideMeasurement],
    samples: &[SampleId],
) -> Vec<(SampleId, f64)> {
    if samples.is_empty() {
        return Vec::new();
    }

    let sample_index = samples
        .iter()
        .enumerate()
        .map(|(index, sample)| (*sample, index))
        .collect::<HashMap<_, _>>();
    let rows = peptide_intensity_rows(measurements, &sample_index, samples.len());
    if rows.is_empty() {
        return Vec::new();
    }

    if samples.len() == 1 {
        let mut values = rows
            .iter()
            .filter_map(|(_, values)| values[0])
            .collect::<Vec<_>>();
        return median(&mut values)
            .map(|value| vec![(samples[0], value)])
            .unwrap_or_default();
    }

    if rows.len() == 1 {
        return rows[0]
            .1
            .iter()
            .enumerate()
            .filter_map(|(index, value)| value.map(|value| (samples[index], value)))
            .collect();
    }

    let original_sum = rows
        .iter()
        .flat_map(|(_, values)| values.iter().filter_map(|value| *value))
        .sum::<f64>();
    if original_sum <= 0.0 {
        return Vec::new();
    }

    let log_rows = rows
        .iter()
        .map(|(_, values)| {
            values
                .iter()
                .map(|value| value.map(f64::log2))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let Some(reference) = reference_peptide_row(&log_rows) else {
        return Vec::new();
    };
    let aligned = align_to_reference(&log_rows, reference);
    let mut intensities = median_aligned_by_sample(&aligned)
        .into_iter()
        .enumerate()
        .map(|(index, value)| {
            value.map(|value| value.exp2()).or_else(|| {
                let mut raw_values = rows
                    .iter()
                    .filter_map(|(_, values)| values[index])
                    .collect::<Vec<_>>();
                median(&mut raw_values)
            })
        })
        .collect::<Vec<_>>();

    let current_sum = intensities.iter().filter_map(|value| *value).sum::<f64>();
    if current_sum > 0.0 {
        let scale = original_sum / current_sum;
        for value in intensities.iter_mut().flatten() {
            *value *= scale;
        }
    }

    intensities
        .into_iter()
        .enumerate()
        .filter_map(|(index, value)| {
            value
                .filter(|value| value.is_finite() && *value > 0.0)
                .map(|value| (samples[index], value))
        })
        .collect()
}

fn peptide_intensity_rows(
    measurements: &[PeptideMeasurement],
    sample_index: &HashMap<SampleId, usize>,
    sample_count: usize,
) -> Vec<(PeptideId, Vec<Option<f64>>)> {
    let mut rows = HashMap::<PeptideId, Vec<Option<f64>>>::new();
    for measurement in measurements {
        if !measurement.intensity.is_finite() || measurement.intensity <= 0.0 {
            continue;
        }
        let Some(sample) = sample_index.get(&measurement.sample).copied() else {
            continue;
        };
        let values = rows
            .entry(measurement.peptide)
            .or_insert_with(|| vec![None; sample_count]);
        // Python's built-in MaxLFQ pivots the (peptide, sample) matrix with
        // aggfunc="sum", so duplicate (peptide, sample) intensities are summed.
        values[sample] = values[sample]
            .map(|current| current + measurement.intensity)
            .or(Some(measurement.intensity));
    }

    let mut rows = rows.into_iter().collect::<Vec<_>>();
    rows.sort_by_key(|(peptide, _)| peptide.get());
    rows
}

fn reference_peptide_row(log_rows: &[Vec<Option<f64>>]) -> Option<usize> {
    log_rows
        .iter()
        .enumerate()
        .max_by_key(|(_, values)| values.iter().filter(|value| value.is_some()).count())
        .map(|(index, _)| index)
}

fn align_to_reference(log_rows: &[Vec<Option<f64>>], reference: usize) -> Vec<Vec<Option<f64>>> {
    let reference_trace = &log_rows[reference];
    log_rows
        .iter()
        .enumerate()
        .map(|(index, trace)| {
            if index == reference {
                return trace.clone();
            }
            let mut shifts = reference_trace
                .iter()
                .zip(trace)
                .filter_map(|(reference_value, value)| Some((*reference_value)? - (*value)?))
                .collect::<Vec<_>>();
            let Some(shift) = median(&mut shifts) else {
                return trace.clone();
            };
            trace
                .iter()
                .map(|value| value.map(|value| value + shift))
                .collect()
        })
        .collect()
}

fn median_aligned_by_sample(aligned: &[Vec<Option<f64>>]) -> Vec<Option<f64>> {
    let Some(sample_count) = aligned.first().map(Vec::len) else {
        return Vec::new();
    };
    (0..sample_count)
        .map(|sample| {
            let mut values = aligned
                .iter()
                .filter_map(|trace| trace[sample])
                .collect::<Vec<_>>();
            median(&mut values)
        })
        .collect()
}
