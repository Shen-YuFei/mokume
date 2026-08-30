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
    let Some(reference) = reference_peptide_row(&rows, &log_rows) else {
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

fn reference_peptide_row(
    rows: &[(PeptideId, Vec<Option<f64>>)],
    log_rows: &[Vec<Option<f64>>],
) -> Option<usize> {
    let max_support = log_rows
        .iter()
        .map(|values| values.iter().filter(|value| value.is_some()).count())
        .max()?;
    let mut candidates = log_rows
        .iter()
        .enumerate()
        .filter_map(|(index, values)| {
            (values.iter().filter(|value| value.is_some()).count() == max_support).then_some(index)
        })
        .collect::<Vec<_>>();
    if candidates.len() == 1 {
        return candidates.pop();
    }

    let trace_totals = candidates
        .iter()
        .map(|&index| {
            let mut values = log_rows[index]
                .iter()
                .filter_map(|value| *value)
                .collect::<Vec<_>>();
            values.sort_by(f64::total_cmp);
            (index, values.into_iter().sum::<f64>())
        })
        .collect::<Vec<_>>();
    let max_total = trace_totals
        .iter()
        .map(|(_, total)| *total)
        .max_by(f64::total_cmp)?;

    trace_totals
        .into_iter()
        .filter(|(_, total)| *total == max_total)
        .min_by_key(|(index, _)| rows[*index].0.get())
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

#[cfg(test)]
mod tests {
    use super::{max_lfq_with_samples, reference_peptide_row};
    use crate::PeptideMeasurement;
    use mokume_core::{PeptideId, SampleId};

    fn measurement(peptide: u32, sample: u32, intensity: f64) -> PeptideMeasurement {
        PeptideMeasurement {
            peptide: PeptideId::new(peptide),
            sample: SampleId::new(sample),
            intensity,
        }
    }

    fn assert_profiles_close(actual: &[(SampleId, f64)], expected: &[(SampleId, f64)]) {
        assert_eq!(actual.len(), expected.len());
        for ((actual_sample, actual_value), (expected_sample, expected_value)) in
            actual.iter().zip(expected)
        {
            assert_eq!(actual_sample, expected_sample);
            let tolerance = expected_value.abs().max(1.0) * 1e-12;
            assert!(
                (actual_value - expected_value).abs() <= tolerance,
                "sample {actual_sample:?}: {actual_value} vs {expected_value}"
            );
        }
    }

    #[test]
    fn reference_selection_matches_python_tie_break_order() {
        let rows = vec![
            (
                PeptideId::new(20),
                vec![Some(128.0), None, Some(2.0), Some(64.0)],
            ),
            (PeptideId::new(10), vec![Some(1e30), None, None, Some(1e30)]),
            (
                PeptideId::new(30),
                vec![Some(2.0), Some(32.0), Some(32.0), None],
            ),
            (
                PeptideId::new(5),
                vec![Some(16.0), Some(32.0), None, Some(32.0)],
            ),
        ];
        let log_rows = vec![
            vec![Some(7.0), None, Some(1.0), Some(6.0)],
            vec![Some(100.0), None, None, Some(100.0)],
            vec![Some(1.0), Some(5.0), Some(5.0), None],
            vec![Some(4.0), Some(5.0), None, Some(5.0)],
        ];

        // Rows 0, 2, and 3 have the largest support (three samples). Rows 0
        // and 3 then tie at a sorted log-total of 14, so the smallest stable
        // peptide id wins the final tie, matching Python's peptide-name key
        // once callers map names to lexical ids.
        assert_eq!(reference_peptide_row(&rows, &log_rows), Some(3));
    }

    #[test]
    fn max_lfq_is_invariant_to_encounter_order_ids() {
        let samples = [
            SampleId::new(0),
            SampleId::new(1),
            SampleId::new(2),
            SampleId::new(3),
        ];
        let forward = vec![
            measurement(0, 0, 128.0),
            measurement(0, 2, 2.0),
            measurement(0, 3, 64.0),
            measurement(1, 0, 64.0),
            measurement(1, 3, 4.0),
            measurement(2, 0, 2.0),
            measurement(2, 1, 32.0),
            measurement(2, 2, 32.0),
        ];
        // Same physical peptide traces in reverse encounter order. The caller's
        // insertion-order registry therefore assigns the first and last
        // peptides opposite numeric ids.
        let reversed = vec![
            measurement(0, 2, 32.0),
            measurement(0, 1, 32.0),
            measurement(0, 0, 2.0),
            measurement(1, 3, 4.0),
            measurement(1, 0, 64.0),
            measurement(2, 3, 64.0),
            measurement(2, 2, 2.0),
            measurement(2, 0, 128.0),
        ];

        let forward = max_lfq_with_samples(&forward, &samples);
        let reversed = max_lfq_with_samples(&reversed, &samples);
        assert_profiles_close(&forward, &reversed);
        let expected = samples
            .into_iter()
            .zip([
                173.941_622_437_385_3,
                86.970_811_218_692_65,
                15.374_412_594_508_161,
                51.713_153_749_413_884,
            ])
            .collect::<Vec<_>>();
        assert_profiles_close(&forward, &expected);
    }
}
