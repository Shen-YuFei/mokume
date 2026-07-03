use std::collections::HashMap;

use mokume_core::stats::median;
use mokume_core::Result;

use crate::math::{unsupported, valid_scale};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleNormalizationMethod {
    GlobalMedian,
    ConditionMedian,
    Quantile,
    MedianCenter,
    MeanCenter,
    Rlr,
    Loess,
    Hierarchical,
    Tmm,
}

pub fn parse_sample_normalization_method(
    method: &str,
) -> Result<Option<SampleNormalizationMethod>> {
    match method.trim().to_ascii_lowercase().as_str() {
        "" | "none" | "no" | "false" => Ok(None),
        "globalmedian" | "global_median" | "median" => {
            Ok(Some(SampleNormalizationMethod::GlobalMedian))
        }
        "conditionmedian" | "condition_median" => {
            Ok(Some(SampleNormalizationMethod::ConditionMedian))
        }
        "quantile" => Ok(Some(SampleNormalizationMethod::Quantile)),
        "mediancenter" | "median_center" => Ok(Some(SampleNormalizationMethod::MedianCenter)),
        "meancenter" | "mean_center" => Ok(Some(SampleNormalizationMethod::MeanCenter)),
        "rlr" => Ok(Some(SampleNormalizationMethod::Rlr)),
        "loess" => Ok(Some(SampleNormalizationMethod::Loess)),
        "hierarchical" => Ok(Some(SampleNormalizationMethod::Hierarchical)),
        "tmm" => Ok(Some(SampleNormalizationMethod::Tmm)),
        _ => unsupported("sample-normalization-method"),
    }
}

pub fn global_median_sample_factors(
    sample_values: HashMap<String, Vec<f64>>,
) -> HashMap<String, f64> {
    let sample_medians = sample_medians(sample_values);
    let mut medians = sample_medians.values().copied().collect::<Vec<_>>();
    let Some(global_median) = median(&mut medians) else {
        return HashMap::new();
    };
    sample_medians
        .into_iter()
        .filter_map(|(sample, sample_median)| {
            valid_scale(sample_median).then_some((sample, global_median / sample_median))
        })
        .collect()
}

pub fn condition_median_sample_factors(
    condition_sample_values: HashMap<String, HashMap<String, Vec<f64>>>,
) -> HashMap<String, f64> {
    let mut factors = HashMap::new();
    for sample_values in condition_sample_values.into_values() {
        let sample_medians = sample_medians(sample_values);
        if sample_medians.is_empty() {
            continue;
        }
        let condition_mean = sample_medians.values().sum::<f64>() / sample_medians.len() as f64;
        if !valid_scale(condition_mean) {
            continue;
        }
        factors.extend(
            sample_medians
                .into_iter()
                .filter_map(|(sample, sample_median)| {
                    valid_scale(sample_median).then_some((sample, condition_mean / sample_median))
                }),
        );
    }
    factors
}

fn sample_medians(sample_values: HashMap<String, Vec<f64>>) -> HashMap<String, f64> {
    sample_values
        .into_iter()
        .filter_map(|(sample, mut values)| median(&mut values).map(|median| (sample, median)))
        .collect()
}

/// TMM (Trimmed Mean of M-values) trimming proportions and reference-selection
/// defaults, matching `TMMNormalizer(m_trim=0.3, a_trim=0.05, ref_sample=None,
/// log_transform=False)` in `python/mokume/normalization/tmm.py`.
const TMM_M_TRIM: f64 = 0.3;
const TMM_A_TRIM: f64 = 0.05;

/// Compute per-sample TMM normalization factors for a raw-intensity wide matrix.
///
/// This reproduces `TMMNormalizer` (`python/mokume/normalization/tmm.py`) with
/// `m_trim=0.3, a_trim=0.05, ref_sample=None, log_transform=False`. The matrix
/// is given column-major: `columns[j]` is the name of sample `j`, and
/// `matrix[j]` holds that sample's per-feature raw intensities in a fixed row
/// order shared across every column (row `i` is the same feature in every
/// column). Missing/zero entries must be encoded as `0.0` or a non-finite value
/// (`NaN`/`inf`), exactly matching Python's `replace(0, NaN)` semantics.
///
/// Returns the rescaled factors keyed by sample name (`factor_s /
/// geometric_mean`, tmm.py:296-299). Samples with no valid data still receive a
/// factor (defaulting to 1.0 before rescaling, tmm.py:159/175/231). Divide each
/// sample column by its returned factor to obtain the normalized intensities
/// (tmm.py:334-343 uses division, not multiplication).
///
/// All math is done in `f64` to match NumPy's float64, even when the caller's
/// intensities originate from `f32` cells.
pub fn tmm_norm_factors(columns: &[String], matrix: &[Vec<f64>]) -> HashMap<String, f64> {
    let n_samples = columns.len();
    if n_samples == 0 || matrix.len() != n_samples {
        return HashMap::new();
    }

    // Library sizes: sum of the non-zero finite values per sample
    // (`X_raw.replace(0, np.nan).sum(skipna=True)`, tmm.py:287). Zeros and
    // non-finite entries are treated as missing.
    let lib_sizes: Vec<f64> = matrix
        .iter()
        .map(|column| {
            column
                .iter()
                .filter(|value| value.is_finite() && **value != 0.0)
                .sum::<f64>()
        })
        .collect();

    // Reference sample selection (`_select_reference_sample`, tmm.py:100-132):
    // the upper quartile (75th percentile, linear interpolation) of the
    // non-zero finite values per sample, then the sample whose UQ is closest to
    // the mean UQ. NaN UQ columns are excluded from the mean but `(uq -
    // mean).abs().idxmin()` over a Series that still contains NaN never selects
    // a NaN (NaN compares as not-minimal), so we mirror that by skipping NaN
    // columns and breaking ties toward the first column in order.
    let upper_quartiles: Vec<f64> = matrix
        .iter()
        .map(|column| {
            let mut vals: Vec<f64> = column
                .iter()
                .copied()
                .filter(|value| value.is_finite() && *value != 0.0)
                .collect();
            if vals.is_empty() {
                f64::NAN
            } else {
                percentile_linear(&mut vals, 75.0)
            }
        })
        .collect();

    let finite_uqs: Vec<f64> = upper_quartiles
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    let ref_index = if finite_uqs.is_empty() {
        0
    } else {
        let mean_uq = finite_uqs.iter().sum::<f64>() / finite_uqs.len() as f64;
        let mut best_index = 0usize;
        let mut best_distance = f64::INFINITY;
        for (index, uq) in upper_quartiles.iter().enumerate() {
            if !uq.is_finite() {
                continue;
            }
            let distance = (uq - mean_uq).abs();
            if distance < best_distance {
                best_distance = distance;
                best_index = index;
            }
        }
        best_index
    };

    // Per-sample factor vs the reference (`_compute_tmm_factor`, tmm.py:134-245).
    let mut raw_factors = vec![1.0f64; n_samples];
    for sample_index in 0..n_samples {
        raw_factors[sample_index] = if sample_index == ref_index {
            1.0
        } else {
            compute_tmm_factor(
                &matrix[sample_index],
                &matrix[ref_index],
                lib_sizes[sample_index],
                lib_sizes[ref_index],
            )
        };
    }

    // Rescale so the geometric mean of the factors is 1
    // (`geometric_mean = exp(mean(log(factor)))`, `norm = factor /
    // geometric_mean`, tmm.py:296-299). `np.log` of a non-positive factor is
    // NaN in NumPy; guard against it so the rescale stays finite.
    let log_sum: f64 = raw_factors
        .iter()
        .map(|factor| if *factor > 0.0 { factor.ln() } else { 0.0 })
        .sum();
    let geometric_mean = (log_sum / n_samples as f64).exp();

    columns
        .iter()
        .zip(raw_factors.iter())
        .map(|(name, factor)| {
            let rescaled = if geometric_mean > 0.0 && geometric_mean.is_finite() {
                factor / geometric_mean
            } else {
                *factor
            };
            (name.clone(), rescaled)
        })
        .collect()
}

/// TMM factor for one sample vs the reference (`_compute_tmm_factor`,
/// tmm.py:134-245). Both slices are per-feature raw intensities in the same row
/// order; `n_sample`/`n_ref` are the library sizes.
/// A feature is usable in the fold-change when it is strictly positive and
/// finite in both the sample and reference columns (tmm.py:171-173).
fn tmm_valid_pair(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
}

fn compute_tmm_factor(y_sample: &[f64], y_ref: &[f64], n_sample: f64, n_ref: f64) -> f64 {
    // Valid features: strictly positive and finite in both columns
    // (tmm.py:171-173). Fewer than 10 -> factor 1.0 (tmm.py:175).
    let mut ys = Vec::<f64>::new();
    let mut yr = Vec::<f64>::new();
    for (a, b) in y_sample.iter().zip(y_ref.iter()) {
        if tmm_valid_pair(*a, *b) {
            ys.push(*a);
            yr.push(*b);
        }
    }
    let n = ys.len();
    if n < 10 {
        return 1.0;
    }

    // M-values (log2 fold-change), A-values (average log expression) and
    // weights (inverse asymptotic variance) (tmm.py:187-199).
    let m_values: Vec<f64> = (0..n)
        .map(|i| ((ys[i] / n_sample) / (yr[i] / n_ref)).log2())
        .collect();
    let a_values: Vec<f64> = (0..n)
        .map(|i| 0.5 * ((ys[i] / n_sample).log2() + (yr[i] / n_ref).log2()))
        .collect();
    let weights: Vec<f64> = (0..n)
        .map(|i| {
            let w = (n_sample - ys[i]) / (n_sample * ys[i]) + (n_ref - yr[i]) / (n_ref * yr[i]);
            if w > 0.0 {
                1.0 / w
            } else {
                0.0
            }
        })
        .collect();

    // Double trimming bounds (tmm.py:202-206). `int(np.floor(...))` and
    // `int(np.ceil(...))` on the feature count.
    let n_f = n as f64;
    let m_lo = (n_f * TMM_M_TRIM).floor() as usize;
    let m_hi = (n_f * (1.0 - TMM_M_TRIM)).ceil() as usize;
    let a_lo = (n_f * TMM_A_TRIM).floor() as usize;
    let a_hi = (n_f * (1.0 - TMM_A_TRIM)).ceil() as usize;

    // keep_m / keep_a via ascending stable argsort of the M and A values
    // (`np.argsort`, tmm.py:209-216). A feature is kept when its rank falls in
    // the half-open [lo, hi) window.
    let m_order = stable_argsort(&m_values);
    let a_order = stable_argsort(&a_values);
    let mut keep_m = vec![false; n];
    for &feature in m_order.iter().take(m_hi).skip(m_lo) {
        keep_m[feature] = true;
    }
    let mut keep_a = vec![false; n];
    for &feature in a_order.iter().take(a_hi).skip(a_lo) {
        keep_a[feature] = true;
    }

    // keep = keep_m & keep_a; fall back to keep_m if fewer than 5 survive; if
    // still fewer than 5 -> factor 1.0 (tmm.py:219-231).
    let mut keep: Vec<bool> = (0..n).map(|i| keep_m[i] && keep_a[i]).collect();
    if keep.iter().filter(|k| **k).count() < 5 {
        keep = keep_m;
    }
    if keep.iter().filter(|k| **k).count() < 5 {
        return 1.0;
    }

    // Weighted mean of M over the kept features; if the kept weights sum to
    // zero fall back to the unweighted mean (tmm.py:233-243).
    let mut weighted_sum = 0.0f64;
    let mut weight_sum = 0.0f64;
    let mut m_sum = 0.0f64;
    let mut kept_count = 0usize;
    for i in 0..n {
        if keep[i] {
            weighted_sum += weights[i] * m_values[i];
            weight_sum += weights[i];
            m_sum += m_values[i];
            kept_count += 1;
        }
    }
    let tmm_log = if weight_sum > 0.0 {
        weighted_sum / weight_sum
    } else {
        m_sum / kept_count as f64
    };

    // Convert from log2 back to a linear factor (tmm.py:243).
    2.0f64.powf(tmm_log)
}

/// Ascending argsort matching `np.argsort` default: a stable sort that breaks
/// ties by original index. Sorts values with `f64::total_cmp` so ordering is
/// deterministic; `np.argsort` on distinct keys is deterministic and this
/// stable tie-break reproduces its behaviour when keys collide.
fn stable_argsort(values: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..values.len()).collect();
    indices.sort_by(|&a, &b| values[a].total_cmp(&values[b]).then(a.cmp(&b)));
    indices
}

/// 75th-percentile-style linear interpolation matching `np.percentile(vals, q)`
/// with the default `interpolation="linear"`. `values` is sorted in place with
/// `f64::total_cmp`; `q` is a percentile in `[0, 100]`.
fn percentile_linear(values: &mut [f64], q: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(f64::total_cmp);
    if values.len() == 1 {
        return values[0];
    }
    // NumPy maps the percentile to a fractional rank over [0, n-1].
    let rank = (q / 100.0) * (values.len() - 1) as f64;
    let lower = rank.floor() as usize;
    let upper = rank.ceil() as usize;
    if lower == upper {
        return values[lower];
    }
    let fraction = rank - lower as f64;
    values[lower] + (values[upper] - values[lower]) * fraction
}

#[cfg(test)]
mod tmm_tests {
    use super::*;

    /// Golden factors computed with the pure-Python `TMMNormalizer` on the same
    /// tiny matrix (see the crate test module's docstring for the recipe). The
    /// matrix has 12 features across 3 samples so the >=10 valid-feature guard
    /// is satisfied; `Sample1` is the reference (UQ closest to the mean UQ).
    #[test]
    fn tmm_factors_match_python() {
        let columns = vec![
            "Sample1".to_string(),
            "Sample2".to_string(),
            "Sample3".to_string(),
        ];
        // Row-major features; transposed into column-major below.
        let features: [[f64; 3]; 12] = [
            [100.0, 210.0, 95.0],
            [200.0, 390.0, 205.0],
            [300.0, 610.0, 295.0],
            [400.0, 820.0, 390.0],
            [500.0, 1010.0, 505.0],
            [600.0, 1180.0, 610.0],
            [700.0, 1420.0, 690.0],
            [800.0, 1590.0, 810.0],
            [900.0, 1810.0, 895.0],
            [1000.0, 2020.0, 990.0],
            [1100.0, 2180.0, 1110.0],
            [1200.0, 2410.0, 1190.0],
        ];
        let matrix: Vec<Vec<f64>> = (0..3)
            .map(|sample| features.iter().map(|row| row[sample]).collect())
            .collect();

        let factors = tmm_norm_factors(&columns, &matrix);
        // Expected values frozen from the Python reference implementation.
        let expected = [
            ("Sample1", 0.999_723_651_114_078_7),
            ("Sample2", 1.002_352_970_935_732_2),
            ("Sample3", 0.997_928_328_921_841_2),
        ];
        for (name, want) in expected {
            let got = factors[name];
            assert!(
                (got - want).abs() < 1e-9,
                "sample {name}: got {got}, want {want}"
            );
        }
    }
}
