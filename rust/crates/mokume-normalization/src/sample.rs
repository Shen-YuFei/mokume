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
