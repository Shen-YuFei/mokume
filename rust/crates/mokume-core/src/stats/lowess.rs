//! LOWESS following statsmodels 0.14.6 `_smoothers_lowess.pyx`.
//!
//! Supports evaluation at the input points with `delta=0`. Missing-value
//! filtering and sorting belong to the caller, as with the upstream Cython core.

use super::median_sorted;

/// Smooth finite `y` at sorted, finite `x`, using local linear regression and
/// `iterations` residual reweightings. The inputs must have the same length and
/// `frac` must lie in `[0, 1]`.
pub fn lowess_fit(x: &[f64], y: &[f64], frac: f64, iterations: usize) -> Vec<f64> {
    assert_eq!(x.len(), y.len(), "LOWESS x/y lengths differ");
    assert!((0.0..=1.0).contains(&frac), "LOWESS frac must be in [0, 1]");
    let n = x.len();
    if n < 2 {
        return y.to_vec();
    }
    let window = ((frac * n as f64 + 1e-10) as usize).clamp(2, n);
    let mut fitted = vec![0.0; n];
    let mut robustness = vec![1.0; n];
    let mut weights = vec![0.0; window];
    for iteration in 0..=iterations {
        let mut left = 0;
        for i in 0..n {
            // Upstream reuses the fit for tied x values, even if regression
            // failed and the first tied point's original y was returned.
            if i > 0 && x[i] == x[i - 1] {
                fitted[i] = fitted[i - 1];
                continue;
            }
            while left + window < n && x[i] > (x[left] + x[left + window]) / 2.0 {
                left += 1;
            }
            let right = left + window;
            let radius = (x[i] - x[left]).max(x[right - 1] - x[i]);
            fitted[i] = local_fit(
                &x[left..right],
                &y[left..right],
                &robustness[left..right],
                &mut weights,
                x[i],
                radius,
                y[i],
            );
        }
        if iteration < iterations {
            update_robustness(y, &fitted, &mut robustness);
        }
    }
    fitted
}

fn local_fit(
    x: &[f64],
    y: &[f64],
    robustness: &[f64],
    weights: &mut [f64],
    xval: f64,
    radius: f64,
    original: f64,
) -> f64 {
    let mut nonzero = 0;
    let mut sum = 0.0;
    for ((weight, &xj), &robust) in weights.iter_mut().zip(x).zip(robustness) {
        let distance = (xj - xval).abs() / radius;
        let tricube = 1.0 - distance * distance * distance;
        *weight = tricube * tricube * tricube * robust;
        nonzero += usize::from(*weight > 1e-12);
        sum += *weight;
    }
    if nonzero < 2 {
        return original;
    }
    let mut center = 0.0;
    for (weight, &xj) in weights.iter_mut().zip(x) {
        *weight /= sum;
        center += *weight * xj;
    }
    let variance = weights
        .iter()
        .zip(x)
        .map(|(&w, &xj)| w * (xj - center).powi(2))
        .sum::<f64>()
        .max(1e-12);
    weights
        .iter()
        .zip(x)
        .zip(y)
        .map(|((&w, &xj), &yj)| w * (1.0 + (xval - center) * (xj - center) / variance) * yj)
        .sum()
}

fn update_robustness(y: &[f64], fitted: &[f64], robustness: &mut [f64]) {
    let mut residuals = y
        .iter()
        .zip(fitted)
        .map(|(y, f)| (y - f).abs())
        .collect::<Vec<_>>();
    let median = median_sorted(&mut residuals).unwrap_or(0.0);
    for ((weight, &yj), &fit) in robustness.iter_mut().zip(y).zip(fitted) {
        let residual = (yj - fit).abs();
        let scaled = if median == 0.0 {
            f64::from(residual > 0.0)
        } else {
            (residual / (6.0 * median)).min(1.0)
        };
        *weight = (1.0 - scaled * scaled).powi(2);
    }
}

#[cfg(test)]
mod tests {
    use super::lowess_fit;
    use std::collections::BTreeMap;

    #[test]
    fn matches_statsmodels_0_14_6() -> Result<(), Box<dyn std::error::Error>> {
        let mut cases = BTreeMap::<_, Vec<_>>::new();
        for line in include_str!("../../tests/data/lowess_statsmodels.tsv").lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            cases
                .entry((fields[0], fields[1], fields[2]))
                .or_default()
                .push((
                    fields[3].parse::<f64>()?,
                    fields[4].parse::<f64>()?,
                    fields[5].parse::<f64>()?,
                ));
        }
        for ((name, frac, iterations), rows) in cases {
            let x = rows.iter().map(|r| r.0).collect::<Vec<_>>();
            let y = rows.iter().map(|r| r.1).collect::<Vec<_>>();
            let actual = lowess_fit(&x, &y, frac.parse()?, iterations.parse()?);
            for (i, (&value, row)) in actual.iter().zip(&rows).enumerate() {
                let tolerance = 1e-10 + 1e-9 * row.2.abs();
                assert!(
                    (value - row.2).abs() <= tolerance,
                    "{name}[{i}]: {value} != {}",
                    row.2
                );
            }
        }
        Ok(())
    }
}
