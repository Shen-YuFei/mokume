//! Shared numeric primitives reused across quantification, normalization, and
//! imputation crates. These are intentionally dependency-free so any crate that
//! depends on `mokume-core` can reuse them without cross-crate coupling.

/// Linear-interpolated quantile of `values`. Sorts `values` in place. Returns
/// `None` for an empty slice or an out-of-range quantile.
pub fn quantile_linear(values: &mut [f64], quantile: f64) -> Option<f64> {
    if values.is_empty() || !quantile.is_finite() || !(0.0..=1.0).contains(&quantile) {
        return None;
    }
    values.sort_by(f64::total_cmp);
    let position = (values.len() - 1) as f64 * quantile;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        return values.get(lower).copied();
    }
    let lower_value = *values.get(lower)?;
    let upper_value = *values.get(upper)?;
    Some(lower_value + (upper_value - lower_value) * (position - lower as f64))
}

/// Median over the finite, strictly positive values of `values` (others are
/// dropped). Mutates `values`.
pub fn median(values: &mut Vec<f64>) -> Option<f64> {
    values.retain(|value| value.is_finite() && *value > 0.0);
    median_sorted(values)
}

/// Median over the finite values of `values` (non-finite are dropped). Mutates
/// `values`.
pub fn median_finite(values: &mut Vec<f64>) -> Option<f64> {
    values.retain(|value| value.is_finite());
    median_sorted(values)
}

/// Mean over the finite, strictly positive values of `values`.
pub fn mean_positive(values: &[f64]) -> Option<f64> {
    let values = values
        .iter()
        .copied()
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect::<Vec<_>>();
    (!values.is_empty()).then(|| values.iter().sum::<f64>() / values.len() as f64)
}

/// Mean over the finite values of `values`.
pub fn mean_finite(values: &[f64]) -> Option<f64> {
    let values = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    (!values.is_empty()).then(|| values.iter().sum::<f64>() / values.len() as f64)
}

/// Population standard deviation over the finite values of `values`.
pub fn finite_sd(values: &[f64]) -> Option<f64> {
    let values = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect::<Vec<_>>();
    if values.is_empty() {
        return None;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / values.len() as f64;
    Some(variance.sqrt())
}

/// LOWESS smoother (Cleveland's locally-weighted scatterplot smoothing),
/// mirroring `statsmodels.nonparametric.smoothers_lowess.lowess(y, x, frac,
/// it=iterations, delta=0)` evaluated at the input points (`return_sorted=False`
/// order is the caller's responsibility — `x` must be pre-sorted ascending).
///
/// For each point a window of the `r = ceil(frac * n)` nearest neighbours (by
/// `x`) is selected by sliding a contiguous interval, weighted by the tricube
/// kernel `w_j = (1 - (|x_j - x_i| / h)^3)^3` (with `h` the max in-window
/// distance) times the current robustness weights. A weighted least-squares line
/// is fit and evaluated at `x_i`. After each of `iterations` passes the
/// robustness weights are recomputed from the residuals via the bisquare
/// function `(1 - (e / (6 * mad))^2)^2`, clamped to `[0, 1]`.
pub fn lowess_fit(x: &[f64], y: &[f64], frac: f64, iterations: usize) -> Vec<f64> {
    let n = x.len();
    let mut fitted = vec![0.0; n];
    if n == 0 {
        return fitted;
    }
    if n == 1 {
        fitted[0] = y[0];
        return fitted;
    }
    let window = ((frac * n as f64).ceil() as usize).clamp(2, n);
    let mut robustness = vec![1.0_f64; n];

    for iteration in 0..=iterations {
        let mut left = 0usize;
        for i in 0..n {
            // Slide the width-`window` interval right while doing so brings its
            // far edge closer to x[i] than its near edge (Cleveland's update).
            left = left.min(n - window);
            while left + window < n && (x[i] - x[left]) > (x[left + window - 1] - x[i]) {
                left += 1;
            }
            let right = left + window;
            let span = (x[i] - x[left]).max(x[right - 1] - x[i]);

            let (mut sw, mut swx, mut swy, mut swxx, mut swxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for j in left..right {
                let weight = if span > 0.0 {
                    let d = (x[j] - x[i]).abs() / span;
                    if d >= 1.0 {
                        0.0
                    } else {
                        let tricube = 1.0 - d * d * d;
                        tricube * tricube * tricube
                    }
                } else {
                    1.0
                } * robustness[j];
                sw += weight;
                swx += weight * x[j];
                swy += weight * y[j];
                swxx += weight * x[j] * x[j];
                swxy += weight * x[j] * y[j];
            }
            if sw <= 0.0 {
                fitted[i] = y[i];
                continue;
            }
            let mean_x = swx / sw;
            let mean_y = swy / sw;
            let var_x = swxx / sw - mean_x * mean_x;
            let cov_xy = swxy / sw - mean_x * mean_y;
            let slope = if var_x > 0.0 { cov_xy / var_x } else { 0.0 };
            let intercept = mean_y - slope * mean_x;
            fitted[i] = intercept + slope * x[i];
        }

        if iteration < iterations {
            let mut abs_residuals = (0..n).map(|i| (y[i] - fitted[i]).abs()).collect::<Vec<_>>();
            match median_finite(&mut abs_residuals) {
                Some(mad) if mad > 0.0 => {
                    for (rw, i) in robustness.iter_mut().zip(0..n) {
                        let u = ((y[i] - fitted[i]) / (6.0 * mad)).clamp(-1.0, 1.0);
                        let bisquare = 1.0 - u * u;
                        *rw = bisquare * bisquare;
                    }
                }
                _ => robustness.iter_mut().for_each(|rw| *rw = 1.0),
            }
        }
    }
    fitted
}

/// Median of an already-filtered slice; sorts in place. No value filtering.
pub fn median_sorted(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(f64::total_cmp);
    let midpoint = values.len() / 2;
    if values.len().is_multiple_of(2) {
        Some((values[midpoint - 1] + values[midpoint]) / 2.0)
    } else {
        Some(values[midpoint])
    }
}
