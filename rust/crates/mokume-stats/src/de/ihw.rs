//! Independent Hypothesis Weighting (IHW) multiple-testing correction.
//!
//! Port of mokume's `_ihw_correction` and its helpers in
//! `mokume/analysis/differential_expression.py` (the IHW block,
//! `differential_expression.py:245-420`). IHW (Ignatiadis et al. 2016) groups
//! the tested p-values into covariate bins, optimises a per-bin weight, then
//! applies weighted Benjamini-Hochberg over the finite-p subset. mokume's
//! variant is fully deterministic, so this port is cell-exact (1e-9) versus the
//! Python oracle.
//!
//! The flow mirrors `_ihw_correction` (`differential_expression.py:389`):
//!   1. drop entries whose p-value or covariate is non-finite -- they get `NaN`
//!      adjusted p-values, matching `_ihw_adjust_valid`'s
//!      `valid = isfinite(p) & isfinite(covariate)`
//!      (`differential_expression.py:376`);
//!   2. if fewer than `2 * n_bins` valid entries remain, fall back to plain BH
//!      over the finite p-values (`differential_expression.py:377-378, 405-406`);
//!   3. bin the valid covariates with `pd.qcut(q=n_bins, duplicates="drop")`
//!      (`differential_expression.py:380`), reproduced cell-exactly in
//!      [`qcut_labels`];
//!   4. optimise per-bin weights ([`optimize_weights`],
//!      `_ihw_optimize_weights`, `differential_expression.py:336`);
//!   5. divide each valid p-value by its bin weight (clipped to >= 0.01), clip
//!      the result to [0, 1], and BH-adjust
//!      (`differential_expression.py:383-386`).
//!
//! `bh_adjust` here is the same finite-only, NaN-passthrough step-up procedure
//! mokume uses everywhere else (`statsmodels` `fdr_bh`), shared from [`super`].

use super::correct::bh_adjust;

/// Number of iterations of the weight-optimisation loop, matching Python's
/// `for _ in range(10)` (`differential_expression.py:345`).
const MAX_ITERATIONS: usize = 10;

/// Lower clip on a bin weight before dividing a p-value, matching Python's
/// `np.clip(..., 0.01, None)` (`differential_expression.py:347, 384`).
const MIN_WEIGHT: f64 = 0.01;

/// Lower/upper clip on the re-estimated per-bin weight, matching Python's
/// `np.clip(rej / total * n_bins, 0.1, 10.0)` (`differential_expression.py:356`).
const WEIGHT_LOWER: f64 = 0.1;
const WEIGHT_UPPER: f64 = 10.0;

/// Convergence tolerances for the weight loop's early stop, matching
/// `np.allclose(weights, nw, atol=0.01)` (default `rtol=1e-5`)
/// (`differential_expression.py:358`).
const ALLCLOSE_ATOL: f64 = 0.01;
const ALLCLOSE_RTOL: f64 = 1e-5;

/// IHW-adjusted p-values aligned to `pvalues`, given a per-entry `covariate`.
///
/// Cell-exact port of `_ihw_correction` (`differential_expression.py:389`) with
/// `alpha` = the FDR threshold and `n_bins` = 5 (mokume's default). Entries with
/// a non-finite p-value or covariate receive `NaN`; if fewer than `2 * n_bins`
/// entries are jointly finite, the result is plain BH over the finite p-values
/// (the Python fallback). `pvalues` and `covariate` must be the same length;
/// a length mismatch yields all-`NaN` defensively rather than panicking.
pub fn ihw_correction(pvalues: &[f64], covariate: &[f64], alpha: f64, n_bins: usize) -> Vec<f64> {
    let mut adjusted = vec![f64::NAN; pvalues.len()];
    if pvalues.is_empty() || pvalues.len() != covariate.len() {
        return adjusted;
    }

    // `valid = isfinite(p) & isfinite(covariate)` (differential_expression.py:376).
    let valid: Vec<usize> = (0..pvalues.len())
        .filter(|&i| pvalues[i].is_finite() && covariate[i].is_finite())
        .collect();

    // Too few valid points: fall back to plain BH over the finite p-values, the
    // `_ihw_adjust_valid` -> None branch (differential_expression.py:377, 405).
    if valid.len() < n_bins.saturating_mul(2) {
        return bh_adjust(pvalues);
    }

    let valid_p: Vec<f64> = valid.iter().map(|&i| pvalues[i]).collect();
    let valid_cov: Vec<f64> = valid.iter().map(|&i| covariate[i]).collect();

    // qcut bin label per valid entry; labels are contiguous 0..=max, so they
    // double as positions into the per-bin weight vector.
    let bins = qcut_labels(&valid_cov, n_bins);
    let n_weight_bins = bins.iter().copied().max().map_or(0, |m| m + 1);
    let weights = optimize_weights(&valid_p, &bins, n_weight_bins, alpha);

    // weighted_p = clip(p / clip(weight, 0.01, None), 0, 1)
    // (differential_expression.py:383-385).
    let weighted_p: Vec<f64> = valid_p
        .iter()
        .zip(&bins)
        .map(|(&p, &bin)| {
            let weight = weights.get(bin).copied().unwrap_or(1.0).max(MIN_WEIGHT);
            (p / weight).clamp(0.0, 1.0)
        })
        .collect();

    let adj_valid = bh_adjust(&weighted_p);
    for (slot, &index) in adj_valid.iter().zip(&valid) {
        if let Some(target) = adjusted.get_mut(index) {
            *target = *slot;
        }
    }
    adjusted
}

/// Per-bin weights from `_ihw_optimize_weights` (`differential_expression.py:336`).
///
/// Up to [`MAX_ITERATIONS`] rounds: divide each p-value by its current bin
/// weight (clipped to >= [`MIN_WEIGHT`]), count BH rejections at `alpha` within
/// each bin, set the new weight proportional to a bin's rejection share scaled
/// by the bin count (clipped to `[0.1, 10.0]`) and renormalised to mean 1, and
/// stop early once the weights stop moving (`np.allclose(atol=0.01)`).
fn optimize_weights(pvalues: &[f64], bins: &[usize], n_bins: usize, alpha: f64) -> Vec<f64> {
    let mut weights = vec![1.0; n_bins];
    if n_bins == 0 {
        return weights;
    }
    for _ in 0..MAX_ITERATIONS {
        // w_p = p / clip(weight_of_bin, 0.01, None) (differential_expression.py:346-347).
        let weighted_p: Vec<f64> = pvalues
            .iter()
            .zip(bins)
            .map(|(&p, &bin)| p / weights.get(bin).copied().unwrap_or(1.0).max(MIN_WEIGHT))
            .collect();

        // rej[bin] = (BH-adjusted within-bin p < alpha).sum()
        // (differential_expression.py:348-353).
        let mut rejections = vec![0.0; n_bins];
        for (bin, rejection) in rejections.iter_mut().enumerate() {
            let bin_p: Vec<f64> = weighted_p
                .iter()
                .zip(bins)
                .filter_map(|(&p, &b)| (b == bin).then_some(p))
                .collect();
            let adjusted = bh_adjust(&bin_p);
            *rejection = adjusted.iter().filter(|&&a| a < alpha).count() as f64;
        }

        let total: f64 = rejections.iter().sum();
        if total <= 0.0 {
            // Python only updates weights when total > 0; otherwise the loop just
            // re-runs with unchanged weights, so it is equivalent to stopping.
            break;
        }

        // nw = clip(rej / total * n_bins, 0.1, 10.0); nw /= nw.mean()
        // (differential_expression.py:356-357).
        let scale = n_bins as f64 / total;
        let mut next: Vec<f64> = rejections
            .iter()
            .map(|&r| (r * scale).clamp(WEIGHT_LOWER, WEIGHT_UPPER))
            .collect();
        let mean = next.iter().sum::<f64>() / n_bins as f64;
        if mean > 0.0 {
            for weight in &mut next {
                *weight /= mean;
            }
        }

        if allclose(&weights, &next) {
            break;
        }
        weights = next;
    }
    weights
}

/// numpy `np.allclose(current, next, atol=0.01)` (default `rtol=1e-5`):
/// `|current - next| <= atol + rtol * |next|` for every element.
fn allclose(current: &[f64], next: &[f64]) -> bool {
    current
        .iter()
        .zip(next)
        .all(|(&a, &b)| (a - b).abs() <= ALLCLOSE_ATOL + ALLCLOSE_RTOL * b.abs())
}

/// Cell-exact port of `pd.qcut(x, q=n_bins, labels=False, duplicates="drop")`
/// (`differential_expression.py:380`).
///
/// Computes the `n_bins + 1` quantile edges at `linspace(0, 1, n_bins + 1)` via
/// numpy's default linear interpolation, drops consecutive duplicate edges
/// (`duplicates="drop"`), nudges the first edge down by `0.001 * |range|` (so
/// the minimum lands in bin 0, as pandas does), then labels each value by the
/// number of edges strictly below it (right-closed intervals), clamped to the
/// last bin. Returns one contiguous label in `0..n_unique_edges-1` per value.
fn qcut_labels(values: &[f64], n_bins: usize) -> Vec<usize> {
    if values.is_empty() || n_bins == 0 {
        return vec![0; values.len()];
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);

    // Quantile edges with numpy "linear" interpolation, then drop consecutive
    // duplicates (the edges are sorted, so consecutive-dedup == pd.unique here).
    let mut edges: Vec<f64> = Vec::with_capacity(n_bins + 1);
    for k in 0..=n_bins {
        let quantile = k as f64 / n_bins as f64;
        let edge = linear_quantile(&sorted, quantile);
        if edges.last().is_none_or(|&last| last != edge) {
            edges.push(edge);
        }
    }

    // Nudge the first edge down so the minimum value falls into bin 0, matching
    // pandas' `bins[0] -= 0.001 * (max - min)` (or 0.001 for a zero range).
    let span = edges.last().copied().unwrap_or(0.0) - edges[0];
    let nudge = if span != 0.0 {
        0.001 * span.abs()
    } else {
        0.001
    };
    if let Some(first) = edges.first_mut() {
        *first -= nudge;
    }

    let last_bin = edges.len().saturating_sub(2);
    values
        .iter()
        .map(|&value| {
            // Count of edges strictly less than `value` (searchsorted side="left")
            // minus one, clamped to the last bin -- right-closed intervals.
            let below = edges.iter().filter(|&&edge| edge < value).count();
            below.saturating_sub(1).min(last_bin)
        })
        .collect()
}

/// numpy `np.quantile(sorted, q, method="linear")` for a pre-sorted slice:
/// position `h = (n - 1) * q`, then linear-interpolate between the two
/// surrounding order statistics.
fn linear_quantile(sorted: &[f64], quantile: f64) -> f64 {
    match sorted.len() {
        0 => f64::NAN,
        1 => sorted[0],
        n => {
            let position = (n - 1) as f64 * quantile;
            let lower = position.floor() as usize;
            let upper = (lower + 1).min(n - 1);
            let frac = position - lower as f64;
            sorted[lower] + frac * (sorted[upper] - sorted[lower])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{ihw_correction, qcut_labels};
    use crate::de::correct::bh_adjust;

    fn assert_close(actual: f64, expected: f64, label: &str) {
        let tol = 1e-9 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tol,
            "{label}: actual={actual} expected={expected}"
        );
    }

    // qcut labels must match pd.qcut(q=5, duplicates="drop") exactly, including
    // the duplicate-edge merge that collapses a 1,1,1,1 prefix into one bin.
    #[test]
    fn qcut_labels_match_pandas() {
        assert_eq!(
            qcut_labels(
                &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
                5
            ),
            vec![0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
        );
        // Duplicated low covariates collapse the first bins (duplicates="drop").
        assert_eq!(
            qcut_labels(
                &[1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 7.0, 7.0, 9.0, 9.0, 11.0, 12.0],
                5
            ),
            vec![0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3]
        );
    }

    // Golden oracle: mokume's `_ihw_correction(p, de_df, alpha=0.05)` with a
    // mean covariate, 5 bins. Captured cell-exact from
    //   conda run -n Bigbio python  # against the IHW block of
    //   mokume.analysis.differential_expression
    // on 12 (p, covariate) pairs (covariate = row-mean of the two mean_ columns,
    // here just the distinct covariate). IHW differs from plain BH on the same
    // p-values (the non-vacuity check below), proving the weighting is active.
    const CASE1_P: [f64; 12] = [
        0.001, 0.002, 0.2, 0.3, 0.0005, 0.04, 0.5, 0.6, 0.008, 0.9, 0.012, 0.45,
    ];
    const CASE1_COV: [f64; 12] = [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    ];
    const CASE1_IHW: [f64; 12] = [
        0.003_06,
        0.004_08,
        0.204,
        0.524_571_428_571_428_6,
        0.003_06,
        0.612,
        1.0,
        0.734_399_999_999_999_9,
        0.024_48,
        1.0,
        0.029_376,
        0.612,
    ];

    #[test]
    fn ihw_matches_python_oracle_case1() {
        let adjusted = ihw_correction(&CASE1_P, &CASE1_COV, 0.05, 5);
        for (index, &expected) in CASE1_IHW.iter().enumerate() {
            assert_close(adjusted[index], expected, &format!("case1[{index}]"));
        }
    }

    // Non-vacuity: IHW must genuinely differ from plain BH on the same inputs.
    // On case1, IHW down-weights bin 1 hits and up-weights others, so several
    // cells move materially (e.g. index 5: BH=0.08 -> IHW=0.612).
    #[test]
    fn ihw_differs_from_plain_bh() {
        let ihw = ihw_correction(&CASE1_P, &CASE1_COV, 0.05, 5);
        let bh = bh_adjust(&CASE1_P);
        let max_delta = ihw
            .iter()
            .zip(&bh)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        // Largest cell delta (index 5): |0.612 - 0.08| = 0.532.
        assert!(
            max_delta > 0.5,
            "IHW must differ from BH (max delta {max_delta})"
        );
        assert_close(ihw[5], 0.612, "ihw[5]");
        assert_close(bh[5], 0.08, "bh[5]");
    }

    // Edge case: a NaN p-value passes through as NaN; a NaN covariate excludes
    // that entry from binning (NaN adjusted p), and the remaining 10 valid
    // points (>= 2 * 5) still run full IHW. Oracle captured the same way.
    const CASE2_P: [f64; 12] = [
        0.001,
        0.002,
        f64::NAN,
        0.3,
        0.0005,
        0.04,
        0.5,
        0.6,
        0.008,
        0.9,
        0.012,
        0.45,
    ];
    const CASE2_COV: [f64; 12] = [
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
        f64::NAN,
        7.0,
        8.0,
        9.0,
        10.0,
        11.0,
        12.0,
    ];
    const CASE2_IHW: [f64; 12] = [
        0.002_55,
        0.003_4,
        f64::NAN,
        0.51,
        0.002_55,
        f64::NAN,
        1.0,
        1.0,
        0.020_4,
        1.0,
        0.024_48,
        0.655_714_285_714_285_8,
    ];

    #[test]
    fn ihw_handles_nan_p_and_nan_covariate() {
        let adjusted = ihw_correction(&CASE2_P, &CASE2_COV, 0.05, 5);
        for (index, &expected) in CASE2_IHW.iter().enumerate() {
            if expected.is_nan() {
                assert!(adjusted[index].is_nan(), "case2[{index}] should be NaN");
            } else {
                assert_close(adjusted[index], expected, &format!("case2[{index}]"));
            }
        }
    }

    // Edge case: fewer than 2 * n_bins valid points -> fall back to plain BH
    // (the `_ihw_adjust_valid` -> None branch). With 3 points and 5 bins the
    // result equals `bh_adjust` exactly (a single "bin").
    #[test]
    fn ihw_falls_back_to_bh_when_too_few_points() {
        let pvalues = [0.01, 0.2, 0.3];
        let covariate = [1.0, 2.0, 3.0];
        let adjusted = ihw_correction(&pvalues, &covariate, 0.05, 5);
        let bh = bh_adjust(&pvalues);
        for (index, (&got, &expected)) in adjusted.iter().zip(&bh).enumerate() {
            assert_close(got, expected, &format!("fallback[{index}]"));
        }
        // Oracle values for this fallback case.
        assert_close(adjusted[0], 0.03, "fallback[0]");
        assert_close(adjusted[1], 0.3, "fallback[1]");
        assert_close(adjusted[2], 0.3, "fallback[2]");
    }

    // Duplicate-covariate case forcing duplicates="drop" to merge edges, so the
    // valid set lands in only 4 bins. Confirms the merged-bin weighting path is
    // cell-exact too.
    const CASE4_P: [f64; 12] = [
        0.001, 0.002, 0.2, 0.3, 0.0005, 0.04, 0.5, 0.6, 0.008, 0.9, 0.012, 0.45,
    ];
    const CASE4_COV: [f64; 12] = [1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 7.0, 7.0, 9.0, 9.0, 11.0, 12.0];
    const CASE4_IHW: [f64; 12] = [
        0.002_306_25,
        0.003_075,
        0.131_785_714_285_714_3,
        0.172_968_75,
        0.002_306_25,
        0.036_9,
        1.0,
        1.0,
        0.036_9,
        1.0,
        0.036_9,
        0.922_5,
    ];

    #[test]
    fn ihw_matches_python_oracle_merged_bins() {
        let adjusted = ihw_correction(&CASE4_P, &CASE4_COV, 0.05, 5);
        for (index, &expected) in CASE4_IHW.iter().enumerate() {
            assert_close(adjusted[index], expected, &format!("case4[{index}]"));
        }
    }
}
