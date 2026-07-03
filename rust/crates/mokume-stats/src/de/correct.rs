//! Multiple-testing correction for differential expression.
//!
//! [`bh_adjust`] reproduces `statsmodels.stats.multitest.multipletests(method =
//! "fdr_bh")` as used by mokume's `_finalize_results`: the Benjamini-Hochberg
//! step-up procedure over the finite p-values, with non-finite entries passed
//! through as `NaN`.

/// Benjamini-Hochberg adjusted p-values.
///
/// Finite p-values are ranked ascending, scaled by `n / rank`, made monotone by
/// a reverse cumulative minimum, and clipped to `1.0`. `NaN` inputs stay `NaN`
/// and are excluded from the count `n`, matching how mokume drops non-finite
/// p-values before correcting.
pub(crate) fn bh_adjust(pvalues: &[f64]) -> Vec<f64> {
    let mut finite = pvalues
        .iter()
        .enumerate()
        .filter(|(_, value)| value.is_finite())
        .map(|(index, value)| (index, *value))
        .collect::<Vec<_>>();
    let n = finite.len();
    let mut adjusted = vec![f64::NAN; pvalues.len()];
    if n == 0 {
        return adjusted;
    }

    finite.sort_by(|(_, a), (_, b)| a.total_cmp(b));
    let scaled = finite
        .iter()
        .enumerate()
        .map(|(rank, (_, value))| value * n as f64 / (rank + 1) as f64)
        .collect::<Vec<_>>();

    // Reverse cumulative minimum, then clip to 1.0.
    let mut running_min = f64::INFINITY;
    let mut monotone = vec![0.0; n];
    for index in (0..n).rev() {
        running_min = running_min.min(scaled[index]);
        monotone[index] = running_min.min(1.0);
    }

    for (rank, (original_index, _)) in finite.iter().enumerate() {
        adjusted[*original_index] = monotone[rank];
    }
    adjusted
}

#[cfg(test)]
mod tests {
    use super::bh_adjust;

    fn assert_close(actual: f64, expected: f64) {
        let tol = 1e-9 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tol,
            "actual={actual} expected={expected}"
        );
    }

    // Reference from statsmodels fdr_bh on the limma oracle p-values (6 tests,
    // including a tie at 3.1734e-05). Adjusted values verified against
    // `multipletests(p, method="fdr_bh")`.
    #[test]
    fn matches_statsmodels_fdr_bh() {
        let p = [
            1.712_569_861_574_927e-7,
            3.1734308960086415e-05,
            3.1734308960086415e-05,
            0.544_763_715_713_143_1,
            1.0,
            0.941_320_895_121_956_2,
        ];
        let adjusted = bh_adjust(&p);
        assert_close(adjusted[0], 1.0275419169449564e-06);
        assert_close(adjusted[1], 6.346861792017283e-05);
        assert_close(adjusted[2], 6.346861792017283e-05);
        assert_close(adjusted[3], 0.817_145_573_569_714_6);
        assert_close(adjusted[4], 1.0);
        assert_close(adjusted[5], 1.0);
    }

    #[test]
    fn passes_through_non_finite() {
        let p = [0.01, f64::NAN, 0.02];
        let adjusted = bh_adjust(&p);
        // n = 2 finite values: 0.01 -> 0.02 (0.01*2/1), 0.02 -> 0.02 (0.02*2/2).
        assert_close(adjusted[0], 0.02);
        assert!(adjusted[1].is_nan());
        assert_close(adjusted[2], 0.02);
    }
}
