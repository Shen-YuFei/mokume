//! DEqMS count-aware moderated t-test for two-group differential expression.
//!
//! This is a port of mokume's `analysis/deqms.py` (itself a pure-Python
//! reimplementation of R `DEqMS`): the limma pipeline (`lmFit` ->
//! `contrasts.fit` -> `eBayes`, reused verbatim from [`super::limma`]) followed
//! by `spectraCounteBayes` count-aware variance moderation (Zhu et al. 2020).
//! The flow per protein mirrors `_spectra_count_ebayes` (deqms.py:43)
//! operation-by-operation:
//!   - `sigma2 = sigma**2`, `log_var = log(sigma2)`, `df_valid = df` with 0
//!     replaced by NaN (deqms.py:51-57);
//!   - `x = log2(counts)` (after `+1` if any count is 0) (deqms.py:59-63);
//!   - the `valid.sum() < 10` and `ptp(x_valid) < 1e-10` FALLBACK to standard
//!     eBayes (deqms.py:65-80), faithful to Python's default-counts behaviour;
//!   - a statsmodels-equivalent LOWESS (`frac=0.75`, `it=3`) of `log_var` on
//!     `x` over the valid points, then linear interpolation onto any unseen
//!     non-valid `x` (deqms.py:82-101) via [`mokume_core::stats::lowess_fit`];
//!   - the `eg`/`egpred`/`myfct`/`mean_myfct` count-trend statistic
//!     (deqms.py:103-107);
//!   - the `_grid_search_d0` prior-df search with its early break
//!     (deqms.py:122-139);
//!   - the moderated `post_var`, `sca_t`, and two-sided `sca_pvalue`
//!     (deqms.py:111-117).
//!
//! TIERED TOLERANCE: this method is NOT cell-exact. The LOWESS smoother matches
//! statsmodels only to ~2e-3 (a window-size and `delta` difference), and that
//! error propagates through `mean_myfct`, the `d0` grid search, and `post_var`
//! into `sca_t`/`sca_pvalue`. The oracle test asserts a relative tolerance of
//! ~1e-2 on `sca_t`/`sca_pvalue` plus an exact protein RANK match; `log2FC` is
//! cell-exact (it is the limma coefficient, untouched by the count moderation).
//! The all-equal-counts fallback test is cell-exact (1e-9) because it bypasses
//! LOWESS entirely and reduces to plain limma `eBayes`.
//!
//! Reference
//! ---------
//! Zhu Y, Orre LM, Zhou Tran Y, et al. DEqMS: a method for accurate variance
//! estimation in differential protein expression analysis. *Mol Cell
//! Proteomics*. 2020;19(6):1047-1057.

use mokume_core::stats::lowess_fit;

use super::correct::bh_adjust;
use super::limma::{run_limma_moderation, LimmaModeration};
use super::special::{digamma, trigamma};
use super::student_t::{student_t_sf, student_t_two_sided_log_pvalue};
use super::{classify, DeResult};

/// Minimum number of valid (finite `log_var`, `x`, `df`) proteins required for
/// the spectraCounteBayes count path; below this the method falls back to plain
/// eBayes (deqms.py:66).
const MIN_VALID_POINTS: usize = 10;

/// Minimum spread of `log2(counts)` over the valid points for the count path;
/// below this all counts are effectively identical and the method falls back to
/// plain eBayes (deqms.py:76).
const MIN_COUNT_SPREAD: f64 = 1e-10;

/// LOWESS span fraction (statsmodels `frac=0.75`) used by spectraCounteBayes.
const LOWESS_FRAC: f64 = 0.75;

/// LOWESS robustifying iterations (statsmodels default `it=3`).
const LOWESS_ITERATIONS: usize = 3;

/// Two-group DEqMS differential expression.
///
/// `proteins[i]` labels `rows[i]`; each row holds log2 intensities with the
/// first `n_a` columns in condition A and the next `n_b` in condition B.
/// `counts[i]` is the per-protein peptide/spectra count aligned to `rows[i]`
/// (Python's `_build_count_vector` default is all-ones, which triggers the
/// all-equal fallback to plain eBayes). Returns one [`DeResult`] per testable
/// protein, sorted by adjusted p-value ascending (BH), reproducing mokume's
/// `DifferentialExpression(method="deqms")` output contract: `pvalue` and
/// `adj_pvalue` carry the spectra-count-moderated `sca_pvalue` / its BH
/// adjustment, and `t_statistic` carries `sca_t`.
pub fn deqms_two_group(
    proteins: &[String],
    rows: &[&[f64]],
    n_a: usize,
    n_b: usize,
    counts: &[f64],
    fdr_threshold: f64,
    log2fc_threshold: f64,
) -> Vec<DeResult> {
    let (moderations, df_prior) = run_limma_moderation(rows, n_a, n_b);
    if moderations.is_empty() {
        return Vec::new();
    }

    // Align counts to the kept proteins by original row index (defaulting to 1,
    // matching `_build_count_vector`'s default fill).
    let kept_counts = moderations
        .iter()
        .map(|(index, _)| counts.get(*index).copied().unwrap_or(1.0))
        .collect::<Vec<_>>();

    let stats = spectra_count_ebayes(&moderations, &kept_counts, df_prior);

    let p_values = stats.iter().map(|stat| stat.sca_pvalue).collect::<Vec<_>>();
    let adjusted = bh_adjust(&p_values);

    let mut results = moderations
        .iter()
        .zip(&stats)
        .zip(adjusted)
        .filter_map(|(((index, moderation), stat), adj_p_value)| {
            let protein = proteins.get(*index)?.clone();
            let significance = classify(
                adj_p_value,
                moderation.log2_fold_change,
                fdr_threshold,
                log2fc_threshold,
            );
            Some(DeResult {
                protein,
                log2_fold_change: moderation.log2_fold_change,
                p_value: stat.sca_pvalue,
                log_p_value: stat.log_pvalue,
                adj_p_value,
                t_statistic: stat.sca_t,
                // DEqMS emits no AveExpr/B columns, so these limma-only quirks
                // are filled with NaN/0.0 and must not be written for deqms.
                ave_expr: f64::NAN,
                b: 0.0,
                mean_a: moderation.mean_a,
                mean_b: moderation.mean_b,
                n_a: moderation.count_a,
                n_b: moderation.count_b,
                significance,
            })
        })
        .collect::<Vec<_>>();

    results.sort_by(|left, right| left.adj_p_value.total_cmp(&right.adj_p_value));
    results
}

/// One protein's spectraCounteBayes output: the count-moderated t-statistic and
/// its two-sided p-value (deqms.py returns `(sca_t, sca_pvalue)` as columns).
struct ScaStat {
    sca_t: f64,
    sca_pvalue: f64,
    log_pvalue: f64,
}

/// Port of `_spectra_count_ebayes` (deqms.py:43). Applies the count-aware Bayes
/// moderation to the kept-protein limma fits; falls back to the standard eBayes
/// statistics (`t_statistic`, `p_value`) when the count path is not exercised.
fn spectra_count_ebayes(
    moderations: &[(usize, LimmaModeration)],
    counts: &[f64],
    df_prior: f64,
) -> Vec<ScaStat> {
    let n_genes = moderations.len();
    // sigma2 = sigma**2; log_var = log(sigma2); df_valid = df (0 -> NaN).
    let sigma2 = moderations
        .iter()
        .map(|(_, fit)| fit.sigma * fit.sigma)
        .collect::<Vec<_>>();
    let log_var = sigma2.iter().map(|value| value.ln()).collect::<Vec<_>>();
    let df = moderations
        .iter()
        .map(|(_, fit)| fit.df_residual)
        .collect::<Vec<_>>();
    let df_valid = df
        .iter()
        .map(|value| if *value == 0.0 { f64::NAN } else { *value })
        .collect::<Vec<_>>();

    // counts_f += 1 if any count is 0 (deqms.py:60-61); x = log2(counts_f).
    let min_count = counts
        .iter()
        .copied()
        .filter(|value| !value.is_nan())
        .fold(f64::INFINITY, f64::min);
    let counts_f = counts
        .iter()
        .map(|value| {
            if min_count == 0.0 {
                value + 1.0
            } else {
                *value
            }
        })
        .collect::<Vec<_>>();
    let x = counts_f
        .iter()
        .map(|value| value.log2())
        .collect::<Vec<_>>();

    // valid = isfinite(log_var) & isfinite(x) & isfinite(df_valid).
    let valid = (0..n_genes)
        .map(|i| log_var[i].is_finite() && x[i].is_finite() && df_valid[i].is_finite())
        .collect::<Vec<_>>();
    let valid_count = valid.iter().filter(|flag| **flag).count();
    if valid_count < MIN_VALID_POINTS {
        return fallback(moderations, df_prior);
    }

    // ptp(x_valid) < 1e-10 -> all counts identical -> fall back.
    let x_valid = (0..n_genes)
        .filter(|i| valid[*i])
        .map(|i| x[i])
        .collect::<Vec<_>>();
    let x_range = ptp(&x_valid);
    if x_range < MIN_COUNT_SPREAD {
        return fallback(moderations, df_prior);
    }

    let y_pred = lowess_predict(&log_var, &x, &valid);

    // eg / egpred / myfct / mean_myfct (deqms.py:103-107).
    let mut myfct_sum = 0.0;
    let mut myfct_count = 0usize;
    let mut egpred = vec![f64::NAN; n_genes];
    for i in 0..n_genes {
        if !df_valid[i].is_finite() || !y_pred[i].is_finite() {
            continue;
        }
        let half_df = df_valid[i] / 2.0;
        let shift = -digamma(half_df) + (half_df).ln();
        let eg_i = log_var[i] + shift;
        let egpred_i = y_pred[i] + shift;
        egpred[i] = egpred_i;
        let residual = eg_i - egpred_i;
        let myfct_i = residual * residual - trigamma(half_df);
        if myfct_i.is_finite() {
            myfct_sum += myfct_i;
            myfct_count += 1;
        }
    }
    let mean_myfct = if myfct_count > 0 {
        myfct_sum / myfct_count as f64
    } else {
        f64::NAN
    };

    let d0 = grid_search_d0(mean_myfct, n_genes);

    // s02 = exp(egpred + digamma(d0/2) - log(d0/2)); post_var; post_df; sca_t.
    let half_d0 = d0 / 2.0;
    let d0_shift = digamma(half_d0) - half_d0.ln();
    moderations
        .iter()
        .enumerate()
        .map(|(i, (_, fit))| {
            let s02 = (egpred[i] + d0_shift).exp();
            let post_var = (d0 * s02 + df[i] * sigma2[i]) / (d0 + df[i]);
            let post_df = d0 + df[i];
            let sca_t = fit.log2_fold_change / (fit.stdev_unscaled * post_var.sqrt());
            let sca_pvalue = 2.0 * student_t_sf(sca_t.abs(), post_df);
            let log_pvalue = student_t_two_sided_log_pvalue(sca_t, post_df);
            ScaStat {
                sca_t,
                sca_pvalue,
                log_pvalue,
            }
        })
        .collect()
}

/// The standard-eBayes fallback (deqms.py:72): `sca_t = t_stat`,
/// `sca_pvalue = p_value`. Used when too few valid points or all-equal counts.
fn fallback(moderations: &[(usize, LimmaModeration)], df_prior: f64) -> Vec<ScaStat> {
    let df_pooled = moderations
        .iter()
        .map(|(_, fit)| fit.df_residual)
        .sum::<f64>();
    moderations
        .iter()
        .map(|(_, fit)| ScaStat {
            sca_t: fit.t_statistic,
            sca_pvalue: fit.p_value,
            log_pvalue: student_t_two_sided_log_pvalue(
                fit.t_statistic,
                (fit.df_residual + df_prior).min(df_pooled),
            ),
        })
        .collect()
}

/// LOWESS prediction of `log_var` on `x` over the valid points, then linear
/// interpolation onto any unseen non-valid `x` (deqms.py:82-101).
///
/// statsmodels `lowess(endog, exog, frac, return_sorted=False)` sorts the pairs
/// internally and maps the fit back to input order; [`lowess_fit`] requires `x`
/// pre-sorted ascending, so the valid points are sorted, fitted, then scattered
/// back. Non-valid points with a finite `x` are linearly interpolated against
/// the sorted valid `(x, fitted)` curve, exactly matching the Python `np.interp`
/// on the unseen counts.
fn lowess_predict(log_var: &[f64], x: &[f64], valid: &[bool]) -> Vec<f64> {
    let n = log_var.len();
    let mut valid_indices = (0..n).filter(|i| valid[*i]).collect::<Vec<_>>();
    valid_indices.sort_by(|left, right| x[*left].total_cmp(&x[*right]));
    let sorted_x = valid_indices.iter().map(|i| x[*i]).collect::<Vec<_>>();
    let sorted_y = valid_indices
        .iter()
        .map(|i| log_var[*i])
        .collect::<Vec<_>>();
    let fitted = lowess_fit(&sorted_x, &sorted_y, LOWESS_FRAC, LOWESS_ITERATIONS);

    let mut y_pred = vec![f64::NAN; n];
    for (position, original) in valid_indices.iter().enumerate() {
        y_pred[*original] = fitted[position];
    }
    for i in 0..n {
        if !valid[i] && x[i].is_finite() {
            y_pred[i] = interp(x[i], &sorted_x, &fitted);
        }
    }
    y_pred
}

/// Linear interpolation matching `numpy.interp(query, xs, ys)` with `xs` sorted
/// ascending: clamps to the endpoints outside the range, linear between.
fn interp(query: f64, xs: &[f64], ys: &[f64]) -> f64 {
    let Some(first_x) = xs.first().copied() else {
        return f64::NAN;
    };
    let Some(last_x) = xs.last().copied() else {
        return f64::NAN;
    };
    if query <= first_x {
        return ys.first().copied().unwrap_or(f64::NAN);
    }
    if query >= last_x {
        return ys.last().copied().unwrap_or(f64::NAN);
    }
    let upper = xs.partition_point(|value| *value <= query);
    let lower = upper - 1;
    let (x_lo, x_hi) = (xs[lower], xs[upper]);
    let (y_lo, y_hi) = (ys[lower], ys[upper]);
    let span = x_hi - x_lo;
    if span == 0.0 {
        return y_lo;
    }
    y_lo + (y_hi - y_lo) * (query - x_lo) / span
}

/// Peak-to-peak (`np.ptp`) of a slice: max minus min, or 0 when empty.
fn ptp(values: &[f64]) -> f64 {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for value in values {
        if *value < min {
            min = *value;
        }
        if *value > max {
            max = *value;
        }
    }
    if min.is_finite() && max.is_finite() {
        max - min
    } else {
        0.0
    }
}

/// Grid search for the prior df `d0` (port of `_grid_search_d0`, deqms.py:122).
/// Scans `test_d0 = i/10` for `i` in `1..=n_genes*10`, tracking the `test_d0`
/// whose `trigamma(test_d0/2)` is closest to `mean_myfct`, with the same
/// early break: stop once the previous-previous diff is below the previous diff
/// (after at least two steps). Defaults to `0.1` when no step improves.
fn grid_search_d0(mean_myfct: f64, n_genes: usize) -> f64 {
    let max_iter = n_genes * 10;
    let mut best_d0 = 0.1;
    let mut best_diff = f64::INFINITY;
    let mut prev_prev = f64::INFINITY;
    let mut prev = f64::INFINITY;
    for i in 1..=max_iter {
        let test_d0 = i as f64 / 10.0;
        let diff = (mean_myfct - trigamma(test_d0 / 2.0)).abs();
        if diff < best_diff {
            best_diff = diff;
            best_d0 = test_d0;
        }
        if i > 2 && prev_prev < prev {
            break;
        }
        prev_prev = prev;
        prev = diff;
    }
    best_d0
}

#[cfg(test)]
mod tests {
    use super::deqms_two_group;
    use crate::de::{
        limma_two_group, Significance, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD,
    };

    // 15-protein / 3-vs-3 log2 fixture shared by the oracle scripts in the
    // scratchpad (deqms_oracle.py, deqms_lowess_probe.py).
    fn fixture_rows() -> Vec<Vec<f64>> {
        vec![
            vec![10.00, 10.20, 9.80, 12.00, 12.10, 11.90],
            vec![15.00, 15.10, 14.90, 13.00, 13.20, 12.80],
            vec![8.00, 8.10, 7.90, 8.05, 7.95, 8.00],
            vec![20.00, 22.00, 18.00, 21.00, 19.00, 23.00],
            vec![5.00, 5.01, 4.99, 6.00, 6.01, 5.99],
            vec![9.00, 9.50, 8.50, 9.20, 10.50, 7.50],
            vec![11.00, 11.30, 10.70, 12.50, 12.80, 12.20],
            vec![3.00, 3.20, 2.80, 4.00, 3.90, 4.20],
            vec![7.00, 7.05, 6.95, 7.50, 7.55, 7.45],
            vec![13.00, 13.40, 12.60, 11.00, 11.20, 10.80],
            vec![6.00, 6.30, 5.70, 6.10, 6.40, 5.80],
            vec![18.00, 18.20, 17.80, 16.00, 16.30, 15.70],
            vec![4.00, 4.40, 3.60, 5.00, 5.30, 4.70],
            vec![9.50, 9.70, 9.30, 9.55, 9.45, 9.60],
            vec![12.00, 12.50, 11.50, 14.00, 14.30, 13.70],
        ]
    }

    fn fixture_proteins() -> Vec<String> {
        (0..15).map(|i| format!("P{i:02}")).collect()
    }

    fn fixture_counts() -> Vec<f64> {
        vec![
            1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 4.0, 7.0, 11.0, 6.0, 9.0, 15.0, 2.0, 18.0,
        ]
    }

    // Golden oracle captured verbatim from `conda run -n Bigbio python` on
    // mokume's `run_deqms` / `_spectra_count_ebayes` (see scratchpad
    // deqms_oracle.py). The spectraCounteBayes path genuinely runs here (d0 = 1.0,
    // valid.sum() == 15 >= 10, ptp(log2 counts) == 4.392 > 1e-10), NOT the
    // fallback. Per-protein (protein, log2FC, sca_t, sca_pvalue).
    const ORACLE: &[(&str, f64, f64, f64)] = &[
        ("P00", -2.0, -16.678289372901432, 1.4156705245728974e-05),
        ("P01", 2.0, 16.52831838925866, 1.4800509999750121e-05),
        ("P02", 0.0, 0.0, 1.0),
        ("P03", -1.0, -0.6842879424608374, 0.5242143386505393),
        ("P04", -1.0, -16.192422308943645, 1.6374181949382853e-05),
        (
            "P05",
            -0.06666666666666643,
            -0.08112534703275041,
            0.9384894938829584,
        ),
        ("P06", -1.5, -6.432268772004748, 0.0013493270988088847),
        (
            "P07",
            -1.0333333333333337,
            -7.522857250144055,
            0.0006568918407631296,
        ),
        (
            "P08",
            -0.5000000000000009,
            -7.358554909913771,
            0.0007278314810727803,
        ),
        ("P09", 2.0, 8.301510632638784, 0.0004142498025629735),
        (
            "P10",
            -0.09999999999999964,
            -0.44402014374156507,
            0.6755855902200026,
        ),
        (
            "P11",
            2.0000000000000018,
            10.142604275847223,
            0.00015971865918412806,
        ),
        ("P12", -1.0, -3.7241137860450872, 0.01365385543271829),
        (
            "P13",
            -0.033333333333333215,
            -0.2865592057962772,
            0.7859443831678532,
        ),
        ("P14", -2.0, -6.435006203233639, 0.0013467251630308547),
    ];

    // Python's adj_pvalue rank as ORDERED TIE GROUPS. Within each group the
    // proteins share a bit-identical BH-adjusted p-value (verified in the
    // oracle), so their internal order is a sort-stability artifact, not a
    // numeric ranking. The groups themselves are strictly ordered and that
    // ordering must match Python exactly. (The only place the LOWESS divergence
    // could otherwise flip the visible order is inside these ties, e.g. P07/P08.)
    const ORACLE_RANK_GROUPS: &[&[&str]] = &[
        &["P00", "P01", "P04"], // adj = 8.187e-05
        &["P11"],               // adj = 5.989e-04
        &["P09"],               // adj = 1.243e-03
        &["P07", "P08"],        // adj = 1.560e-03
        &["P06", "P14"],        // adj = 2.249e-03
        &["P12"],               // adj = 2.048e-02
        &["P03"],               // adj = 7.148e-01
        &["P10"],               // adj = 8.445e-01
        &["P13"],               // adj = 9.069e-01
        &["P02", "P05"],        // adj = 1.0
    ];

    // ACHIEVED tolerance (measured against the Python oracle on this fixture):
    //   - log2FC: cell-exact (1e-9) -- it is the limma coefficient, untouched by
    //     the count moderation;
    //   - sca_t: max relative error 5.06e-3 (P08); asserted at relative 6e-3;
    //   - sca_pvalue: max relative error 2.31e-2 (P08), max absolute 8.13e-5;
    //     asserted at relative 3e-2 with an absolute 1e-4 floor for tiny tails.
    // The divergence comes entirely from the LOWESS smoother differing from
    // statsmodels by ~2e-3 (window-size / `delta` differences), which propagates
    // through mean_myfct -> the d0 grid search -> post_var into sca_t, and the
    // steep t-tail amplifies the residual t-error into a larger relative
    // p-error. The protein RANK by adj_pvalue matches Python's strict
    // between-group ordering exactly; order within a bit-identical BH tie group
    // (e.g. P07/P08) is undefined and checked set-wise.
    #[test]
    fn deqms_two_group_matches_python_oracle() {
        let rows = fixture_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = fixture_proteins();
        let counts = fixture_counts();
        let results = deqms_two_group(
            &proteins,
            &refs,
            3,
            3,
            &counts,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
        );
        assert_eq!(results.len(), ORACLE.len());

        for &(protein, log2fc, sca_t, sca_pvalue) in ORACLE {
            let Some(row) = results.iter().find(|r| r.protein == protein) else {
                panic!("{protein} missing from deqms results");
            };
            // log2FC is the untouched limma coefficient: cell-exact.
            assert!(
                (row.log2_fold_change - log2fc).abs() <= 1e-9,
                "{protein} log2FC actual={} expected={log2fc}",
                row.log2_fold_change
            );
            // sca_t: relative 6e-3 (LOWESS-bounded; measured max 5.06e-3). The
            // absolute floor of 1.0 handles the exact-zero P02.
            let t_tol = sca_t.abs().max(1.0) * 6e-3;
            assert!(
                (row.t_statistic - sca_t).abs() <= t_tol,
                "{protein} sca_t actual={} expected={sca_t} tol={t_tol}",
                row.t_statistic
            );
            // sca_pvalue: relative 3e-2 (measured max 2.31e-2) with an absolute
            // 1e-4 floor (measured max abs 8.13e-5) for the steep small-p tail.
            let p_tol = (sca_pvalue.abs() * 3e-2).max(1e-4);
            assert!(
                (row.p_value - sca_pvalue).abs() <= p_tol,
                "{protein} sca_pvalue actual={} expected={sca_pvalue} tol={p_tol}",
                row.p_value
            );
        }

        // RANK match by adj_pvalue, tie-aware: walk the Rust ranking group by
        // group; each emitted group must contain exactly the Python tie group's
        // proteins (order within a bit-identical-adj tie is undefined). This
        // proves the strict between-group ordering matches Python exactly.
        let rank = results
            .iter()
            .map(|r| r.protein.as_str())
            .collect::<Vec<_>>();
        let mut cursor = 0;
        for group in ORACLE_RANK_GROUPS {
            let actual_slice = &rank[cursor..cursor + group.len()];
            let mut actual_sorted = actual_slice.to_vec();
            actual_sorted.sort_unstable();
            let mut expected_sorted = group.to_vec();
            expected_sorted.sort_unstable();
            assert_eq!(
                actual_sorted, expected_sorted,
                "deqms rank tie-group mismatch at position {cursor}: got {actual_slice:?}, \
                 expected group {group:?} (order within a BH tie is undefined)"
            );
            cursor += group.len();
        }
        assert_eq!(cursor, rank.len(), "deqms rank length mismatch");
    }

    // Non-vacuity guard: the oracle sca_t must differ materially from plain limma
    // t on the SAME matrix, proving the spectraCounteBayes count path is
    // exercised (not silently equal to eBayes). The largest divergence in the
    // Python oracle is on P04 (|sca_t - limma_t| ~= 6.85).
    #[test]
    fn deqms_oracle_differs_from_plain_limma() {
        let rows = fixture_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = fixture_proteins();
        let counts = fixture_counts();

        let deqms = deqms_two_group(
            &proteins,
            &refs,
            3,
            3,
            &counts,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
        );
        let limma = limma_two_group(
            &proteins,
            &refs,
            3,
            3,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
        );

        let mut max_diff = 0.0_f64;
        for d in &deqms {
            let Some(l) = limma.iter().find(|r| r.protein == d.protein) else {
                panic!("{} missing from limma results", d.protein);
            };
            max_diff = max_diff.max((d.t_statistic - l.t_statistic).abs());
        }
        // Python oracle: max |sca_t - limma_t| ~= 6.85 (P04). Require a large
        // divergence so the count path is unmistakably exercised.
        assert!(
            max_diff > 5.0,
            "deqms vs limma max |t| divergence too small ({max_diff}); count path not exercised"
        );
    }

    // Fallback proof: all-equal counts make `ptp(log2 counts) < 1e-10`, so DEqMS
    // falls back to plain eBayes and is cell-exact (1e-9) against limma_two_group
    // on the same matrix (deqms.py:76-80).
    #[test]
    fn deqms_all_equal_counts_falls_back_to_ebayes() {
        let rows = fixture_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = fixture_proteins();
        let counts = vec![1.0; rows.len()];

        let deqms = deqms_two_group(
            &proteins,
            &refs,
            3,
            3,
            &counts,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
        );
        let limma = limma_two_group(
            &proteins,
            &refs,
            3,
            3,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
        );
        assert_eq!(deqms.len(), limma.len());

        for d in &deqms {
            let Some(l) = limma.iter().find(|r| r.protein == d.protein) else {
                panic!("{} missing from limma results", d.protein);
            };
            assert!(
                (d.log2_fold_change - l.log2_fold_change).abs() <= 1e-9,
                "{} log2FC deqms={} limma={}",
                d.protein,
                d.log2_fold_change,
                l.log2_fold_change
            );
            assert!(
                (d.t_statistic - l.t_statistic).abs() <= 1e-9,
                "{} t deqms={} limma={}",
                d.protein,
                d.t_statistic,
                l.t_statistic
            );
            assert!(
                (d.p_value - l.p_value).abs() <= 1e-9,
                "{} pvalue deqms={} limma={}",
                d.protein,
                d.p_value,
                l.p_value
            );
            assert!(
                (d.adj_p_value - l.adj_p_value).abs() <= 1e-9,
                "{} adj_pvalue deqms={} limma={}",
                d.protein,
                d.adj_p_value,
                l.adj_p_value
            );
            assert_eq!(
                d.significance, l.significance,
                "{} significance deqms vs limma mismatch",
                d.protein
            );
        }

        // Sanity: the fallback still produces a non-trivial significance call.
        assert!(deqms
            .iter()
            .any(|r| r.significance != Significance::Unchanged));
    }
}
