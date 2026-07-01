//! LimROTS (limma empirical-Bayes + ROTS reproducibility optimisation) two-group
//! differential expression.
//!
//! This is a FAITHFUL ALGORITHM PORT of mokume's `analysis/limrots.py` (a pure
//! numpy/scipy reimplementation of the LimROTS hybrid method), NOT a cell-for-
//! cell match. LimROTS computes the moderated d-statistic
//! `d = beta / max(a1 + a2*s_post*u, 1e-10)` (limrots.py:54), where `beta`,
//! `s_post = sqrt(s2_post)`, and `u = stdev_unscaled` all come from a limma
//! `lmFit -> contrasts.fit -> eBayes` run (limrots.py:170 `_boot_ebayes`). The
//! `(a1, a2)` pair is chosen by a ROTS-style bootstrap/permutation reproducibility
//! search (limrots.py:70 `_limrots_bootstrap_optimize`), and the per-protein
//! p-value comes from a permutation null (limrots.py:197 `_pvalue_permutation`).
//!
//! Because the optimisation and the p-value both consume random draws, Python's
//! own output is unstable across seeds. Bit-matching numpy's PCG64 stream is
//! therefore neither possible (the workspace has no `rand`/numpy) nor meaningful,
//! so it is intentionally NOT attempted -- the same policy ROTS' port uses.
//!
//! WHAT IS MATCHED EXACTLY (cell-exact, 1e-9, RNG-independent):
//!   - `_limrots_d_stat` (limrots.py:54), including the `1e-10` denominator
//!     clamp.
//!   - `_boot_ebayes` (limrots.py:170): `beta`, `s_post`, `u` for a fixed
//!     bootstrap `sample_idx` AND a fixed permuted design (the reorder path);
//!     both reuse the cell-exact `run_limma_moderation` engine (limma is
//!     order-invariant within a group, so reordering columns to `[A.., B..]`
//!     reproduces Python's per-protein fit exactly -- proven in
//!     `boot_ebayes_*_matches_python_oracle`).
//!   - The p-value counting given a fixed set of permuted d_p arrays
//!     (limrots.py:227): `count` STARTS AT 1, accumulates the POOLED
//!     `sum_g(|d_p|_g >= |d_obs|_j)` across all genes for each gene `j`, and
//!     divides by `n_genes * (n_perm + 1)` clamped to 1.0 (limrots.py:229).
//!   - `log2FC` (= `beta`, the full-data eBayes contrast coefficient
//!     `mean_A - mean_B`, limrots.py:287/309): deterministic, RNG-independent.
//!
//! WHAT IS NOT MATCHED: the exact stochastic stream. The bootstrap/permutation
//! here uses the shared self-contained splitmix64 PRNG from `rots.rs` (no
//! external `rand` dependency), seeded by a fixed constant. The result matches
//! Python in ALGORITHM and DISTRIBUTION (within LimROTS' own seed-noise), not
//! cell-for-cell. The stochastic wrapper is covered by property tests (p in
//! [0,1], finite `d_stat`, the selected `(a1, a2)` lands on a valid grid point),
//! and `log2FC` is asserted cell-exact since it is RNG-independent.
//!
//! Reference
//! ---------
//! Anwar AM, Jeba A, Lahti L, Coffey E. LimROTS: a hybrid method integrating
//! empirical Bayes and reproducibility-optimized statistics for robust
//! differential expression analysis. *Bioinformatics*. 2025.

use super::correct::bh_adjust;
use super::limma::{run_limma_moderation, LimmaModeration};
use super::rots::{build_n_grid, build_ssq_grid, rank_abs, reproducibility_score, SplitMix64};
use super::{classify, DeResult};

/// Default bootstrap iterations, matching `run_limrots(n_boot=100)`
/// (limrots.py:42).
const DEFAULT_N_BOOT: usize = 100;

/// Fixed PRNG seed for the bootstrap/permutation draws. LimROTS is seed-unstable,
/// so this is a deterministic constant rather than `run_limrots`'s `seed=42`
/// (which only makes sense against numpy's PCG64 stream). A fixed seed keeps the
/// Rust output reproducible run-to-run, matching the ROTS port's policy.
const PRNG_SEED: u64 = 0x9E37_79B9_7F4A_7C15;

/// Lower clamp on the d-statistic denominator (limrots.py:66).
const DENOM_FLOOR: f64 = 1e-10;

/// Minimum standard deviation in the reproducibility table, matching
/// `std_real[std_real < 1e-10] = 1e-10` (limrots.py:163).
const SD_FLOOR: f64 = 1e-10;

/// One protein's per-fit `(beta, s_post, u)` from `_boot_ebayes`
/// (limrots.py:170).
struct BootFit {
    beta: Vec<f64>,
    s_post: Vec<f64>,
    u: Vec<f64>,
}

/// Port of `_limrots_d_stat` (limrots.py:54):
/// `d = beta / max(a1 + a2 * s_post * u, 1e-10)`, element-wise.
fn limrots_d_stat(beta: &[f64], s_post: &[f64], u: &[f64], a1: f64, a2: f64) -> Vec<f64> {
    beta.iter()
        .zip(s_post)
        .zip(u)
        .map(|((b, s), uu)| {
            let denom = a1 + a2 * s * uu;
            let denom = if denom < DENOM_FLOOR {
                DENOM_FLOOR
            } else {
                denom
            };
            b / denom
        })
        .collect()
}

/// `_boot_ebayes` for a bootstrap `sample_idx` (limrots.py:170): build
/// `boot_mat = mat[:, sample_idx]` with `boot_design = design[sample_idx]`, run
/// limma `lmFit -> contrasts.fit -> eBayes`, return `(beta, s_post, u)`.
///
/// In the bootstrap path the indices resample WITHIN groups, so the column layout
/// stays `[A.. , B..]` (group A block of length `n_a`, then group B). The limma
/// cell-means fit depends only on the per-group cell values, not their order
/// within a group, so reusing `run_limma_moderation(rows, n_a, n_b)` on the
/// reindexed columns reproduces Python's `_boot_ebayes` cell-exact.
fn boot_ebayes_bootstrap(mat: &[Vec<f64>], sample_idx: &[usize], n_a: usize) -> BootFit {
    let n_b = sample_idx.len() - n_a;
    let boot_rows = mat
        .iter()
        .map(|row| sample_idx.iter().map(|&col| row[col]).collect::<Vec<f64>>())
        .collect::<Vec<_>>();
    boot_ebayes_from_rows(&boot_rows, n_a, n_b)
}

/// `_boot_ebayes` for a PERMUTED DESIGN (limrots.py:138/218): Python runs eBayes
/// on the ORIGINAL data with `design[perm_idx]` (the group LABELS shuffled) and
/// `sample_idx = arange`. limma's two-group cell-means fit depends only on the
/// A/B partition (not the column order within a group), so we apply the
/// permutation to the LABELS, then reorder each protein row's columns to put all
/// group-A-labeled columns first (preserving relative order) then group-B,
/// and call `run_limma_moderation`. `perm_idx` permutes the row indices of the
/// `[A=1, B=0]` design column, so the A-label set is exactly `{i : perm_idx[i] <
/// n_a_original}` -- a permutation of the original labels, so the A count is
/// preserved.
fn boot_ebayes_permuted(mat: &[Vec<f64>], perm_idx: &[usize], n_a: usize) -> BootFit {
    // After `design[perm_idx]`, sample position `i` carries the label of original
    // design row `perm_idx[i]`: group A iff `perm_idx[i] < n_a` (design col0 == 1
    // for rows 0..n_a). Collect A-labeled positions then B-labeled positions,
    // each preserving ascending order.
    let mut a_cols = Vec::new();
    let mut b_cols = Vec::new();
    for (position, &source) in perm_idx.iter().enumerate() {
        if source < n_a {
            a_cols.push(position);
        } else {
            b_cols.push(position);
        }
    }
    let new_n_a = a_cols.len();
    let new_n_b = b_cols.len();
    let reorder = a_cols.into_iter().chain(b_cols).collect::<Vec<_>>();
    let reordered = mat
        .iter()
        .map(|row| reorder.iter().map(|&col| row[col]).collect::<Vec<f64>>())
        .collect::<Vec<_>>();
    boot_ebayes_from_rows(&reordered, new_n_a, new_n_b)
}

/// Shared limma `lmFit -> contrasts.fit -> eBayes` for an already `[A.., B..]`
/// laid-out matrix, returning `(beta, s_post, u)` aligned to `rows`. Proteins
/// `run_limma_moderation` drops (fewer than two finite obs in a group) get
/// `beta = 0, s_post = 1, u = 0` so the d-statistic is finite and the gene index
/// alignment is preserved (Python's `_lm_fit_with_na` leaves such genes as NaN,
/// but they cannot pass `filter_testable` upstream, so this only guards the
/// degenerate bootstrap-resample-into-NaN edge).
fn boot_ebayes_from_rows(rows: &[Vec<f64>], n_a: usize, n_b: usize) -> BootFit {
    let n_genes = rows.len();
    let mut beta = vec![0.0; n_genes];
    let mut s_post = vec![1.0; n_genes];
    let mut u = vec![0.0; n_genes];
    let row_refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let (moderations, _df_prior) = run_limma_moderation(&row_refs, n_a, n_b);
    for (index, moderation) in &moderations {
        fill_boot_fit(&mut beta, &mut s_post, &mut u, *index, moderation);
    }
    BootFit { beta, s_post, u }
}

/// Scatter one protein's moderation into the `(beta, s_post, u)` vectors:
/// `beta = coefficients[:,0]`, `s_post = sqrt(s2_post)`, `u = stdev_unscaled[:,0]`
/// (limrots.py:186-188).
fn fill_boot_fit(
    beta: &mut [f64],
    s_post: &mut [f64],
    u: &mut [f64],
    index: usize,
    moderation: &LimmaModeration,
) {
    if let (Some(b), Some(s), Some(uu)) =
        (beta.get_mut(index), s_post.get_mut(index), u.get_mut(index))
    {
        *b = moderation.log2_fold_change;
        *s = moderation.s2_post.sqrt();
        *uu = moderation.stdev_unscaled;
    }
}

/// Port of `_limrots_bootstrap_optimize` (limrots.py:70): find the optimal
/// `(a1, a2)` by bootstrap reproducibility. `a1_grid` is the ROTS ssq grid,
/// `a2_grid` is all ones (so `a2 == 1` always; there is NO raw-D last-row branch,
/// unlike ROTS). Returns `(a1, a2)`; only the a1-axis index of the argmax is
/// used.
fn bootstrap_optimize(
    mat: &[Vec<f64>],
    n_a: usize,
    n_total: usize,
    n_boot: usize,
    rng: &mut SplitMix64,
) -> (f64, f64) {
    let a1_grid = build_ssq_grid();
    let n_params = a1_grid.len();
    let n_b = n_total - n_a;

    // k_grid == ROTS n_grid kept `< k_max`, k_max = max(1, n_genes // 4)
    // (limrots.py:92-104). `build_n_grid` already applies the `< k_max` filter.
    let k_max = (mat.len() / 4).max(1);
    let mut k_values = build_n_grid(k_max);
    if k_values.is_empty() {
        k_values = vec![k_max.max(1)];
    }
    let n_k = k_values.len();

    // repro_real[b][pi][ki], repro_perm[b][pi][ki].
    let mut repro_real = vec![vec![vec![0.0; n_k]; n_params]; n_boot];
    let mut repro_perm = vec![vec![vec![0.0; n_k]; n_params]; n_boot];

    for b in 0..n_boot {
        // idx1/idx2: two INDEPENDENT bootstrap resamples, each = choice(n_a, n_a)
        // ++ choice(n_b, n_b)+n_a (limrots.py:113-124).
        let idx1 = bootstrap_indices(rng, n_a, n_b);
        let idx2 = bootstrap_indices(rng, n_a, n_b);
        let fit1 = boot_ebayes_bootstrap(mat, &idx1, n_a);
        let fit2 = boot_ebayes_bootstrap(mat, &idx2, n_a);

        // perm_idx applied to the DESIGN (limrots.py:135).
        let perm_idx = rng.permutation(n_total);
        let fit_p = boot_ebayes_permuted(mat, &perm_idx, n_a);

        for (pi, &a1) in a1_grid.iter().enumerate() {
            // a2_grid is all ones (limrots.py:89).
            let d1 = limrots_d_stat(&fit1.beta, &fit1.s_post, &fit1.u, a1, 1.0);
            let d2 = limrots_d_stat(&fit2.beta, &fit2.s_post, &fit2.u, a1, 1.0);
            let dp = limrots_d_stat(&fit_p.beta, &fit_p.s_post, &fit_p.u, a1, 1.0);
            let r1 = rank_abs(&d1);
            let r2 = rank_abs(&d2);
            let rp = rank_abs(&dp);
            for (ki, &k) in k_values.iter().enumerate() {
                repro_real[b][pi][ki] = reproducibility_score(&r1, &r2, k);
                repro_perm[b][pi][ki] = reproducibility_score(&r1, &rp, k);
            }
        }
    }

    // z = (mean_real - mean_perm) / max(std_real(ddof=1), 1e-10); best =
    // nanargmax(z) unraveled; return (a1_grid[best_pi], a2_grid[best_pi]==1.0)
    // (limrots.py:160-167).
    let mut best_value = f64::NEG_INFINITY;
    let mut best_pi = 0usize;
    for pi in 0..n_params {
        for ki in 0..n_k {
            let real_col = (0..n_boot)
                .map(|b| repro_real[b][pi][ki])
                .collect::<Vec<_>>();
            let perm_col = (0..n_boot)
                .map(|b| repro_perm[b][pi][ki])
                .collect::<Vec<_>>();
            let mean_real = mean(&real_col);
            let mean_perm = mean(&perm_col);
            let mut std_real = std_ddof1(&real_col);
            if std_real < SD_FLOOR {
                std_real = SD_FLOOR;
            }
            let z = (mean_real - mean_perm) / std_real;
            // nanargmax: skip NaN, first finite max wins (row-major over pi, ki).
            if z.is_finite() && z > best_value {
                best_value = z;
                best_pi = pi;
            }
        }
    }

    (a1_grid[best_pi], 1.0)
}

/// One bootstrap index vector for LimROTS (limrots.py:113-118): `choice(n_a, n_a,
/// replace=True)` for the group-A block, then `choice(n_b, n_b, replace=True) +
/// n_a` for the group-B block. The layout stays `[A.., B..]`.
fn bootstrap_indices(rng: &mut SplitMix64, n_a: usize, n_b: usize) -> Vec<usize> {
    let mut indices = rng.choice_with_replacement(n_a, n_a);
    for value in rng.choice_with_replacement(n_b, n_b) {
        indices.push(value + n_a);
    }
    indices
}

/// Port of `_pvalue_permutation` (limrots.py:197): permutation null of the
/// LimROTS d-statistic. `count` STARTS AT 1 (limrots.py:212); for each
/// permutation it adds the POOLED count `sum_g(|d_p|_g >= |d_obs|_j)` over all
/// genes `g` for each observed gene `j` (limrots.py:227); the p-value is
/// `min(count / (n_genes * (n_perm + 1)), 1.0)` (limrots.py:229).
fn pvalue_permutation(
    d_obs: &[f64],
    mat: &[Vec<f64>],
    n_a: usize,
    n_total: usize,
    alpha: (f64, f64),
    n_perm: usize,
    rng: &mut SplitMix64,
) -> Vec<f64> {
    let (a1, a2) = alpha;
    let n_genes = d_obs.len();
    let abs_d = d_obs.iter().map(|value| value.abs()).collect::<Vec<_>>();
    let mut count = vec![1.0; n_genes];

    for _ in 0..n_perm {
        let perm_idx = rng.permutation(n_total);
        let fit = boot_ebayes_permuted(mat, &perm_idx, n_a);
        let d_p = limrots_d_stat(&fit.beta, &fit.s_post, &fit.u, a1, a2);
        accumulate_counts(&mut count, &d_p, &abs_d);
    }

    let denom = (n_genes * (n_perm + 1)) as f64;
    count.iter().map(|value| (value / denom).min(1.0)).collect()
}

/// Add, for every gene `j`, the number of permuted `|d_p|_g` (over all genes `g`)
/// that are `>=` `|d_obs|_j` (the pooled `np.sum(...)` of limrots.py:227).
fn accumulate_counts(count: &mut [f64], d_p: &[f64], abs_d_obs: &[f64]) {
    for (slot, &threshold) in count.iter_mut().zip(abs_d_obs) {
        let hits = d_p.iter().filter(|value| value.abs() >= threshold).count();
        *slot += hits as f64;
    }
}

/// Arithmetic mean (`np.mean`), 0.0 for an empty slice.
fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Sample standard deviation with `ddof=1` (`np.std(ddof=1)`); 0.0 when fewer
/// than two points (numpy yields NaN there, but the caller floors it to 1e-10).
fn std_ddof1(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return 0.0;
    }
    let m = mean(values);
    let var = values.iter().map(|value| (value - m).powi(2)).sum::<f64>() / (n - 1) as f64;
    var.sqrt()
}

/// `filter_testable` predicate (\_helpers.py:62-64): at least two finite
/// observations in BOTH groups.
fn testable(row: &[f64], n_a: usize, n_b: usize) -> bool {
    const MIN_PER_GROUP: usize = 2;
    let finite_a = row[..n_a].iter().filter(|value| value.is_finite()).count();
    let finite_b = row[n_a..n_a + n_b]
        .iter()
        .filter(|value| value.is_finite())
        .count();
    finite_a >= MIN_PER_GROUP && finite_b >= MIN_PER_GROUP
}

/// `np.nanmean` plus the finite count, for the per-group `mean_*` / `n_*`
/// summary columns (\_helpers.py:47-50). An all-nan slice yields a NaN mean and
/// a zero count.
fn nan_mean_count(values: &[f64]) -> (f64, usize) {
    let mut sum = 0.0;
    let mut count = 0usize;
    for value in values {
        if value.is_finite() {
            sum += value;
            count += 1;
        }
    }
    if count == 0 {
        (f64::NAN, 0)
    } else {
        (sum / count as f64, count)
    }
}

/// Two-group LimROTS differential expression.
///
/// `proteins[i]` labels `rows[i]`; each row holds log2 intensities with the
/// first `n_a` columns in condition A and the next `n_b` in condition B. Proteins
/// without at least two finite observations in BOTH groups are dropped before
/// testing (matching `filter_testable`). Returns one [`DeResult`] per testable
/// protein, sorted by adjusted p-value ascending (BH), reproducing mokume's
/// `run_limrots` output contract: `log2FC` is the full-data eBayes contrast
/// coefficient (`beta`, deterministic), `pvalue` the permutation p-value,
/// `adj_pvalue` its BH adjustment, and `t_statistic` carries the LimROTS `d_stat`
/// (the method's extra column).
///
/// NOTE: this is a faithful algorithm port; with the internal fixed-seed PRNG the
/// p-values match Python in distribution (within LimROTS' seed-noise), not
/// cell-for-cell. `log2_fold_change` IS deterministic and matches Python exactly.
pub fn limrots_two_group(
    proteins: &[String],
    rows: &[&[f64]],
    n_a: usize,
    n_b: usize,
    fdr_threshold: f64,
    log2fc_threshold: f64,
) -> Vec<DeResult> {
    limrots_two_group_seeded(
        proteins,
        rows,
        n_a,
        n_b,
        fdr_threshold,
        log2fc_threshold,
        DEFAULT_N_BOOT,
        PRNG_SEED,
    )
}

/// Seedable core of [`limrots_two_group`], exposed for deterministic-property
/// tests that want a fixed `n_boot`/`seed`. The public entry point hardcodes the
/// Python defaults (`n_boot=100`) and the fixed module PRNG seed.
#[allow(clippy::too_many_arguments)]
fn limrots_two_group_seeded(
    proteins: &[String],
    rows: &[&[f64]],
    n_a: usize,
    n_b: usize,
    fdr_threshold: f64,
    log2fc_threshold: f64,
    n_boot: usize,
    seed: u64,
) -> Vec<DeResult> {
    // filter_testable: keep proteins with >= 2 finite obs per group.
    let kept = rows
        .iter()
        .enumerate()
        .filter(|(_, row)| testable(row, n_a, n_b))
        .map(|(index, row)| (index, row.to_vec()))
        .collect::<Vec<_>>();
    if kept.is_empty() {
        return Vec::new();
    }

    let mat = kept.iter().map(|(_, row)| row.clone()).collect::<Vec<_>>();
    let n_total = n_a + n_b;
    let mut rng = SplitMix64::new(seed);

    // Optimise (a1, a2) over the bootstrap reproducibility grid.
    let (best_a1, best_a2) = bootstrap_optimize(&mat, n_a, n_total, n_boot, &mut rng);

    // Observed d-statistic from the full-data eBayes fit (limrots.py:283-290).
    let full_fit = boot_ebayes_from_rows(&mat, n_a, n_b);
    let d_stat = limrots_d_stat(
        &full_fit.beta,
        &full_fit.s_post,
        &full_fit.u,
        best_a1,
        best_a2,
    );

    // n_perm = min(n_boot, 100) (limrots.py:292).
    let n_perm = n_boot.min(100);
    let p_values = pvalue_permutation(
        &d_stat,
        &mat,
        n_a,
        n_total,
        (best_a1, best_a2),
        n_perm,
        &mut rng,
    );
    let adjusted = bh_adjust(&p_values);

    let mut results = Vec::with_capacity(kept.len());
    for (position, (original_index, row)) in kept.iter().enumerate() {
        let Some(protein) = proteins.get(*original_index) else {
            continue;
        };
        let (mean_a, count_a) = nan_mean_count(&row[..n_a]);
        let (mean_b, count_b) = nan_mean_count(&row[n_a..]);
        // log2FC = beta (the full-data eBayes contrast coefficient = mean_A -
        // mean_B), deterministic (limrots.py:287/309).
        let log2_fold_change = full_fit.beta[position];
        let p_value = p_values[position];
        let adj_p_value = adjusted[position];
        let significance = classify(
            adj_p_value,
            log2_fold_change,
            fdr_threshold,
            log2fc_threshold,
        );
        results.push(DeResult {
            protein: protein.clone(),
            log2_fold_change,
            p_value,
            adj_p_value,
            t_statistic: d_stat[position],
            // LimROTS emits no AveExpr/B columns (its extra column is t_stat ==
            // d_stat, carried in t_statistic), so these limma-only quirks are
            // NaN/0.0 and must not be written for limrots.
            ave_expr: f64::NAN,
            b: 0.0,
            mean_a,
            mean_b,
            n_a: count_a,
            n_b: count_b,
            significance,
        });
    }

    results.sort_by(|left, right| left.adj_p_value.total_cmp(&right.adj_p_value));
    results
}

#[cfg(test)]
mod tests {
    use super::super::rots::{build_ssq_grid, SplitMix64};
    use super::{
        boot_ebayes_bootstrap, boot_ebayes_permuted, bootstrap_optimize, limrots_d_stat,
        limrots_two_group, limrots_two_group_seeded, pvalue_permutation,
    };
    use crate::de::{Significance, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD};

    const TOL: f64 = 1e-9;

    fn assert_close(actual: f64, expected: f64, label: &str) {
        assert!(
            (actual - expected).abs() <= TOL,
            "{label}: actual={actual} expected={expected} diff={}",
            (actual - expected).abs()
        );
    }

    // 6-gene, 3-vs-3 fixture shared with the Python oracle (limrots_oracle.py).
    fn fixture() -> Vec<Vec<f64>> {
        vec![
            vec![10.0, 10.2, 9.8, 12.0, 12.1, 11.9],
            vec![15.0, 15.1, 14.9, 13.0, 13.2, 12.8],
            vec![8.0, 8.1, 7.9, 8.05, 7.95, 8.0],
            vec![20.0, 22.0, 18.0, 21.0, 19.0, 23.0],
            vec![5.0, 5.01, 4.99, 6.0, 6.01, 5.99],
            vec![9.0, 9.5, 8.5, 9.2, 10.5, 7.5],
        ]
    }

    // ---- _limrots_d_stat incl. clamp (oracle: limrots_oracle.py) ----
    #[test]
    fn limrots_d_stat_matches_python_oracle() {
        let beta = [2.0, -1.5, 0.0, 3.3, -4.0];
        let s_post = [0.5, 1.0, 2.0, 0.0, 100.0];
        let u = [0.8, 1.2, 0.5, 0.0, 0.9];
        let d = limrots_d_stat(&beta, &s_post, &u, 0.05, 1.0);
        let expected = [
            4.444444444444445,
            -1.2,
            0.0,
            65.99999999999999,
            -0.04441976679622432,
        ];
        for (i, (a, e)) in d.iter().zip(expected).enumerate() {
            assert_close(*a, e, &format!("d_stat[{i}] a1=0.05 a2=1.0"));
        }
        // a1=1, a2=0 -> d == beta.
        let raw = limrots_d_stat(&beta, &s_post, &u, 1.0, 0.0);
        for (i, (a, e)) in raw.iter().zip(beta).enumerate() {
            assert_close(*a, e, &format!("d_stat[{i}] a1=1 a2=0"));
        }
        // denom clamp: denom = 0 + 0*0 = 0 -> max(.,1e-10) -> 1/1e-10 = 1e10.
        let clamped = limrots_d_stat(&[1.0], &[0.0], &[0.0], 0.0, 0.0);
        assert_close(clamped[0], 1e10, "d_stat clamp");
    }

    // ---- _boot_ebayes FIXED bootstrap sample_idx (oracle: limrots_oracle.py).
    // Proves the bootstrap reindex reuses the cell-exact limma engine. ----
    #[test]
    fn boot_ebayes_bootstrap_matches_python_oracle() {
        let mat = fixture();
        // sample_idx [0,2,1,4,3,5]: within-group reorder of the [A.., B..] layout.
        let fit = boot_ebayes_bootstrap(&mat, &[0, 2, 1, 4, 3, 5], 3);
        let expected_beta = [-2.0, 2.0, 0.0, -1.0, -1.0, -0.06666666666666643];
        let expected_s_post = [
            0.14967183943258613,
            0.14967183943258613,
            0.07716474268929527,
            1.8732915428578778,
            0.023657886142284618,
            1.0501447642495674,
        ];
        let expected_u = [0.816496580927726; 6];
        for (i, (a, e)) in fit.beta.iter().zip(expected_beta).enumerate() {
            assert_close(*a, e, &format!("boot beta[{i}]"));
        }
        for (i, (a, e)) in fit.s_post.iter().zip(expected_s_post).enumerate() {
            assert_close(*a, e, &format!("boot s_post[{i}]"));
        }
        for (i, (a, e)) in fit.u.iter().zip(expected_u).enumerate() {
            assert_close(*a, e, &format!("boot u[{i}]"));
        }
    }

    // ---- _boot_ebayes FIXED permuted design, the reorder path (oracle:
    // limrots_oracle.py). PROVES the permuted-design reorder equivalence: Python
    // runs eBayes on `design[perm_idx]` (labels shuffled), we reorder columns to
    // `[A.., B..]` by label and reuse the same engine -- cell-exact. ----
    #[test]
    fn boot_ebayes_permuted_matches_python_oracle() {
        let mat = fixture();
        // perm_idx [3,0,4,1,5,2] -> design col0 = [0,1,0,1,0,1]: positions 1,3,5
        // are group A, 0,2,4 are group B.
        let fit = boot_ebayes_permuted(&mat, &[3, 0, 4, 1, 5, 2], 3);
        let expected_beta = [
            0.7333333333333343,
            -0.7333333333333343,
            0.09999999999999964,
            3.0,
            0.33333333333333304,
            -0.5999999999999996,
        ];
        let expected_s_post = [
            1.0452868618399844,
            1.0452868618399846,
            0.15546172151340235,
            0.9116314605384261,
            0.5402489123490589,
            0.9648002104031316,
        ];
        for (i, (a, e)) in fit.beta.iter().zip(expected_beta).enumerate() {
            assert_close(*a, e, &format!("perm beta[{i}]"));
        }
        for (i, (a, e)) in fit.s_post.iter().zip(expected_s_post).enumerate() {
            assert_close(*a, e, &format!("perm s_post[{i}]"));
        }

        // A second, more thoroughly mixed permutation [5,2,0,4,1,3] -> col0 =
        // [0,1,1,0,1,0]: positions 1,2,4 group A, 0,3,5 group B.
        let fit2 = boot_ebayes_permuted(&mat, &[5, 2, 0, 4, 1, 3], 3);
        let expected_beta2 = [
            -0.6000000000000014,
            0.8000000000000007,
            -0.033333333333333215,
            -1.6666666666666679,
            -0.3266666666666671,
            0.9333333333333336,
        ];
        let expected_s_post2 = [
            1.0673891644315912,
            1.0277066314557972,
            0.2081833587939854,
            1.636473512619515,
            0.5513978103642335,
            0.8811427634509565,
        ];
        for (i, (a, e)) in fit2.beta.iter().zip(expected_beta2).enumerate() {
            assert_close(*a, e, &format!("perm2 beta[{i}]"));
        }
        for (i, (a, e)) in fit2.s_post.iter().zip(expected_s_post2).enumerate() {
            assert_close(*a, e, &format!("perm2 s_post[{i}]"));
        }
    }

    // ---- s2_post confirmation: sqrt(LimmaModeration.s2_post) == Python
    // sqrt(fit_eb.s2_post) on the fixed full-data fit (oracle: limrots_oracle.py).
    // boot_ebayes_from_rows on the unchanged matrix is the full-data eBayes path,
    // so its s_post == sqrt(s2_post). ----
    #[test]
    fn s2_post_matches_python_full_data_oracle() {
        let mat = fixture();
        let fit = boot_ebayes_bootstrap(&mat, &[0, 1, 2, 3, 4, 5], 3);
        let expected_s_post = [
            0.14967183943258613,
            0.14967183943258613,
            0.07716474268929527,
            1.8732915428578778,
            0.023657886142284618,
            1.0501447642495674,
        ];
        for (i, (a, e)) in fit.s_post.iter().zip(expected_s_post).enumerate() {
            assert_close(*a, e, &format!("s_post == sqrt(s2_post)[{i}]"));
        }
    }

    // ---- p-value counting given a FIXED set of permuted d_p arrays, isolated
    // from the RNG (oracle: limrots_oracle.py). We hand a stub matrix and use the
    // count formula directly: count starts at 1, pools |d_p| >= |d_obs|, divides
    // by n_genes*(n_perm+1). Because pvalue_permutation regenerates d_p from the
    // matrix via the RNG, we replicate the counting in a small local helper that
    // matches accumulate_counts and assert against Python's pvalues. ----
    #[test]
    fn pvalue_counting_matches_python_oracle() {
        // d_obs and the three fixed d_p arrays from the oracle.
        let d_obs: [f64; 6] = [3.0, -1.0, 0.5, -2.5, 0.0, 1.2];
        let fixed_dps: [[f64; 6]; 3] = [
            [2.0, -1.0, 0.3, 1.0, 0.0, -0.5],
            [0.5, 1.5, -2.0, 3.5, -0.1, 0.2],
            [1.0, 4.0, -0.5, 0.0, 0.05, -1.3],
        ];
        let n_genes = d_obs.len();
        let n_perm = fixed_dps.len();
        let abs_d = d_obs.iter().map(|v| v.abs()).collect::<Vec<_>>();
        // Replicate the production counting (count starts at 1; pooled >=).
        let mut count = vec![1.0; n_genes];
        for dp in &fixed_dps {
            super::accumulate_counts(&mut count, dp, &abs_d);
        }
        let expected_count = [3.0, 10.0, 13.0, 3.0, 19.0, 7.0];
        for (i, (a, e)) in count.iter().zip(expected_count).enumerate() {
            assert_close(*a, e, &format!("count[{i}]"));
        }
        let denom = (n_genes * (n_perm + 1)) as f64;
        let pvalues = count
            .iter()
            .map(|v| (v / denom).min(1.0))
            .collect::<Vec<_>>();
        let expected_p = [
            0.125,
            0.4166666666666667,
            0.5416666666666666,
            0.125,
            0.7916666666666666,
            0.2916666666666667,
        ];
        for (i, (a, e)) in pvalues.iter().zip(expected_p).enumerate() {
            assert_close(*a, e, &format!("pvalue[{i}]"));
        }
    }

    // ---- pvalue_permutation end-to-end on the fixture: every p in [0,1], the
    // count-starts-at-1 / (n_perm+1) denominator and the min(.,1.0) clamp hold
    // (RNG-tolerant, structural). ----
    #[test]
    fn pvalue_permutation_is_in_range_and_clamped() {
        let mat = fixture();
        let mut rng = SplitMix64::new(super::PRNG_SEED);
        let d_obs = [2.0, -2.0, 0.0, -1.0, -3.0, 0.1];
        let p = pvalue_permutation(&d_obs, &mat, 3, 6, (0.05, 1.0), 10, &mut rng);
        assert_eq!(p.len(), 6);
        for (i, value) in p.iter().enumerate() {
            assert!(
                (0.0..=1.0).contains(value),
                "pvalue[{i}] out of range: {value}"
            );
            assert!(value.is_finite(), "pvalue[{i}] not finite");
        }
        // The denominator is n_genes*(n_perm+1) and count starts at 1, so the
        // smallest possible p is 1/(6*11) > 0.
        let floor = 1.0 / (6.0 * 11.0);
        for value in &p {
            assert!(*value + 1e-12 >= floor, "pvalue below 1/(n*(n_perm+1))");
        }
    }

    // Deterministic log2FC oracle for the 30x8, 4-vs-4 e2e fixture: limrots'
    // full-data eBayes beta == NEGATION of the rots E2E_LOG2FC (limrots uses the
    // [1,-1] contrast = mean_A - mean_B; rots uses nanmean(B)-nanmean(A)).
    // Captured from limrots_oracle.py.
    const E2E_LOG2FC: &[(usize, f64)] = &[
        (0, -2.3194578965825734),
        (1, -2.931353656903349),
        (2, -2.0977123988471753),
        (3, -0.9602528723815524),
        (4, -3.1187316602549),
        (5, -1.3732610141404162),
        (6, 0.6534347030345558),
        (7, 0.18229245477504286),
        (8, -0.08872443593928736),
        (9, 1.233392283025486),
        (10, 0.9518595362316393),
        (11, 0.32522117747048007),
        (12, 0.18440244273792494),
        (13, -0.044205566044979605),
        (14, -0.5677719771180847),
        (15, 0.022735918751076056),
        (16, 0.4190151627149028),
        (17, 0.8166528743519947),
        (18, 0.8351899914488765),
        (19, -0.6444003491502173),
        (20, -1.3867569584436854),
        (21, -0.8770610325497419),
        (22, 0.9413520718320889),
        (23, -1.0399964650885156),
        (24, 0.6658888175338671),
        (25, -1.4855932509659233),
        (26, 0.019778429193756608),
        (27, -1.6099247963736314),
        (28, 0.4258634750602255),
        (29, 0.44572816568048523),
    ];

    // Stochastic-wrapper property test. NOT bit-exact vs Python: LimROTS is RNG-
    // driven and Python's own output is seed-unstable (the selected a1 and the
    // adj_pvalues move between seeds), so we deliberately do not assert cell
    // equality against numpy's PCG64 stream. We DO assert:
    //   - log2FC cell-exact (it is the deterministic full-data eBayes coefficient,
    //     RNG-independent);
    //   - every p-value in [0,1] and finite;
    //   - d_stat finite for every protein;
    //   - the selected (a1, a2) lands on a valid grid point (a1 in the ssq grid,
    //     a2 == 1).
    #[test]
    fn limrots_two_group_property_holds_and_log2fc_is_cell_exact() {
        let rows = e2e_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = (0..30).map(|i| format!("P{i:02}")).collect::<Vec<_>>();
        let results = limrots_two_group_seeded(
            &proteins,
            &refs,
            4,
            4,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
            20,
            super::PRNG_SEED,
        );
        assert_eq!(results.len(), 30, "all proteins testable");

        // log2FC cell-exact vs Python (deterministic full-data eBayes beta).
        for &(idx, expected) in E2E_LOG2FC {
            let name = format!("P{idx:02}");
            let Some(row) = results.iter().find(|r| r.protein == name) else {
                panic!("{name} missing");
            };
            assert_close(row.log2_fold_change, expected, &format!("{name} log2FC"));
        }

        // p in [0,1], finite d_stat for every protein.
        for row in &results {
            assert!(
                (0.0..=1.0).contains(&row.p_value),
                "{} pvalue out of range: {}",
                row.protein,
                row.p_value
            );
            assert!(row.p_value.is_finite(), "{} pvalue not finite", row.protein);
            assert!(
                row.t_statistic.is_finite(),
                "{} d_stat not finite",
                row.protein
            );
            assert!(
                row.adj_p_value.is_finite(),
                "{} adj_pvalue not finite",
                row.protein
            );
        }
    }

    // The selected (a1, a2) must land on a valid grid point: a1 in the ssq grid
    // and a2 == 1 (there is NO raw-D branch in limrots, unlike rots). RNG-
    // tolerant: it checks the optimiser only ever returns grid coordinates.
    #[test]
    fn bootstrap_optimize_returns_grid_point() {
        let rows = e2e_rows();
        let ssq = build_ssq_grid();
        for seed in [super::PRNG_SEED, 1, 7, 2024] {
            let mut rng = SplitMix64::new(seed);
            let (a1, a2) = bootstrap_optimize(&rows, 4, 8, 20, &mut rng);
            assert!(
                (a2 - 1.0).abs() <= 1e-12,
                "seed {seed}: a2 must be 1.0 (no raw-D branch), got {a2}"
            );
            assert!(
                ssq.iter().any(|value| (value - a1).abs() <= 1e-12),
                "seed {seed}: a1={a1} not on the ssq grid"
            );
        }
    }

    // End-to-end public entry point runs, BH-sorts, and emits valid p-values.
    #[test]
    fn limrots_two_group_public_entry_runs() {
        let rows = e2e_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = (0..30).map(|i| format!("P{i:02}")).collect::<Vec<_>>();
        let results = limrots_two_group(
            &proteins,
            &refs,
            4,
            4,
            DEFAULT_FDR_THRESHOLD,
            DEFAULT_LOG2FC_THRESHOLD,
        );
        assert_eq!(results.len(), 30);
        // Adjusted p-values are sorted ascending.
        for window in results.windows(2) {
            assert!(
                window[0].adj_p_value <= window[1].adj_p_value || window[1].adj_p_value.is_nan(),
                "results not sorted by adj_pvalue"
            );
        }
        // Significance calls are well-formed.
        for row in &results {
            assert!(matches!(
                row.significance,
                Significance::Up | Significance::Down | Significance::Unchanged
            ));
        }
    }

    fn e2e_rows() -> Vec<Vec<f64>> {
        E2E_MATRIX.iter().map(|row| row.to_vec()).collect()
    }

    // The same seed-123 numpy matrix the rots e2e oracle uses (4-vs-4), so the
    // limrots full-data beta oracle (E2E_LOG2FC) feeds identical inputs.
    #[rustfmt::skip]
    const E2E_MATRIX: &[[f64; 8]] = &[
        [9.010878649652149, 9.632213348532117, 11.28792526128925, 10.193974419132614, 12.920230899639858, 12.577103791257251, 11.36353635362902, 12.541952220410293],
        [9.683404548834185, 9.677610883841039, 10.097167318670458, 8.47406959348105, 13.19216610410166, 11.328910324825891, 13.00026941965946, 12.136321123853117],
        [11.532033079628796, 9.34003058620818, 9.688205143530082, 10.337769126558825, 9.792528901800196, 12.827921441558736, 13.541630394690618, 13.126806793265029],
        [10.754769644312251, 9.854022106884775, 11.281902227059712, 11.074030621971943, 12.392620844577271, 12.005114312828983, 11.638233127839078, 10.769767804509556],
        [11.22622929282115, 7.827956113314818, 9.629852654147685, 10.164380069674667, 12.859881184612737, 13.761661236511811, 12.993323775951811, 11.70847857390156],
        [10.728127557889144, 8.738399683080303, 11.429938526688707, 9.843524675170595, 11.326240850012942, 11.360939899567795, 11.938638672379627, 11.607215077430057],
        [12.289909947314579, 9.281818852119404, 10.032607743156971, 10.028049895585639, 10.028272122739738, 10.05534586195271, 9.51843714181005, 9.41659249953588],
        [9.137839497928717, 8.51182538674841, 10.216306833109211, 10.984376350695877, 9.456915858987372, 9.441384960956215, 9.683517170677353, 9.539360258761096],
        [8.563730250975055, 11.365108032036918, 10.438999887557467, 9.288304972801967, 10.297171761539852, 9.561542727384268, 9.78836256662532, 10.363963831579115],
        [10.952964491974534, 11.519524129241278, 11.703909448949037, 9.751141292569057, 9.50025140886649, 10.099597501922018, 10.128343212288312, 9.26577810755514],
        [9.379524711765237, 10.813273720420849, 11.64180101374076, 9.773499151620827, 9.35203478900466, 9.71662879337581, 9.004868640029956, 9.727128230210692],
        [10.422444141467754, 9.918657038410645, 11.234577597061605, 10.150888032220434, 10.481119527334224, 9.851242467590508, 11.315665706560356, 8.777654397793428],
        [9.69640865969977, 8.826311324319809, 10.826273507011688, 10.850322289623024, 9.484232406986733, 11.658113318303489, 9.702737404360747, 8.616622880051628],
        [9.718795495163823, 10.360020509941075, 9.765607985072931, 12.265520599867209, 10.855386650779941, 11.73127943872283, 11.385885443090944, 8.314215321631243],
        [9.62224139333485, 7.271514313127154, 9.35360326409304, 11.115104255783713, 9.156788904060573, 9.363311562033822, 10.326134167526362, 10.787316501190338],
        [9.643686615867946, 9.744187824559853, 10.808981702502642, 10.254058684864173, 9.714713228886096, 10.314503492254888, 10.076121508150356, 10.254632923498972],
        [12.001231181802662, 9.694792076324612, 9.460237142053542, 11.41363079410976, 9.29429935592141, 11.719889475677249, 9.805802193788042, 10.073839518044268],
        [10.733628281140872, 11.2131096954131, 10.996998651667894, 9.79636464474424, 9.633571164502154, 10.347472646054985, 10.13471907223136, 9.357726892769634],
        [10.410728074840504, 10.994119620161474, 10.166506704973866, 11.563999717050974, 10.410302129274417, 9.844186999193733, 10.100214889931419, 9.439890132831742],
        [8.968702031873006, 10.476528511333012, 9.146749682439074, 9.641708027670198, 9.91079994629505, 10.385104174182407, 10.46116198661093, 10.054223542827762],
        [8.06869234595192, 9.79340730165992, 9.22446059009289, 9.88865244793397, 10.991326162670685, 10.80482429127482, 10.220330268055328, 10.50575979741261],
        [8.51159179926529, 10.194867465692276, 10.624510778948961, 9.293872713913213, 10.404987601909117, 9.619173748690965, 11.000506751356722, 11.108418786061897],
        [9.465216906342956, 11.462922589144863, 9.809919703234534, 12.061944558167134, 9.508445785135853, 10.063756314040551, 10.213189418250687, 9.249203952134037],
        [8.274620471728936, 9.200396053287452, 11.079605669439184, 9.671333392473215, 9.673293393034323, 11.516801903442548, 10.493537792192598, 10.702308358613383],
        [10.849448659386821, 9.092582030833373, 10.12472619119045, 10.150216891456425, 9.841133386939504, 8.906968190360423, 10.463779755206073, 8.341537170225596],
        [9.062933234818388, 9.190661855274174, 9.587868312306712, 10.841092596271773, 11.656680431754527, 11.72229532456346, 10.80625937846645, 10.439693867750302],
        [7.663451133175238, 11.130119646001527, 12.628946566875129, 10.452876432658636, 10.234039314763036, 10.698695826852546, 10.728628585472741, 10.134916334847182],
        [9.626184133953055, 10.361759209669309, 8.628210790454695, 8.11757345670472, 11.360938640721507, 9.541938669664201, 10.715955476440959, 11.554593989449632],
        [11.220280349607807, 9.446607291342193, 10.324754029510578, 8.83391302200727, 9.895426898582674, 11.04601695100833, 10.032026240362665, 7.1486307022732785],
        [7.748736592956661, 9.198414429430384, 11.333241659298887, 10.199085261163354, 10.141273124012217, 8.684032816429761, 10.146431651468244, 7.724827688217122],
    ];
}
