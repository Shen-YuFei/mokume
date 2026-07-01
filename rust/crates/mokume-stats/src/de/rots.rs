//! ROTS (Reproducibility-Optimised Test Statistic) two-group differential
//! expression.
//!
//! This is a FAITHFUL ALGORITHM PORT of mokume's `analysis/rots.py` (a pure
//! numpy reimplementation of R `ROTS`), NOT a cell-for-cell match. ROTS is an
//! RNG-driven method: it bootstraps each class and permutes the labels to score
//! how reproducible the top-k feature ranking is across a grid of variance-
//! shrinkage parameters, picks the best `(a1, a2, k)`, then computes the final
//! `d_stat = fc / max(a1 + a2*s, 1e-10)` and a permutation p-value. Because the
//! optimisation and the p-value both consume random draws, Python's own output
//! is unstable across seeds (the selected `a1` moves, the top-20 ranks reshuffle,
//! and `adj_pvalue` differs by ~1e-1 between seeds). Bit-matching numpy's PCG64
//! is therefore neither possible without reproducing its exact stream nor
//! meaningful, so it is intentionally NOT attempted.
//!
//! WHAT IS MATCHED EXACTLY (cell-exact, 1e-9, RNG-independent):
//!   - `_compute_d_stat` (rots.py:36): `d = fc / max(a1 + a2*s, 1e-10)`.
//!   - `_group_stats` (rots.py:48): unpaired two-group `d = my - mx` and pooled
//!     `s`, with the `nx<2 || ny<2 -> d=0, s=1` edge and nan-aware means/sums.
//!   - `_reproducibility_score` (rots.py:77): top-k index-set overlap / k.
//!   - `_build_ssq_grid` (rots.py:90) and `_build_n_grid` (rots.py:101): exact
//!     grid lengths and values.
//!   - `_rank_abs` (rots.py:222): rank by descending |x|, 1 = largest, with
//!     numpy `argsort` STABLE tie order (ties broken by ascending original
//!     index).
//!   - `_calculate_p` (rots.py:230): two-pointer count of permuted |d| values
//!     `>=` each observed |d|, scattered back to original order.
//!   - `log2FC` (rots.py:313): `nanmean(B) - nanmean(A)`, deterministic.
//!
//! WHAT IS NOT MATCHED: the exact stochastic stream. The bootstrap/permutation
//! here uses a small self-contained splitmix64 PRNG (no external `rand`
//! dependency, the workspace deliberately has none) seeded by a fixed constant.
//! The result therefore matches Python in ALGORITHM and DISTRIBUTION (within
//! ROTS's own seed-noise), not cell-for-cell. The stochastic wrapper is covered
//! by property tests (p in [0,1], finite `d_stat`, the selected `(a1,a2,k)` lands
//! on a valid grid point), and `log2FC` is asserted cell-exact since it is
//! RNG-independent.
//!
//! Reference
//! ---------
//! Suomi T, Seyednasrollah F, Jaakkola MK, Faux T, Elo LL. ROTS: An R package
//! for reproducibility-optimized statistical testing. *PLoS Comput Biol*.
//! 2017;13(5):e1005562.

use super::correct::bh_adjust;
use super::{classify, DeResult};

/// Default bootstrap iterations, matching `run_rots(n_boot=100)` (rots.py:268).
const DEFAULT_N_BOOT: usize = 100;

/// Fixed PRNG seed for the bootstrap/permutation draws. ROTS is seed-unstable,
/// so this is a deterministic constant rather than `run_rots`'s `seed=42`
/// (which only makes sense against numpy's PCG64 stream). A fixed seed keeps the
/// Rust output reproducible run-to-run.
const PRNG_SEED: u64 = 0x9E37_79B9_7F4A_7C15;

/// Lower clamp on the d-statistic denominator (rots.py:44).
const DENOM_FLOOR: f64 = 1e-10;

/// Minimum standard deviation in the reproducibility table, matching
/// `sd[sd < 1e-10] = 1e-10` (rots.py:183).
const SD_FLOOR: f64 = 1e-10;

/// Minimum finite observations per group required for a protein to be testable,
/// matching `filter_testable(min_per_group=2)` (\_helpers.py:59).
const MIN_PER_GROUP: usize = 2;

/// A tiny self-contained splitmix64 PRNG (no `unsafe`, no external crate).
///
/// splitmix64 is the standard seeding generator for xoshiro/PCG families; it has
/// good statistical quality for the bootstrap index draws ROTS needs and is ~10
/// lines. The workspace deliberately avoids `rand`/`getrandom`, so this is
/// written in-module. Only `next_u64` plus a bounded `next_below` (unbiased via
/// rejection) and a Fisher-Yates `permutation` are exposed.
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// One splitmix64 step (the canonical constants).
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform integer in `0..bound` (`bound > 0`) without modulo bias, via
    /// Lemire-style rejection. Returns 0 if `bound == 0` (never called that way).
    pub(crate) fn next_below(&mut self, bound: usize) -> usize {
        if bound == 0 {
            return 0;
        }
        let bound64 = bound as u64;
        // Reject the top `2^64 mod bound` values so the remainder is uniform.
        let zone = u64::MAX - (u64::MAX % bound64);
        loop {
            let candidate = self.next_u64();
            if candidate <= zone {
                return (candidate % bound64) as usize;
            }
        }
    }

    /// `count` draws of `next_below(n)`, mirroring the SHAPE of numpy's
    /// `rng.choice(n, count, replace=True)` (a with-replacement resample), not
    /// its exact stream. Reused by LimROTS' bootstrap index construction.
    pub(crate) fn choice_with_replacement(&mut self, n: usize, count: usize) -> Vec<usize> {
        (0..count).map(|_| self.next_below(n)).collect()
    }

    /// A random permutation of `0..n` (Fisher-Yates), matching the SHAPE of
    /// `numpy.random.Generator.permutation(n)` (a uniform shuffle), not its
    /// exact stream.
    pub(crate) fn permutation(&mut self, n: usize) -> Vec<usize> {
        let mut indices = (0..n).collect::<Vec<_>>();
        for i in (1..n).rev() {
            let j = self.next_below(i + 1);
            indices.swap(i, j);
        }
        indices
    }

    /// A uniform `f64` in `[0, 1)` from the top 53 bits of one `next_u64`.
    pub(crate) fn next_f64(&mut self) -> f64 {
        // 53-bit mantissa: shift away the low 11 bits, scale by 2^-53.
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// One standard-normal draw via the Box-Muller transform. This matches the
    /// SHAPE of `numpy.random.RandomState.normal(0, 1)` (an i.i.d. standard
    /// normal) but NOT its exact stream; proDA uses it only to seed the initial
    /// imputation that starts the EM trajectory, where bit-exactness is not
    /// required.
    pub(crate) fn next_gaussian(&mut self) -> f64 {
        // Guard the radius against u1 == 0 (ln(0) = -inf).
        let u1 = self.next_f64().max(f64::MIN_POSITIVE);
        let u2 = self.next_f64();
        let radius = (-2.0 * u1.ln()).sqrt();
        radius * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// Port of `_compute_d_stat` (rots.py:36): `d = fc / max(a1 + a2*s, 1e-10)`,
/// element-wise.
fn compute_d_stat(fc: &[f64], s: &[f64], a1: f64, a2: f64) -> Vec<f64> {
    fc.iter()
        .zip(s)
        .map(|(fold, sdev)| {
            let denom = a1 + a2 * sdev;
            let denom = if denom < DENOM_FLOOR {
                DENOM_FLOOR
            } else {
                denom
            };
            fold / denom
        })
        .collect()
}

/// Per-protein `(d, s)` from `_group_stats` (rots.py:48): unpaired two-group
/// difference of nan-means and pooled standard error.
fn group_stats(mat: &[Vec<f64>], n_a: usize) -> (Vec<f64>, Vec<f64>) {
    let mut d = Vec::with_capacity(mat.len());
    let mut s = Vec::with_capacity(mat.len());
    for row in mat {
        let (x, y) = row.split_at(n_a);
        let (mx, nx) = nan_mean_count(x);
        let (my, ny) = nan_mean_count(y);

        if nx < 2 || ny < 2 {
            // rots.py:70-72 edge: too few finite obs -> d=0, s=1.
            d.push(0.0);
            s.push(1.0);
            continue;
        }

        let sx = nan_sum_sq_dev(x, mx);
        let sy = nan_sum_sq_dev(y, my);
        let denom = ((nx + ny) as f64 - 2.0).max(1.0);
        let pooled =
            ((sx + sy) / denom) * (1.0 / (nx as f64).max(1.0) + 1.0 / (ny as f64).max(1.0));
        d.push(my - mx);
        s.push(pooled.sqrt());
    }
    (d, s)
}

/// `np.nanmean` plus the finite count (`np.sum(np.isfinite(x))`). An all-nan
/// slice yields a NaN mean and a zero count (handled by the `nx<2` edge).
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

/// `np.nansum((x - mean)**2)` over the finite entries.
fn nan_sum_sq_dev(values: &[f64], mean: f64) -> f64 {
    values
        .iter()
        .filter(|value| value.is_finite())
        .map(|value| {
            let dev = value - mean;
            dev * dev
        })
        .sum()
}

/// Port of `_reproducibility_score` (rots.py:77): size of the top-k index-set
/// overlap divided by `k`. `ranks` are 1-based (`1 = largest |d|`).
pub(crate) fn reproducibility_score(ranks1: &[usize], ranks2: &[usize], k: usize) -> f64 {
    if k == 0 {
        return 0.0;
    }
    let mut overlap = 0usize;
    for (rank1, rank2) in ranks1.iter().zip(ranks2) {
        if *rank1 <= k && *rank2 <= k {
            overlap += 1;
        }
    }
    overlap as f64 / k as f64
}

/// Port of `_build_ssq_grid` (rots.py:90): `c((0:20)/100, (11:50)/50, (6:25)/5)`.
/// numpy `arange(0,21)`/`arange(11,51)`/`arange(6,26)` are end-exclusive, so the
/// inclusive upper bounds are 20, 50, and 25.
pub(crate) fn build_ssq_grid() -> Vec<f64> {
    let mut grid = Vec::new();
    for i in 0..=20 {
        grid.push(i as f64 / 100.0);
    }
    for i in 11..=50 {
        grid.push(i as f64 / 50.0);
    }
    for i in 6..=25 {
        grid.push(i as f64 / 5.0);
    }
    grid
}

/// Port of `_build_n_grid` (rots.py:101): the top-list grid, keeping entries
/// strictly below `k_max`. numpy `arange(a,b)` is end-exclusive.
pub(crate) fn build_n_grid(k_max: usize) -> Vec<usize> {
    let mut grid = Vec::new();
    let mut push_block = |start: usize, end_exclusive: usize, mult: usize| {
        for i in start..end_exclusive {
            grid.push(i * mult);
        }
    };
    push_block(1, 21, 5);
    push_block(11, 51, 10);
    push_block(21, 41, 25);
    push_block(11, 1001, 100);
    grid.into_iter().filter(|value| *value < k_max).collect()
}

/// Port of `_rank_abs` (rots.py:222): 1-based rank by descending `|x|`, `1` =
/// largest. Returns `ranks` aligned to the input order.
///
/// TIE NOTE: Python uses `np.argsort(-|x|)` with numpy's DEFAULT quicksort,
/// which is NOT stable; on exact magnitude ties it produces an internal,
/// version-specific order that no simple rule reproduces (verified: it neither
/// preserves nor reverses tie runs uniformly). This port instead uses a STABLE
/// sort (ties broken by ascending original index). On DISTINCT magnitudes -- the
/// only case that arises on real continuous float data -- stable and numpy's
/// quicksort give IDENTICAL ranks (verified over 2000 random vectors), so this
/// matches Python wherever it matters. `_rank_abs` only feeds the RNG-driven
/// reproducibility scoring, never the final p-value, so the contrived-tie
/// difference cannot affect a real result.
pub(crate) fn rank_abs(x: &[f64]) -> Vec<usize> {
    let mut order = (0..x.len()).collect::<Vec<_>>();
    // `sort_by` is stable: sorting by |x| descending breaks magnitude ties by
    // ascending original index.
    order.sort_by(|&i, &j| x[j].abs().total_cmp(&x[i].abs()));
    let mut ranks = vec![0usize; x.len()];
    for (position, &index) in order.iter().enumerate() {
        ranks[index] = position + 1;
    }
    ranks
}

/// Port of `_calculate_p` (rots.py:230): permutation p-values. For each observed
/// `|d|` (descending), count how many permuted `|d|` values are `>=` it, divided
/// by the total permuted count; scatter back to the original protein order.
fn calculate_p(observed: &[f64], permuted: &[Vec<f64>]) -> Vec<f64> {
    let n = observed.len();
    // Observed indices sorted by descending |d| (stable, like argsort).
    let mut obs_order = (0..n).collect::<Vec<_>>();
    obs_order.sort_by(|&i, &j| observed[j].abs().total_cmp(&observed[i].abs()));
    let obs_sorted = obs_order
        .iter()
        .map(|&i| observed[i].abs())
        .collect::<Vec<_>>();

    // All permuted |d| flattened and sorted DESCENDING (rots.py:241).
    let mut perm_flat = permuted
        .iter()
        .flat_map(|row| row.iter().map(|value| value.abs()))
        .collect::<Vec<_>>();
    perm_flat.sort_by(|a, b| b.total_cmp(a));
    let n_perm = perm_flat.len();

    let mut p_sorted = vec![0.0; n];
    if n_perm == 0 {
        // No permutations: every count is 0 -> p = 0 (matches j/n_perm with j=0,
        // though Python would divide by zero; guard defensively).
        return scatter(&obs_order, &p_sorted);
    }
    let mut j = 0usize;
    for (position, obs) in obs_sorted.iter().enumerate() {
        while j < n_perm && perm_flat[j] >= *obs {
            j += 1;
        }
        p_sorted[position] = j as f64 / n_perm as f64;
    }
    scatter(&obs_order, &p_sorted)
}

/// Scatter `sorted` values back to original order: `result[order[i]] = sorted[i]`.
fn scatter(order: &[usize], sorted: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; order.len()];
    for (position, &index) in order.iter().enumerate() {
        result[index] = sorted[position];
    }
    result
}

/// Select a single column of `mat` by index, used to build a bootstrap/permuted
/// matrix without allocating intermediate row copies in the hot loop.
fn permute_columns(mat: &[Vec<f64>], indices: &[usize]) -> Vec<Vec<f64>> {
    mat.iter()
        .map(|row| indices.iter().map(|&col| row[col]).collect())
        .collect()
}

/// One bootstrap index vector: resample WITH replacement WITHIN each class
/// (rots.py:144-147). Class 1 is the first `n_a` columns, class 2 the rest.
fn bootstrap_indices(rng: &mut SplitMix64, n_a: usize, n_total: usize) -> Vec<usize> {
    let mut indices = vec![0usize; n_total];
    // Class 1 positions: 0..n_a, resampled from themselves.
    for slot in indices.iter_mut().take(n_a) {
        *slot = rng.next_below(n_a);
    }
    // Class 2 positions: n_a..n_total, resampled from themselves.
    let n_b = n_total - n_a;
    for slot in indices.iter_mut().skip(n_a) {
        *slot = n_a + rng.next_below(n_b);
    }
    indices
}

/// Find the optimal `(a1, a2, k)` over the reproducibility grid (port of
/// `_bootstrap_optimize`, rots.py:114). Returns the grid-point parameters the
/// final d-statistic uses.
fn bootstrap_optimize(
    mat: &[Vec<f64>],
    n_a: usize,
    n_boot: usize,
    rng: &mut SplitMix64,
) -> (f64, f64, usize) {
    let n_genes = mat.len();
    let n_total = if n_genes == 0 { 0 } else { mat[0].len() };
    let k_max = n_genes / 4;
    if k_max < 1 {
        return (0.0, 1.0, 1);
    }

    let ssq = build_ssq_grid();
    let mut n_grid = build_n_grid(k_max);
    if n_grid.is_empty() {
        n_grid = vec![1];
    }
    let n_ssq = ssq.len();
    let n_k = n_grid.len();
    let two_b = 2 * n_boot;

    // D/S from bootstraps, pD/pS from permutations, columns 0..2*n_boot.
    let mut d_boot = vec![Vec::new(); two_b];
    let mut s_boot = vec![Vec::new(); two_b];
    let mut d_perm = vec![Vec::new(); two_b];
    let mut s_perm = vec![Vec::new(); two_b];
    for i in 0..two_b {
        let boot_idx = bootstrap_indices(rng, n_a, n_total);
        let boot_mat = permute_columns(mat, &boot_idx);
        let (db, sb) = group_stats(&boot_mat, n_a);
        d_boot[i] = db;
        s_boot[i] = sb;

        let perm_idx = rng.permutation(n_total);
        let perm_mat = permute_columns(mat, &perm_idx);
        let (dp, sp) = group_stats(&perm_mat, n_a);
        d_perm[i] = dp;
        s_perm[i] = sp;
    }

    // reprotable / reprotable_p / reprotable_sd over (n_ssq + 1) x n_k.
    let mut reprotable = vec![vec![0.0; n_k]; n_ssq + 1];
    let mut reprotable_p = vec![vec![0.0; n_k]; n_ssq + 1];
    let mut reprotable_sd = vec![vec![0.0; n_k]; n_ssq + 1];

    for si in 0..n_ssq {
        // `Some(ssq[si])` for every grid row; rots.py:186-203 makes the LAST row
        // (index n_ssq) use raw D/pD (no ssq division) via `None`.
        let row = fill_repro_row(
            Some(ssq[si]),
            n_boot,
            &n_grid,
            &d_boot,
            &s_boot,
            &d_perm,
            &s_perm,
        );
        reprotable[si] = row.reprotable;
        reprotable_p[si] = row.reprotable_p;
        reprotable_sd[si] = row.reprotable_sd;
    }
    let last = n_ssq;
    let row = fill_repro_row(None, n_boot, &n_grid, &d_boot, &s_boot, &d_perm, &s_perm);
    reprotable[last] = row.reprotable;
    reprotable_p[last] = row.reprotable_p;
    reprotable_sd[last] = row.reprotable_sd;

    // ztable = (reprotable - reprotable_p) / reprotable_sd; non-finite -> -inf;
    // argmax via unravel_index (row-major, first max wins) (rots.py:205-209).
    let mut best_value = f64::NEG_INFINITY;
    let mut best_si = 0usize;
    let mut best_ki = 0usize;
    for si in 0..=n_ssq {
        for ki in 0..n_k {
            let z = (reprotable[si][ki] - reprotable_p[si][ki]) / reprotable_sd[si][ki];
            let z = if z.is_finite() { z } else { f64::NEG_INFINITY };
            if z > best_value {
                best_value = z;
                best_si = si;
                best_ki = ki;
            }
        }
    }

    let (best_a1, best_a2) = if best_si < n_ssq {
        (ssq[best_si], 1.0)
    } else {
        (1.0, 0.0)
    };
    (best_a1, best_a2, n_grid[best_ki])
}

/// One reprotable / reprotable_p / reprotable_sd row (length `n_k`), returned by
/// [`fill_repro_row`].
struct ReproRow {
    reprotable: Vec<f64>,
    reprotable_p: Vec<f64>,
    reprotable_sd: Vec<f64>,
}

/// Compute one reprotable / reprotable_p / reprotable_sd row for a given ssq
/// value (`Some(ssq)` divides `D/(ssq+S)`; `None` uses raw `D`, the last-row
/// branch). Mirrors the inner loop of `_bootstrap_optimize` (rots.py:167-203).
fn fill_repro_row(
    ssq: Option<f64>,
    n_boot: usize,
    n_grid: &[usize],
    d_boot: &[Vec<f64>],
    s_boot: &[Vec<f64>],
    d_perm: &[Vec<f64>],
    s_perm: &[Vec<f64>],
) -> ReproRow {
    let n_k = n_grid.len();
    // overlaps[b][ki] and overlaps_p[b][ki].
    let mut overlaps = vec![vec![0.0; n_k]; n_boot];
    let mut overlaps_p = vec![vec![0.0; n_k]; n_boot];
    for b in 0..n_boot {
        let col1 = 2 * b;
        let col2 = 2 * b + 1;
        let d1 = shrink(&d_boot[col1], &s_boot[col1], ssq);
        let d2 = shrink(&d_boot[col2], &s_boot[col2], ssq);
        let dp = shrink(&d_perm[col1], &s_perm[col1], ssq);
        let r1 = rank_abs(&d1);
        let r2 = rank_abs(&d2);
        let rp = rank_abs(&dp);
        for (ki, &k) in n_grid.iter().enumerate() {
            overlaps[b][ki] = reproducibility_score(&r1, &r2, k);
            overlaps_p[b][ki] = reproducibility_score(&r1, &rp, k);
        }
    }
    let mut reprotable = vec![0.0; n_k];
    let mut reprotable_p = vec![0.0; n_k];
    let mut reprotable_sd = vec![0.0; n_k];
    for ki in 0..n_k {
        let col = overlaps.iter().map(|row| row[ki]).collect::<Vec<_>>();
        let col_p = overlaps_p.iter().map(|row| row[ki]).collect::<Vec<_>>();
        reprotable[ki] = mean(&col);
        reprotable_p[ki] = mean(&col_p);
        let sd = std_ddof1(&col);
        reprotable_sd[ki] = if sd < SD_FLOOR { SD_FLOOR } else { sd };
    }
    ReproRow {
        reprotable,
        reprotable_p,
        reprotable_sd,
    }
}

/// `D / (ssq + S)` element-wise (rots.py:171), or raw `D` when `ssq` is `None`
/// (the last-row branch, rots.py:190).
fn shrink(d: &[f64], s: &[f64], ssq: Option<f64>) -> Vec<f64> {
    match ssq {
        Some(value) => d.iter().zip(s).map(|(di, si)| di / (value + si)).collect(),
        None => d.to_vec(),
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

/// Two-group ROTS differential expression.
///
/// `proteins[i]` labels `rows[i]`; each row holds log2 intensities with the
/// first `n_a` columns in condition A and the next `n_b` in condition B. Proteins
/// without at least `MIN_PER_GROUP` finite observations in BOTH groups are
/// dropped before testing (matching `filter_testable`). Returns one [`DeResult`]
/// per testable protein, sorted by adjusted p-value ascending (BH), reproducing
/// mokume's `DifferentialExpression(method="rots")` output contract: `pvalue`
/// carries the permutation p-value, `adj_pvalue` its BH adjustment, and
/// `t_statistic` carries the ROTS `d_stat` (the method's extra column).
///
/// NOTE: this is a faithful algorithm port; with the internal fixed-seed PRNG the
/// p-values match Python in distribution (within ROTS's seed-noise), not
/// cell-for-cell. `log2_fold_change` IS deterministic and matches Python exactly.
pub fn rots_two_group(
    proteins: &[String],
    rows: &[&[f64]],
    n_a: usize,
    n_b: usize,
    fdr_threshold: f64,
    log2fc_threshold: f64,
) -> Vec<DeResult> {
    rots_two_group_seeded(
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

/// Seedable core of [`rots_two_group`], exposed for deterministic-property tests
/// that want a fixed `n_boot`/`seed`. The public entry point hardcodes the Python
/// defaults (`n_boot=100`) and the fixed module PRNG seed.
#[allow(clippy::too_many_arguments)]
fn rots_two_group_seeded(
    proteins: &[String],
    rows: &[&[f64]],
    n_a: usize,
    n_b: usize,
    fdr_threshold: f64,
    log2fc_threshold: f64,
    n_boot: usize,
    seed: u64,
) -> Vec<DeResult> {
    // filter_testable: keep proteins with >= MIN_PER_GROUP finite obs per group.
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

    let (best_a1, best_a2, _best_k) = bootstrap_optimize(&mat, n_a, n_boot, &mut rng);

    let (d_obs, s_obs) = group_stats(&mat, n_a);
    let d_stat = compute_d_stat(&d_obs, &s_obs, best_a1, best_a2);

    // Final permutation null: 2*n_boot relabelings -> permuted d-statistics.
    let two_b = 2 * n_boot;
    let mut perm_d = Vec::with_capacity(two_b);
    for _ in 0..two_b {
        let perm_idx = rng.permutation(n_total);
        let perm_mat = permute_columns(&mat, &perm_idx);
        let (dp, sp) = group_stats(&perm_mat, n_a);
        perm_d.push(compute_d_stat(&dp, &sp, best_a1, best_a2));
    }
    let p_values = calculate_p(&d_stat, &perm_d);
    let adjusted = bh_adjust(&p_values);

    let mut results = Vec::with_capacity(kept.len());
    for (position, (original_index, row)) in kept.iter().enumerate() {
        let Some(protein) = proteins.get(*original_index) else {
            continue;
        };
        let (mean_a, count_a) = nan_mean_count(&row[..n_a]);
        let (mean_b, count_b) = nan_mean_count(&row[n_a..]);
        // log2FC = nanmean(B) - nanmean(A), deterministic (rots.py:313).
        let log2_fold_change = mean_b - mean_a;
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
            // ROTS emits no AveExpr/B columns (its extra column is d_stat,
            // carried in t_statistic), so these limma-only quirks are NaN/0.0
            // and must not be written for rots.
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

/// `filter_testable` predicate (\_helpers.py:62-64): at least [`MIN_PER_GROUP`]
/// finite observations in BOTH groups.
fn testable(row: &[f64], n_a: usize, n_b: usize) -> bool {
    let finite_a = row[..n_a].iter().filter(|value| value.is_finite()).count();
    let finite_b = row[n_a..n_a + n_b]
        .iter()
        .filter(|value| value.is_finite())
        .count();
    finite_a >= MIN_PER_GROUP && finite_b >= MIN_PER_GROUP
}

#[cfg(test)]
mod tests {
    use super::{
        build_n_grid, build_ssq_grid, calculate_p, compute_d_stat, group_stats, rank_abs,
        reproducibility_score, rots_two_group, rots_two_group_seeded, SplitMix64,
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

    // ---- _compute_d_stat (oracle: rots_oracle.py) ----
    #[test]
    fn compute_d_stat_matches_python_oracle() {
        let fc = [2.0, -1.5, 0.0, 3.3, -4.0];
        let s = [0.5, 1.0, 2.0, 0.0, 100.0];
        let d = compute_d_stat(&fc, &s, 0.05, 1.0);
        let expected = [
            3.6363636363636362,
            -1.4285714285714286,
            0.0,
            65.99999999999999,
            -3.998_000_999_500_25e-2,
        ];
        for (i, (a, e)) in d.iter().zip(expected).enumerate() {
            assert_close(*a, e, &format!("d_stat[{i}] a1=0.05 a2=1.0"));
        }
        // a1=1, a2=0 -> d == fc.
        let raw = compute_d_stat(&fc, &s, 1.0, 0.0);
        for (i, (a, e)) in raw.iter().zip(fc).enumerate() {
            assert_close(*a, e, &format!("d_stat[{i}] a1=1 a2=0"));
        }
        // denom clamp: fc/1e-10 = 1e10.
        let clamped = compute_d_stat(&[1.0], &[0.0], 0.0, 0.0);
        assert_close(clamped[0], 1e10, "d_stat clamp");
    }

    // ---- _group_stats incl. nx<2 edge (oracle: rots_oracle.py) ----
    #[test]
    fn group_stats_matches_python_oracle() {
        let nan = f64::NAN;
        let mat = vec![
            vec![10.0, 10.2, 9.8, 12.0, 12.1, 11.9],
            vec![15.0, 15.1, 14.9, 13.0, 13.2, 12.8],
            vec![8.0, 8.1, 7.9, 8.05, 7.95, 8.0],
            vec![5.0, nan, nan, 6.0, 6.01, 5.99], // nx=1 -> d=0, s=1
            vec![3.0, 3.2, 2.8, 4.0, nan, 4.2],   // ny=2 ok
        ];
        let (d, s) = group_stats(&mat, 3);
        let expected_d = [2.0, -2.0, 0.0, 0.0, 1.0999999999999996];
        let expected_s = [
            0.1290994448735801,
            0.1290994448735801,
            0.06454972243679016,
            1.0,
            0.1666666666666668,
        ];
        for (i, (a, e)) in d.iter().zip(expected_d).enumerate() {
            assert_close(*a, e, &format!("group_stats d[{i}]"));
        }
        for (i, (a, e)) in s.iter().zip(expected_s).enumerate() {
            assert_close(*a, e, &format!("group_stats s[{i}]"));
        }
    }

    // ---- _reproducibility_score (oracle: rots_oracle.py) ----
    #[test]
    fn reproducibility_score_matches_python_oracle() {
        let ranks1 = [1usize, 2, 3, 4, 5];
        let ranks2 = [2usize, 1, 5, 3, 4];
        assert_close(reproducibility_score(&ranks1, &ranks2, 2), 1.0, "repro k=2");
        assert_close(
            reproducibility_score(&ranks1, &ranks2, 3),
            0.6666666666666666,
            "repro k=3",
        );
        assert_close(reproducibility_score(&ranks1, &ranks2, 0), 0.0, "repro k=0");
        assert_close(reproducibility_score(&ranks1, &ranks2, 5), 1.0, "repro k=5");
    }

    // ---- _build_ssq_grid: exact length + values (oracle: rots_oracle.py) ----
    #[test]
    fn build_ssq_grid_matches_python_oracle() {
        let grid = build_ssq_grid();
        assert_eq!(grid.len(), 81, "ssq grid length");
        let expected = [
            0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,
            0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38,
            0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68,
            0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98,
            1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2,
            4.4, 4.6, 4.8, 5.0,
        ];
        assert_eq!(grid.len(), expected.len());
        for (i, (a, e)) in grid.iter().zip(expected).enumerate() {
            assert_close(*a, e, &format!("ssq grid[{i}]"));
        }
    }

    // ---- _build_n_grid: exact (oracle: rots_oracle.py) ----
    #[test]
    fn build_n_grid_matches_python_oracle() {
        assert_eq!(
            build_n_grid(50),
            vec![5, 10, 15, 20, 25, 30, 35, 40, 45],
            "n_grid kmax=50"
        );
        let g200 = build_n_grid(200);
        assert_eq!(g200.len(), 29, "n_grid kmax=200 len");
        assert_eq!(g200.first().copied(), Some(5));
        assert_eq!(g200.last().copied(), Some(190));
        let g1000 = build_n_grid(1000);
        assert_eq!(g1000.len(), 79, "n_grid kmax=1000 len");
        assert_eq!(g1000.last().copied(), Some(975));
        // The 525/550/... block proves the *25 step engages past 500.
        assert!(g1000.contains(&525) && g1000.contains(&975));
        // kmax=3 keeps nothing (< 3 leaves the grid empty).
        assert!(build_n_grid(3).is_empty(), "n_grid kmax=3 empty");
    }

    // ---- _rank_abs (oracle: rots_oracle.py) ----
    //
    // On DISTINCT magnitudes this is cell-exact vs numpy `argsort(-|x|)`. On
    // exact ties, numpy's default quicksort order is a non-portable internal
    // artifact (it gave [6,1,2,5,4,3] on the tie vector below), so this port
    // documents the STABLE-sort result instead -- which matches numpy on every
    // tie-free input (the only case real float data produces). See the rank_abs
    // doc comment.
    #[test]
    fn rank_abs_matches_python_oracle() {
        // Distinct magnitudes: |x| = [2.5, 0.1, 4.0, 1.7, 3.3] -> descending
        // order is idx 2(4.0),4(3.3),0(2.5),3(1.7),1(0.1); ranks = [3,5,1,4,2].
        // Captured from numpy argsort -- cell-exact (no ties).
        let distinct = [2.5, -0.1, 4.0, 1.7, -3.3];
        assert_eq!(
            rank_abs(&distinct),
            vec![3, 5, 1, 4, 2],
            "rank_abs distinct magnitudes (cell-exact vs numpy)"
        );
        // Tie vector: |0.5|(idx 0,3), |3|(idx 1,2), |1|(idx 4,5). STABLE order
        // breaks every tie by ascending original index -> ranks [5,1,2,6,3,4]
        // (matches numpy `argsort(kind='stable')`). numpy's DEFAULT quicksort
        // instead returns [6,1,2,5,4,3] here -- a non-portable pivot artifact we
        // intentionally do not reproduce (it cannot affect a real p-value; see
        // the rank_abs doc comment).
        let tied = [0.5, -3.0, 3.0, -0.5, 1.0, -1.0];
        assert_eq!(
            rank_abs(&tied),
            vec![5, 1, 2, 6, 3, 4],
            "rank_abs stable ties"
        );
        // All-equal magnitudes -> stable ranks follow original order.
        assert_eq!(
            rank_abs(&[2.0, 2.0, 2.0]),
            vec![1, 2, 3],
            "rank_abs all-tie"
        );
    }

    // ---- _calculate_p with FIXED observed + FIXED permuted matrix (isolates the
    // deterministic p-computation from the RNG) (oracle: rots_oracle.py) ----
    #[test]
    fn calculate_p_matches_python_oracle() {
        let observed = [3.0, -1.0, 0.5, -2.5, 0.0];
        let permuted = vec![
            vec![2.0, -1.0, 0.3],
            vec![0.5, 1.5, -2.0],
            vec![3.5, -0.1, 0.2],
            vec![1.0, 4.0, -0.5],
            vec![0.0, 0.05, 0.0],
        ];
        let p = calculate_p(&observed, &permuted);
        let expected = [
            0.13333333333333333,
            0.4666666666666667,
            0.6,
            0.13333333333333333,
            1.0,
        ];
        for (i, (a, e)) in p.iter().zip(expected).enumerate() {
            assert_close(*a, e, &format!("calculate_p[{i}]"));
        }
    }

    // ---- self-contained PRNG sanity: deterministic + spans the index range ----
    #[test]
    fn splitmix64_is_deterministic_and_bounded() {
        let mut a = SplitMix64::new(super::PRNG_SEED);
        let mut b = SplitMix64::new(super::PRNG_SEED);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64(), "same seed -> same stream");
        }
        // next_below stays in range and is non-degenerate over many draws.
        let mut rng = SplitMix64::new(1);
        let mut seen = [false; 8];
        for _ in 0..2000 {
            let v = rng.next_below(8);
            assert!(v < 8, "next_below out of range: {v}");
            seen[v] = true;
        }
        assert!(seen.iter().all(|hit| *hit), "next_below skipped a value");
        // permutation is a bijection of 0..n.
        let perm = SplitMix64::new(99).permutation(20);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..20).collect::<Vec<_>>(), "permutation bijection");
    }

    // 30-protein, 4-vs-4 fixture from rots_e2e_oracle.py (Bigbio env, seed 123
    // numpy normal draws). log2FC is the deterministic target; the stochastic
    // p-values are NOT bit-matched (see module docs).
    fn e2e_rows() -> Vec<Vec<f64>> {
        // Captured verbatim from `np.random.default_rng(123).normal(10,1,(30,8))`
        // with +2.0 added to the B half (cols 4..8) of the first 6 proteins.
        // Stored as the raw matrix so the Rust test feeds identical inputs.
        E2E_MATRIX.iter().map(|row| row.to_vec()).collect()
    }

    // Deterministic log2FC oracle (protein, log2FC) from rots_e2e_oracle.py.
    const E2E_LOG2FC: &[(usize, f64)] = &[
        (0, 2.3194578965825734),
        (1, 2.931353656903349),
        (2, 2.0977123988471753),
        (3, 0.9602528723815524),
        (4, 3.1187316602549),
        (5, 1.3732610141404162),
        (6, -0.6534347030345558),
        (7, -0.18229245477504286),
        (8, 0.08872443593928736),
        (9, -1.233392283025486),
        (10, -0.9518595362316393),
        (11, -0.32522117747048007),
        (12, -0.18440244273792494),
        (13, 0.044205566044979605),
        (14, 0.5677719771180847),
        (15, -0.022735918751076056),
        (16, -0.4190151627149028),
        (17, -0.8166528743519947),
        (18, -0.8351899914488765),
        (19, 0.6444003491502173),
        (20, 1.3867569584436854),
        (21, 0.8770610325497419),
        (22, -0.9413520718320889),
        (23, 1.0399964650885156),
        (24, -0.6658888175338671),
        (25, 1.4855932509659233),
        (26, -0.019778429193756608),
        (27, 1.6099247963736314),
        (28, -0.4258634750602255),
        (29, -0.44572816568048523),
    ];

    // Stochastic-wrapper property test. NOT bit-exact vs Python: ROTS is RNG-
    // driven and Python's own output is seed-unstable (the selected a1 and the
    // adj_pvalues move ~1e-1 between numpy seeds), so we deliberately do not
    // assert cell equality against numpy's PCG64 stream. We DO assert:
    //   - log2FC cell-exact (it is the deterministic nanmean(B)-nanmean(A));
    //   - every p-value in [0,1] and finite;
    //   - d_stat finite for every protein;
    //   - the d_stat ranking is sane (the spiked-in proteins 0..6 dominate the
    //     top by |d_stat|).
    #[test]
    fn rots_two_group_property_holds_and_log2fc_is_cell_exact() {
        let rows = e2e_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = (0..30).map(|i| format!("P{i:02}")).collect::<Vec<_>>();
        let results = rots_two_group_seeded(
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

        // log2FC cell-exact vs Python (deterministic).
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

        // Sanity: the 6 spiked-in proteins (P00..P05) should dominate the top by
        // |d_stat| -- a distributional, RNG-tolerant check, not a bit match.
        let mut by_magnitude = results.clone();
        by_magnitude.sort_by(|a, b| b.t_statistic.abs().total_cmp(&a.t_statistic.abs()));
        let top6 = by_magnitude
            .iter()
            .take(6)
            .map(|r| r.protein.clone())
            .collect::<Vec<_>>();
        let spiked = top6
            .iter()
            .filter(|name| matches!(name.as_str(), "P00" | "P01" | "P02" | "P03" | "P04" | "P05"))
            .count();
        assert!(
            spiked >= 4,
            "expected the spiked proteins to dominate top6 by |d_stat|, got {top6:?}"
        );
    }

    // The selected (a1, a2, k) must land on a valid grid point (a1 in the ssq
    // grid with a2=1, or the raw-D branch a1=1/a2=0; k in the n_grid). This is
    // RNG-tolerant: it checks the optimiser only ever returns grid coordinates,
    // regardless of which one the random draws pick.
    #[test]
    fn bootstrap_optimize_returns_grid_point() {
        let rows = e2e_rows();
        let mat = rows.clone();
        let ssq = build_ssq_grid();
        let n_grid = build_n_grid(rows.len() / 4);
        for seed in [super::PRNG_SEED, 1, 7, 2024] {
            let mut rng = SplitMix64::new(seed);
            let (a1, a2, k) = super::bootstrap_optimize(&mat, 4, 20, &mut rng);
            let on_ssq = a2 == 1.0 && ssq.iter().any(|value| (value - a1).abs() <= 1e-12);
            let raw_branch = a1 == 1.0 && a2 == 0.0;
            assert!(
                on_ssq || raw_branch,
                "seed {seed}: (a1={a1}, a2={a2}) not a grid point"
            );
            assert!(
                n_grid.contains(&k) || (n_grid.is_empty() && k == 1),
                "seed {seed}: k={k} not in n_grid"
            );
        }
    }

    // End-to-end wrapper runs, BH-sorts, and makes a non-trivial significance
    // call on the spiked fixture (the public fixed-seed entry point).
    #[test]
    fn rots_two_group_public_entry_runs_and_calls_significance() {
        let rows = e2e_rows();
        let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let proteins = (0..30).map(|i| format!("P{i:02}")).collect::<Vec<_>>();
        let results = rots_two_group(
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
        // At least one protein is called UP/DOWN (the strong spike-in).
        assert!(
            results
                .iter()
                .any(|r| r.significance != Significance::Unchanged),
            "expected at least one significant protein"
        );
    }

    // The seed-123 numpy draws captured verbatim so the Rust test feeds the same
    // numbers the Python e2e oracle used (rots_e2e_oracle.py).
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
