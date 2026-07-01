//! Robust sequential imputation (impSeqRob) ported from the pure-Python
//! `mokume.imputation.impseqrob` (a reimplementation of `rrcovNA::impSeqRob`).
//!
//! The matrix is laid out with proteins as rows and samples as columns, the same
//! orientation the Python pipeline feeds to `impute_impseqrob`. Values arrive
//! already in log2 space, so no extra log/exp is applied here. Only the
//! originally-missing cells are returned.
//!
//! Determinism note: the algorithm is DETERMINISTIC. The Python reference draws
//! random difference directions only when the number of complete-row pairs
//! exceeds `MAX_EXTRADIR_DIRECTIONS` (500), and even then it uses a fixed seed
//! (`numpy.random.default_rng(0)`); for any matrix with at most 32 complete rows
//! the full pairwise set is used and no RNG fires. This port reproduces the
//! deterministic full-pairwise branch exactly. The capped-and-seeded branch for
//! very tall matrices is not reproduced bit-for-bit (it would require matching
//! NumPy's PCG64 stream); instead this port keeps generating the full pairwise
//! set, which is the same set the Python code would have used below the cap and
//! a superset above it. Outside that tall-matrix regime the routine is
//! deterministic and matches Python cell-for-cell up to floating-point
//! accumulation order.
//!
//! Algorithm (matching the Python reference `_impseqrob_core`):
//! 1. Sort rows by ascending missing-count (stable). Rows with zero missing
//!    values are "complete".
//! 2. If there are fewer complete rows than `ceil(p / alpha)`, pre-impute the
//!    rows with the fewest missing values using a conditional multivariate
//!    normal estimated from pairwise-complete covariance (`_impnorm_preimpute`).
//!    The covariance is regularised to be positive-definite by shifting its
//!    diagonal when the smallest eigenvalue is below `EIG_FLOOR`.
//! 3. If pre-imputation removes all missing cells the result is returned.
//! 4. Otherwise estimate a robust (Stahel-Donoho / MCD) covariance from the
//!    complete rows, then visit the remaining incomplete rows in order. Each
//!    missing entry is filled from the conditional normal mean computed via the
//!    inverse covariance; the inverse is maintained with Sherman-Morrison rank
//!    one updates and the mean/covariance grow as non-outlying imputed rows are
//!    accepted into the support set.
//!
//! Reference: Beguin C, Hulliger B. The BACON-EEM Algorithm for Multivariate
//! Outlier Detection in Incomplete Survey Data. Survey Methodology. 2008.

use mokume_core::{ProteinId, SampleId};

use crate::linalg::solve_linear_system;

/// MCD support fraction, matching the Python `alpha` default.
const ALPHA: f64 = 0.9;
/// Cap on the number of difference directions; below this every row pair is
/// used (deterministic). See the module determinism note for the tall-matrix
/// case.
const MAX_EXTRADIR_DIRECTIONS: usize = 500;
/// Numerical floor reused throughout for "effectively zero" magnitudes, matching
/// the `1e-12` guards in the Python reference.
const TINY: f64 = 1e-12;
/// Consistency constant turning a MAD into a robust standard-deviation estimate
/// (NumPy `* 1.4826`).
const MAD_TO_SD: f64 = 1.4826;
/// Eigenvalue floor below which the pairwise-complete covariance is shifted to
/// become positive-definite (`eigvals.min() < 1e-10` in Python).
const EIG_FLOOR: f64 = 1e-10;
/// Diagonal shift added on top of `|min_eig|` when regularising (`+ 1e-6`).
const EIG_SHIFT: f64 = 1e-6;
/// Maximum sweeps for the symmetric Jacobi eigenvalue routine.
const JACOBI_MAX_SWEEPS: usize = 100;

pub(crate) fn impseqrob_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    value_at: &mut F,
) -> Vec<(ProteinId, SampleId, f64)>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let n_rows = proteins.len();
    let n_cols = samples.len();
    if n_rows == 0 || n_cols == 0 {
        return Vec::new();
    }

    // Read the matrix as `Option<f64>` cells (None = missing / non-finite).
    let matrix = proteins
        .iter()
        .map(|protein| {
            samples
                .iter()
                .map(|sample| value_at(*protein, *sample).filter(|value| value.is_finite()))
                .collect::<Vec<Option<f64>>>()
        })
        .collect::<Vec<_>>();

    let original_missing = matrix
        .iter()
        .map(|row| row.iter().map(Option::is_none).collect::<Vec<bool>>())
        .collect::<Vec<_>>();

    // Nothing to do when there are no missing cells.
    if original_missing.iter().all(|row| row.iter().all(|m| !*m)) {
        return Vec::new();
    }

    let filled = impseqrob_core(&matrix, n_cols);

    // Emit only the originally-missing cells whose fill is finite.
    let mut imputed = Vec::new();
    for (row_index, protein) in proteins.iter().enumerate() {
        for (col_index, sample) in samples.iter().enumerate() {
            if !original_missing[row_index][col_index] {
                continue;
            }
            let value = filled[row_index][col_index];
            if value.is_finite() {
                imputed.push((*protein, *sample, value));
            }
        }
    }
    imputed
}

/// Robust sequential imputation core, mirroring Python `_impseqrob_core`.
///
/// Returns a fully filled `n x p` matrix in the original row order. `matrix`
/// cells are `Some(value)` when observed and `None` when missing.
fn impseqrob_core(matrix: &[Vec<Option<f64>>], p: usize) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let miss_per_row = matrix
        .iter()
        .map(|row| row.iter().filter(|cell| cell.is_none()).count())
        .collect::<Vec<usize>>();

    // Stable argsort by ascending missing-count, then the inverse permutation
    // so the result can be restored to the original row order.
    let mut sort_order = (0..n).collect::<Vec<usize>>();
    sort_order.sort_by(|&left, &right| {
        miss_per_row[left]
            .cmp(&miss_per_row[right])
            .then(left.cmp(&right))
    });
    let mut unsort = vec![0usize; n];
    for (sorted_position, &original_index) in sort_order.iter().enumerate() {
        unsort[original_index] = sorted_position;
    }

    // `x` holds the working matrix in sorted-row order; `None` marks a hole.
    let mut x = sort_order
        .iter()
        .map(|&original| matrix[original].clone())
        .collect::<Vec<Vec<Option<f64>>>>();

    let mut complete_rows = sorted_complete_rows(&x);
    let n_complete = complete_rows.len();
    let mut incomplete_rows = sorted_incomplete_rows(&x);

    if incomplete_rows.is_empty() {
        return restore_order(&materialize(&x), &unsort);
    }

    // ceil(p / alpha): smallest number of complete rows the robust core needs.
    let min_complete = (p as f64 / ALPHA).ceil() as usize;
    if n_complete < min_complete {
        impnorm_preimpute(&mut x, p, min_complete - n_complete);
        complete_rows = sorted_complete_rows(&x);
        incomplete_rows = sorted_incomplete_rows(&x);
        if incomplete_rows.is_empty() {
            return restore_order(&materialize(&x), &unsort);
        }
    }

    robust_sequential_core(&mut x, p, &complete_rows, &incomplete_rows);
    restore_order(&materialize(&x), &unsort)
}

/// Indices (in sorted order) of rows that have no missing cell.
fn sorted_complete_rows(x: &[Vec<Option<f64>>]) -> Vec<usize> {
    x.iter()
        .enumerate()
        .filter(|(_, row)| row.iter().all(Option::is_some))
        .map(|(index, _)| index)
        .collect()
}

/// Indices (in sorted order) of rows that still have a missing cell.
fn sorted_incomplete_rows(x: &[Vec<Option<f64>>]) -> Vec<usize> {
    x.iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(Option::is_none))
        .map(|(index, _)| index)
        .collect()
}

/// Replace every cell with a concrete `f64` (`NaN` for any residual hole).
fn materialize(x: &[Vec<Option<f64>>]) -> Vec<Vec<f64>> {
    x.iter()
        .map(|row| {
            row.iter()
                .map(|cell| cell.unwrap_or(f64::NAN))
                .collect::<Vec<f64>>()
        })
        .collect()
}

/// Undo the missing-count sort, returning rows to their original positions.
fn restore_order(sorted: &[Vec<f64>], unsort: &[usize]) -> Vec<Vec<f64>> {
    unsort
        .iter()
        .map(|&sorted_position| sorted[sorted_position].clone())
        .collect()
}

/// Pre-impute the `n_needed` rows with the fewest missing values using a
/// conditional multivariate normal, mirroring Python `_impnorm_preimpute`.
///
/// `mu` is the column-wise mean over observed entries (`np.nanmean`). `sigma` is
/// the pairwise-complete-case covariance (off-diagonal from pairs observed in
/// both columns with `ddof = 1`; diagonal from the single-column variance, or
/// `1.0` when a column has fewer than two observed values). The covariance is
/// shifted to positive-definiteness when its smallest eigenvalue is below the
/// floor. Each selected row's missing block is filled with
/// `mu_m + Sigma_mo Sigma_oo^-1 (x_o - mu_o)`.
fn impnorm_preimpute(x: &mut [Vec<Option<f64>>], p: usize, n_needed: usize) {
    let mu = column_observed_means(x, p);
    let mut sigma = pairwise_complete_covariance(x, p);
    regularize_to_positive_definite(&mut sigma, p);

    // Candidate rows (those with at least one hole) ordered by ascending
    // missing-count; the first `n_needed` are imputed. Ties keep ascending row
    // index (NumPy's argsort is not guaranteed stable, but the candidate counts
    // here are processed first-to-last and the relative order of equal counts
    // does not change the per-row conditional fill).
    let mut candidates = x
        .iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(Option::is_none))
        .map(|(index, row)| (index, row.iter().filter(|cell| cell.is_none()).count()))
        .collect::<Vec<(usize, usize)>>();
    candidates.sort_by(|left, right| left.1.cmp(&right.1).then(left.0.cmp(&right.0)));

    for &(row_index, _) in candidates.iter().take(n_needed) {
        let observed = observed_columns(&x[row_index]);
        let missing = missing_columns(&x[row_index]);
        if observed.is_empty() {
            for &col in &missing {
                x[row_index][col] = Some(mu[col]);
            }
            continue;
        }

        // residual = x_o - mu_o.
        let residual = observed
            .iter()
            .map(|&col| x[row_index][col].unwrap_or(mu[col]) - mu[col])
            .collect::<Vec<f64>>();

        // Solve Sigma_oo z = residual instead of forming the explicit inverse.
        let sigma_oo = submatrix(&sigma, &observed, &observed);
        let Some(z) = solve_linear_system(&sigma_oo, &residual) else {
            for &col in &missing {
                x[row_index][col] = Some(mu[col]);
            }
            continue;
        };

        // x_m = mu_m + Sigma_mo z.
        for &m in &missing {
            let mut acc = mu[m];
            for (k, &o) in observed.iter().enumerate() {
                acc += sigma[m][o] * z[k];
            }
            x[row_index][m] = Some(acc);
        }
    }
}

/// `np.nanmean(x, axis=0)`: per-column mean over observed entries (`0.0` for a
/// fully missing column, which the pipeline does not produce).
fn column_observed_means(x: &[Vec<Option<f64>>], p: usize) -> Vec<f64> {
    (0..p)
        .map(|col| {
            let mut sum = 0.0;
            let mut count = 0usize;
            for row in x {
                if let Some(value) = row[col] {
                    sum += value;
                    count += 1;
                }
            }
            if count > 0 {
                sum / count as f64
            } else {
                0.0
            }
        })
        .collect()
}

/// Pairwise-complete covariance matching the Python `_impnorm_preimpute` build:
/// off-diagonal `(j, k)` from the rows observed in both columns (`ddof = 1`),
/// the diagonal from the single-column variance (or `1.0` when a column has
/// fewer than two observed values).
fn pairwise_complete_covariance(x: &[Vec<Option<f64>>], p: usize) -> Vec<Vec<f64>> {
    let mut sigma = vec![vec![0.0; p]; p];
    for j in 0..p {
        for k in j..p {
            let pairs = x
                .iter()
                .filter_map(|row| match (row[j], row[k]) {
                    (Some(value_j), Some(value_k)) => Some((value_j, value_k)),
                    _ => None,
                })
                .collect::<Vec<(f64, f64)>>();
            if pairs.len() > 1 {
                let value = pair_covariance(&pairs);
                sigma[j][k] = value;
                sigma[k][j] = value;
            } else if j == k {
                let column = x.iter().filter_map(|row| row[j]).collect::<Vec<f64>>();
                sigma[j][j] = if column.len() > 1 {
                    sample_variance(&column)
                } else {
                    1.0
                };
            }
        }
    }
    sigma
}

/// Sample covariance of paired values with `ddof = 1` (NumPy `np.cov(..)[0, 1]`).
fn pair_covariance(pairs: &[(f64, f64)]) -> f64 {
    let count = pairs.len() as f64;
    let mean_x = pairs.iter().map(|(x, _)| x).sum::<f64>() / count;
    let mean_y = pairs.iter().map(|(_, y)| y).sum::<f64>() / count;
    pairs
        .iter()
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>()
        / (count - 1.0)
}

/// Sample variance with `ddof = 1` (NumPy `np.var(.., ddof=1)`).
fn sample_variance(values: &[f64]) -> f64 {
    let count = values.len() as f64;
    let mean = values.iter().sum::<f64>() / count;
    values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (count - 1.0)
}

/// Shift the diagonal of a symmetric covariance until it is positive-definite,
/// mirroring `sigma += (|min_eig| + 1e-6) I` when `min_eig < 1e-10`.
fn regularize_to_positive_definite(sigma: &mut [Vec<f64>], p: usize) {
    let min_eig = smallest_symmetric_eigenvalue(sigma, p);
    if min_eig < EIG_FLOOR {
        let shift = min_eig.abs() + EIG_SHIFT;
        for (index, row) in sigma.iter_mut().enumerate().take(p) {
            row[index] += shift;
        }
    }
}

/// Smallest eigenvalue of a symmetric matrix via the cyclic Jacobi method
/// (matches NumPy `np.linalg.eigvalsh(..).min()` up to float accumulation). For
/// the small covariance matrices used here (one row per sample) this converges
/// in a handful of sweeps.
fn smallest_symmetric_eigenvalue(matrix: &[Vec<f64>], p: usize) -> f64 {
    if p == 0 {
        return 0.0;
    }
    if p == 1 {
        return matrix[0][0];
    }

    let mut a = matrix.iter().map(Clone::clone).collect::<Vec<Vec<f64>>>();
    for _ in 0..JACOBI_MAX_SWEEPS {
        let mut off_diagonal = 0.0;
        for (row_index, row) in a.iter().enumerate().take(p) {
            for &value in row.iter().take(p).skip(row_index + 1) {
                off_diagonal += value * value;
            }
        }
        if off_diagonal <= TINY * TINY {
            break;
        }
        for q in 1..p {
            for r in 0..q {
                jacobi_rotate(&mut a, p, r, q);
            }
        }
    }

    (0..p)
        .map(|index| a[index][index])
        .fold(f64::INFINITY, f64::min)
}

/// Apply one Jacobi rotation zeroing the symmetric `(r, q)` off-diagonal entry.
fn jacobi_rotate(a: &mut [Vec<f64>], p: usize, r: usize, q: usize) {
    let a_rq = a[r][q];
    if a_rq.abs() <= TINY {
        return;
    }
    let a_rr = a[r][r];
    let a_qq = a[q][q];
    let theta = (a_qq - a_rr) / (2.0 * a_rq);
    let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
    let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
    let c = 1.0 / (t * t + 1.0).sqrt();
    let s = t * c;

    // Compute the rotated `(k, r)` / `(k, q)` entries into a temporary buffer
    // first. Each iteration depends only on the pre-rotation row `k`, which is
    // never touched by another iteration (`k != r` and `k != q`), so the writes
    // can be deferred without changing any value.
    let rotated: Vec<(usize, f64, f64)> = a
        .iter()
        .enumerate()
        .take(p)
        .filter(|&(k, _)| k != r && k != q)
        .map(|(k, row_k)| {
            let a_kr = row_k[r];
            let a_kq = row_k[q];
            (k, c * a_kr - s * a_kq, s * a_kr + c * a_kq)
        })
        .collect();
    for (k, new_kr, new_kq) in rotated {
        a[k][r] = new_kr;
        a[r][k] = new_kr;
        a[k][q] = new_kq;
        a[q][k] = new_kq;
    }
    let new_rr = c * c * a_rr - 2.0 * s * c * a_rq + s * s * a_qq;
    let new_qq = s * s * a_rr + 2.0 * s * c * a_rq + c * c * a_qq;
    a[r][r] = new_rr;
    a[q][q] = new_qq;
    a[r][q] = 0.0;
    a[q][r] = 0.0;
}

/// Robust sequential imputation over the incomplete rows, mirroring the
/// main loop of Python `_impseqrob_core` once a robust covariance is available.
fn robust_sequential_core(
    x: &mut [Vec<Option<f64>>],
    p: usize,
    complete_rows: &[usize],
    incomplete_rows: &[usize],
) {
    let n_complete = complete_rows.len();
    let complete_data = complete_rows
        .iter()
        .map(|&index| {
            x[index]
                .iter()
                .map(|cell| cell.unwrap_or(f64::NAN))
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>();

    let h0 = h_alpha_n(n_complete, p);
    let robust = match cov_stahel_donoho(&complete_data, p, h0) {
        Some(robust) => robust,
        // Degenerate complete set (no usable directions): fall back to filling
        // every remaining hole with its column mean.
        None => {
            fill_remaining_with_column_means(x, p, incomplete_rows);
            return;
        }
    };

    let mut covx = robust.cov;
    let mut mx = robust.mean;
    let directions = robust.directions;
    let med_y = robust.med_y;
    let mad_y = robust.mad_y;

    // Sorted outlyingness list (`ds_vals`); grows as imputed rows are added.
    let mut ds_vals = robust.outlyingness.clone();
    ds_vals.sort_by(f64::total_cmp);

    let mut h = h0;
    let mut n_update = h0;

    // `flag == false` marks an accepted (non-outlying) row used to grow the
    // estimates; the MCD support of the complete rows starts accepted.
    let mut accepted_prev = false;
    let mut icovx: Option<Vec<Vec<f64>>> = None;
    let mut mxo = mx.clone();

    for (loop_index, &row_index) in incomplete_rows.iter().enumerate() {
        let missing = missing_columns(&x[row_index]);
        let observed = observed_columns(&x[row_index]);

        if loop_index == 0 {
            icovx = invert_matrix(&covx, p);
        } else if accepted_prev {
            if let Some(current) = icovx.as_ref() {
                let denom = (n_update as f64 - 1.0).max(1.0).sqrt();
                let prev_row = &x_row_values(x, incomplete_rows[loop_index - 1]);
                let f = prev_row
                    .iter()
                    .zip(mxo.iter())
                    .map(|(value, mean)| (value - mean) / denom)
                    .collect::<Vec<f64>>();
                let scale = (n_update as f64 - 1.0) / (n_update as f64 - 2.0).max(1.0);
                let scaled = scale_matrix(current, scale);
                icovx = Some(update_invcov(&scaled, &f));
            }
        }

        let Some(inverse) = icovx.as_ref() else {
            for &col in &missing {
                x[row_index][col] = Some(mx[col]);
            }
            accepted_prev = false;
            continue;
        };

        // Conditional mean via the inverse covariance blocks:
        // x_m = mu_m - icov_mm^-1 icov_mo (x_o - mu_o).
        if observed.is_empty() {
            for &col in &missing {
                x[row_index][col] = Some(mx[col]);
            }
        } else {
            let icov_mm = submatrix(inverse, &missing, &missing);
            let residual = observed
                .iter()
                .map(|&col| x[row_index][col].unwrap_or(mx[col]) - mx[col])
                .collect::<Vec<f64>>();
            // rhs = icov_mo (x_o - mu_o).
            let rhs = missing
                .iter()
                .map(|&m| {
                    observed
                        .iter()
                        .enumerate()
                        .map(|(k, &o)| inverse[m][o] * residual[k])
                        .sum::<f64>()
                })
                .collect::<Vec<f64>>();
            match solve_linear_system(&icov_mm, &rhs) {
                Some(adjustment) => {
                    for (k, &m) in missing.iter().enumerate() {
                        x[row_index][m] = Some(mx[m] - adjustment[k]);
                    }
                }
                None => {
                    for &m in &missing {
                        x[row_index][m] = Some(mx[m]);
                    }
                }
            }
        }

        // Outlyingness of the now-complete row against the robust directions.
        let row_values = x_row_values(x, row_index);
        let dx = outlyingness(&row_values, &directions, &med_y, &mad_y);

        // Accept the row when its outlyingness is within the current threshold.
        let threshold_index = h.min(ds_vals.len().saturating_sub(1));
        let threshold = ds_vals
            .get(threshold_index)
            .copied()
            .unwrap_or(f64::INFINITY);
        if dx <= threshold {
            accepted_prev = true;
            n_update += 1;
            h += 1;
            mxo = mx.clone();
            // mx = ((n_update - 1) mxo + row) / n_update.
            mx = mxo
                .iter()
                .zip(row_values.iter())
                .map(|(old, value)| ((n_update as f64 - 1.0) * old + value) / n_update as f64)
                .collect();
            // covx = (n_update - 2)/(n_update - 1) covx + outer(row - mxo) / (n_update - 1).
            let cov_scale = (n_update as f64 - 2.0) / (n_update as f64 - 1.0);
            let outer_scale = 1.0 / (n_update as f64 - 1.0);
            let delta = row_values
                .iter()
                .zip(mxo.iter())
                .map(|(value, old)| value - old)
                .collect::<Vec<f64>>();
            for a in 0..p {
                for b in 0..p {
                    covx[a][b] = cov_scale * covx[a][b] + outer_scale * delta[a] * delta[b];
                }
            }
        } else {
            accepted_prev = false;
            mxo = mx.clone();
        }

        // Append dx to the sorted outlyingness list.
        let insert_at = ds_vals.partition_point(|value| *value < dx);
        ds_vals.insert(insert_at, dx);
    }
}

/// Snapshot a row's current values (holes already filled by this point).
fn x_row_values(x: &[Vec<Option<f64>>], row_index: usize) -> Vec<f64> {
    x[row_index]
        .iter()
        .map(|cell| cell.unwrap_or(f64::NAN))
        .collect()
}

/// Fall back to filling each remaining hole with its column mean.
fn fill_remaining_with_column_means(
    x: &mut [Vec<Option<f64>>],
    p: usize,
    incomplete_rows: &[usize],
) {
    let means = column_observed_means(x, p);
    for &row_index in incomplete_rows {
        for col in 0..p {
            if x[row_index][col].is_none() {
                x[row_index][col] = Some(means[col]);
            }
        }
    }
}

/// `h` for the MCD support, matching `robustbase:::h.alpha.n`:
/// `floor(2 * h2 - n + 2 (n - h2) alpha)` with `h2 = (n + p + 1) / 2`.
fn h_alpha_n(n: usize, p: usize) -> usize {
    let h2 = (n + p).div_ceil(2);
    let value = 2.0 * h2 as f64 - n as f64 + 2.0 * (n as f64 - h2 as f64) * ALPHA;
    value.floor().max(1.0) as usize
}

/// Robust covariance summary returned by the Stahel-Donoho estimator.
struct RobustCov {
    cov: Vec<Vec<f64>>,
    mean: Vec<f64>,
    directions: Vec<Vec<f64>>,
    med_y: Vec<f64>,
    mad_y: Vec<f64>,
    outlyingness: Vec<f64>,
}

/// Stahel-Donoho robust covariance over the complete rows, mirroring Python
/// `_covSD`. Returns `None` when no usable direction survives normalisation.
fn cov_stahel_donoho(data: &[Vec<f64>], p: usize, h: usize) -> Option<RobustCov> {
    let raw_directions = extradir(data, p);
    let mut directions = Vec::new();
    for direction in &raw_directions {
        let norm = direction
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        if norm > TINY {
            directions.push(
                direction
                    .iter()
                    .map(|value| value / norm)
                    .collect::<Vec<f64>>(),
            );
        }
    }
    if directions.is_empty() {
        return None;
    }

    // Projections Y = data @ A^T, robust centre/scale per direction.
    let projections = data
        .iter()
        .map(|row| {
            directions
                .iter()
                .map(|direction| {
                    row.iter()
                        .zip(direction.iter())
                        .map(|(value, weight)| value * weight)
                        .sum::<f64>()
                })
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<Vec<f64>>>();

    let direction_count = directions.len();
    let mut med_y = vec![0.0; direction_count];
    let mut mad_y = vec![0.0; direction_count];
    for column in 0..direction_count {
        let mut values = projections
            .iter()
            .map(|row| row[column])
            .collect::<Vec<f64>>();
        let median = median_of(&mut values);
        med_y[column] = median;
        let mut deviations = values
            .iter()
            .map(|value| (value - median).abs())
            .collect::<Vec<f64>>();
        let mad = median_of(&mut deviations) * MAD_TO_SD;
        mad_y[column] = if mad < TINY { TINY } else { mad };
    }

    // Per-row outlyingness d = max over directions of |Y - medY| / madY.
    let outlyingness = projections
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(column, value)| (value - med_y[column]).abs() / mad_y[column])
                .fold(0.0_f64, f64::max)
        })
        .collect::<Vec<f64>>();

    // MCD support: the `h` rows with smallest outlyingness (stable argsort).
    let mut order = (0..data.len()).collect::<Vec<usize>>();
    order.sort_by(|&left, &right| {
        outlyingness[left]
            .total_cmp(&outlyingness[right])
            .then(left.cmp(&right))
    });
    let support = &order[..h.min(order.len())];
    let support_rows = support
        .iter()
        .map(|&index| data[index].clone())
        .collect::<Vec<_>>();

    let mean = column_means(&support_rows, p);
    let cov = sample_covariance(&support_rows, &mean, p);

    Some(RobustCov {
        cov,
        mean,
        directions,
        med_y,
        mad_y,
        outlyingness,
    })
}

/// Difference directions for Stahel-Donoho outlyingness, mirroring Python
/// `_extradir`. For at most `MAX_EXTRADIR_DIRECTIONS` row pairs every ordered
/// pair `(i, j), i < j` contributes `data[i] - data[j]` (deterministic). Above
/// the cap the Python reference samples with a fixed seed; this port keeps the
/// full pairwise set (see the module determinism note).
fn extradir(data: &[Vec<f64>], _p: usize) -> Vec<Vec<f64>> {
    let n = data.len();
    let mut directions = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let diff = data[i]
                .iter()
                .zip(data[j].iter())
                .map(|(a, b)| a - b)
                .collect::<Vec<f64>>();
            directions.push(diff);
        }
    }
    let _ = MAX_EXTRADIR_DIRECTIONS;
    directions
}

/// Sherman-Morrison rank-one inverse update, mirroring `_update_invcov`:
/// `A - (A f)(A f)^T / (1 + f^T A f)`. Returns `A` unchanged when the
/// denominator is effectively zero.
fn update_invcov(a: &[Vec<f64>], f: &[f64]) -> Vec<Vec<f64>> {
    let p = f.len();
    let af = (0..p)
        .map(|row| {
            a[row]
                .iter()
                .zip(f.iter())
                .map(|(value, weight)| value * weight)
                .sum::<f64>()
        })
        .collect::<Vec<f64>>();
    let denom = 1.0
        + f.iter()
            .zip(af.iter())
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
    if denom.abs() <= TINY {
        return a.iter().map(Clone::clone).collect();
    }
    let mut result = a.iter().map(Clone::clone).collect::<Vec<Vec<f64>>>();
    for row in 0..p {
        for column in 0..p {
            result[row][column] -= af[row] * af[column] / denom;
        }
    }
    result
}

/// Multiply every entry of a square matrix by a scalar.
fn scale_matrix(matrix: &[Vec<f64>], scale: f64) -> Vec<Vec<f64>> {
    matrix
        .iter()
        .map(|row| row.iter().map(|value| value * scale).collect())
        .collect()
}

/// Outlyingness of a single row against the robust directions:
/// `max_k |row . A_k - medY_k| / madY_k`.
fn outlyingness(row: &[f64], directions: &[Vec<f64>], med_y: &[f64], mad_y: &[f64]) -> f64 {
    directions
        .iter()
        .enumerate()
        .map(|(column, direction)| {
            let projection = row
                .iter()
                .zip(direction.iter())
                .map(|(value, weight)| value * weight)
                .sum::<f64>();
            (projection - med_y[column]).abs() / mad_y[column]
        })
        .fold(0.0_f64, f64::max)
}

/// Invert a square matrix by solving against the identity columns; `None` when
/// singular. Reuses the shared Gauss-Jordan solver one column at a time.
fn invert_matrix(matrix: &[Vec<f64>], p: usize) -> Option<Vec<Vec<f64>>> {
    let mut columns = Vec::with_capacity(p);
    for column in 0..p {
        let mut unit = vec![0.0; p];
        unit[column] = 1.0;
        columns.push(solve_linear_system(matrix, &unit)?);
    }
    // `columns[c]` is the c-th column of the inverse; transpose into rows.
    let mut inverse = vec![vec![0.0; p]; p];
    for (column, solution) in columns.iter().enumerate() {
        for (row, &value) in solution.iter().enumerate() {
            inverse[row][column] = value;
        }
    }
    Some(inverse)
}

/// Column indices of the observed cells of a row.
fn observed_columns(row: &[Option<f64>]) -> Vec<usize> {
    row.iter()
        .enumerate()
        .filter(|(_, cell)| cell.is_some())
        .map(|(index, _)| index)
        .collect()
}

/// Column indices of the missing cells of a row.
fn missing_columns(row: &[Option<f64>]) -> Vec<usize> {
    row.iter()
        .enumerate()
        .filter(|(_, cell)| cell.is_none())
        .map(|(index, _)| index)
        .collect()
}

/// Extract the submatrix selected by `rows` and `cols`.
fn submatrix(matrix: &[Vec<f64>], rows: &[usize], cols: &[usize]) -> Vec<Vec<f64>> {
    rows.iter()
        .map(|&row| cols.iter().map(|&col| matrix[row][col]).collect())
        .collect()
}

/// Plain per-column mean over the supplied rows.
fn column_means(rows: &[Vec<f64>], p: usize) -> Vec<f64> {
    let count = rows.len();
    (0..p)
        .map(|col| {
            if count == 0 {
                0.0
            } else {
                rows.iter().map(|row| row[col]).sum::<f64>() / count as f64
            }
        })
        .collect()
}

/// Sample covariance (`ddof = 1`) over the supplied rows, matching NumPy
/// `np.cov(rows, rowvar=False, ddof=1)`. A single row yields all zeros.
fn sample_covariance(rows: &[Vec<f64>], mean: &[f64], p: usize) -> Vec<Vec<f64>> {
    let count = rows.len();
    let mut cov = vec![vec![0.0; p]; p];
    if count < 2 {
        return cov;
    }
    let denom = (count - 1) as f64;
    for a in 0..p {
        for b in 0..p {
            let sum = rows
                .iter()
                .map(|row| (row[a] - mean[a]) * (row[b] - mean[b]))
                .sum::<f64>();
            cov[a][b] = sum / denom;
        }
    }
    cov
}

/// Median of a slice (NumPy `np.median`: lower-of-middle-two averaged). Sorts in
/// place with `total_cmp`; returns `0.0` for an empty slice.
fn median_of(values: &mut [f64]) -> f64 {
    let count = values.len();
    if count == 0 {
        return 0.0;
    }
    values.sort_by(f64::total_cmp);
    if count % 2 == 1 {
        values[count / 2]
    } else {
        (values[count / 2 - 1] + values[count / 2]) / 2.0
    }
}
