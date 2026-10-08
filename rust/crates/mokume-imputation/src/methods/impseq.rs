//! Sequential conditional-normal imputation matching rrcovNA::impSeq.
//! Protein rows are processed by increasing missing count. When the current
//! complete set has at most p rows, the official covariance gets 0.01 I.
//! Inputs and returned fills use the supplied scale (normally log2).

use mokume_core::{MokumeError, ProteinId, Result, SampleId};

use crate::linalg::solve_linear_system;

pub(crate) fn impseq_imputed_values<F>(
    proteins: &[ProteinId],
    samples: &[SampleId],
    value_at: &mut F,
) -> Result<Vec<(ProteinId, SampleId, f64)>>
where
    F: FnMut(ProteinId, SampleId) -> Option<f64>,
{
    let n_rows = proteins.len();
    let n_cols = samples.len();
    if n_rows == 0 || n_cols == 0 {
        return Ok(Vec::new());
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
        return Ok(Vec::new());
    }

    // Working copy holding the imputed values; complete cells keep their value.
    let mut result = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|cell| cell.unwrap_or(f64::NAN))
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<_>>();

    let complete_mask = matrix
        .iter()
        .map(|row| row.iter().all(Option::is_some))
        .collect::<Vec<bool>>();
    let n_complete = complete_mask.iter().filter(|c| **c).count();

    if n_complete < 2 {
        return Err(MokumeError::InvalidInput {
            message: "impSeq requires at least two complete protein rows to estimate covariance"
                .to_owned(),
        });
    }
    impseq_conditional_normal(&matrix, &complete_mask, n_cols, &mut result)?;

    // Emit only the originally-missing cells whose fill is finite.
    let mut imputed = Vec::new();
    for (row_index, protein) in proteins.iter().enumerate() {
        for (col_index, sample) in samples.iter().enumerate() {
            if !original_missing[row_index][col_index] {
                continue;
            }
            let value = result[row_index][col_index];
            if value.is_finite() {
                imputed.push((*protein, *sample, value));
            }
        }
    }
    Ok(imputed)
}

/// Conditional-normal sequential imputation over the incomplete rows.
fn impseq_conditional_normal(
    matrix: &[Vec<Option<f64>>],
    complete_mask: &[bool],
    n_cols: usize,
    result: &mut [Vec<f64>],
) -> Result<()> {
    let complete_rows = complete_mask
        .iter()
        .zip(result.iter())
        .filter(|(complete, _)| **complete)
        .map(|(_, row)| row.clone())
        .collect::<Vec<_>>();

    let (mut mu, mut cov) = mean_and_cov(&complete_rows, n_cols);

    // Incomplete rows ordered by ascending missing-count, ties by row index
    // (stable), matching the Python `incomplete_idx[np.argsort(n_miss)]`.
    let mut order = complete_mask
        .iter()
        .enumerate()
        .filter_map(|(index, complete)| (!*complete).then_some(index))
        .collect::<Vec<usize>>();
    order.sort_by(|left, right| {
        let left_miss = matrix[*left].iter().filter(|cell| cell.is_none()).count();
        let right_miss = matrix[*right].iter().filter(|cell| cell.is_none()).count();
        left_miss.cmp(&right_miss).then_with(|| left.cmp(right))
    });

    let mut used_rows = complete_rows;

    for row_index in order {
        let observed = matrix[row_index]
            .iter()
            .enumerate()
            .filter_map(|(col, cell)| cell.is_some().then_some(col))
            .collect::<Vec<usize>>();
        let missing = matrix[row_index]
            .iter()
            .enumerate()
            .filter_map(|(col, cell)| cell.is_none().then_some(col))
            .collect::<Vec<usize>>();

        if observed.is_empty() {
            for &col in &missing {
                if let Some(value) = mu.get(col) {
                    result[row_index][col] = *value;
                }
            }
        } else {
            let mut estimation_cov = cov.clone();
            if used_rows.len() <= n_cols {
                for (i, row) in estimation_cov.iter_mut().enumerate() {
                    row[i] += 0.01;
                }
            }
            let imputed = conditional_normal_impute(
                &result[row_index],
                &observed,
                &missing,
                &mu,
                &estimation_cov,
            )?;
            for (slot, &col) in imputed.into_iter().zip(missing.iter()) {
                result[row_index][col] = slot;
            }
        }

        used_rows.push(result[row_index].clone());
        let (next_mu, next_cov) = mean_and_cov(&used_rows, n_cols);
        mu = next_mu;
        if used_rows.len() > 1 {
            cov = next_cov;
        }
    }
    Ok(())
}

/// Impute the missing entries of `row` from the conditional normal mean:
/// `x_m = mu_m + Sigma_mo * Sigma_oo^-1 * (x_o - mu_o)`.
/// An unsolvable covariance produces an error rather than a different imputer.
fn conditional_normal_impute(
    row: &[f64],
    observed: &[usize],
    missing: &[usize],
    mu: &[f64],
    cov: &[Vec<f64>],
) -> Result<Vec<f64>> {
    // Sigma_oo (k x k) and the right-hand side (x_o - mu_o).
    let k = observed.len();
    let mut sigma_oo = vec![vec![0.0; k]; k];
    for (a, &oa) in observed.iter().enumerate() {
        for (b, &ob) in observed.iter().enumerate() {
            sigma_oo[a][b] = cov[oa][ob];
        }
    }
    let mut diff = vec![0.0; k];
    for (a, &oa) in observed.iter().enumerate() {
        diff[a] = row[oa] - mu[oa];
    }

    // Solve Sigma_oo * z = (x_o - mu_o) instead of forming the explicit
    // inverse; mathematically identical to `Sigma_oo^-1 (x_o - mu_o)`.
    let z = solve_linear_system(&sigma_oo, &diff).ok_or_else(|| MokumeError::InvalidInput {
        message: "impSeq observed covariance is singular".to_owned(),
    })?;

    // x_m = mu_m + Sigma_mo * z.
    Ok(missing
        .iter()
        .map(|&m| {
            let mut acc = mu[m];
            for (a, &oa) in observed.iter().enumerate() {
                acc += cov[m][oa] * z[a];
            }
            acc
        })
        .collect())
}

/// Column means and sample covariance (ddof = 1) of `rows` (each length
/// `n_cols`). The covariance matches NumPy `np.cov(rows, rowvar=False,
/// ddof=1)`; with a single row every covariance entry is zero.
fn mean_and_cov(rows: &[Vec<f64>], n_cols: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut mu = vec![0.0; n_cols];
    let mut cov = vec![vec![0.0; n_cols]; n_cols];
    let m = rows.len();
    if m == 0 {
        return (mu, cov);
    }
    for col in 0..n_cols {
        let sum: f64 = rows.iter().map(|row| row[col]).sum();
        mu[col] = sum / m as f64;
    }
    if m < 2 {
        return (mu, cov);
    }
    let denom = (m - 1) as f64;
    for a in 0..n_cols {
        for b in 0..n_cols {
            let sum: f64 = rows
                .iter()
                .map(|row| (row[a] - mu[a]) * (row[b] - mu[b]))
                .sum();
            cov[a][b] = sum / denom;
        }
    }
    (mu, cov)
}
