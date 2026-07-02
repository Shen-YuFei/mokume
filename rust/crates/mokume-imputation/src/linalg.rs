//! Small dense linear-algebra helpers shared by the covariance-based imputers
//! (impSeq, impSeqRob, BPCA). Pure Rust, no external crates.

/// Solve the square system `a * x = b` by Gauss-Jordan elimination with partial
/// pivoting. `a` is `n x n` and `b` is length `n`. Returns `None` when `a` is
/// singular (mirroring `np.linalg.solve` raising `LinAlgError`). The
/// factorization matches LAPACK's up to floating-point accumulation order.
pub(crate) fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }

    // Augmented matrix [a | b].
    let mut augmented = a
        .iter()
        .zip(b)
        .map(|(row, &rhs)| {
            let mut extended = row.clone();
            extended.push(rhs);
            extended
        })
        .collect::<Vec<Vec<f64>>>();

    for column in 0..n {
        // Partial pivot: the row at or below `column` with the largest
        // magnitude in this column.
        let pivot = (column..n).max_by(|&left, &right| {
            augmented[left][column]
                .abs()
                .total_cmp(&augmented[right][column].abs())
        })?;
        if augmented[pivot][column].abs() < 1e-300 {
            return None;
        }
        augmented.swap(column, pivot);

        // Normalize the pivot row so the pivot becomes 1.
        let pivot_value = augmented[column][column];
        for value in &mut augmented[column] {
            *value /= pivot_value;
        }

        // Eliminate this column from every other row.
        let (before, rest) = augmented.split_at_mut(column);
        let (pivot_row, after) = rest.split_first_mut()?;
        for row in before.iter_mut().chain(after.iter_mut()) {
            let factor = row[column];
            if factor != 0.0 {
                for (target, &source) in row.iter_mut().zip(pivot_row.iter()) {
                    *target -= factor * source;
                }
            }
        }
    }

    Some(augmented.iter().map(|row| row[n]).collect())
}
