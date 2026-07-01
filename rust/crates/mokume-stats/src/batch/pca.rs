//! Deterministic PCA for the batch-correction outlier pass.
//!
//! Port of `sklearn.decomposition.PCA` as `mokume.postprocessing.batch_correction`
//! `compute_pca` invokes it (`PCA(n_components).fit_transform`). The Python path
//! runs PCA on the `samples x features` matrix (`df.T`), so a point is a sample
//! and the returned coordinates are the per-sample principal-component scores fed
//! to HDBSCAN.
//!
//! sklearn centers each feature, runs an exact SVD `X_centered = U S V^T`, applies
//! the `svd_flip` sign convention, and returns the scores `U S = X_centered V`.
//! `PCA._fit_full` calls `svd_flip(U, Vt, u_based_decision=False)`: the sign of
//! each component `k` is chosen so the largest-magnitude entry of the *right*
//! singular vector `V[:, k]` (a row of `components_`) is positive. This module
//! reproduces that result through the symmetric eigendecomposition of the
//! (unscaled) covariance `X_centered^T X_centered`, whose eigenvectors are the
//! right singular vectors `V` and whose eigenvalues are `S^2`. The scores are
//! `X_centered V`, with the same `u_based_decision=False` sign rule applied from
//! the eigenvector columns `V[:, k]`.
//!
//! The eigensolver is the cyclic-Jacobi routine ([`symmetric_eigen_descending`]),
//! which is deterministic and order-independent, so no RNG or randomized SVD is
//! involved. Outlier detection consumes only Euclidean distances between score
//! rows, which are invariant to per-component sign, so the sign convention only
//! matters for reproducing the raw coordinates in tests.

/// Number of cyclic Jacobi sweeps; ample for the small covariance matrices
/// (`features x features`, capped by component count) the outlier pass builds.
const JACOBI_SWEEPS: usize = 100;

/// Below-this off-diagonal magnitude the Jacobi sweep is treated as converged.
const JACOBI_EPS: f64 = 1e-300;

/// Principal-component scores for `data` (rows are points/samples, columns are
/// features), keeping the leading `n_components` components, mirroring
/// `compute_pca`. The number of components is clamped to
/// `min(n_components, n_samples, n_features)` exactly as sklearn does. Returns one
/// score row per input row, each of length equal to the retained component count.
pub fn pca_scores(data: &[Vec<f64>], n_components: usize) -> Vec<Vec<f64>> {
    let n_samples = data.len();
    if n_samples == 0 {
        return Vec::new();
    }
    let n_features = data[0].len();
    if n_features == 0 {
        return vec![Vec::new(); n_samples];
    }
    let keep = n_components.min(n_samples).min(n_features);
    if keep == 0 {
        return vec![Vec::new(); n_samples];
    }

    // Center each feature column (sklearn subtracts the per-feature mean).
    let mut means = vec![0.0_f64; n_features];
    for row in data {
        for (feature, &value) in row.iter().enumerate() {
            means[feature] += value;
        }
    }
    for mean in &mut means {
        *mean /= n_samples as f64;
    }
    let centered: Vec<Vec<f64>> = data
        .iter()
        .map(|row| {
            row.iter()
                .zip(means.iter())
                .map(|(&value, &mean)| value - mean)
                .collect()
        })
        .collect();

    // Covariance-style Gram matrix C = X_centered^T X_centered (features x features).
    // Its eigenvectors are the right singular vectors V; eigenvalues are S^2.
    let mut gram = vec![vec![0.0_f64; n_features]; n_features];
    for row in &centered {
        for (i, &xi) in row.iter().enumerate() {
            for (j, &xj) in row.iter().enumerate() {
                gram[i][j] += xi * xj;
            }
        }
    }
    // Only the eigenvectors (right singular vectors V) are needed: the scores are
    // X_centered V, and the component count is fixed by `keep`, independent of the
    // eigenvalues / singular values themselves.
    let (_eigenvalues, eigenvectors) = symmetric_eigen_descending(&gram, n_features);

    // svd_flip(u_based_decision=False): for each component k, choose the sign so
    // the largest-magnitude entry of the eigenvector column V[:, k] (a row of
    // sklearn's components_) is positive. numpy `argmax` takes the first index on
    // ties, mirrored by the strict `>` comparison.
    let mut signs = vec![1.0_f64; keep];
    for (k, sign) in signs.iter_mut().enumerate() {
        let mut best_feature = 0usize;
        let mut best_abs = -1.0_f64;
        for (feature, row) in eigenvectors.iter().enumerate() {
            let magnitude = row[k].abs();
            if magnitude > best_abs {
                best_abs = magnitude;
                best_feature = feature;
            }
        }
        // np.sign(0.0) == 0.0, but an all-zero eigenvector cannot be the argmax of
        // a non-degenerate component; default to +1 in that pathological case.
        if eigenvectors[best_feature][k] < 0.0 {
            *sign = -1.0;
        }
    }

    // Scores T = X_centered V (eigenvectors are V's columns), with the sign flip
    // folded in. Component k of sample s is sign_k * sum_f centered[s][f] * V[f][k].
    let mut scores = vec![vec![0.0_f64; keep]; n_samples];
    for (s, centered_row) in centered.iter().enumerate() {
        for (k, score) in scores[s].iter_mut().enumerate() {
            let mut acc = 0.0;
            for (f, &value) in centered_row.iter().enumerate() {
                acc += value * eigenvectors[f][k];
            }
            *score = acc * signs[k];
        }
    }

    scores
}

/// Eigendecomposition of a symmetric matrix by cyclic Jacobi rotations, returning
/// `(eigenvalues, eigenvectors)` sorted by descending eigenvalue. For a symmetric
/// PSD Gram matrix this matches `numpy.linalg.svd`'s right singular vectors and
/// squared singular values. Deterministic and order-independent (the pivot is the
/// largest-magnitude off-diagonal entry). Mirrors the routine
/// `mokume-imputation` uses for BPCA so the two stay numerically consistent.
fn symmetric_eigen_descending(matrix: &[Vec<f64>], dimension: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut a = matrix.to_vec();
    let mut vectors = identity(dimension);

    for _ in 0..JACOBI_SWEEPS {
        let mut off = 0.0;
        let mut pivot_p = 0usize;
        let mut pivot_q = if dimension > 1 { 1 } else { 0 };
        for (i, row_i) in a.iter().enumerate() {
            for (j, &entry) in row_i.iter().enumerate().skip(i + 1) {
                let value = entry.abs();
                if value > off {
                    off = value;
                    pivot_p = i;
                    pivot_q = j;
                }
            }
        }
        if off < JACOBI_EPS {
            break;
        }

        let app = a[pivot_p][pivot_p];
        let aqq = a[pivot_q][pivot_q];
        let apq = a[pivot_p][pivot_q];
        if apq.abs() < JACOBI_EPS {
            break;
        }
        let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
        let (sin_phi, cos_phi) = phi.sin_cos();

        for row in a.iter_mut() {
            let aip = row[pivot_p];
            let aiq = row[pivot_q];
            row[pivot_p] = cos_phi * aip - sin_phi * aiq;
            row[pivot_q] = sin_phi * aip + cos_phi * aiq;
        }
        let (lower_rows, upper_rows) = a.split_at_mut(pivot_q);
        let row_p = &mut lower_rows[pivot_p];
        let row_q = &mut upper_rows[0];
        for column in 0..dimension {
            let api = row_p[column];
            let aqi = row_q[column];
            row_p[column] = cos_phi * api - sin_phi * aqi;
            row_q[column] = sin_phi * api + cos_phi * aqi;
        }
        for row in vectors.iter_mut() {
            let vip = row[pivot_p];
            let viq = row[pivot_q];
            row[pivot_p] = cos_phi * vip - sin_phi * viq;
            row[pivot_q] = sin_phi * vip + cos_phi * viq;
        }
    }

    let eigenvalues = (0..dimension).map(|i| a[i][i]).collect::<Vec<f64>>();
    let mut order = (0..dimension).collect::<Vec<usize>>();
    order.sort_by(|&left, &right| {
        eigenvalues[right]
            .partial_cmp(&eigenvalues[left])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let sorted_values = order.iter().map(|&i| eigenvalues[i]).collect::<Vec<f64>>();
    let sorted_vectors = (0..dimension)
        .map(|row| {
            order
                .iter()
                .map(|&col| vectors[row][col])
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<_>>();
    (sorted_values, sorted_vectors)
}

/// `dimension x dimension` identity used to accumulate the Jacobi rotations.
fn identity(dimension: usize) -> Vec<Vec<f64>> {
    (0..dimension)
        .map(|i| {
            (0..dimension)
                .map(|j| if i == j { 1.0 } else { 0.0 })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pca_clamps_components_and_centers() {
        // Two features, three samples; rank is at most 2.
        let data = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 9.0]];
        let scores = pca_scores(&data, 5);
        // keep = min(5, 3 samples, 2 features) = 2.
        assert_eq!(scores.len(), 3);
        assert!(scores.iter().all(|row| row.len() == 2));
        // Centered scores sum to (near) zero per component.
        for k in 0..2 {
            let sum: f64 = scores.iter().map(|row| row[k]).sum();
            assert!(sum.abs() < 1e-9, "component {k} not centered: {sum}");
        }
    }

    #[test]
    fn pca_preserves_pairwise_distances_on_full_rank() {
        // With n_components >= rank, PCA is a rigid rotation: pairwise Euclidean
        // distances in score space equal those in centered feature space.
        let data = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 3.0],
        ];
        let scores = pca_scores(&data, 3);
        let dist = |a: &[f64], b: &[f64]| -> f64 {
            a.iter()
                .zip(b.iter())
                .map(|(x, y)| (x - y).powi(2))
                .sum::<f64>()
                .sqrt()
        };
        for i in 0..data.len() {
            for j in 0..data.len() {
                let original = dist(&data[i], &data[j]);
                let projected = dist(&scores[i], &scores[j]);
                assert!(
                    (original - projected).abs() < 1e-9,
                    "distance {i},{j} changed: {original} vs {projected}"
                );
            }
        }
    }
}
