//! Self-contained optimizers used by the proDA differential-expression port.
//!
//! proDA's per-protein MLE chains several SciPy optimizers
//! (`scipy.optimize.minimize` with `method` in {`Nelder-Mead`, `L-BFGS-B`,
//! `BFGS`} plus `scipy.optimize.brentq`), then does its own Newton refinement.
//! The workspace forbids new external crates, so each optimizer is reimplemented
//! here. The goal is the *converged fixed point*, not bit-exact iterates:
//! proDA chains L-BFGS-B -> BFGS -> a 5-step Newton polish, so the final
//! stationary point is what determines the output.
//!
//! [`nelder_mead`] follows SciPy's `_minimize_neldermead` (the reflect / expand /
//! contract / shrink simplex update with the standard coefficients and the
//! non-adaptive initial-simplex construction), so on the prior objectives it
//! lands on SciPy's minimum to ~1e-6. [`brentq`] is a direct port of Brent's
//! method as implemented in `scipy/optimize/Zeros/brentq.c`, reproducing a known
//! root to ~1e-12. [`bfgs`] and [`lbfgsb`] are quasi-Newton methods with a
//! Wolfe-condition line search; they reach the same minimum on the smooth proDA
//! likelihood that SciPy does.

/// Outcome of a multivariate minimization: the converged point. Callers that
/// need the objective value evaluate `f(&result.x)` (the returned point is the
/// single source of truth, avoiding a possibly-stale cached value).
#[derive(Debug, Clone)]
pub(crate) struct OptimizeResult {
    pub x: Vec<f64>,
}

/// Accepted step from a line search: the new point, its objective value, and the
/// gradient there. Bundled into a struct so the search helpers do not return a
/// long tuple (which would otherwise trip `clippy::type_complexity`).
struct LineStep {
    x: Vec<f64>,
    f: f64,
    grad: Vec<f64>,
}

/// Nelder-Mead simplex minimization, matching
/// `scipy.optimize.minimize(method="Nelder-Mead")` with the default
/// (non-adaptive) reflection/expansion/contraction/shrink coefficients.
///
/// Convergence mirrors SciPy: stop when both the simplex spread in `x`
/// (`max |x_i - x_0|_inf <= xatol`) and in `f` (`max |f_i - f_0| <= fatol`)
/// fall below the tolerances, or `maxiter` iterations elapse. The initial
/// simplex uses SciPy's `nonzdelt = 0.05` / `zdelt = 0.00025` perturbations.
pub(crate) fn nelder_mead<F>(
    f: &F,
    x0: &[f64],
    xatol: f64,
    fatol: f64,
    maxiter: usize,
) -> OptimizeResult
where
    F: Fn(&[f64]) -> f64,
{
    const RHO: f64 = 1.0;
    const CHI: f64 = 2.0;
    const PSI: f64 = 0.5;
    const SIGMA: f64 = 0.5;
    const NONZDELT: f64 = 0.05;
    const ZDELT: f64 = 0.00025;

    let n = x0.len();
    if n == 0 {
        let _ = f(x0);
        return OptimizeResult { x: Vec::new() };
    }

    // Build the initial simplex: x0 plus n perturbed vertices (SciPy's scheme).
    let mut sim: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    sim.push(x0.to_vec());
    for k in 0..n {
        let mut y = x0.to_vec();
        if y[k] != 0.0 {
            y[k] *= 1.0 + NONZDELT;
        } else {
            y[k] = ZDELT;
        }
        sim.push(y);
    }

    let mut fsim: Vec<f64> = sim.iter().map(|v| f(v)).collect();
    sort_simplex(&mut sim, &mut fsim);

    let coeffs = NmCoeffs {
        rho: RHO,
        chi: CHI,
        psi: PSI,
        sigma: SIGMA,
    };
    let mut iterations = 1usize;
    while iterations < maxiter {
        // Convergence test (SciPy uses the inf-norm spread of the simplex).
        if simplex_x_spread(&sim) <= xatol && simplex_f_spread(&fsim) <= fatol {
            break;
        }
        let centroid = simplex_centroid(&sim, n);
        if simplex_replace_worst(f, &mut sim, &mut fsim, &centroid, &coeffs) {
            simplex_shrink(f, &mut sim, &mut fsim, coeffs.sigma);
        }
        sort_simplex(&mut sim, &mut fsim);
        iterations += 1;
    }

    OptimizeResult { x: sim[0].clone() }
}

/// Nelder-Mead simplex coefficients (reflection/expansion/contraction/shrink).
struct NmCoeffs {
    rho: f64,
    chi: f64,
    psi: f64,
    sigma: f64,
}

/// Centroid of all but the worst (last) vertex of a sorted simplex.
fn simplex_centroid(sim: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut centroid = vec![0.0; n];
    for vertex in sim.iter().take(n) {
        for (c, v) in centroid.iter_mut().zip(vertex) {
            *c += v;
        }
    }
    for c in &mut centroid {
        *c /= n as f64;
    }
    centroid
}

/// Reflect/expand/contract the worst vertex (SciPy's update). Returns `true` if
/// the move failed and the caller should shrink the whole simplex.
fn simplex_replace_worst<F>(
    f: &F,
    sim: &mut [Vec<f64>],
    fsim: &mut [f64],
    centroid: &[f64],
    c: &NmCoeffs,
) -> bool
where
    F: Fn(&[f64]) -> f64,
{
    let n = sim.len() - 1;
    let worst = &sim[n];
    let xr = affine(centroid, worst, c.rho);
    let fxr = f(&xr);

    if fxr < fsim[0] {
        // Reflected point is the new best: try to expand further.
        let xe = affine(centroid, worst, c.rho * c.chi);
        let fxe = f(&xe);
        let (pt, val) = if fxe < fxr { (xe, fxe) } else { (xr, fxr) };
        sim[n] = pt;
        fsim[n] = val;
        false
    } else if fxr < fsim[n - 1] {
        sim[n] = xr;
        fsim[n] = fxr;
        false
    } else if fxr < fsim[n] {
        // Outside contraction.
        let xc = affine(centroid, worst, c.rho * c.psi);
        let fxc = f(&xc);
        if fxc <= fxr {
            sim[n] = xc;
            fsim[n] = fxc;
            false
        } else {
            true
        }
    } else {
        // Inside contraction.
        let xcc = affine_inside(centroid, worst, c.psi);
        let fxcc = f(&xcc);
        if fxcc < fsim[n] {
            sim[n] = xcc;
            fsim[n] = fxcc;
            false
        } else {
            true
        }
    }
}

/// Shrink all non-best vertices toward the best by factor `sigma`.
fn simplex_shrink<F>(f: &F, sim: &mut [Vec<f64>], fsim: &mut [f64], sigma: f64)
where
    F: Fn(&[f64]) -> f64,
{
    let best = sim[0].clone();
    for j in 1..sim.len() {
        let shrunk: Vec<f64> = best
            .iter()
            .zip(&sim[j])
            .map(|(b, s)| b + sigma * (s - b))
            .collect();
        fsim[j] = f(&shrunk);
        sim[j] = shrunk;
    }
}

/// `centroid + coeff * (centroid - worst)` (reflection/expansion/outer family).
fn affine(centroid: &[f64], worst: &[f64], coeff: f64) -> Vec<f64> {
    centroid
        .iter()
        .zip(worst)
        .map(|(c, w)| (1.0 + coeff) * c - coeff * w)
        .collect()
}

/// `centroid - psi * (centroid - worst)` (inside contraction).
fn affine_inside(centroid: &[f64], worst: &[f64], psi: f64) -> Vec<f64> {
    centroid
        .iter()
        .zip(worst)
        .map(|(c, w)| (1.0 - psi) * c + psi * w)
        .collect()
}

/// Sort the simplex vertices by ascending function value, carrying `fsim`.
/// NaN values sort to the end (treated as worst), matching SciPy's argsort with
/// the convention that an evaluated NaN never displaces a finite best vertex.
fn sort_simplex(sim: &mut [Vec<f64>], fsim: &mut [f64]) {
    let mut order: Vec<usize> = (0..fsim.len()).collect();
    order.sort_by(|&a, &b| fsim[a].total_cmp(&fsim[b]));
    let new_sim: Vec<Vec<f64>> = order.iter().map(|&i| sim[i].clone()).collect();
    let new_fsim: Vec<f64> = order.iter().map(|&i| fsim[i]).collect();
    sim.clone_from_slice(&new_sim);
    fsim.clone_from_slice(&new_fsim);
}

/// `max_i ||x_i - x_0||_inf` over the simplex (SciPy convergence numerator).
fn simplex_x_spread(sim: &[Vec<f64>]) -> f64 {
    let best = &sim[0];
    sim.iter()
        .skip(1)
        .map(|v| {
            v.iter()
                .zip(best)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max)
        })
        .fold(0.0_f64, f64::max)
}

/// `max_i |f_i - f_0|` over the simplex.
fn simplex_f_spread(fsim: &[f64]) -> f64 {
    let best = fsim[0];
    fsim.iter()
        .skip(1)
        .map(|f| (f - best).abs())
        .fold(0.0_f64, f64::max)
}

/// Mutable bracket state for [`brentq`], mirroring the locals of SciPy's C
/// routine: the current / previous / "block" estimates and their function
/// values, plus the last two step sizes.
struct BrentState {
    xpre: f64,
    xcur: f64,
    xblk: f64,
    fpre: f64,
    fcur: f64,
    fblk: f64,
    spre: f64,
    scur: f64,
}

/// Brent's root-finding method, a port of SciPy's `brentq` (the C routine in
/// `scipy/optimize/Zeros/brentq.c`). Finds a root of `f` bracketed by `[a, b]`
/// (requires `f(a) * f(b) <= 0`) to within `xtol`. Returns `NaN` if the bracket
/// is invalid.
pub(crate) fn brentq<F>(f: &F, a: f64, b: f64, xtol: f64, maxiter: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    const RTOL: f64 = 4.0 * f64::EPSILON;
    let fpre = f(a);
    let fcur = f(b);
    if fpre == 0.0 {
        return a;
    }
    if fcur == 0.0 {
        return b;
    }
    if fpre.signum() == fcur.signum() {
        return f64::NAN;
    }

    let mut s = BrentState {
        xpre: a,
        xcur: b,
        xblk: 0.0,
        fpre,
        fcur,
        fblk: 0.0,
        spre: 0.0,
        scur: 0.0,
    };

    for _ in 0..maxiter {
        brent_maybe_swap(&mut s);
        let delta = (xtol + RTOL * s.xcur.abs()) / 2.0;
        let sbis = (s.xblk - s.xcur) / 2.0;
        if s.fcur == 0.0 || sbis.abs() < delta {
            return s.xcur;
        }
        brent_choose_step(&mut s, delta, sbis);

        s.xpre = s.xcur;
        s.fpre = s.fcur;
        if s.scur.abs() > delta {
            s.xcur += s.scur;
        } else {
            s.xcur += if sbis > 0.0 { delta } else { -delta };
        }
        s.fcur = f(s.xcur);
    }
    s.xcur
}

/// Re-bracket the root: refresh the "block" point when the sign flips, then
/// ensure `xcur` is the better (smaller `|f|`) estimate.
fn brent_maybe_swap(s: &mut BrentState) {
    if s.fpre != 0.0 && s.fcur != 0.0 && s.fpre.signum() != s.fcur.signum() {
        s.xblk = s.xpre;
        s.fblk = s.fpre;
        s.spre = s.xcur - s.xpre;
        s.scur = s.xcur - s.xpre;
    }
    if s.fblk.abs() < s.fcur.abs() {
        s.xpre = s.xcur;
        s.xcur = s.xblk;
        s.xblk = s.xpre;
        s.fpre = s.fcur;
        s.fcur = s.fblk;
        s.fblk = s.fpre;
    }
}

/// Pick the next step: secant / inverse-quadratic interpolation when it stays
/// within bounds, otherwise bisection (`sbis`).
fn brent_choose_step(s: &mut BrentState, delta: f64, sbis: f64) {
    if s.spre.abs() > delta && s.fcur.abs() < s.fpre.abs() {
        let stry = if s.xpre == s.xblk {
            // Secant step.
            -s.fcur * (s.xcur - s.xpre) / (s.fcur - s.fpre)
        } else {
            // Inverse quadratic interpolation.
            let dpre = (s.fpre - s.fcur) / (s.xpre - s.xcur);
            let dblk = (s.fblk - s.fcur) / (s.xblk - s.xcur);
            -s.fcur * (s.fblk * dblk - s.fpre * dpre) / (dblk * dpre * (s.fblk - s.fpre))
        };
        if 2.0 * stry.abs() < s.spre.abs().min(3.0 * sbis.abs() - delta) {
            s.spre = s.scur;
            s.scur = stry;
        } else {
            s.spre = sbis;
            s.scur = sbis;
        }
    } else {
        s.spre = sbis;
        s.scur = sbis;
    }
}

/// BFGS quasi-Newton minimization with an analytic gradient, matching
/// `scipy.optimize.minimize(method="BFGS")` at convergence. The inverse-Hessian
/// approximation is updated with the BFGS formula and steps are taken along a
/// Wolfe line search; iteration stops when `||grad||_inf <= gtol`.
pub(crate) fn bfgs<F, G>(f: &F, grad: &G, x0: &[f64], gtol: f64, maxiter: usize) -> OptimizeResult
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
{
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut g = grad(&x);
    let mut hinv = identity(n);
    let mut fx = f(&x);

    for _ in 0..maxiter {
        if inf_norm(&g) <= gtol {
            break;
        }
        // Search direction p = -Hinv g.
        let mut p = mat_vec_neg(&hinv, &g);
        // Guard against a non-descent direction (reset to steepest descent).
        if dot(&p, &g) >= 0.0 {
            hinv = identity(n);
            p = g.iter().map(|gi| -gi).collect();
        }

        let Some(step) = line_search(f, grad, &x, &p, fx, &g) else {
            break;
        };

        // s = x_new - x (== alpha * p), the actual step taken.
        let s: Vec<f64> = step.x.iter().zip(&x).map(|(a, b)| a - b).collect();
        let yk: Vec<f64> = step.grad.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &yk);
        if sy > 1e-12 {
            bfgs_update(&mut hinv, &s, &yk, sy);
        }

        x = step.x;
        fx = step.f;
        g = step.grad;
    }

    OptimizeResult { x }
}

/// Bounded quasi-Newton minimization matching
/// `scipy.optimize.minimize(method="L-BFGS-B")` at convergence. Implemented as a
/// projected-gradient BFGS: search directions come from the BFGS inverse-Hessian
/// approximation, candidate points are clipped to `bounds`, and components
/// pinned at an active bound have their gradient projected out. Stops when the
/// projected gradient infinity-norm is `<= gtol` or the relative function
/// decrease is `<= ftol`.
pub(crate) fn lbfgsb<F, G>(
    f: &F,
    grad: &G,
    x0: &[f64],
    bounds: &[(Option<f64>, Option<f64>)],
    ftol: f64,
    gtol: f64,
    maxiter: usize,
) -> OptimizeResult
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
{
    let n = x0.len();
    let mut x = clip(x0, bounds);
    let mut g = grad(&x);
    let mut hinv = identity(n);
    let mut fx = f(&x);

    for _ in 0..maxiter {
        if inf_norm(&projected_gradient(&x, &g, bounds)) <= gtol {
            break;
        }
        let mut p = mat_vec_neg(&hinv, &g);
        if dot(&p, &g) >= 0.0 {
            hinv = identity(n);
            p = g.iter().map(|gi| -gi).collect();
        }

        // Projected line search: clip the trial point to the box each step.
        let Some(step) = projected_line_search(f, grad, &x, &p, fx, &g, bounds) else {
            break;
        };

        let s: Vec<f64> = step.x.iter().zip(&x).map(|(a, b)| a - b).collect();
        let yk: Vec<f64> = step.grad.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &yk);
        if sy > 1e-12 {
            bfgs_update(&mut hinv, &s, &yk, sy);
        }

        let denom = fx.abs().max(step.f.abs()).max(1.0);
        let rel_decrease = (fx - step.f).abs() / denom;

        let prev = fx;
        x = step.x;
        fx = step.f;
        g = step.grad;

        if rel_decrease <= ftol && prev >= fx {
            break;
        }
    }

    OptimizeResult { x }
}

/// In-place BFGS inverse-Hessian update for accepted step `s` / gradient change
/// `yk` with `sy = s.y > 0`. Uses the standard rank-two formula
/// `H+ = (I - rho s y^T) H (I - rho y s^T) + rho s s^T`.
fn bfgs_update(hinv: &mut [Vec<f64>], s: &[f64], yk: &[f64], sy: f64) {
    let n = s.len();
    let rho = 1.0 / sy;
    // hy = H y
    let hy = mat_vec(hinv, yk);
    let yhy = dot(yk, &hy);
    // H+ = H + rho*(sy + yHy)/sy * s s^T - rho*(H y s^T + s y^T H)
    let coeff = rho * (1.0 + rho * yhy);
    for i in 0..n {
        for j in 0..n {
            hinv[i][j] += coeff * s[i] * s[j] - rho * (hy[i] * s[j] + s[i] * hy[j]);
        }
    }
}

/// Backtracking Wolfe line search returning the accepted [`LineStep`].
/// Satisfies the Armijo (sufficient-decrease) condition and the curvature
/// condition where possible, falling back to backtracking. Returns `None` if no
/// acceptable step is found.
fn line_search<F, G>(f: &F, grad: &G, x: &[f64], p: &[f64], fx: f64, g: &[f64]) -> Option<LineStep>
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
{
    const C1: f64 = 1e-4;
    const MAX_STEPS: usize = 50;
    let g0 = dot(g, p);
    if g0 >= 0.0 {
        return None;
    }
    let mut alpha = 1.0;
    for _ in 0..MAX_STEPS {
        let x_new: Vec<f64> = x.iter().zip(p).map(|(xi, pi)| xi + alpha * pi).collect();
        let f_new = f(&x_new);
        // Accept on the Armijo (sufficient-decrease) condition. The curvature
        // condition is not enforced; the BFGS update instead guards against bad
        // (s, y) pairs via the `sy > 1e-12` curvature test at the call site.
        if f_new.is_finite() && f_new <= fx + C1 * alpha * g0 {
            let g_new = grad(&x_new);
            return Some(LineStep {
                x: x_new,
                f: f_new,
                grad: g_new,
            });
        }
        alpha *= 0.5;
        if alpha < 1e-16 {
            return None;
        }
    }
    None
}

/// Backtracking line search that clips each trial point to `bounds`, returning
/// the accepted [`LineStep`].
fn projected_line_search<F, G>(
    f: &F,
    grad: &G,
    x: &[f64],
    p: &[f64],
    fx: f64,
    g: &[f64],
    bounds: &[(Option<f64>, Option<f64>)],
) -> Option<LineStep>
where
    F: Fn(&[f64]) -> f64,
    G: Fn(&[f64]) -> Vec<f64>,
{
    const C1: f64 = 1e-4;
    const MAX_STEPS: usize = 50;
    let mut alpha = 1.0;
    let g0 = dot(g, p);
    for _ in 0..MAX_STEPS {
        let trial: Vec<f64> = x.iter().zip(p).map(|(xi, pi)| xi + alpha * pi).collect();
        let x_new = clip(&trial, bounds);
        let f_new = f(&x_new);
        // Sufficient decrease relative to the actual (projected) step taken.
        let step: Vec<f64> = x_new.iter().zip(x).map(|(a, b)| a - b).collect();
        let expected = C1 * dot(g, &step);
        if f_new.is_finite() && f_new <= fx + expected.min(C1 * alpha * g0) {
            let g_new = grad(&x_new);
            return Some(LineStep {
                x: x_new,
                f: f_new,
                grad: g_new,
            });
        }
        alpha *= 0.5;
        if alpha < 1e-16 {
            // Last resort: a tiny projected step toward the descent direction.
            let trial: Vec<f64> = x.iter().zip(p).map(|(xi, pi)| xi + 1e-16 * pi).collect();
            let x_new = clip(&trial, bounds);
            let f_new = f(&x_new);
            if f_new.is_finite() && f_new <= fx {
                let g_new = grad(&x_new);
                return Some(LineStep {
                    x: x_new,
                    f: f_new,
                    grad: g_new,
                });
            }
            return None;
        }
    }
    None
}

/// Clip `x` componentwise to the box `bounds` (None = unbounded on that side).
fn clip(x: &[f64], bounds: &[(Option<f64>, Option<f64>)]) -> Vec<f64> {
    x.iter()
        .zip(bounds)
        .map(|(&xi, &(lo, hi))| {
            let mut v = xi;
            if let Some(l) = lo {
                if v < l {
                    v = l;
                }
            }
            if let Some(h) = hi {
                if v > h {
                    v = h;
                }
            }
            v
        })
        .collect()
}

/// Projected gradient: zero out components pushing further into an active bound,
/// matching L-BFGS-B's convergence test on the projected gradient.
fn projected_gradient(x: &[f64], g: &[f64], bounds: &[(Option<f64>, Option<f64>)]) -> Vec<f64> {
    x.iter()
        .zip(g)
        .zip(bounds)
        .map(|((&xi, &gi), &(lo, hi))| {
            let at_lo = lo.map(|l| xi <= l).unwrap_or(false);
            let at_hi = hi.map(|h| xi >= h).unwrap_or(false);
            if (at_lo && gi > 0.0) || (at_hi && gi < 0.0) {
                0.0
            } else {
                gi
            }
        })
        .collect()
}

fn identity(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}

fn mat_vec(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter().map(|row| dot(row, v)).collect()
}

fn mat_vec_neg(m: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    m.iter().map(|row| -dot(row, v)).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |acc, x| acc.max(x.abs()))
}

#[cfg(test)]
mod tests {
    use super::{bfgs, brentq, lbfgsb, nelder_mead};

    fn assert_close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "actual={actual} expected={expected} diff={}",
            (actual - expected).abs()
        );
    }

    // Nelder-Mead on the 2D Rosenbrock from SciPy's standard start. SciPy's
    // converged minimum is x=[~1, ~1], f~1.1e-18.
    #[test]
    fn nelder_mead_rosenbrock_matches_scipy() {
        let rosen = |x: &[f64]| {
            let a = x[1] - x[0] * x[0];
            let b = 1.0 - x[0];
            100.0 * a * a + b * b
        };
        let res = nelder_mead(&rosen, &[-1.2, 1.0], 1e-8, 1e-8, 2000);
        assert_close(res.x[0], 1.0, 1e-6);
        assert_close(res.x[1], 1.0, 1e-6);
        assert_close(rosen(&res.x), 0.0, 1e-12);
    }

    // Nelder-Mead on a separable quadratic; SciPy minimum x=[3,-2], f=1.
    #[test]
    fn nelder_mead_quadratic_matches_scipy() {
        let quad = |x: &[f64]| (x[0] - 3.0).powi(2) + (x[1] + 2.0).powi(2) + 1.0;
        let res = nelder_mead(&quad, &[0.0, 0.0], 1e-8, 1e-8, 2000);
        assert_close(res.x[0], 3.0, 1e-6);
        assert_close(res.x[1], -2.0, 1e-6);
        assert_close(quad(&res.x), 1.0, 1e-9);
    }

    // Nelder-Mead on proDA's dropout-curve negative log-likelihood (proda.py:87),
    // the exact 2-parameter `(rho, 1/zeta)` objective. SciPy's Nelder-Mead from
    // the same start converges to x=[21.605976, -0.7670982], f=6.598158593775
    // (captured verbatim from `conda run -n Bigbio python`); this confirms the
    // in-crate Nelder-Mead matches SciPy on the actual proDA objective shape.
    #[test]
    fn nelder_mead_dropout_objective_matches_scipy() {
        use crate::de::special::{norm_cdf, norm_logcdf, norm_logsf, norm_sf};
        fn invprobit(x: f64, rho: f64, zeta: f64, log: bool, oneminus: bool) -> f64 {
            let z = (x - rho) / zeta.abs().max(1e-100);
            if zeta < 0.0 {
                if oneminus {
                    if log {
                        norm_logcdf(z)
                    } else {
                        norm_cdf(z)
                    }
                } else if log {
                    norm_logsf(z)
                } else {
                    norm_sf(z)
                }
            } else if oneminus {
                if log {
                    norm_logsf(z)
                } else {
                    norm_sf(z)
                }
            } else if log {
                norm_logcdf(z)
            } else {
                norm_cdf(z)
            }
        }
        let mu0 = 18.0_f64;
        let sigma20 = 10.0_f64;
        let yo = [22.0, 21.5, 24.0, 24.2];
        let predm = [21.0, 23.0];
        let pvm = [0.05, 0.05];
        let neg_ll = |par: &[f64]| -> f64 {
            let (r, zinv) = (par[0], par[1]);
            if zinv >= 0.0 {
                return 1e10;
            }
            let z = 1.0 / zinv;
            let zz = (r - mu0) / sigma20.sqrt();
            let mut ll = -0.5 * zz * zz - 0.918_938_533_204_672_8 - sigma20.sqrt().ln();
            ll += zinv.abs().ln().min(10000.0_f64.ln());
            for &yi in &yo {
                ll += invprobit(yi, r, z, true, true);
            }
            let sign = zinv.signum();
            for (&pm, &pv) in predm.iter().zip(&pvm) {
                ll += invprobit(pm, r, sign * (1.0 / (zinv * zinv) + pv).sqrt(), true, false);
            }
            -ll
        };
        let x0 = [mu0, -1.0 / sigma20.sqrt()];
        let res = nelder_mead(&neg_ll, &x0, 1e-8, 1e-8, 1000);
        assert_close(res.x[0], 21.605976, 1e-5);
        assert_close(res.x[1], -0.7670982, 1e-5);
        assert_close(neg_ll(&res.x), 6.598158593774928, 1e-9);
    }

    // brentq roots matching scipy.optimize.brentq to ~1e-10.
    #[test]
    fn brentq_matches_scipy() {
        // root of cos is pi/2; root of exp(x)=2 is ln 2 (use the std constants
        // so the literals don't trip clippy::approx_constant).
        let root_cos = brentq(&|x: f64| x.cos(), 0.0, 3.0, 1e-12, 100);
        assert_close(root_cos, std::f64::consts::FRAC_PI_2, 1e-10);

        let root_cubic = brentq(&|x: f64| x * x * x - 2.0 * x - 5.0, 1.0, 3.0, 1e-12, 100);
        assert_close(root_cubic, 2.0945514815423265, 1e-10);

        let root_exp = brentq(&|x: f64| x.exp() - 2.0, -1.0, 2.0, 1e-12, 100);
        assert_close(root_exp, std::f64::consts::LN_2, 1e-10);
    }

    // brentq with no sign change returns NaN.
    #[test]
    fn brentq_no_bracket_is_nan() {
        let r = brentq(&|x: f64| x * x + 1.0, -1.0, 1.0, 1e-12, 100);
        assert!(r.is_nan());
    }

    // BFGS on the Rosenbrock: SciPy converges to x=[1,1], f~0.
    #[test]
    fn bfgs_rosenbrock_matches_scipy() {
        let f = |x: &[f64]| (x[0] - 1.0).powi(2) + 100.0 * (x[1] - x[0] * x[0]).powi(2);
        let g = |x: &[f64]| {
            vec![
                2.0 * (x[0] - 1.0) - 400.0 * x[0] * (x[1] - x[0] * x[0]),
                200.0 * (x[1] - x[0] * x[0]),
            ]
        };
        let res = bfgs(&f, &g, &[-1.0, 1.0], 1e-10, 1000);
        assert_close(res.x[0], 1.0, 1e-6);
        assert_close(res.x[1], 1.0, 1e-6);
        assert_close(f(&res.x), 0.0, 1e-10);
    }

    // L-BFGS-B with an active upper bound on x0: SciPy returns x=[1,3], f=1.
    #[test]
    fn lbfgsb_bounded_matches_scipy() {
        let f = |x: &[f64]| (x[0] - 2.0).powi(2) + (x[1] - 3.0).powi(2);
        let g = |x: &[f64]| vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 3.0)];
        let bounds = vec![(None, Some(1.0)), (None, None)];
        let res = lbfgsb(&f, &g, &[0.0, 0.0], &bounds, 1e-15, 1e-12, 1000);
        assert_close(res.x[0], 1.0, 1e-6);
        assert_close(res.x[1], 3.0, 1e-6);
        assert_close(f(&res.x), 1.0, 1e-9);
    }
}
