"""
Adaptive-pi0 false discovery rate: Storey q-values and pi0 reliability.

Estimates the proportion of true null hypotheses (pi0) from the p-value
distribution and turns it into Storey q-values, which relax the conservative
Benjamini-Hochberg procedure when pi0 < 1 (many true alternatives). BH is the
special case pi0 = 1; when pi0 cannot be trusted the caller falls back to BH.

This is a native re-implementation of the R ``qvalue`` package's ``pi0est`` and
``qvalue``. It deliberately keeps only the pieces mokume needs and exposes the
pi0 reliability signals that gate whether the adaptive procedure is used at all.

Port alignment: this port was compared against R ``qvalue`` 2.42.0 through rpy2
on real spike-in data (13 PXD datasets, 77 differential-expression tables, 72 of
them comparable — R errors out on the other 5). Holding pi0 fixed at the R
estimate, :func:`qvalues` is a line-for-line equivalent of ``qvalue::qvalue``:
max |diff| = 2.2e-16 (double rounding) and the q <= 0.05 discovery set is
identical on 72/72 tables. :func:`estimate_pi0` with ``method="smoother"`` agrees
with ``qvalue::pi0est`` to max |diff| = 7.2e-05 / median 4.8e-06; that residual
comes from how each side pins the spline's effective df to 3 (R's Reinsch solver
vs the bisection in :func:`_smooth_spline_df`), i.e. a convergence tolerance, and
it changes no discovery. The environment, commands and per-table raw numbers are
recorded in ``docs/internal/port-alignment.md`` (sections 0-2).

References
----------
Storey JD. A direct approach to false discovery rates. *J. R. Stat. Soc. B*.
2002;64(3):479-498.
Storey JD, Taylor JE, Siegmund D. Strong control, conservative point estimation
and simultaneous conservative consistency of false discovery rates.
*J. R. Stat. Soc. B*. 2004;66(1):187-205.
Storey JD, Tibshirani R. Statistical significance for genome-wide studies.
*PNAS*. 2003;100(16):9440-9445.
"""

from __future__ import annotations

import numpy as np

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.adaptive_fdr")

# qvalue's default lambda grid for pi0 estimation.
_DEFAULT_LAMBDAS = np.arange(0.05, 0.96, 0.05)


def _multipletests():
    """Lazily import statsmodels' ``multipletests`` (needs ``mokume-py[analysis]``)."""
    try:
        from statsmodels.stats.multitest import (  # pylint: disable=import-outside-toplevel
            multipletests,
        )
    except ImportError as exc:
        raise ImportError(
            "statsmodels is required for differential-expression FDR correction. "
            "Install it with: pip install mokume-py[analysis]"
        ) from exc
    return multipletests


def _finite_pvalues(pvalues: np.ndarray) -> np.ndarray:
    """Return the finite p-values in [0, 1] as a float array."""
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if p.size and (p.min() < 0 or p.max() > 1):
        raise ValueError("p-values must lie in [0, 1]")
    return p


def _pi0_curve(p: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    """pi0(lambda) = mean(p >= lambda) / (1 - lambda) for each lambda."""
    return np.array([np.mean(p >= lam) / (1.0 - lam) for lam in lambdas])


def estimate_pi0(
    pvalues: np.ndarray,
    lambdas: np.ndarray | None = None,
    method: str = "smoother",
    smooth_df: int = 3,
    conservative_bound: bool = True,
) -> float:
    """Estimate the true-null proportion pi0 from p-values.

    Parameters
    ----------
    method : {"smoother", "bootstrap", "fixed"}
        ``smoother`` (qvalue default) fits a cubic smoothing spline to
        pi0(lambda) and evaluates it at the largest lambda; ``bootstrap`` picks
        the lambda minimising the bootstrap MSE (Storey 2004); ``fixed`` uses a
        single lambda (0.5) — the exact, trivially-alignable estimator.
    conservative_bound : bool
        Raise the estimate to :func:`pi0_lower_bound` when it falls below it (on
        by default). This can only make pi0 larger, i.e. the resulting FDR
        control weaker-but-safer, never more aggressive. See
        :func:`pi0_lower_bound` for the derivation.

    Returns
    -------
    float
        pi0 clamped to (0, 1].
    """
    return _estimate_pi0_detail(
        pvalues,
        lambdas=lambdas,
        method=method,
        smooth_df=smooth_df,
        conservative_bound=conservative_bound,
    )["pi0"]


def _estimate_pi0_detail(
    pvalues: np.ndarray,
    lambdas: np.ndarray | None = None,
    method: str = "smoother",
    smooth_df: int = 3,
    conservative_bound: bool = True,
) -> dict:
    """pi0 estimation with the raw estimate, the lower bound and which one won."""
    p = _finite_pvalues(pvalues)
    if p.size == 0:
        return {
            "pi0": 1.0,
            "pi0_raw": 1.0,
            "pi0_lower_bound": 1.0,
            "bound_applied": False,
        }
    lambdas = _DEFAULT_LAMBDAS if lambdas is None else np.asarray(lambdas, dtype=float)

    if method == "fixed" or lambdas.size == 1:
        lam = float(lambdas[-1]) if lambdas.size == 1 else 0.5
        pi0 = float(np.mean(p >= lam) / (1.0 - lam))
    elif method == "bootstrap":
        pi0 = _pi0_bootstrap(p, lambdas)
    elif method == "smoother":
        pi0 = _pi0_smoother(p, lambdas, smooth_df)
    else:
        raise ValueError(f"unknown pi0 method {method!r}")

    pi0_raw = float(min(max(pi0, np.finfo(float).tiny), 1.0))
    if not conservative_bound:
        return {
            "pi0": pi0_raw,
            "pi0_raw": pi0_raw,
            "pi0_lower_bound": float("nan"),
            "bound_applied": False,
        }

    bound = pi0_lower_bound(p, lambdas=lambdas)
    applied = bound > pi0_raw
    if applied:
        logger.info(
            "pi0 estimate %.4f below its conservative lower bound %.4f -> using the bound",
            pi0_raw,
            bound,
        )
    return {
        "pi0": float(max(pi0_raw, bound)),
        "pi0_raw": pi0_raw,
        "pi0_lower_bound": bound,
        "bound_applied": bool(applied),
    }


def pi0_lower_bound(pvalues: np.ndarray, lambdas: np.ndarray | None = None) -> float:
    """Data-derived conservative lower bound for pi0: min of the pi0(lambda) curve.

    Every point of the empirical curve ``pi0(lambda) = #{p >= lambda} / ((1 -
    lambda) m)`` is an upward-biased estimator of pi0: null p-values are
    (super-)uniform, so ``E[#{p >= lambda}] >= (1 - lambda) m0`` for every lambda
    (Storey 2002, Thm 1; Storey, Taylor & Siegmund 2004, Thm 1). The minimum over
    the lambda grid is therefore the lowest value the data itself supports, and
    is used here as the bound. It is not a magic constant: it is recomputed from
    each dataset, and it tracks genuinely low pi0 (strong signal keeps the whole
    curve low, so the bound moves down with it).

    The bound only ever binds against the ``smoother`` estimator, which is the
    only one that can leave the curve: it extrapolates the spline to lambda -> 1,
    where the curve is noisiest, and can land below every observed point — an
    artefact of the fit rather than evidence, which would hand the adaptive FDR
    an unsupported rejection budget. The ``fixed`` and ``bootstrap`` estimators
    return curve points, so they are bounded by construction and unaffected.

    Selection over the grid biases the minimum downward relative to any single
    point, so the bound is deliberately weak — but not rare. On real spike-in
    data (13 PXD datasets, 77 differential-expression tables; see
    ``docs/internal/port-alignment.md`` section 2.2) it binds on 31% of tables
    (24/77) and then raises pi0 by a median of +0.0098 (max +0.0902) in absolute
    terms, i.e. a median factor of 1.054 (max 3.258). The bound only ever makes
    FDR control more conservative: it removes the part of the estimate that no
    lambda supports and never overrides an estimate the data does support.

    Returns
    -------
    float
        The bound in (0, 1]. Callers apply it as ``max(pi0_hat, bound)``, which
        can only make the procedure more conservative, never more aggressive.
    """
    p = _finite_pvalues(pvalues)
    if p.size == 0:
        return 1.0
    lambdas = _DEFAULT_LAMBDAS if lambdas is None else np.asarray(lambdas, dtype=float)
    bound = float(_pi0_curve(p, lambdas).min())
    return float(min(max(bound, np.finfo(float).tiny), 1.0))


def _pi0_smoother(p: np.ndarray, lambdas: np.ndarray, smooth_df: int) -> float:
    """qvalue 'smoother': cubic smoothing spline (df) of pi0(lambda), value at max lambda.

    Matches R ``qvalue``'s default ``smooth.spline(lambda, pi0, df=3)`` by fitting
    a natural cubic smoothing spline whose effective degrees of freedom is pinned
    to ``smooth_df`` (Reinsch algorithm), then evaluating at the largest lambda.
    """
    pi0_lam = _pi0_curve(p, lambdas)
    fitted = _smooth_spline_df(lambdas, pi0_lam, float(smooth_df))
    return float(fitted[-1])


def _smooth_spline_df(x: np.ndarray, y: np.ndarray, df: float) -> np.ndarray:
    """Natural cubic smoothing spline with effective df pinned to ``df`` (Reinsch).

    Returns the fitted values at ``x``. The smoothing parameter is found by
    bisection so that trace((I + lam*K)^-1) == df, matching R's
    ``smooth.spline(..., df=df)``.
    """
    n = x.size
    if n <= 3 or df >= n:
        # Not enough points to smooth to this df; return the raw curve.
        return np.asarray(y, dtype=float)
    h = np.diff(x)
    # Q: n x (n-2); R: (n-2) x (n-2) tridiagonal (Green & Silverman 1994).
    q = np.zeros((n, n - 2))
    r = np.zeros((n - 2, n - 2))
    for j in range(1, n - 1):
        c = j - 1
        q[j - 1, c] = 1.0 / h[j - 1]
        q[j, c] = -(1.0 / h[j - 1] + 1.0 / h[j])
        q[j + 1, c] = 1.0 / h[j]
        r[c, c] = (h[j - 1] + h[j]) / 3.0
        if c + 1 < n - 2:
            r[c, c + 1] = r[c + 1, c] = h[j] / 6.0
    k = q @ np.linalg.solve(r, q.T)  # n x n penalty operator

    def fitted_and_df(lam: float):
        a = np.eye(n) + lam * k
        smoother = np.linalg.solve(a, np.eye(n))
        return smoother @ y, np.trace(smoother)

    # Bisection on log(lam): large lam -> df->2 (linear), lam->0 -> df->n.
    lo, hi = 1e-8, 1e8
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if fitted_and_df(mid)[1] > df:  # too wiggly -> more smoothing
            lo = mid
        else:
            hi = mid
    return fitted_and_df(np.sqrt(lo * hi))[0]


def _pi0_bootstrap(p: np.ndarray, lambdas: np.ndarray, n_boot: int = 100) -> float:
    """Storey (2004) bootstrap pi0: lambda minimising bootstrap MSE to min pi0."""
    pi0_lam = _pi0_curve(p, lambdas)
    min_pi0 = np.quantile(pi0_lam, 0.1)
    m = p.size
    rng = np.random.default_rng(0)
    mse = np.zeros(lambdas.size)
    for _ in range(n_boot):
        boot = p[rng.integers(0, m, m)]
        pi0_boot = _pi0_curve(boot, lambdas)
        mse += (pi0_boot - min_pi0) ** 2
    return float(pi0_lam[int(np.argmin(mse))])


def qvalues(pvalues: np.ndarray, pi0: float | None = None) -> np.ndarray:
    """Storey q-values, aligned to non-finite p-values as NaN.

    q(p_(m)) = pi0 * p_(m); q(p_(i)) = min(pi0 * m * p_(i) / i, q(p_(i+1))).
    """
    p_all = np.asarray(pvalues, dtype=float)
    finite = np.isfinite(p_all)
    q = np.full(p_all.shape, np.nan)
    p = p_all[finite]
    if p.size == 0:
        return q
    if pi0 is None:
        pi0 = estimate_pi0(p)

    m = p.size
    order = np.argsort(p)
    ranks = np.arange(1, m + 1)
    q_sorted = pi0 * m * p[order] / ranks
    # Enforce monotonicity from the largest p-value downward.
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.minimum(q_sorted, 1.0)

    q_finite = np.empty(m)
    q_finite[order] = q_sorted
    q[finite] = q_finite
    return q


def pi0_reliability(
    pvalues: np.ndarray,
    pi0: float | None = None,
    boundary_tol: float = 1e-3,
) -> dict:
    """Signals that gate whether adaptive-pi0 FDR should be used, else fall back to BH.

    Adaptive-pi0 degenerates to BH when pi0 = 1 (a safe default), but is
    untrustworthy when pi0 hits the boundary, the null p-value distribution is
    not anti-conservative, or there are too few hypotheses.

    Returns a dict with ``reliable`` (bool) plus the individual diagnostics, so
    the caller can log why it fell back. ``pi0`` is the value to use downstream
    and already carries the conservative lower bound of :func:`pi0_lower_bound`;
    ``pi0_raw``, ``pi0_lower_bound`` and ``bound_applied`` expose whether that
    bound actually bound. When ``pi0`` is supplied by the caller it is taken as
    given: the bound is not recomputed and those three keys report it unbounded
    (``pi0_lower_bound`` is NaN).
    """
    p = _finite_pvalues(pvalues)
    m = p.size
    if pi0 is None:
        detail = (
            _estimate_pi0_detail(p)
            if m
            else {
                "pi0": 1.0,
                "pi0_raw": 1.0,
                "pi0_lower_bound": 1.0,
                "bound_applied": False,
            }
        )
    else:
        detail = {
            "pi0": float(pi0),
            "pi0_raw": float(pi0),
            "pi0_lower_bound": float("nan"),
            "bound_applied": False,
        }
    pi0 = detail["pi0"]

    at_boundary = pi0 >= 1.0 - boundary_tol
    too_few = m < 100  # pi0 estimation is unstable with few hypotheses
    anticonservative = _is_anticonservative(p) if m else False

    reliable = (m > 0) and (not at_boundary) and (not too_few) and anticonservative
    return {
        "reliable": bool(reliable),
        "pi0": float(pi0),
        "pi0_raw": float(detail["pi0_raw"]),
        "pi0_lower_bound": float(detail["pi0_lower_bound"]),
        "bound_applied": bool(detail["bound_applied"]),
        "n": int(m),
        "at_boundary": bool(at_boundary),
        "too_few": bool(too_few),
        "anticonservative": bool(anticonservative),
    }


def _is_anticonservative(p: np.ndarray) -> bool:
    """True if the p-value histogram is anti-conservative (a spike near 0).

    A crude, robust check: the density in the first bin [0, 0.05) exceeds the
    density in the flat right tail [0.5, 1.0). Only then is pi0 trustworthy.
    """
    left = np.mean(p < 0.05) / 0.05
    right = np.mean(p >= 0.5) / 0.5
    return left > right * 1.1


def adjust_pvalues(
    pvalues: np.ndarray,
    method: str = "bh",
    alpha: float = 0.05,
) -> tuple[np.ndarray, str]:
    """Return (adjusted, method_used): dispatch bh/bky/storey with pi0-reliability fallback.

    Parameters
    ----------
    method : {"bh", "bky", "storey"}
        ``bh`` Benjamini-Hochberg; ``bky`` Benjamini-Krieger-Yekutieli two-stage;
        ``storey`` Storey q-values. The two adaptive-pi0 procedures recover the
        FDR budget BH wastes when pi0 < 1, but are untrustworthy when pi0 is
        unreliable (see :func:`pi0_reliability`) — there they fall back to BH.
        Any other method (e.g. ``ihw``, which this function does not implement)
        also falls back to BH, so callers with a wider method space can route
        everything through here.
    alpha : float
        Target FDR level. Only ``bky`` consumes it; BH-adjusted p-values and
        Storey q-values do not depend on the level.

    Returns
    -------
    (numpy.ndarray, str)
        The adjusted p-values (Storey q-values for ``storey``), aligned to the
        input with NaN wherever the input p-value was not finite, and the method
        actually applied — which differs from ``method`` on any fallback, so the
        caller can report what it really got.
    """
    p_all = np.asarray(pvalues, dtype=float)
    finite = np.isfinite(p_all)
    adjusted = np.full(p_all.shape, np.nan)
    method = method.lower()
    if not finite.any():
        return adjusted, method

    p = p_all[finite]
    if method not in ("bh", "bky", "storey"):
        logger.warning("Unknown FDR method %r -> falling back to BH", method)
        method = "bh"

    if method in ("bky", "storey"):
        rel = pi0_reliability(p)
        if rel["reliable"]:
            if method == "storey":
                adjusted[finite] = qvalues(p, pi0=rel["pi0"])
            else:
                adjusted[finite] = _multipletests()(p, method="fdr_tsbky", alpha=alpha)[
                    1
                ]
            return adjusted, method
        logger.info(
            "Adaptive FDR (%s) not trustworthy (pi0=%.3f n=%d boundary=%s "
            "too_few=%s anticons=%s) -> falling back to BH",
            method,
            rel["pi0"],
            rel["n"],
            rel["at_boundary"],
            rel["too_few"],
            rel["anticonservative"],
        )
        method = "bh"

    adjusted[finite] = _multipletests()(p, method="fdr_bh")[1]
    return adjusted, method
