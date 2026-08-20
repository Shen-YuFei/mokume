"""
Data-driven effect-size gate for differential expression.

Estimates a fold-change threshold from the data's own null distribution instead
of hard-coding it. The classic fixed ``|log2FC| > 0.5`` gate is calibrated for
label-free 2x designs; it is miscalibrated for isobaric labelling (TMT/iTRAQ),
where reporter-ion co-isolation compresses observed ratios toward 1:1 and pushes
true positives below the fixed gate. Estimating the gate from the null width lets
it track the data (a tighter, compressed null -> a lower gate).

The estimator fits a two-component mixture (null + signal) to the log2 fold
changes and returns the crossover, i.e. the ``|log2FC|`` at which a protein is
more likely to belong to the signal component than the null. This is the
local-FDR = 0.5 boundary of the mixture model, in the spirit of the R
``fdrtool`` / ``locfdr`` local-fdr thresholding, but native and dependency-free.

The gate is a *recommendation* with strong reasoning; it is not backed by a
single citation, and the local-fdr target level it implies is still a soft
user preference. Specificity is ultimately controlled by the FDR procedure, not
this gate, so the estimate is used as-is, with no floor: clamping it from below
would silently override the null width the data reports, which is the very
miscalibration this module exists to remove (a null narrow enough to put the
crossover at 0.1 is what compressed isobaric ratios look like, not an error).
The fixed ``fallback`` applies only when the null cannot be characterised at
all, and says so in the log.

References (methodological lineage, not a direct citation for the compound claim)
--------------------------------------------------------------------------------
Efron B. Local False Discovery Rates. 2005.
Strimmer K. fdrtool: local and tail-area FDR. Bioinformatics. 2008.
"""

from __future__ import annotations

import numpy as np

from mokume.core.logger import get_logger

logger = get_logger("mokume.analysis.effect_size_gate")


def estimate_effect_size_gate(
    log2fc: np.ndarray,
    method: str = "mixture",
    fallback: float = 0.5,
) -> float:
    """Estimate a |log2FC| gate from the fold-change distribution.

    The estimate is returned unclamped: there is no lower bound on the gate,
    because a narrow null (compressed isobaric ratios) legitimately implies a
    small gate, and specificity is controlled by the FDR procedure rather than
    here. ``fallback`` is not an estimate but a last-resort default, used only
    when the null cannot be characterised; every use of it is logged as a
    warning, so a caller that sees no warning knows the gate is data-derived.

    Parameters
    ----------
    method : {"mixture", "null_quantile"}
        ``mixture`` fits a two-component Gaussian mixture and returns the
        null/signal crossover (local-fdr = 0.5). ``null_quantile`` returns a high
        quantile of the null-centred |log2FC| (a robust, assumption-light
        alternative).
    fallback : float
        Returned when the gate cannot be estimated at all — too few points, a
        mixture fit that fails or finds no separable signal component, or a
        degenerate (non-finite or non-positive) estimate. Not data-derived: it
        is the conventional label-free default, and is miscalibrated for
        isobaric data by exactly the reason this module exists.

    Returns
    -------
    float
        The estimated |log2FC| gate, or ``fallback`` if it is not estimable.
    """
    x = np.asarray(log2fc, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 50:
        logger.warning(
            "Effect-size gate: only %d finite log2FC values (<50), the null is not "
            "estimable -> using the fixed fallback gate %.3f",
            x.size,
            fallback,
        )
        return fallback

    if method == "null_quantile":
        gate = _null_quantile_gate(x)
    elif method == "mixture":
        gate = _mixture_gate(x)
    else:
        raise ValueError(f"unknown gate method {method!r}")

    if gate is None or not np.isfinite(gate) or gate <= 0:
        logger.warning(
            "Effect-size gate: the '%s' estimator returned a degenerate gate (%s); "
            "the fold-change null has no separable signal component or is not "
            "estimable -> using the fixed fallback gate %.3f",
            method,
            gate,
            fallback,
        )
        return fallback
    return float(gate)


def _null_quantile_gate(x: np.ndarray, quantile: float = 0.95) -> float:
    """Gate = high quantile of |log2FC| within the robust central null."""
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    if mad <= 0:
        return float(np.quantile(np.abs(x - med), quantile))
    # Central (null-like) proteins: within ~2 robust SD of the null centre.
    central = x[np.abs(x - med) <= 2.0 * mad]
    if central.size < 20:
        central = x
    return float(np.quantile(np.abs(central - med), quantile))


def _mixture_gate(x: np.ndarray) -> float | None:
    """Two-component mixture crossover on |log2FC| (null vs signal)."""
    try:
        from sklearn.mixture import GaussianMixture  # pylint: disable=import-outside-toplevel
    except ImportError:  # pragma: no cover
        return _null_quantile_gate(x)

    # Work on |log2FC| centred at the null median: component 0 = null (near 0),
    # component 1 = signal (larger magnitude).
    med = np.median(x)
    a = np.abs(x - med).reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=0)
    try:
        gmm.fit(a)
    except ValueError as exc:  # degenerate input: the mixture cannot be fitted
        logger.warning("Effect-size gate: mixture fit failed (%s)", exc)
        return None
    means = gmm.means_.ravel()
    null_c = int(np.argmin(means))
    sig_c = 1 - null_c
    if means[sig_c] <= means[null_c]:
        return None  # no separable signal component

    # Scan |log2FC| and find where the signal posterior first exceeds the null.
    grid = np.linspace(0, float(np.quantile(a, 0.999)), 400).reshape(-1, 1)
    post = gmm.predict_proba(grid)
    crossing = np.where(post[:, sig_c] >= post[:, null_c])[0]
    if crossing.size == 0:
        return None
    return float(grid[crossing[0], 0])
