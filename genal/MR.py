import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.stats import norm, chi2, binomtest, t
from scipy.special import log_ndtr
from scipy.optimize import minimize_scalar
from tqdm import tqdm

from .constants import MR_METHODS_NAMES

# Maximum number of (bootstrap replicate x SNP) values held in one array: the bootstrap replicates are processed in chunks
BOOTSTRAP_CHUNK_VALUES = 2**15
# Grid points per block of the kernel density evaluated on the grid of the mode methods, and per part of a block for
# the tighter upper bounds (see grid_density)
KDE_GRID_BLOCK = 64
KDE_GRID_SUB_BLOCK = 8


def _null_result(method_name):
    """Generate a null result dict for early returns when insufficient data."""
    return {
        "method": method_name,
        "b": np.nan,
        "se": np.nan,
        "pval": np.nan,
        "nSNP": np.nan,
    }


"""
Mode methods
"""

def mr_simple_mode(BETA_e, SE_e, BETA_o, SE_o, phi, nboot, cpus, show_progress=True, seed=None):
    """
    Perform a Mendelian Randomization analysis using the simple mode method.

    Args:
    BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
    SE_e (numpy array): Standard errors corresponding to `BETA_e`.
    BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
    SE_o (numpy array): Standard errors corresponding to `BETA_o`.
    phi (int): Factor for the bandwidth parameter used in the kernel density estimation
    nboot (int): Number of boostrap iterations to obtain the standard error and p-value
    cpus (int): Not used (the bootstrap runs on a single core).
    seed (int, optional): Random seed for the bootstrap draws. Default None (non-reproducible).

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Mode of the ratio estimates, representing the causal estimate.
            - "se": Bootstrap standard error (1.4826 * median absolute deviation of the bootstrap estimates).
            - "pval": P-value for the causal estimate (t distribution with nSNP - 1 degrees of freedom).
            - "nSNP": Number of genetic variants used in the analysis.
    """
    l = len(BETA_e)
    if l < 3:
        return [_null_result(MR_METHODS_NAMES["Simple-mode"])]

    BETA_IV   = BETA_o/BETA_e
    SE_IV = np.sqrt((SE_o**2) / (BETA_e**2) + ((BETA_o**2) * (SE_e**2)) / (BETA_e**4))

    BETA_mode = beta_mode(BETA_IV, np.ones(l), phi, method="Simple") #Point estimate calculated using unweighted KDE
    BETA_boot = bootstrap_mode(BETA_IV, SE_IV, phi, nboot, cpus, method="Simple", show_progress=show_progress, seed=seed) #Generating bootstrapped point estimates with normal sampling and unweighted KDE
    SE_mode = stats.median_abs_deviation(BETA_boot, scale=1/1.4826)
    P_mode = 2 * t.sf(np.abs(BETA_mode / SE_mode), df=l - 1)
    
    return [
    {
        "method": MR_METHODS_NAMES["Simple-mode"],
        "b": BETA_mode[0],
        "se": SE_mode[0],
        "pval": P_mode[0],
        "nSNP": l,
    }
]

def mr_weighted_mode(BETA_e, SE_e, BETA_o, SE_o, phi, nboot, cpus, show_progress=True, seed=None):
    """
    Perform a Mendelian Randomization analysis using the weighted mode method.

    Args:
    BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
    SE_e (numpy array): Standard errors corresponding to `BETA_e`.
    BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
    SE_o (numpy array): Standard errors corresponding to `BETA_o`.
    phi (int): Factor for the bandwidth parameter used in the kernel density estimation
    nboot (int): Number of boostrap iterations to obtain the standard error and p-value
    cpus (int): Not used (the bootstrap runs on a single core).
    seed (int, optional): Random seed for the bootstrap draws. Default None (non-reproducible).

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Weighted mode of the ratio estimates, representing the causal estimate.
            - "se": Bootstrap standard error (1.4826 * median absolute deviation of the bootstrap estimates).
            - "pval": P-value for the causal estimate (t distribution with nSNP - 1 degrees of freedom).
            - "nSNP": Number of genetic variants used in the analysis.
    """
    l = len(BETA_e)
    if l < 3:
        return [_null_result(MR_METHODS_NAMES["Weighted-mode"])]

    BETA_IV   = BETA_o/BETA_e
    SE_IV = np.sqrt((SE_o**2) / (BETA_e**2) + ((BETA_o**2) * (SE_e**2)) / (BETA_e**4))

    BETA_mode = beta_mode(BETA_IV, SE_IV, phi, method="Weighted") #Point estimate calculated using KDE weighted with standard errors
    BETA_boot = bootstrap_mode(BETA_IV, SE_IV, phi, nboot, cpus, method="Weighted", show_progress=show_progress, seed=seed) #Generating bootstrapped point estimates with normal sampling and weighted KDE
    SE_mode = stats.median_abs_deviation(BETA_boot, scale=1/1.4826)
    P_mode = 2 * t.sf(np.abs(BETA_mode / SE_mode), df=l - 1)
    
    return [
    {
        "method": MR_METHODS_NAMES["Weighted-mode"],
        "b": BETA_mode[0],
        "se": SE_mode[0],
        "pval": P_mode[0],
        "nSNP": l,
    }
]

def _kde(x, b_sorted, w_sorted, h, block=512):
    """
    Unnormalised (weighted) Gaussian kernel density with bandwidth h of the sorted values b_sorted (weights w_sorted),
    evaluated at the sorted points x, by blocks of `block` points. Kernels more than 10 bandwidths away from a block
    (relative contribution < 2e-22) are skipped.
    """
    out = np.empty(len(x))
    for start in range(0, len(x), block):
        xs = x[start:start + block]
        lo = np.searchsorted(b_sorted, xs[0] - 10 * h, side="left")
        hi = np.searchsorted(b_sorted, xs[-1] + 10 * h, side="right")
        out[start:start + block] = np.exp(-0.5 * ((xs[:, None] - b_sorted[None, lo:hi]) / h) ** 2) @ w_sorted[lo:hi]
    return out


def _kde_at(x, b_sorted, w_sorted, h):
    """Kernel density of :func:`_kde` at the single point x (same arithmetic as _kde(np.array([x]), ...)[0])."""
    lo = np.searchsorted(b_sorted, x - 10 * h, side="left")
    hi = np.searchsorted(b_sorted, x + 10 * h, side="right")
    return (np.exp(-0.5 * ((np.array([[x]]) - b_sorted[None, lo:hi]) / h) ** 2) @ w_sorted[lo:hi])[0]


def _density_bound(first, last, b, w, h):
    """
    Upper bound of the kernel density of :func:`_kde` on each interval [first, last]: each kernel within 10 bandwidths
    of the interval taken at its maximum over the interval (one exponential per kernel and interval).
    """
    lo = np.searchsorted(b, first - 10 * h, side="left")
    hi = np.searchsorted(b, last + 10 * h, side="right")
    counts = hi - lo
    interval = np.repeat(np.arange(len(first)), counts)
    kernel = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts) + np.repeat(lo, counts)
    distance = np.maximum(0, np.maximum(first[interval] - b[kernel], b[kernel] - last[interval]))
    return np.bincount(interval, weights=w[kernel] * np.exp(-0.5 * (distance / h) ** 2), minlength=len(first)) * (1 + 1e-9)


def grid_density(grid, b, w, h):
    """
    Kernel density (see :func:`_kde`) on the sorted grid of :func:`kde_mode`, by blocks of KDE_GRID_BLOCK points, where
    it can reach (1 - 1/400) of its maximum. Upper bounds of the density (see :func:`_density_bound`) on each block, then
    on the KDE_GRID_SUB_BLOCK-point parts of the blocks that are not excluded, identify the blocks that cannot: their
    points are set to -inf. The points that are evaluated get the values of the full evaluation, so the local maxima
    within (1 - 1/400) of the maximum are the same.
    """
    block, sub_block = KDE_GRID_BLOCK, KDE_GRID_SUB_BLOCK
    starts = np.arange(0, len(grid), block)
    bound = _density_bound(grid[starts], grid[np.minimum(starts + block, len(grid)) - 1], b, w, h)

    # Block with the highest bound first, then the blocks whose bound reaches (1 - 1/400) of the density found
    density = np.full(len(grid), -np.inf)
    point_block = np.arange(len(grid)) // block
    top = point_block == np.argmax(bound)
    density[top] = _kde(grid[top], b, w, h, block=block)
    threshold = (1 - 1 / 400) * density[top].max()
    candidates = np.flatnonzero(bound >= threshold)
    candidates = candidates[candidates != np.argmax(bound)]
    if len(candidates):
        # Tighter bound of each candidate block: the highest bound of its parts
        sub_starts = (starts[candidates][:, None] + np.arange(0, block, sub_block)).ravel()
        in_grid = sub_starts < len(grid)
        sub_bound = _density_bound(grid[sub_starts[in_grid]], grid[np.minimum(sub_starts[in_grid] + sub_block, len(grid)) - 1], b, w, h)
        block_bound = np.full(len(candidates), -np.inf)
        np.maximum.at(block_bound, np.repeat(np.arange(len(candidates)), block // sub_block)[in_grid], sub_bound)
        rest = np.isin(point_block, candidates[block_bound >= threshold])
        density[rest] = _kde(grid[rest], b, w, h, block=block)
    return density


def mode_bandwidth(BETA_IV, phi):
    """
    Bandwidth of the kernel density of the mode methods, one per row of BETA_IV: modified Silverman's rule proposed by
    Bickel (2002), 0.9 * min(SD, MAD) / n^(1/5) with MAD = 1.4826 * median absolute deviation, times phi (at least 1e-8).
    """
    std = np.std(BETA_IV, axis=-1, ddof=1)
    mad = np.median(np.abs(BETA_IV - np.median(BETA_IV, axis=-1, keepdims=True)), axis=-1) / (1 / 1.4826)
    s = 0.9 * np.where(mad < std, mad, std) / BETA_IV.shape[-1] ** (1 / 5)
    return np.where(s * phi > 0.00000001, s * phi, 0.00000001)


def mode_weights(SE_IV, method, n):
    """Standardised kernel weights of the mode methods: equal ("Simple") or proportional to SE_IV^-2 ("Weighted")."""
    if method == "Weighted":  # Using weighted KDE
        weights = np.asarray(SE_IV, dtype=float) ** -2
    elif method == "Simple":  # Using unweighted KDE
        weights = np.ones(n)
    else:
        raise ValueError("Method should be either 'Simple' or 'Weighted'.")
    return weights / np.sum(weights)


def kde_mode(b, w, h):
    """
    Global maximum of the Gaussian kernel density with bandwidth h of the sorted values b (standardised weights w).
    The density is evaluated on a grid of spacing h/10 covering the union of the intervals [b - 4h, b + 4h] (see
    :func:`grid_density`); every grid local maximum within 0.25% of the highest one (the largest possible grid error)
    is refined by a bounded scalar maximisation, and the highest refined maximum is returned. Distinct maxima of equal
    density (up to rounding errors, e.g. with symmetric data) are ties: the smallest of them is returned.

    Returns:
        tuple: (mode, other maxima tied with it)
    """
    # Grid of spacing h/10 on the union of the intervals [b - 4h, b + 4h] (one segment per group of close estimates)
    step = h / 10
    breaks = np.flatnonzero(np.diff(b) > 8 * h) + 1
    segments = [
        np.linspace(seg[0] - 4 * h, seg[-1] + 4 * h, int(np.ceil((seg[-1] - seg[0] + 8 * h) / step)) + 1)
        for seg in np.split(b, breaks)
    ]
    grid = np.concatenate(segments)
    density = grid_density(grid, b, w, h)

    # Local maxima of the grid density within each segment
    bounds = np.cumsum([0] + [len(seg) for seg in segments])
    peaks = []
    for start, stop in zip(bounds[:-1], bounds[1:]):
        d = np.concatenate([[-np.inf], density[start:stop], [-np.inf]])
        peaks.extend(start + np.flatnonzero((d[1:-1] >= d[:-2]) & (d[1:-1] > d[2:])))
    peaks = np.array(peaks)
    candidates = peaks[density[peaks] >= (1 - 1 / 400) * density.max()]
    candidates = candidates[np.argsort(-density[candidates])]

    # Bounded maximisation around each candidate; global maximum among the refined candidates
    mode, best_density = grid[candidates[0]], density[candidates[0]]
    refined = []
    for i in candidates:
        opt = minimize_scalar(
            lambda x: -_kde_at(x, b, w, h),
            bounds=(grid[i] - step, grid[i] + step),
            method="bounded",
            options={"xatol": 1e-9 * h},
        )
        refined.append((opt.x, -opt.fun))
        if -opt.fun > best_density:
            mode, best_density = opt.x, -opt.fun
    # Ties: other maxima (more than a grid step away) whose density equals the highest one up to rounding errors
    maxima = [mode]
    for x, value in refined:
        if value >= best_density * (1 - 1e-12) and all(abs(x - m) > step for m in maxima):
            maxima.append(x)
    maxima.sort()
    return maxima[0], maxima[1:]


def beta_mode(BETA_IV, SE_IV, phi, method):
    """
    Function to compute the point estimate of mode methods: the global maximum of the Gaussian kernel density of the
    ratio estimates (unweighted for "Simple", weighted by SE_IV^-2 for "Weighted"), see :func:`kde_mode`.
    """
    BETA_IV = np.asarray(BETA_IV, dtype=float)
    weights = mode_weights(SE_IV, method, len(BETA_IV))
    order = np.argsort(BETA_IV)
    mode, tied = kde_mode(BETA_IV[order], weights[order], float(mode_bandwidth(BETA_IV, phi)))
    if tied:
        warnings.warn(
            f"The {method.lower()} mode estimate is ambiguous: the kernel density has maxima of equal height at "
            f"{', '.join(f'{x:.6g}' for x in [mode] + tied)}. The smallest one ({mode:.6g}) is reported."
        )
    return np.array([mode])


def bootstrap_mode(BETA_IV, SE_IV, phi, nboot, cpus, method, show_progress=True, seed=None):
    """
    Function to obtain arrays of bootstrapped beta estimates (weighted or unweighted) to estimate SEs of mode methods.
    The bootstrap replicates of the ratio estimates are drawn from a single seeded generator, in chunks (bandwidths and
    sorting computed for a whole chunk at once). cpus is not used (single core).
    """
    BETA_IV, SE_IV = np.asarray(BETA_IV, dtype=float), np.asarray(SE_IV, dtype=float)
    weights = mode_weights(SE_IV, method, len(BETA_IV))
    rng = np.random.default_rng(seed)
    modes = np.empty((nboot, 1))
    with tqdm(total=nboot, desc=f"{method} mode bootstrapping", ncols=100, disable=not show_progress) as progress:
        for start, stop in bootstrap_chunks(nboot, len(BETA_IV)):
            # Replicates drawn as rng.normal(BETA_IV, SE_IV, size=(nboot, nSNP)), chunk by chunk
            BETA_IV_boot = rng.normal(BETA_IV, SE_IV, size=(stop - start, len(BETA_IV)))
            h = mode_bandwidth(BETA_IV_boot, phi)
            order = np.argsort(BETA_IV_boot, axis=1)
            b_sorted, w_sorted = np.take_along_axis(BETA_IV_boot, order, axis=1), weights[order]
            for row in range(stop - start):
                # Fresh arrays (aligned in memory like any new array): some BLAS libraries round the dot products of
                # the kernel density differently depending on the memory alignment of their inputs
                modes[start + row, 0] = kde_mode(b_sorted[row].copy(), w_sorted[row].copy(), float(h[row]))[0]
            progress.update(stop - start)
    return modes


"""
Egger regression methods
"""

def mr_egger_regression(BETA_e, SE_e, BETA_o, SE_o):
    """
    Perform a Mendelian Randomization analysis using the egger regression method. See :func:`mr_egger_regression_bootstrap` for a version with bootstrapping.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.

    Returns:
        list of dict: A list containing two dictionaries with the results for the egger regression estimate and the egger regression intercept (horizontal pleiotropy estimate):
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate or the intercept.
            - "se": Adjusted standard error of the coefficient or intercept.
            - "pval": P-value for the causal estimate or intercept.
            - "nSNP": Number of genetic variants used in the analysis.
    """
    null_result = [_null_result(name) for name in MR_METHODS_NAMES["Egger"]]

    l = len(BETA_e)
    # Early return if insufficient data
    if l < 3:
        return null_result

    BETA_e, SE_o, BETA_o = (np.asarray(x, dtype=float) for x in (BETA_e, SE_o, BETA_o))
    sign0 = np.where(BETA_e < 0, -1.0, 1.0)  # Sign of the exposure effect (0 counted as positive)

    X = sm.add_constant(np.abs(BETA_e))
    y = BETA_o * sign0
    weights = 1 / SE_o**2

    mod = sm.WLS(y, X, weights=weights).fit()

    if len(mod.params) > 1:
        # Standard errors sqrt(diag((X'WX)^-1)) * max(1, sigma): computed analytically so that a perfect fit (sigma = 0)
        # gives the fixed-effect standard errors
        b_i, b = np.asarray(mod.params)
        se_i, se = np.sqrt(np.diag(np.asarray(mod.normalized_cov_params))) * max(1, np.sqrt(mod.mse_resid))
        pval = 2 * t.sf(abs(b / se), l - 2)
        pval_i = 2 * t.sf(abs(b_i / se_i), l - 2)

        Q = mod.mse_resid * (l - 2)
        Q_df = l - 2
        Q_pval = chi2.sf(Q, Q_df)

        return [
            {
                "method": MR_METHODS_NAMES["Egger"][0],
                "b": b,
                "se": se,
                "pval": pval,
                "nSNP": l,
                "Q": Q,
                "Q_df": Q_df,
                "Q_pval": Q_pval,
            },
            {
                "method": MR_METHODS_NAMES["Egger"][1],
                "b": b_i,
                "se": se_i,
                "pval": pval_i,
                "nSNP": l,
                "Q": Q,
                "Q_df": Q_df,
                "Q_pval": Q_pval,
            },
        ]
    else:
        print(
            "Warning: Collinearities in MR Egger, try LD pruning the exposure (can be done with .clump())."
        )
        return null_result

def mr_egger_regression_bootstrap(BETA_e, SE_e, BETA_o, SE_o, nboot, seed=None):
    """
    Perform a Mendelian Randomization analysis using the egger regression method with boostrapped standard errors. See :func:`mr_egger_regression` for a version without bootstrapping.

    The estimates are the regular MR Egger slope and intercept on the observed data (:func:`mr_egger_regression`).
    Their standard errors are the standard deviations of nboot bootstrap replicates: each replicate draws new exposure
    and outcome effects from N(BETA, SE) and fits the same model (SNPs oriented to a positive exposure effect, outcome
    effects regressed on the exposure effects with an intercept, weighted by 1/SE_o^2).
    The bootstrap standard error reflects the sampling error of the SNP effects only, not the heterogeneity between
    variants: it can be smaller than the regular MR Egger standard error when the heterogeneity (Q) is large.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value
        seed (int, optional): Random seed for the bootstrap draws. Default None (non-reproducible).

    Returns:
        list of dict: A list containing two dictionaries with the results for the egger regression estimate and the egger regression intercept (horizontal pleiotropy estimate):
            - "method": Name of the analysis method.
            - "b": MR Egger slope (causal estimate) or intercept on the observed data.
            - "se": Standard deviation of the bootstrap estimates.
            - "pval": Two-sided p-value of b / se (normal distribution).
            - "nSNP": Number of genetic variants used in the analysis.
    """

    l = len(BETA_e)
    if l < 3:
        return [_null_result(name) for name in MR_METHODS_NAMES["Egger-boot"]]

    BETA_e, SE_e, BETA_o, SE_o = (np.asarray(x, dtype=float) for x in (BETA_e, SE_e, BETA_o, SE_o))
    weights = 1 / SE_o**2

    # Estimates on the observed data (regular MR Egger)
    egger = mr_egger_regression(pd.Series(BETA_e), pd.Series(SE_e), pd.Series(BETA_o), pd.Series(SE_o))
    b, b_i = egger[0]["b"], egger[1]["b"]

    rng = np.random.default_rng(seed)
    intercepts, slopes = np.empty(nboot), np.empty(nboot)
    for start, stop in bootstrap_chunks(nboot, l):
        xs, ys = draw_bootstrap_effects(rng, BETA_e, SE_e, BETA_o, SE_o, stop - start)
        intercepts[start:stop], slopes[start:stop] = egger_wls(xs, ys, weights)

    se, se_i = np.std(slopes, ddof=1), np.std(intercepts, ddof=1)

    return [
        {
            "method": MR_METHODS_NAMES["Egger-boot"][0],
            "b": b,
            "se": se,
            "pval": 2 * norm.sf(abs(b / se)),
            "nSNP": l,
        },
        {
            "method": MR_METHODS_NAMES["Egger-boot"][1],
            "b": b_i,
            "se": se_i,
            "pval": 2 * norm.sf(abs(b_i / se_i)),
            "nSNP": l,
        },
    ]


def egger_wls(BETA_e, BETA_o, weights):
    """
    Helper function returning the (intercept, slope) of the MR Egger model fitted by :func:`mr_egger_regression`:
    SNPs oriented to a positive exposure effect, weighted least squares of BETA_o on BETA_e with an intercept.
    BETA_e and BETA_o can hold one dataset per row (last axis: SNPs), with one intercept and slope per row.
    """
    sign0 = np.where(BETA_e < 0, -1.0, 1.0)  # Sign of the exposure effect (0 counted as positive)
    x = np.abs(BETA_e)
    y = BETA_o * sign0
    x_mean = np.sum(weights * x, axis=-1, keepdims=True) / np.sum(weights)
    y_mean = np.sum(weights * y, axis=-1, keepdims=True) / np.sum(weights)
    slope = np.sum(weights * (x - x_mean) * (y - y_mean), axis=-1, keepdims=True) / np.sum(
        weights * (x - x_mean) ** 2, axis=-1, keepdims=True
    )
    intercept = y_mean - slope * x_mean
    if np.ndim(BETA_e) == 1:
        return intercept[0], slope[0]
    return intercept[:, 0], slope[:, 0]


def bootstrap_chunks(nboot, n_snps):
    """(start, stop) ranges of bootstrap replicates processed together (at most BOOTSTRAP_CHUNK_VALUES values per array)."""
    rows = max(1, BOOTSTRAP_CHUNK_VALUES // max(n_snps, 1))
    return [(start, min(start + rows, nboot)) for start in range(0, nboot, rows)]


def draw_bootstrap_effects(rng, BETA_e, SE_e, BETA_o, SE_o, n_rows):
    """
    Bootstrap replicates of the SNP effects, one per row: exposure effects ~ N(BETA_e, SE_e) and outcome effects
    ~ N(BETA_o, SE_o). Each replicate uses the standard normal draws of rng.normal(BETA_e, SE_e) followed by
    rng.normal(BETA_o, SE_o), so drawing the replicates in chunks gives the same values as drawing them one by one.
    """
    if (SE_e < 0).any() or (SE_o < 0).any():
        raise ValueError("Standard errors (SE_e and SE_o) must be non-negative.")
    z = rng.standard_normal(size=(n_rows, 2, len(BETA_e)))
    return BETA_e + SE_e * z[:, 0], BETA_o + SE_o * z[:, 1]

"""
Median methods
"""

def weighted_median(b_iv, weights):
    """Helper function to compute the weighted median estimate (NaN if the weights do not have a finite positive sum)."""
    b_iv = np.asarray(b_iv, dtype=float)
    if b_iv.size == 0:
        return np.nan
    return weighted_median_rows(b_iv[np.newaxis], np.asarray(weights, dtype=float))[0]


def weighted_median_rows(b_iv, weights):
    """
    Weighted medians of the rows of b_iv (one set of ratio estimates per row) with the same weights for every row:
    linear interpolation of the sorted estimates at half of the cumulative weight (each estimate placed at the middle
    of its weight), NaN if the weights do not have a finite positive sum.
    """
    n = b_iv.shape[1]
    if n > 1 and np.all(weights == weights[0]):
        # Equal weights: the cumulative weights do not depend on the order, only two order statistics are needed
        total = np.sum(weights)
        with np.errstate(divide="ignore", invalid="ignore"):
            weights_sum = (np.cumsum(weights) - 0.5 * weights) / total
        below = np.flatnonzero(weights_sum < 0.5)
        if not np.isfinite(total) or total <= 0:
            return np.full(len(b_iv), np.nan)
        if len(below) == 0:
            return np.partition(b_iv, 0, axis=1)[:, 0]
        if below[-1] == n - 1:
            return np.partition(b_iv, n - 1, axis=1)[:, n - 1]
        k = below[-1]
        part = np.partition(b_iv, [k, k + 1], axis=1)
        return part[:, k] + (part[:, k + 1] - part[:, k]) * (0.5 - weights_sum[k]) / (weights_sum[k + 1] - weights_sum[k])

    order = np.argsort(b_iv, axis=1)
    b_sorted = np.take_along_axis(b_iv, order, axis=1)
    w_sorted = weights[order]
    total = np.sum(w_sorted, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        weights_sum = (np.cumsum(w_sorted, axis=1) - 0.5 * w_sorted) / total[:, np.newaxis]
        below_half = weights_sum < 0.5
        below = n - 1 - np.argmax(below_half[:, ::-1], axis=1)  # Last estimate below half of the weight
        left = np.minimum(below, n - 2)[:, np.newaxis]
        b_left, b_right = np.take_along_axis(b_sorted, left, axis=1)[:, 0], np.take_along_axis(b_sorted, left + 1, axis=1)[:, 0]
        s_left, s_right = np.take_along_axis(weights_sum, left, axis=1)[:, 0], np.take_along_axis(weights_sum, left + 1, axis=1)[:, 0]
        b = b_left + (b_right - b_left) * (0.5 - s_left) / (s_right - s_left)
    b = np.where(below == n - 1, b_sorted[:, -1], b)  # All the weight on the largest estimate
    b = np.where(below_half.any(axis=1), b, b_sorted[:, 0])  # All the weight on the smallest estimate
    return np.where(np.isfinite(total) & (total > 0), b, np.nan)


def weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot, seed=None):
    """Helper function returning the bootstrap standard error (standard deviation of nboot bootstrap replicates) of the weighted median."""
    BETA_e, SE_e, BETA_o, SE_o, weights = (np.asarray(x, dtype=float) for x in (BETA_e, SE_e, BETA_o, SE_o, weights))
    rng = np.random.default_rng(seed)
    med = np.empty(nboot)
    for start, stop in bootstrap_chunks(nboot, len(BETA_e)):
        BETA_e_boot, BETA_o_boot = draw_bootstrap_effects(rng, BETA_e, SE_e, BETA_o, SE_o, stop - start)
        med[start:stop] = weighted_median_rows(BETA_o_boot / BETA_e_boot, weights)
    return np.std(med, ddof=1)


def mr_weighted_median(BETA_e, SE_e, BETA_o, SE_o, nboot, seed=None):
    """
    Perform a Mendelian Randomization analysis using the weighted median method.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value
        seed (int, optional): Random seed for the bootstrap draws. Default None (non-reproducible).

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient representing the causal estimate.
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.

    Notes:
        The standard error is obtained with bootstrapping.
    """
    l = len(BETA_e)
    if l < 3:
        return [_null_result(MR_METHODS_NAMES["WM"])]

    b_iv = BETA_o / BETA_e
    VBj = (SE_o**2) / (BETA_e**2) + ((BETA_o**2) * (SE_e**2)) / (BETA_e**4)

    b = weighted_median(b_iv, 1 / VBj)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, 1 / VBj, nboot, seed)
    pval = 2 * norm.sf(abs(b / se))
    return [{"method": MR_METHODS_NAMES["WM"], "nSNP": l, "b": b, "se": se, "pval": pval}]


def mr_pen_wm(BETA_e, SE_e, BETA_o, SE_o, nboot, penk, seed=None):
    """
    Perform a Mendelian Randomization analysis using the penalised weighted median method. See https://arxiv.org/abs/1606.03729.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value.
        penk (float): Constant factor used to penalise the weights.
        seed (int, optional): Random seed for the bootstrap draws. Default None (non-reproducible).

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient representing the causal estimate.
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.
    """
    l = len(BETA_e)
    if l < 3:
        return [_null_result(MR_METHODS_NAMES["WM-pen"])]

    betaIV = BETA_o / BETA_e
    #betaIVW = np.sum(BETA_o * BETA_e / SE_o**2) / np.sum(BETA_e**2 / SE_o**2)
    VBj = (SE_o**2) / (BETA_e**2) + (BETA_o**2) * (SE_e**2) / (BETA_e**4)
    weights = 1 / VBj

    bwm = weighted_median(betaIV, weights)

    # Penalised weights w * min(1, penk * penalty), computed in log space (the penalties of distant SNPs can underflow)
    # and rescaled by their maximum: the weighted median only depends on the relative weights.
    # log of the chi-squared (1 df) survival function: log(2) + log(Phi(-sqrt(x))), finite for large x
    with np.errstate(divide="ignore", invalid="ignore"):
        log_penalty = np.log(2) + log_ndtr(-np.sqrt(np.asarray(weights * (betaIV - bwm) ** 2, dtype=float)))
        log_pen_weights = np.asarray(np.log(weights) + np.minimum(0, np.log(penk) + log_penalty), dtype=float)
    log_pen_weights = np.where(np.isfinite(log_pen_weights), log_pen_weights, -np.inf)
    if not np.isfinite(log_pen_weights).any():
        print("Penalised weighted median: no valid penalised weight (the SNP weights are zero or not finite). Returning NaN.")
        return [_null_result(MR_METHODS_NAMES["WM-pen"])]
    pen_weights = np.exp(log_pen_weights - log_pen_weights.max())

    b = weighted_median(betaIV, pen_weights)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, pen_weights, nboot, seed)
    pval = 2 * norm.sf(abs(b / se))

    return [
        {
            "method": MR_METHODS_NAMES["WM-pen"],
            "b": b,
            "se": se,
            "pval": pval,
            "nSNP": l,
        }
    ]


def mr_simple_median(BETA_e, SE_e, BETA_o, SE_o, nboot, seed=None):
    """
    Perform a Mendelian Randomization analysis using the simple median method.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value
        seed (int, optional): Random seed for the bootstrap draws. Default None (non-reproducible).

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient representing the causal estimate.
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.

    Notes:
        The standard error is obtained with bootstrapping.
    """
    l = len(BETA_e)
    if l < 3:
        return [_null_result(MR_METHODS_NAMES["Simple-median"])]

    b_iv = BETA_o / BETA_e
    weights = np.repeat(1 / len(BETA_e), len(BETA_e))
    b = weighted_median(b_iv, weights)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot, seed)
    pval = 2 * norm.sf(abs(b / se))
    return [{"method": MR_METHODS_NAMES["Simple-median"], "b": b, "se": se, "pval": pval, "nSNP": l}]


"""
Regression methods
"""

def _ivw_fit(BETA_e, BETA_o, SE_o):
    """
    Helper function fitting the weighted regression of BETA_o on BETA_e through the origin (weights 1/SE_o^2).
    Returns the slope, its fixed-effect standard error 1/sqrt(sum(w * BETA_e^2)) and the residual standard error sigma.
    """
    BETA_e, BETA_o, SE_o = (np.asarray(x, dtype=float) for x in (BETA_e, BETA_o, SE_o))
    w = 1 / SE_o**2
    sxx = np.sum(w * BETA_e**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        b = np.sum(w * BETA_e * BETA_o) / sxx
        sigma = np.sqrt(np.sum(w * (BETA_o - b * BETA_e) ** 2) / (len(BETA_e) - 1))
        se_fe = 1 / np.sqrt(sxx)
    return b, se_fe, sigma


def mr_ivw(BETA_e, SE_e, BETA_o, SE_o):
    """
    Perform a Mendelian Randomization analysis using the Inverse Variance Weighted (IVW) method with random effects. Standard Error is corrected for under dispersion (as opposed to :func:`mr_ivw_re`).

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate.
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.
            - "Q": Cochran"s Q statistic for heterogeneity.
            - "Q_df": Degrees of freedom for the Q statistic.
            - "Q_pval": P-value for the Q statistic.

    Notes:
        The function uses weighted least squares regression (WLS) to estimate the causal effect size,
        weighting by the inverse of the variance of the outcome"s effect sizes.
        Cochran"s Q statistics also computed to assess the heterogeneity across the instrumental variables.
    """
    # If less than 2 valid rows, return NA values
    l = len(BETA_e)
    if l < 2:
        return [_null_result(MR_METHODS_NAMES["IVW"])]

    # Weighted regression through the origin; SE = fixed-effect SE * max(1, sigma) (under-dispersion correction)
    b, se_fe, sigma = _ivw_fit(BETA_e, BETA_o, SE_o)
    se = se_fe * max(1, sigma)
    pval = 2 * norm.sf(abs(b / se))

    Q_df = l - 1
    Q = sigma**2 * Q_df
    Q_pval = chi2.sf(Q, Q_df)

    return [
        {
            "method": MR_METHODS_NAMES["IVW"],
            "nSNP": l,
            "b": b,
            "se": se,
            "pval": pval,
            "Q_df": Q_df,
            "Q": Q,
            "Q_pval": Q_pval,
        }
    ]


def mr_ivw_re(BETA_e, SE_e, BETA_o, SE_o):
    """
    Perform a Mendelian Randomization analysis using the Inverse Variance Weighted (IVW) method with random effects. Standard Error is not corrected for under dispersion (as opposed to :func:`mr_ivw`).

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate.
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.
            - "Q": Cochran"s Q statistic for heterogeneity.
            - "Q_df": Degrees of freedom for the Q statistic.
            - "Q_pval": P-value for the Q statistic.

    Notes:
        The function uses weighted least squares regression (WLS) to estimate the causal effect size,
        weighting by the inverse of the variance of the outcome"s effect sizes.
        Cochran"s Q statistics also computed to assess the heterogeneity across the instrumental variables.
    """
    # If less than 2 valid rows, return NA values
    l = len(BETA_e)
    if l < 2:
        return [_null_result(MR_METHODS_NAMES["IVW-RE"])]

    # Weighted regression through the origin; SE = fixed-effect SE * sigma (multiplicative random effects)
    b, se_fe, sigma = _ivw_fit(BETA_e, BETA_o, SE_o)
    se = se_fe * sigma
    with np.errstate(divide="ignore", invalid="ignore"):
        pval = 2 * norm.sf(abs(b / se))
    Q_df = l - 1
    Q = sigma**2 * Q_df
    Q_pval = chi2.sf(Q, Q_df)

    return [
        {
            "method": MR_METHODS_NAMES["IVW-RE"],
            "nSNP": l,
            "b": b,
            "se": se,
            "pval": pval,
            "Q_df": Q_df,
            "Q": Q,
            "Q_pval": Q_pval,
        }
    ]


def mr_ivw_fe(BETA_e, SE_e, BETA_o, SE_o):
    """
    Perform a Mendelian Randomization analysis using the Inverse Variance Weighted (IVW) method with fixed effects.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate.
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.
            - "Q": Cochran"s Q statistic for heterogeneity.
            - "Q_df": Degrees of freedom for the Q statistic.
            - "Q_pval": P-value for the Q statistic.

    Notes:
        The function uses weighted least squares regression (WLS) to estimate the causal effect size,
        weighting by the inverse of the variance of the outcome"s effect sizes.
        Cochran"s Q statistics also computed to assess the heterogeneity across the instrumental variables.
    """
    l = len(BETA_e)
    if l < 2:
        return [_null_result(MR_METHODS_NAMES["IVW-FE"])]

    # Weighted regression through the origin; fixed-effect SE 1/sqrt(sum(w * BETA_e^2))
    b, se, sigma = _ivw_fit(BETA_e, BETA_o, SE_o)
    pval = 2 * norm.sf(abs(b / se))
    Q_df = l - 1
    Q = sigma**2 * Q_df
    Q_pval = chi2.sf(Q, Q_df)
    return [
        {
            "method": MR_METHODS_NAMES["IVW-FE"],
            "nSNP": l,
            "b": b,
            "se": se,
            "pval": pval,
            "Q_df": Q_df,
            "Q": Q,
            "Q_pval": Q_pval,
        }
    ]


def mr_uwr(BETA_e, SE_e, BETA_o, SE_o):
    """
    Performs an unweighted regression Mendelian Randomization analysis.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate.
            - "se": Standard error of the coefficient (residual-based OLS standard error).
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.
            - "Q", "Q_df", "Q_pval": NaN (the residual sum of squares of the unweighted regression is not a heterogeneity test).

    Notes:
        The returned causal estimate is not weighted by the inverse variance.
        The standard error is the ordinary least squares standard error, based on the residual variance of the unweighted regression.
    """

    l = len(BETA_e)
    if l < 2:
        return [_null_result(MR_METHODS_NAMES["UWR"])]

    # Perform regression without weights
    model = sm.OLS(BETA_o, BETA_e).fit()

    # Extract coefficients and compute statistics
    b = model.params.iloc[0]
    se = model.bse.iloc[0]
    pval = 2 * norm.sf(np.abs(b / se))

    return [
        {
            "method": MR_METHODS_NAMES["UWR"],
            "b": b,
            "se": se,
            "pval": pval,
            "nSNP": l,
            "Q": np.nan,
            "Q_df": np.nan,
            "Q_pval": np.nan,
        }
    ]


def mr_wald_ratio(BETA_e, SE_e, BETA_o, SE_o):
    """
    Compute the Wald ratio causal estimate for a single genetic variant.

    Args:
        BETA_e (numpy array): Effect size of the genetic variant on the exposure (length 1).
        SE_e (numpy array): Standard error corresponding to `BETA_e`.
        BETA_o (numpy array): Effect size of the same genetic variant on the outcome (length 1).
        SE_o (numpy array): Standard error corresponding to `BETA_o`.

    Returns:
        list of dict: A list containing a dictionary with the results:
            - "method": Name of the analysis method.
            - "b": Ratio BETA_o / BETA_e, representing the causal estimate.
            - "se": First-order standard error SE_o / |BETA_e|.
            - "pval": P-value for the causal estimate (normal distribution).
            - "nSNP": Number of genetic variants used in the analysis (1).
    """
    if len(BETA_e) != 1:
        return [_null_result(MR_METHODS_NAMES["Wald"])]

    BETA_e, SE_o, BETA_o = (float(np.asarray(x, dtype=float)[0]) for x in (BETA_e, SE_o, BETA_o))
    b = BETA_o / BETA_e
    se = SE_o / abs(BETA_e)
    pval = 2 * norm.sf(abs(b / se))
    return [{"method": MR_METHODS_NAMES["Wald"], "b": b, "se": se, "pval": pval, "nSNP": 1}]


""" 
Sign method
"""

def mr_sign(BETA_e, BETA_o):
    """
    Performs the sign concordance test.

    The sign concordance test is used to determine whether there is a consistent direction
    of effect (i.e., sign) between the exposure and outcome across the variants.
    The consistent directonality is an assumption of Mendelian Randomization.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.

    Returns:
        list of dict: A list containing dictionaries with the following keys:
            - "method": Name of the method.
            - "nSNP": Number of genetic variants used in the analysis.
            - "b": Proportion of concordant signs between exposure and outcome effect
                   sizes minus 0.5, multiplied by 2.
            - "se": Not applicable for this method (returns NaN).
            - "pval": P-value for the sign concordance test based on a binomial distribution.

    Notes:
        Effect sizes that are exactly zero are replaced with NaN and are not included in the analysis. A binomial test is then performed to evaluate the probability of observing the given number of concordant signs by chance alone, assuming a null expectation of 50% concordance.
    """

    # Replace zeros with NaNs
    BETA_e = np.where(BETA_e == 0, np.nan, BETA_e)
    BETA_o = np.where(BETA_o == 0, np.nan, BETA_o)

    # Check for enough non-missing values
    valid_data = (~np.isnan(BETA_e)) & (~np.isnan(BETA_o))
    if np.sum(valid_data) < 6:
        return [_null_result(MR_METHODS_NAMES["Sign"])]

    # Count the number of consistent signs
    x = np.sum(np.sign(BETA_e[valid_data]) == np.sign(BETA_o[valid_data]))
    n = np.sum(valid_data)

    # Binomial test
    pval = binomtest(x, n, p=0.5).pvalue
    b = (x / n - 0.5) * 2

    return [
        {
            "method": MR_METHODS_NAMES["Sign"],
            "nSNP": n,
            "b": b,
            "se": np.nan,
            "pval": pval,
        }
    ]
