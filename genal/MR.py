import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.stats import norm, chi2, binomtest, t
from scipy.special import log_ndtr
from scipy.optimize import minimize_scalar
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from functools import partial

from .constants import MR_METHODS_NAMES


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
    cpus (int): Number of cpu cores to use in parallel for the boostrapping iterations.
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
    cpus (int): Number of cpu cores to use in parallel for the boostrapping iterations.
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

def _kde(x, b_sorted, w_sorted, h):
    """
    Unnormalised (weighted) Gaussian kernel density with bandwidth h of the sorted values b_sorted (weights w_sorted),
    evaluated at the sorted points x. Kernels more than 10 bandwidths away (relative contribution < 2e-22) are skipped.
    """
    out = np.empty(len(x))
    for start in range(0, len(x), 512):
        xs = x[start:start + 512]
        lo = np.searchsorted(b_sorted, xs[0] - 10 * h, side="left")
        hi = np.searchsorted(b_sorted, xs[-1] + 10 * h, side="right")
        out[start:start + 512] = np.exp(-0.5 * ((xs[:, None] - b_sorted[None, lo:hi]) / h) ** 2) @ w_sorted[lo:hi]
    return out


def beta_mode(BETA_IV, SE_IV, phi, method):
    """
    Function to compute the point estimate of mode methods: the global maximum of the Gaussian kernel density of the
    ratio estimates (unweighted for "Simple", weighted by SE_IV^-2 for "Weighted").
    The density is evaluated on a grid of spacing h/10 covering the union of the intervals [b - 4h, b + 4h] around the
    ratio estimates b; every grid local maximum within 0.25% of the highest one (the largest possible grid error) is
    refined by a bounded scalar maximisation, and the highest refined maximum is returned.
    """
    BETA_IV = np.asarray(BETA_IV, dtype=float)

    # Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
    s = 0.9 * min(np.std(BETA_IV, ddof=1), stats.median_abs_deviation(BETA_IV, scale=1/1.4826)) / len(BETA_IV) ** (1/5)
        
    # Define the actual bandwidth
    h = max(0.00000001, s * phi)

    if method == "Weighted": #Using weighted KDE
        weights = np.asarray(SE_IV, dtype=float) ** -2
    elif method == "Simple": #Using unweighted KDE
        weights = np.ones(len(BETA_IV))
    else:
        raise ValueError("Method should be either 'Simple' or 'Weighted'.")
    weights = weights / np.sum(weights) # Standardised weights
    order = np.argsort(BETA_IV)
    b, w = BETA_IV[order], weights[order]

    # Grid of spacing h/10 on the union of the intervals [b - 4h, b + 4h] (one segment per group of close estimates)
    step = h / 10
    breaks = np.flatnonzero(np.diff(b) > 8 * h) + 1
    segments = [
        np.linspace(seg[0] - 4 * h, seg[-1] + 4 * h, int(np.ceil((seg[-1] - seg[0] + 8 * h) / step)) + 1)
        for seg in np.split(b, breaks)
    ]
    grid = np.concatenate(segments)
    density = _kde(grid, b, w, h)

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
    for i in candidates:
        opt = minimize_scalar(
            lambda x: -_kde(np.array([x]), b, w, h)[0],
            bounds=(grid[i] - step, grid[i] + step),
            method="bounded",
            options={"xatol": 1e-9 * h},
        )
        if -opt.fun > best_density:
            mode, best_density = opt.x, -opt.fun

    return np.array([mode])

def bootstrap_mode_iteration(BETA_IV_boot, SE_IV, phi, method):
    """Helper function computing the mode estimate of one bootstrap replicate (pre-drawn ratio estimates)."""
    if method=="Simple": #If Simple mode: the kernel density estimation is not weighted (SE vector of 1)
        return beta_mode(BETA_IV_boot, np.ones(len(BETA_IV_boot)), phi, "Simple")
    elif method=="Weighted": #If Weighted mode: the kernel density estimation is weighted using the standard errors
        return beta_mode(BETA_IV_boot, SE_IV, phi, "Weighted")
    else:
        raise ValueError("Method should be either 'Simple' or 'Weighted'.")

def bootstrap_mode(BETA_IV, SE_IV, phi, nboot, cpus, method, show_progress=True, seed=None):
    """
    Function to obtain arrays of bootstrapped beta estimates (weighted or unweighted) to estimate SEs of mode methods.
    All bootstrap draws are generated here from a single seeded generator; the worker processes only compute
    the mode of each pre-drawn replicate, so the results do not depend on the number of cpus.
    """
    SE_IV = np.asarray(SE_IV, dtype=float)

    # Draw all bootstrap replicates of the ratio estimates (nboot x nSNP)
    rng = np.random.default_rng(seed)
    BETA_IV_boot = rng.normal(np.asarray(BETA_IV, dtype=float), SE_IV, size=(nboot, len(SE_IV)))

    #Wrapper function freezing the bootstrap_iteration for parallel execution
    partial_bootstrap_iteration = partial(bootstrap_mode_iteration, SE_IV=SE_IV, phi=phi, method=method)
    progress = dict(total=nboot, desc=f"{method} mode bootstrapping", ncols=100, disable=not show_progress)

    if cpus == 1:
        results = [partial_bootstrap_iteration(row) for row in tqdm(BETA_IV_boot, **progress)]
    else:
        #Parallel execution with progress bar
        with ProcessPoolExecutor(max_workers=cpus) as executor:
            results = list(
                tqdm(
                    executor.map(partial_bootstrap_iteration, BETA_IV_boot, chunksize=max(1, nboot // (10 * cpus))),
                    **progress,
                )
            )

    return np.array(results)


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
    res = np.zeros((nboot, 2))
    for i in range(nboot):
        xs = rng.normal(BETA_e, SE_e)
        ys = rng.normal(BETA_o, SE_o)
        res[i] = egger_wls(xs, ys, weights)

    se, se_i = np.std(res[:, 1], ddof=1), np.std(res[:, 0], ddof=1)

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
    """
    sign0 = np.where(BETA_e < 0, -1.0, 1.0)  # Sign of the exposure effect (0 counted as positive)
    x = np.abs(BETA_e)
    y = BETA_o * sign0
    x_mean = np.sum(weights * x) / np.sum(weights)
    y_mean = np.sum(weights * y) / np.sum(weights)
    slope = np.sum(weights * (x - x_mean) * (y - y_mean)) / np.sum(weights * (x - x_mean) ** 2)
    return y_mean - slope * x_mean, slope

"""
Median methods
"""

def weighted_median(b_iv, weights):
    """Helper function to compute the weighted median estimate (NaN if the weights do not have a finite positive sum)."""
    order_indices = np.argsort(b_iv)
    betaIV_order = np.asarray(b_iv, dtype=float)[order_indices]
    weights_order = np.asarray(weights, dtype=float)[order_indices]
    total = np.sum(weights_order)
    if len(weights_order) == 0 or not np.isfinite(total) or total <= 0:
        return np.nan
    weights_sum = (np.cumsum(weights_order) - 0.5 * weights_order) / total
    below_indices = np.flatnonzero(weights_sum < 0.5)
    if len(below_indices) == 0:  # All the weight on the smallest estimate
        return betaIV_order[0]
    below = below_indices[-1]
    if below == len(betaIV_order) - 1:  # All the weight on the largest estimate
        return betaIV_order[-1]
    b = betaIV_order[below] + (betaIV_order[below + 1] - betaIV_order[below]) * (
        0.5 - weights_sum[below]
    ) / (weights_sum[below + 1] - weights_sum[below])
    return b


def weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot, seed=None):
    """Helper function returning the bootstrap standard error (standard deviation of nboot bootstrap replicates) of the weighted median."""
    BETA_e, SE_e, BETA_o, SE_o = (np.asarray(x, dtype=float) for x in (BETA_e, SE_e, BETA_o, SE_o))
    rng = np.random.default_rng(seed)
    med = np.zeros(nboot)
    for i in range(nboot):
        BETA_e_boot = rng.normal(loc=BETA_e, scale=SE_e)
        BETA_o_boot = rng.normal(loc=BETA_o, scale=SE_o)
        betaIV_boot = BETA_o_boot / BETA_e_boot
        med[i] = weighted_median(betaIV_boot, weights)
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
