import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.stats import norm, chi2, binomtest, t
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.neighbors import KernelDensity
from tqdm import tqdm
from functools import partial

from .constants import MR_METHODS_NAMES

"""
Mode methods
"""

def mr_simple_mode(BETA_e, SE_e, BETA_o, SE_o, phi, nboot, cpus):
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

    Returns:
        list of dict: A list containing two dictionaries with the results for the egger regression estimate and the egger regression intercept (horizontal pleiotropy estimate):
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate or the intercept.
            - "se": Adjusted standard error of the coefficient or intercept.
            - "pval": P-value for the causal estimate or intercept.
            - "nSNP": Number of genetic variants used in the analysis.
    """
    l = len(BETA_e)
    if l < 3:
        return [
            {
                "method": MR_METHODS_NAMES["Simple-mode"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]
    
    BETA_IV   = BETA_o/BETA_e
    SE_IV = np.sqrt((SE_o**2) / (BETA_e**2) + ((BETA_o**2) * (SE_e**2)) / (BETA_e**4))
    
    BETA_mode = beta_mode(BETA_IV, np.ones(l), phi, method="Simple") #Point estimate calculated using unweighted KDE
    BETA_boot = bootstrap_mode (BETA_IV, SE_IV, phi, nboot, cpus, method="Simple") #Generating bootstrapped point estimates with normal sampling and unweighted KDE
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

def mr_weighted_mode(BETA_e, SE_e, BETA_o, SE_o, phi, nboot, cpus):
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

    Returns:
        list of dict: A list containing two dictionaries with the results for the egger regression estimate and the egger regression intercept (horizontal pleiotropy estimate):
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate or the intercept.
            - "se": Adjusted standard error of the coefficient or intercept.
            - "pval": P-value for the causal estimate or intercept.
            - "nSNP": Number of genetic variants used in the analysis.
    """
    l = len(BETA_e)
    if l < 3:
        return [
            {
                "method": MR_METHODS_NAMES["Weighted-mode"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]
    
    BETA_IV   = BETA_o/BETA_e
    SE_IV = np.sqrt((SE_o**2) / (BETA_e**2) + ((BETA_o**2) * (SE_e**2)) / (BETA_e**4))
    
    BETA_mode = beta_mode(BETA_IV, SE_IV, phi, method="Weighted") #Point estimate calculated using KDE weighted with standard errors
    BETA_boot = bootstrap_mode (BETA_IV, SE_IV, phi, nboot, cpus, method="Weighted") #Generating bootstrapped point estimates with normal sampling and weighted KDE
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

def beta_mode(BETA_IV, SE_IV, phi, method):
    """Function to compute the point estimate of mode methods."""
    # Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
    s = 0.9 * min(np.std(BETA_IV, ddof=1), stats.median_abs_deviation(BETA_IV, scale=1/1.4826)) / len(BETA_IV) ** (1/5)
        
    # Define the actual bandwidth
    h = max(0.00000001, s * phi)

    # Kernel density estimation (gaussian kernel)
    kde = KernelDensity(bandwidth=h, kernel='gaussian')

    if method == "Weighted": #Using weighted KDE 
        weights = SE_IV ** -2 / np.sum(SE_IV ** -2) # Standardised weights
        kde.fit(np.array(BETA_IV).reshape(-1, 1), sample_weight=weights)
    elif method == "Simple": #Using unweighted KDE
        kde.fit(np.array(BETA_IV).reshape(-1, 1))

    # Create points to evaluate the density on (Using 512 points and 3*h as padding to be equivalent to the stas::density function in R)
    x_range = np.linspace(BETA_IV.min()-3*h, BETA_IV.max()+3*h, 512)[:, np.newaxis]

    # Calculate log density and then convert to actual density
    log_density = kde.score_samples(x_range)
    density = np.exp(log_density)

    # Find mode
    mode = x_range[np.argmax(density)]

    return mode

def bootstrap_mode_iteration(i, BETA_IV, SE_IV, phi, method):
    """Helper function to run the simple and weighted mode boostrapping in parallel."""
    BETA_IV_boot = np.random.normal(BETA_IV, SE_IV) #Normally sample BETA_IV
    
    if method=="Simple": #If Simple mode: the kernel density estimation is not weighted (SE vector of 1)
        return beta_mode(BETA_IV_boot, np.ones(len(BETA_IV)), phi, "Simple")
    elif method=="Weighted": #If Weighted mode: the kernel density estimation is weighted using the standard errors
        return beta_mode(BETA_IV_boot, SE_IV, phi, "Weighted")
    else:
        raise ValueError("Method should be either 'Simple' or 'Weighted'.")

def bootstrap_mode (BETA_IV, SE_IV, phi, nboot, cpus, method):
    """Function to obtain arrays of bootstrapped beta estimates (weighted or unweighted) to estimate SEs of mode methods."""
    
    # Declare the array to store results of the bootstrapping
    BETA_boot = np.ones((nboot, 1))
    
    #Wrapper function freezing the bootstrap_iteration for parallel execution
    partial_bootstrap_iteration = partial(bootstrap_mode_iteration, BETA_IV=BETA_IV, SE_IV=SE_IV, phi=phi, method=method) 
    
    #Parallel execution with progress bar
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        results = list(
            tqdm(
                executor.map(partial_bootstrap_iteration, range(nboot)),
                total=nboot,
                desc=f"{method} mode bootstrapping",
                ncols=100,
            )
        )
    BETA_boot = np.array(results)
        
    return BETA_boot


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
    # Initialize null result
    null_result = [
        {
            "method": MR_METHODS_NAMES["Egger"][0],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        },
        {
            "method": MR_METHODS_NAMES["Egger"][1],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        },
    ]

    l = len(BETA_e)
    # Early return if insufficient data
    if l < 3:
        return null_result

    sign0 = np.sign(BETA_e.replace(0, 1))

    BETA_o = BETA_o.copy()
    BETA_e = BETA_e.copy()
    BETA_o *= sign0
    BETA_e = np.abs(BETA_e)
    flipped = sign0 == -1

    X = BETA_e
    X = sm.add_constant(X)
    y = BETA_o
    weights = 1 / SE_o**2

    mod = sm.WLS(y, X, weights=weights).fit()

    if len(mod.params) > 1:
        b = mod.params.iloc[1]
        se = mod.bse.iloc[1] / min(1, np.sqrt(mod.mse_resid))
        pval = 2 * t.sf(abs(b / se), l - 2)

        b_i = mod.params.iloc[0]
        se_i = mod.bse.iloc[0] / min(1, np.sqrt(mod.mse_resid))
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

def mr_egger_regression_bootstrap(BETA_e, SE_e, BETA_o, SE_o, nboot, cpus=4):
    """
    Perform a Mendelian Randomization analysis using the egger regression method with boostrapped standard errors. See :func:`mr_egger_regression` for a version without bootstrapping.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value
        cpus (int): Number of cpu cores to use in parallel for the boostrapping iterations.

    Returns:
        list of dict: A list containing two dictionaries with the results for the egger regression estimate and the egger regression intercept (horizontal pleiotropy estimate):
            - "method": Name of the analysis method.
            - "b": Coefficient of the regression, representing the causal estimate or the intercept.
            - "se": Adjusted standard error of the coefficient or intercept.
            - "pval": P-value for the causal estimate or intercept.
            - "nSNP": Number of genetic variants used in the analysis.
    """

    l = len(BETA_e)
    if l < 3:
        return [
            {
                "method": MR_METHODS_NAMES["Egger-boot"][0],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            },
            {
                "method": MR_METHODS_NAMES["Egger-boot"][1],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            },
        ]

    res = np.zeros((nboot + 1, 2))

    with ProcessPoolExecutor(max_workers=cpus) as executor:
        # Start the load operations and mark each future with its URL
        futures = {
            executor.submit(parallel_bootstrap_func, i, BETA_e, SE_e, BETA_o, SE_o): i
            for i in range(nboot)
        }
        for future in tqdm(
            as_completed(futures), total=nboot, desc="MR Egger bootstrapping", ncols=100
        ):
            ahat, bhat = future.result()
            i = futures[future]  # get the original index/counter
            res[i, 0] = ahat
            res[i, 1] = bhat

    return [
        {
            "method": MR_METHODS_NAMES["Egger-boot"][0],
            "b": np.nanmean(res[:, 1]),
            "se": np.nanstd(res[:, 1]),
            "pval": np.sum(np.sign(np.nanmean(res[:, 1])) * res[:, 1] < 0) / nboot,
            "nSNP": l,
        },
        {
            "method": MR_METHODS_NAMES["Egger-boot"][1],
            "b": np.nanmean(res[:, 0]),
            "se": np.nanstd(res[:, 0]),
            "pval": np.sum(np.sign(np.nanmean(res[:, 0])) * res[:, 0] < 0) / nboot,
            "nSNP": l,
        },
    ]

def parallel_bootstrap_func(i, BETA_e, SE_e, BETA_o, SE_o):
    """Helper function to run the egger regression bootstrapping in parallel."""
    xs = np.random.normal(BETA_e, SE_e)
    ys = np.random.normal(BETA_o, SE_o)

    # Use absolute values for Egger reg
    ys *= np.sign(xs)
    xs = np.abs(xs)

    r = linreg(xs, ys, 1 / SE_o**2)
    return r["ahat"], r["bhat"]


def linreg(x, y, w=None):
    """Helper function to run linear regressions for the parallel egger bootstrapping."""
    if w is None:
        w = np.ones_like(x)

    xp = w * x
    yp = w * y

    bhat = np.cov(xp, yp)[0, 1] / np.var(xp)
    ahat = np.nanmean(y) - np.nanmean(x) * bhat
    yhat = ahat + bhat * x

    residuals = y - yhat
    se = np.sqrt(
        sum(w * residuals**2) / (np.sum(~np.isnan(yhat)) - 2) / np.sum(w * x**2)
    )
    pval = 2 * norm.sf(abs(bhat / se))

    return {"ahat": ahat, "bhat": bhat, "se": se, "pval": pval}

"""
Median methods
"""

def weighted_median(b_iv, weights):
    """Helper function to compute the weighted median estimate."""
    order_indices = np.argsort(b_iv)
    betaIV_order = np.array(b_iv)[order_indices]
    weights_order = np.array(weights)[order_indices]
    weights_sum = np.cumsum(weights_order) - 0.5 * weights_order
    weights_sum /= np.sum(weights_order)
    below = np.max(np.where(weights_sum < 0.5))
    b = betaIV_order[below] + (betaIV_order[below + 1] - betaIV_order[below]) * (
        0.5 - weights_sum[below]
    ) / (weights_sum[below + 1] - weights_sum[below])
    return b


def weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot):
    """Helper function to generate boostrapped replications."""
    med = np.zeros(nboot)
    #for i in tqdm(
    #    range(nboot), total=nboot, desc="Weighted median bootstrapping", ncols=100
    #):
    for i in range(nboot):
        BETA_e_boot = np.random.normal(loc=BETA_e, scale=SE_e)
        BETA_o_boot = np.random.normal(loc=BETA_o, scale=SE_o)
        betaIV_boot = BETA_o_boot / BETA_e_boot
        med[i] = weighted_median(betaIV_boot, weights)
    return np.std(med)


def mr_weighted_median(BETA_e, SE_e, BETA_o, SE_o, nboot):
    """
    Perform a Mendelian Randomization analysis using the weighted median method.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value

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
        return [
            {
                "method": MR_METHODS_NAMES["WM"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]

    b_iv = BETA_o / BETA_e
    VBj = (SE_o**2) / (BETA_e**2) + ((BETA_o**2) * (SE_e**2)) / (BETA_e**4)

    b = weighted_median(b_iv, 1 / VBj)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, 1 / VBj, nboot)
    pval = 2 * norm.sf(abs(b / se))
    return [{"method": MR_METHODS_NAMES["WM"], "nSNP": l, "b": b, "se": se, "pval": pval}]


def mr_pen_wm(BETA_e, SE_e, BETA_o, SE_o, nboot, penk):
    """
    Perform a Mendelian Randomization analysis using the penalised weighted median method. See https://arxiv.org/abs/1606.03729.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value.
        penk (float): Constant factor used to penalise the weights.

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
        return [
            {
                "method": MR_METHODS_NAMES["WM-pen"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]

    betaIV = BETA_o / BETA_e
    #betaIVW = np.sum(BETA_o * BETA_e / SE_o**2) / np.sum(BETA_e**2 / SE_o**2)
    VBj = (SE_o**2) / (BETA_e**2) + (BETA_o**2) * (SE_e**2) / (BETA_e**4)
    weights = 1 / VBj

    bwm = mr_weighted_median(BETA_e, SE_e, BETA_o, SE_o, nboot)
    penalty = chi2.sf(weights * (betaIV - bwm[0]["b"]) ** 2, df=1)
    pen_weights = weights * np.minimum(1, penalty * penk)

    b = weighted_median(betaIV, pen_weights)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, pen_weights, nboot)
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


def mr_simple_median(BETA_e, SE_e, BETA_o, SE_o, nboot):
    """
    Perform a Mendelian Randomization analysis using the simple median method.

    Args:
        BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
        SE_e (numpy array): Standard errors corresponding to `BETA_e`.
        BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
        SE_o (numpy array): Standard errors corresponding to `BETA_o`.
        nboot (int): Number of boostrap iterations to obtain the standard error and p-value

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
        return [
            {
                "method": MR_METHODS_NAMES["Simple-median"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]

    b_iv = BETA_o / BETA_e
    weights = np.repeat(1 / len(BETA_e), len(BETA_e))
    b = weighted_median(b_iv, weights)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot)
    pval = 2 * norm.sf(abs(b / se))
    return [{"method": MR_METHODS_NAMES["Simple-median"], "b": b, "se": se, "pval": pval, "nSNP": l}]


"""
Regression methods
"""

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
        return {
            "method": MR_METHODS_NAMES["IVW"],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        }

    # Create weights and perform weighted regression
    weights = 1 / (SE_o**2)
    try:
        model = sm.WLS(BETA_o, BETA_e, weights=weights).fit()
    except:
        return {
            "method": MR_METHODS_NAMES["IVW"],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        }

    # Extract coefficients
    b = model.params.iloc[0]
    se = model.bse.iloc[0] / min(1, np.sqrt(model.mse_resid))
    pval = 2 * norm.sf(abs(b / se))

    Q_df = l - 1
    Q = model.scale * Q_df
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
        return {
            "method": MR_METHODS_NAMES["IVW-RE"],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        }

    # Create weights and perform weighted regression
    weights = 1 / (SE_o**2)
    model = sm.WLS(BETA_o, BETA_e, weights=weights).fit()

    # Extract coefficients
    b = model.params.iloc[0]
    se = model.bse.iloc[0]
    pval = 2 * norm.sf(abs(b / se))
    Q_df = l - 1
    Q = model.scale * Q_df
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
        return [
            {
                "method": MR_METHODS_NAMES["IVW-FE"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]

    # Create weights and perform weighted regression
    weights = 1 / SE_o**2
    try:
        model = sm.WLS(BETA_o, BETA_e, weights=weights).fit()
    except:
        return {
            "method": MR_METHODS_NAMES["IVW-FE"],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        }

    # Extract coefficients
    b = model.params.iloc[0]
    se = model.bse.iloc[0] / model.mse_resid**0.5
    pval = 2 * norm.sf(abs(b / se))
    Q_df = l - 1
    Q = model.scale * Q_df
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
            - "se": Adjusted standard error of the coefficient.
            - "pval": P-value for the causal estimate.
            - "nSNP": Number of genetic variants used in the analysis.
            - "Q": Cochran"s Q statistic for heterogeneity.
            - "Q_df": Degrees of freedom for the Q statistic.
            - "Q_pval": P-value for the Q statistic.

    Notes:
        The returned causal estimate is not weighted by the inverse variance.
        The standard error is corrected for under dispersion.
        Cochran"s Q statistics also computed to assess the heterogeneity across the instrumental variables.
    """

    l = len(BETA_e)
    if l < 2:
        return {
            "method": MR_METHODS_NAMES["UWR"],
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
            "nSNP": np.nan,
        }

    # Perform regression without weights
    model = sm.OLS(BETA_o, BETA_e).fit()

    # Extract coefficients and compute statistics
    b = model.params.iloc[0]
    se = model.bse.iloc[0] / min(1, model.mse_resid**0.5)  # Adjusted standard error
    pval = 2 * norm.sf(np.abs(b / se))
    Q_df = l - 1
    Q = model.scale * Q_df
    Q_pval = chi2.sf(Q, Q_df)

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
        return [
            {
                "method": MR_METHODS_NAMES["Sign"],
                "b": np.nan,
                "se": np.nan,
                "pval": np.nan,
                "nSNP": np.nan,
            }
        ]

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
