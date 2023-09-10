import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from scipy.stats import norm, chi2, binom_test
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
from numpy.random import default_rng
from functools import partial


"""
Egger regression methods
"""

def mr_egger_regression(BETA_e, SE_e, BETA_o, SE_o):
    # Initialize null result
    null_result = [{"method": "MR Egger", 'b': np.nan, 'se': np.nan, 'pval': np.nan,'nSNP': np.nan}, 
                   {"method": "Egger Intercept", 'b': np.nan, 'se': np.nan, 'pval': np.nan,'nSNP': np.nan}]
    
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
    weights = 1 / SE_o ** 2
    
    mod = sm.WLS(y, X, weights=weights).fit()

    if len(mod.params) > 1:
        b = mod.params[1]
        se = mod.bse[1] / min(1, np.sqrt(mod.mse_resid))
        pval = 2 * (1 - stats.t.cdf(np.abs(b / se), df=l - 2))
        
        b_i = mod.params[0]
        se_i = mod.bse[0] / min(1, np.sqrt(mod.mse_resid))
        pval_i = 2 * (1 - stats.t.cdf(np.abs(b_i / se_i), df=l - 2))
        
        Q = mod.mse_resid * (l - 2)
        Q_df = l - 2
        Q_pval = 1 - chi2.cdf(Q, Q_df)
        
        return [{"method": "MR Egger", 'b': b, 'se': se, 'pval': pval,'nSNP': l, 'Q': Q, 'Q_df': Q_df, 'Q_pval': Q_pval}, 
                   {"method": "Egger Intercept", 'b': b_i, 'se': se_i, 'pval': pval_i,'nSNP': l, 'Q': Q, 'Q_df': Q_df, 'Q_pval': Q_pval}]
    else:
        print("Warning: Collinearities in MR Egger, try LD pruning the exposure (can be done with .clump()).")
        return null_result

def parallel_bootstrap_func(i, BETA_e, SE_e, BETA_o, SE_o):
    xs = np.random.normal(BETA_e, SE_e)
    ys = np.random.normal(BETA_o, SE_o)

    # Use absolute values for Egger reg
    ys *= np.sign(xs)
    xs = np.abs(xs)

    r = linreg(xs, ys, 1/SE_o**2)
    return r["ahat"], r["bhat"]

def mr_egger_regression_bootstrap(BETA_e, SE_e, BETA_o, SE_o, nboot, cpus = 4):
    l = len(BETA_e)
    if l < 3:
        return [{"method": "MR Egger bootstrap", 'b': np.nan, 'se': np.nan, 'pval': np.nan,'nSNP': np.nan}, 
                {"method": "Egger Intercept bootstrap", 'b': np.nan, 'se': np.nan, 'pval': np.nan,'nSNP': np.nan}]
    
    res = np.zeros((nboot+1, 2))

    with ProcessPoolExecutor(max_workers=cpus) as executor:
        # Start the load operations and mark each future with its URL
        futures = {executor.submit(parallel_bootstrap_func, i, BETA_e, SE_e, BETA_o, SE_o): i for i in range(nboot)}
        for future in tqdm(as_completed(futures), total=nboot, desc="MR Egger bootstrapping", ncols=100):
            ahat, bhat = future.result()
            i = futures[future]  # get the original index/counter
            res[i, 0] = ahat
            res[i, 1] = bhat

    return [{"method": "MR Egger bootstrap", 'b': np.nanmean(res[:, 1]), 'se': np.nanstd(res[:, 1]), 'pval': np.sum(np.sign(np.nanmean(res[:, 1])) * res[:, 1] < 0) / nboot, 'nSNP': l}, 
            {"method": "Egger Intercept bootstrap", 'b': np.nanmean(res[:, 0]), 'se': np.nanstd(res[:, 0]), 'pval': np.sum(np.sign(np.nanmean(res[:, 0])) * res[:, 0] < 0) / nboot,'nSNP': l}]


def linreg(x, y, w=None):
    if w is None:
        w = np.ones_like(x)

    xp = w * x
    yp = w * y

    bhat = np.cov(xp, yp)[0, 1] / np.var(xp)
    ahat = np.nanmean(y) - np.nanmean(x) * bhat
    yhat = ahat + bhat * x

    residuals = y - yhat
    se = np.sqrt(sum(w * residuals ** 2) / (np.sum(~np.isnan(yhat)) - 2) / np.sum(w * x**2))
    pval = 2 * (1 - norm.cdf(abs(bhat / se)))

    return {"ahat": ahat, "bhat": bhat, "se": se, "pval": pval}


"""
Median methods
"""

def weighted_median(b_iv, weights):
    order_indices = np.argsort(b_iv)
    betaIV_order = np.array(b_iv)[order_indices]
    weights_order = np.array(weights)[order_indices]
    weights_sum = np.cumsum(weights_order) - 0.5 * weights_order
    weights_sum /= np.sum(weights_order)
    below = np.max(np.where(weights_sum < 0.5))
    b = betaIV_order[below] + (betaIV_order[below + 1] - betaIV_order[below]) * \
        (0.5 - weights_sum[below]) / (weights_sum[below + 1] - weights_sum[below])
    return b

def weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot):
    med = np.zeros(nboot)
    for i in tqdm(range(nboot), total=nboot, desc="Weighted median bootstrapping", ncols=100):
        BETA_e_boot = np.random.normal(loc=BETA_e, scale=SE_e)
        BETA_o_boot = np.random.normal(loc=BETA_o, scale=SE_o)
        betaIV_boot = BETA_o_boot / BETA_e_boot
        med[i] = weighted_median(betaIV_boot, weights)
    return np.std(med)

def mr_weighted_median(BETA_e, SE_e, BETA_o, SE_o, nboot):
    l = len(BETA_e)
    if l < 3:
        return [{"method": "Weighted median",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nSNP': np.nan}]

    b_iv = BETA_o / BETA_e
    VBj = (SE_o ** 2) / (BETA_e ** 2) + ((BETA_o ** 2) * (SE_e ** 2)) / (BETA_e ** 4)

    b = weighted_median(b_iv, 1 / VBj)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, 1 / VBj, nboot)
    pval = 2 * (1 - norm.cdf(abs(b / se)))
    return [{"method": "Weighted Median", 'nSNP': l,'b': b, 'se': se, 'pval': pval}]

def mr_pen_wm(BETA_e, SE_e, BETA_o, SE_o, nboot, penk):
    l = len(BETA_e)
    if l < 3:
        return [{"method": "Penalised weighted median",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nSNP': np.nan}]
    
    betaIV = BETA_o / BETA_e  # Ratio estimates
    betaIVW = np.sum(BETA_o * BETA_e / SE_o**2) / np.sum(BETA_e**2 / SE_o**2)  # IVW estimate
    VBj = (SE_o**2) / (BETA_e**2) + (BETA_o**2) * (SE_e**2) / (BETA_e**4)
    weights = 1 / VBj
    
    bwm = mr_weighted_median(BETA_e, SE_e, BETA_o, SE_o, nboot)
    penalty = chi2.sf(weights * (betaIV - bwm[0]['b'])**2, df=1)
    pen_weights = weights * np.minimum(1, penalty * penk)  # Penalized weights

    b = weighted_median(betaIV, pen_weights)  # Penalized weighted median estimate
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, pen_weights, nboot)
    pval = 2 * (1 - norm.cdf(abs(b / se)))
    
    return [{"method": "Penalised weighted median",'b': b, 'se': se, 'pval': pval, 'nSNP': l}]

def mr_simple_median(BETA_e, SE_e, BETA_o, SE_o, nboot):
    l = len(BETA_e)
    if l < 3:
        return [{"method": "Simple median",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nSNP': np.nan}]

    b_iv = BETA_o / BETA_e
    weights = np.repeat(1/len(BETA_e), len(BETA_e))
    b = weighted_median(b_iv, weights)
    se = weighted_median_bootstrap(BETA_e, SE_e, BETA_o, SE_o, weights, nboot)
    pval = 2 * (1 - norm.cdf(abs(b/se)))
    return [{"method": "Simple median",'b': b, 'se': se, 'pval': pval, 'nSNP': l}]


"""
Regression methods
"""

def mr_ivw(BETA_e, SE_e, BETA_o, SE_o):
    """
    Inverse-variance weighted regression with multiplicative random effects.
    Standard Error is corrected for under dispersion (as opposed to mr_ivw_re()).
    """
    # If less than 2 valid rows, return NA values
    l = len(BETA_e)
    if l < 2:
        return {"method": "Inverse-Variance Weighted", 'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}
    
    # Create weights and perform weighted regression
    weights = 1 / (SE_o ** 2)
    model = sm.WLS(BETA_o, BETA_e, weights=weights).fit()
    
    # Extract coefficients
    b = model.params[0]
    se = model.bse[0] / min(1, np.sqrt(model.mse_resid))
    pval = 2 * (1 - norm.cdf(abs(b / se)))
    
    Q_df = l - 1
    Q = model.scale * Q_df
    Q_pval = 1 - chi2.cdf(Q, Q_df)
    
    return [{"method": "Inverse-Variance Weighted", 'nSNP': l,'b': b, 'se': se, 'pval': pval, "Q_df": Q_df, "Q": Q, "Q_pval": Q_pval}]

def mr_ivw_re(BETA_e, SE_e, BETA_o, SE_o):
    """
    Inverse-variance weighted regression with multiplicative random effects.
    Standard Error is not corrected for under dispersion (as opposed to mr_ivw()).
    """
    # If less than 2 valid rows, return NA values
    l = len(BETA_e)
    if l < 2:
        return {"method": "Inverse-Variance Weighted (Random effects)",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}
    
    # Create weights and perform weighted regression
    weights = 1 / (SE_o ** 2)
    model = sm.WLS(BETA_o, BETA_e, weights=weights).fit()
    
    # Extract coefficients
    b = model.params[0]
    se = model.bse[0]
    pval = 2 * (1 - norm.cdf(abs(b / se)))
    Q_df = l - 1
    Q = model.scale * Q_df
    Q_pval = chi2.sf(Q, Q_df)
    
    return [{"method": "Inverse-Variance Weighted (Random effects)", 'nSNP': l,'b': b, 'se': se, 'pval': pval, "Q_df": Q_df, "Q": Q, "Q_pval": Q_pval}]

def mr_ivw_fe(BETA_e, SE_e, BETA_o, SE_o):
    """
    Perform a Mendelian Randomization analysis using the Inverse Variance Weighted (IVW) method with fixed effects.    
    Parameters:
    - BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
    - SE_e (numpy array): Standard errors associated with BETA_e.
    - BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
    - SE_o (numpy array): Standard errors associated with BETA_o.    
    Returns:
    - list of dict: A list containing a dictionary with the results:
        - 'method': Name of the method.
        - 'nSNP': Number of genetic variants used in the analysis.
        - 'b': Estimated causal effect size.
        - 'se': Standard error of the estimated effect size.
        - 'pval': P-value for the estimated effect size.
        - 'Q_df': Degrees of freedom for Cochran's Q heterogeneity statistic.
        - 'Q': Cochran's Q heterogeneity statistic.
        - 'Q_pval': P-value for the heterogeneity statistic.

    Notes:
    The function uses weighted least squares regression (WLS) to estimate the causal effect size,
    weighting by the inverse of the variance of the outcome's effect sizes. It also calculates the
    heterogeneity of the instrumental variable estimates using Cochran's Q statistic.
    """    
    l = len(BETA_e)
    if l < 2:
        return [{"method": "Inverse Variance weighted (Fixed effects)",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nSNP': np.nan}]
    # Create weights and perform weighted regression
    weights = 1 / SE_o ** 2   
    model = sm.WLS(BETA_o, BETA_e, weights=weights).fit()
    
    # Extract coefficients
    b = model.params[0]
    se = model.bse[0] / model.mse_resid**0.5
    pval = 2 * norm.sf(np.abs(b/se))
    Q_df = l - 1
    Q = model.scale * Q_df
    Q_pval = chi2.sf(Q, Q_df)
    return [{"method": "Inverse Variance weighted (Fixed effects)", 'nSNP': l,'b': b, 'se': se, 'pval': pval, "Q_df": Q_df, "Q": Q, "Q_pval": Q_pval}]    

def mr_uwr(BETA_e, SE_e, BETA_o, SE_o):
    
    l = len(BETA_e)
    if l < 2:
        return {"method": "Unweighted regression", 'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}

    # Perform regression without weights
    model = sm.OLS(BETA_o, BETA_e).fit()

    # Extract coefficients and compute statistics
    b = model.params[0]
    se = model.bse[0] / min(1, model.mse_resid**0.5)  # Adjusted standard error
    pval = 2 * norm.sf(np.abs(b/se))
    Q_df = l - 1
    Q = model.scale * Q_df
    Q_pval = chi2.sf(Q, Q_df)

    return [{"method": "Unweighted regression",'b': b, 'se': se, 'pval': pval, 'nSNP': l, 'Q': Q, 'Q_df': Q_df, 'Q_pval': Q_pval}]


""" 
Sign method
"""

def mr_sign(BETA_e, BETA_o):
    """
    Performs the sign concordance test. This test is used to determine whether there is a consistent direction of effect (i.e., sign) between the exposure and outcome across the variants. b takes value between -1 (all SNPs have opposite signs between exposure and outcome effects) and 1 (all SNPs have the same sign between exposure and outcome effects). The consistent directonality is an assumption of Mendelian Randomization.
    Parameters:
    - BETA_e (numpy array): Effect sizes of genetic variants on the exposure.
    - BETA_o (numpy array): Effect sizes of the same genetic variants on the outcome.
    Returns:
    - list of dict: A list containing a dictionary with the results:
        - 'method': Name of the method.
        - 'nSNP': Number of genetic variants used in the analysis.
        - 'b': Proportion of concordant signs between exposure and outcome effect sizes minus 0.5, multiplied by 2.
        - 'se': Not applicable for this method (returns NaN).
        - 'pval': P-value for the sign concordance test based on a binomial distribution.
    
    Notes:
    The function tests the hypothesis that the signs of the genetic associations with the exposure and the 
    outcome are more concordant than expected by chance. Effect sizes that are exactly zero are replaced with NaN and 
    are not included in the analysis. A binomial test is then performed to evaluate the probability of observing 
    the given number of concordant signs by chance alone, assuming a null expectation of 50% concordance.
    """

    # Replace zeros with NaNs
    BETA_e = np.where(BETA_e == 0, np.nan, BETA_e)
    BETA_o = np.where(BETA_o == 0, np.nan, BETA_o)
    
    # Check for enough non-missing values
    valid_data = (~np.isnan(BETA_e)) & (~np.isnan(BETA_o))
    if np.sum(valid_data) < 6:
        return [{"method": "Sign concordance test",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nSNP': np.nan}]
    
    # Count the number of consistent signs
    x = np.sum(np.sign(BETA_e[valid_data]) == np.sign(BETA_o[valid_data]))
    n = np.sum(valid_data)
    
    # Binomial test
    pval = binom_test(x, n, p=0.5)
    b = (x / n - 0.5) * 2
    
    return [{"method": "Sign concordance test", 'nSNP': n, 'b': b, 'se': np.nan, 'pval': pval}]
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







    