import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.stats import norm, chi2

from .config import *





def mr_egger_regression(df):
    # Initialize null result
    null_result = [{"method": "MR Egger", 'b': np.nan, 'se': np.nan, 'pval': np.nan,'nSNP': np.nan}, 
                   {"method": "Egger Intercept", 'b': np.nan, 'se': np.nan, 'pval': np.nan,'nSNP': np.nan}]
    
    # Early return if insufficient data
    if len(df) < 3:
        return null_result
    
    sign0 = np.sign(df['BETA_e'].replace(0, 1))
    
    df['BETA_o'] *= sign0
    df['BETA_e'] = np.abs(df['BETA_e'])
    df['flipped'] = sign0 == -1
    
    X = df['BETA_e']
    X = sm.add_constant(X)
    y = df['BETA_o']
    weights = 1 / df['SE_o'] ** 2
    
    mod = sm.WLS(y, X, weights=weights).fit()
    
    if len(mod.params) > 1:
        b = mod.params['BETA_e']
        se = mod.bse['BETA_e'] / min(1, np.sqrt(mod.mse_resid))
        pval = 2 * (1 - stats.t.cdf(np.abs(b / se), df=len(df) - 2))
        
        b_i = mod.params['const']
        se_i = mod.bse['const'] / min(1, np.sqrt(mod.mse_resid))
        pval_i = 2 * (1 - stats.t.cdf(np.abs(b_i / se_i), df=len(df) - 2))
        
        return [{"method": "MR Egger", 'b': b, 'se': se, 'pval': pval,'nSNP': len(df)}, 
                   {"method": "Egger Intercept", 'b': b_i, 'se': se_i, 'pval': pval_i,'nSNP': len(df)}]
    else:
        print("Warning: Collinearities in MR Egger, try LD pruning the exposure (can be done with .clump()).")
        return null_result



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

def weighted_median_bootstrap(BETA_e, BETA_o, SE_e, SE_o, weights, nboot):
    med = np.zeros(nboot)
    for i in range(nboot):
        BETA_e_boot = np.random.normal(loc=BETA_e, scale=SE_e)
        BETA_o_boot = np.random.normal(loc=BETA_o, scale=SE_o)
        betaIV_boot = BETA_o_boot / BETA_e_boot
        med[i] = weighted_median(betaIV_boot, weights)
    return np.std(med)

def mr_weighted_median(df, nboot):
    if len(df) < 3:
        return {"method": "Weighted median",'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nSNP': np.nan}

    BETA_e = df['BETA_e']
    BETA_o = df['BETA_o']
    SE_e = df['SE_e']
    SE_o = df['SE_o']

    b_iv = BETA_o / BETA_e
    VBj = (SE_o ** 2) / (BETA_e ** 2) + ((BETA_o ** 2) * (SE_e ** 2)) / (BETA_e ** 4)

    b = weighted_median(b_iv, 1 / VBj)
    se = weighted_median_bootstrap(BETA_e, BETA_o, SE_e, SE_o, 1 / VBj, nboot)

    pval = 2 * (1 - norm.cdf(abs(b / se)))

    return [{"method": "Weighted median", 'nSNP': len(df),'b': b, 'se': se, 'pval': pval}]

def mr_ivw(df):
    # If less than 2 valid rows, return NA values
    if len(df) < 2:
        return {'b': np.nan, 'se': np.nan, 'pval': np.nan, 'nsnp': np.nan}
    
    # Prepare data for linear model
    X = df['BETA_e']
    y = df['BETA_o']
    weights = 1 / (df['SE_o'] ** 2)
    
    # Fit weighted linear model
    model = sm.WLS(y, X, weights=weights).fit()
    
    # Extract coefficients
    b = model.params['BETA_e']
    se = model.bse['BETA_e'] / min(1, np.sqrt(model.mse_resid))
    pval = 2 * (1 - norm.cdf(abs(b / se)))
    
    return [{"method": "Inverse-Variance Weighted", 'nSNP': len(df),'b': b, 'se': se, 'pval': pval}]

    

    
    
    
    





def format_outcome(file):
    """
    Format a variant file to be used as MR outcome. Can read .vcf .vcf.gz/.tbi and .txt (if .txt, the column names should be in genal format: CHR, POS, EA, NEA, 
    file: path to the file to be used as outcome
    return the outcome dataframe in a dask array
    """
    if not os.file.exists(file):
        raise ValueError("The path does not lead to a file")
        

    
    
    return df







    