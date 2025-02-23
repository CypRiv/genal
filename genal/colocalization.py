import numpy as np
import pandas as pd
from numpy import exp, log




def coloc_abf_func(data, trait1_type="quant", trait2_type="quant", 
                   sdY1=None, sdY2=None, n1=None, n2=None,
                   p1=1e-4, p2=1e-4, p12=1e-5):
    """
    Perform colocalization analysis between two GWAS datasets using approximate Bayes factors.
    
    Args:
        data: DataFrame containing merged GWAS results
        trait1_type: Type of trait 1 ("quant" or "cc")
        trait2_type: Type of trait 2 ("quant" or "cc")
        sdY1: Standard deviation of trait 1 (required for quantitative traits)
        sdY2: Standard deviation of trait 2 (required for quantitative traits)
        n1: Sample size for trait 1 (used to estimate sdY if not provided)
        n2: Sample size for trait 2 (used to estimate sdY if not provided)
        p1: Prior probability SNP associated with trait 1
        p2: Prior probability SNP associated with trait 2
        p12: Prior probability SNP associated with both traits
    """
    # Estimate sdY if not provided for quantitative traits
    if trait1_type == "quant" and sdY1 is None:
        if 'EAF_1' not in data.columns or n1 is None:
            print("Neither sdY1 nor EAF and n1 are provided for trait 1. Assuming sdY1 = 1.")
            sdY1 = 1
        else:
            sdY1 = sdY_est(data['SE_1']**2, data['EAF_1'], n1)
            print(f"Using EAF and n1 to estimate sdY1: {sdY1:.2f}")
        
    if trait2_type == "quant" and sdY2 is None:
        if 'EAF_2' not in data.columns or n2 is None:
            print("Neither sdY2 nor EAF and n2 are provided for trait 2. Assuming sdY2 = 1.")
            sdY2 = 1
        else:
            sdY2 = sdY_est(data['SE_2']**2, data['EAF_2'], n2)
            print(f"Using EAF and n2 to estimate sdY2: {sdY2:.2f}")
    # Calculate Bayes factors for each dataset
    lABF_1 = approx_bf_estimates(data['BETA_1'], data['SE_1']**2, 
                                trait_type=trait1_type, sdY=sdY1)
    lABF_2 = approx_bf_estimates(data['BETA_2'], data['SE_2']**2, 
                                trait_type=trait2_type, sdY=sdY2)
    
    # Adjust priors based on number of SNPs
    n_snps = len(data)
    if n_snps * p1 >= 1:
        p1 = 1 / (n_snps + 1)
    if n_snps * p2 >= 1:
        p2 = 1 / (n_snps + 1)
    if n_snps * p12 >= 1:
        p12 = 1 / (n_snps + 1)
    
    # Calculate posterior probabilities
    pp = combine_abf(lABF_1, lABF_2, p1, p2, p12)
    
    # Add SNP-specific results
    results_df = data.copy()
    results_df['lABF_1'] = lABF_1
    results_df['lABF_2'] = lABF_2
    results_df['internal.sum.lABF'] = lABF_1 + lABF_2
    
    # Calculate SNP-specific PP for H4
    my_denom_log_abf = logsum(results_df['internal.sum.lABF'])
    results_df['SNP.PP.H4'] = np.exp(results_df['internal.sum.lABF'] - my_denom_log_abf)
    
    return {
        'summary': {
            'nsnps': n_snps,
            **pp
        },
        'results': results_df,
        'priors': {
            'p1': p1,
            'p2': p2,
            'p12': p12
        }
    }

def approx_bf_estimates(beta, varbeta, trait_type="quant", sdY=1, effect_priors={'quant': 0.15, 'cc': 0.2}):
    """
    Calculate approximate Bayes factors using regression estimates.
    
    Args:
        beta: effect size estimate
        varbeta: variance of the effect size estimate
        trait_type: either "quant" for quantitative trait or "cc" for case-control
        sdY: standard deviation of the trait (for quantitative traits)
        effect_priors: dictionary with prior effect sizes for quantitative and case-control traits
        
    Returns:
        array: log approximate Bayes factors
    """
    z = beta / np.sqrt(varbeta)
    
    # Set prior standard deviation based on trait type
    if trait_type == "quant":
        sd_prior = effect_priors['quant'] * sdY
    else:  # case-control
        sd_prior = effect_priors['cc']
        
    r = sd_prior**2 / (sd_prior**2 + varbeta)
    lABF = 0.5 * (np.log(1 - r) + (r * z**2))
    return lABF

def logsum(x):
    """Calculate log of sum of exponentials"""
    my_max = np.max(x)
    return my_max + np.log(np.sum(np.exp(x - my_max)))

def logdiff(x, y):
    """Calculate log of difference of exponentials"""
    my_max = max(x, y)
    return my_max + np.log(exp(x - my_max) - np.exp(y - my_max))

def combine_abf(l1, l2, p1, p2, p12):
    """Calculate posterior probabilities for different hypotheses"""
    lsum = l1 + l2
    
    lH0_abf = 0
    lH1_abf = np.log(p1) + logsum(l1)
    lH2_abf = np.log(p2) + logsum(l2)
    lH3_abf = np.log(p1) + np.log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
    lH4_abf = np.log(p12) + logsum(lsum)
    
    all_abf = np.array([lH0_abf, lH1_abf, lH2_abf, lH3_abf, lH4_abf])
    denom_log_abf = logsum(all_abf)
    pp_abf = np.exp(all_abf - denom_log_abf)
    
    return {
        'PP.H0.abf': pp_abf[0],
        'PP.H1.abf': pp_abf[1],
        'PP.H2.abf': pp_abf[2],
        'PP.H3.abf': pp_abf[3],
        'PP.H4.abf': pp_abf[4]
    }

def sdY_est(vbeta, maf, n):
    """
    Estimate trait standard deviation given vectors of variance of coefficients, MAF and sample size.
    
    Args:
        vbeta: vector of variance of coefficients
        maf: vector of MAF (same length as vbeta)
        n: sample size
        
    Returns:
        float: estimated standard deviation of Y
    """
    oneover = 1/vbeta
    nvx = 2 * n * maf * (1-maf)
    # Fit linear regression through origin
    coef = np.sum(nvx * oneover) / np.sum(oneover**2)
    if coef < 0:
        raise ValueError("Estimated sdY is negative - this can happen with small datasets, or those with errors. A reasonable estimate of sdY is required to continue.")
    return np.sqrt(coef)