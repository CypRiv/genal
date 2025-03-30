import numpy as np
import pandas as pd
from numpy import exp, log
from genal.geno_tools import check_beta_column, check_allele_column, check_snp_column, check_int_column

# Currently does not support multi-allelic SNPs

def coloc_abf_func(data1, data2, trait1_type="quant", trait2_type="quant", 
                   sdY1=None, sdY2=None, n1=None, n2=None,
                   p1=1e-4, p2=1e-4, p12=1e-5, merge_on_snp=False):
    """
    Perform colocalization analysis between two GWAS datasets using approximate Bayes factors.
    Corresponds to the :meth:`Geno.colocalize` method.
    
    Args:
        data1: DataFrame containing GWAS results for trait 1
        data2: DataFrame containing GWAS results for trait 2
        trait1_type: Type of trait 1 ("quant" for quantitative traits or "cc" for case-control traits)
        trait2_type: Type of trait 2 ("quant" for quantitative traits or "cc" for case-control traits)
        sdY1: Standard deviation of trait 1 (required for quantitative traits)
        sdY2: Standard deviation of trait 2 (required for quantitative traits)
        n1: Sample size for trait 1 (used to estimate sdY if not provided)
        n2: Sample size for trait 2 (used to estimate sdY if not provided)
        p1: Prior probability SNP associated with trait 1
        p2: Prior probability SNP associated with trait 2
        p12: Prior probability SNP associated with both traits
        merge_on_snp: If True, merge the datasets on SNP column. If False, first attempt to merge on CHR and POS columns.

    """

    # Ensure that the BETA columns are preprocessed
    check_beta_column(data1, 'BETA', 'Fill')
    check_beta_column(data2, 'BETA', 'Fill')

    # Adjust EAF column names before merging in case one of the datasets does not have it
    if 'EAF' in data1.columns:
        data1.rename(columns={'EAF': 'EAF_1'}, inplace=True)
    if 'EAF' in data2.columns:
        data2.rename(columns={'EAF': 'EAF_2'}, inplace=True)

    # First determine if we can merge on position, otherwise try SNP
    if all(col in data1.columns for col in ['CHR', 'POS']) and \
       all(col in data2.columns for col in ['CHR', 'POS']) and not merge_on_snp:
        
        print("Merging datasets using genomic positions (CHR, POS)")
        
        # Ensure that the CHR and POS columns are preprocessed
        check_int_column(data1, "CHR")
        check_int_column(data2, "CHR")
        check_int_column(data1, "POS")
        check_int_column(data2, "POS")
        
        # Merge using position
        merged_data = pd.merge(
            data1,
            data2,
            on=['CHR', 'POS'],
            how='left',
            suffixes=('_1', '_2')
        )
        
    elif 'SNP' in data1.columns and 'SNP' in data2.columns:
        print("Position columns (CHR, POS) not present in both datasets. Merging datasets using SNP IDs.")
        
        # Ensure that the SNP column is preprocessed
        check_snp_column(data1)
        check_snp_column(data2)
        
        # Merge using SNP
        merged_data = pd.merge(
            data1,
            data2,
            on='SNP',
            suffixes=('_1', '_2')
        )
    
    else:
        raise ValueError("At least CHR/POS or SNP columns must be present in both datasets for colocalization analysis")

    # After merging, check if we can align alleles
    if all(col in merged_data.columns for col in ['EA_1', 'NEA_1', 'EA_2', 'NEA_2']):
        print("Aligning effect alleles between datasets")
        
        # Ensure allele columns are preprocessed
        check_allele_column(data1, "EA", keep_indel=False)
        check_allele_column(data1, "NEA", keep_indel=False)
        check_allele_column(data2, "EA", keep_indel=False)
        check_allele_column(data2, "NEA", keep_indel=False)
        
        # Adjust BETA from trait 2 to correspond to the same effect allele as trait 1
        conditions = [
            merged_data["EA_1"] == merged_data["EA_2"],
            merged_data["EA_1"] == merged_data["NEA_2"],
            True,
        ]
        choices = [
            merged_data["BETA_2"],
            -merged_data["BETA_2"],
            np.nan,
        ]
        merged_data["BETA_2"] = np.select(conditions, choices)
    else:
        print("Allele columns (EA, NEA) not present in both datasets. "
              "This might lead to incorrect results if the effect estimates (BETA) were not obtained with the same reference allele in both datasets.")

    # Clean up columns
    merged_data.drop(columns=["EA_2", "NEA_2", "SNP_2", "CHR_2", "POS_2"], inplace=True, errors='ignore')
    merged_data.rename(columns={"SNP_1": "SNP", "CHR_1": "CHR", "POS_1": "POS"}, inplace=True, errors='ignore')

    # Drop any rows with duplicate values
    if "SNP" in merged_data.columns:    
        merged_data.drop_duplicates(subset=['SNP'], keep='first', inplace=True)
    if "CHR" in merged_data.columns and "POS" in merged_data.columns:
        merged_data.drop_duplicates(subset=["CHR", "POS"], keep='first', inplace=True)

    # Drop any rows with missing values
    merged_data = merged_data.dropna()
    if merged_data.empty:
        raise ValueError("No overlapping variants found between the datasets")
    
    print(f"Using {len(merged_data)} overlapping variants for colocalization analysis")
    
    # Estimate sdY if not provided for quantitative traits
    if trait1_type == "quant" and sdY1 is None:
        if 'EAF_1' not in merged_data.columns or n1 is None:
            print("Neither sdY1 nor EAF and n1 are provided for trait 1. Assuming sdY1 = 1.")
            sdY1 = 1
        else:
            sdY1 = sdY_est(merged_data['SE_1']**2, merged_data['EAF_1'], n1)
            print(f"Using EAF and n1 to estimate sdY1: {sdY1:.2f}")
        
    if trait2_type == "quant" and sdY2 is None:
        if 'EAF_2' not in merged_data.columns or n2 is None:
            print("Neither sdY2 nor EAF and n2 are provided for trait 2. Assuming sdY2 = 1.")
            sdY2 = 1
        else:
            sdY2 = sdY_est(merged_data['SE_2']**2, merged_data['EAF_2'], n2)
            print(f"Using EAF and n2 to estimate sdY2: {sdY2:.2f}")
    
    # Calculate Bayes factors for each dataset
    lABF_1 = approx_bf_estimates(merged_data['BETA_1'], merged_data['SE_1']**2, 
                                trait_type=trait1_type, sdY=sdY1)
    lABF_2 = approx_bf_estimates(merged_data['BETA_2'], merged_data['SE_2']**2, 
                                trait_type=trait2_type, sdY=sdY2)
    
    # Adjust priors based on number of SNPs
    n_snps = len(merged_data)
    if n_snps * p1 >= 1:
        p1 = 1 / (n_snps + 1)
    if n_snps * p2 >= 1:
        p2 = 1 / (n_snps + 1)
    if n_snps * p12 >= 1:
        p12 = 1 / (n_snps + 1)
    
    # Calculate posterior probabilities
    pp = combine_abf(lABF_1, lABF_2, p1, p2, p12)
    
    # Add SNP-specific results
    results_df = merged_data.copy()
    results_df['lABF_1'] = lABF_1
    results_df['lABF_2'] = lABF_2
    results_df['internal.sum.lABF'] = lABF_1 + lABF_2
    
    # Calculate SNP-specific PP for H4
    my_denom_log_abf = logsum(results_df['internal.sum.lABF'])
    results_df['SNP.PP.H4'] = np.exp(results_df['internal.sum.lABF'] - my_denom_log_abf)
    
    return {
            'nsnps': n_snps,
            **pp
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