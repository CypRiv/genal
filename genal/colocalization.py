import numpy as np
import pandas as pd
from genal.geno_tools import check_beta_column, check_int_column, missing_snp_ids

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _subset_coloc_columns(data, merge_keys):
    """Keep only the columns needed for coloc merging and ABF calculations."""
    required_columns = list(merge_keys) + ["BETA", "SE"]
    optional_columns = [
        col
        for col in ["SNP", "CHR", "POS", "EA", "NEA", "EAF"]
        if col in data.columns and col not in required_columns
    ]
    columns_to_keep = required_columns + optional_columns
    return data.loc[:, columns_to_keep].copy()


def _prepare_coloc_dataset(data, merge_keys, label):
    """
    Check one dataset before merging (as coloc's check_dataset): remove variants with missing or non-finite BETA/SE
    or SE <= 0, variants with a missing merge key (CHR/POS, or a missing or placeholder SNP ID), and duplicated
    variants (same merge key and same alleles, the row with the lowest P-value is kept). Allele columns are upper-cased and
    values not made of A, C, G, T letters are set to NaN.
    """
    nrows = data.shape[0]
    data["SE"] = pd.to_numeric(data["SE"], errors="coerce")
    valid = np.isfinite(data["BETA"]) & np.isfinite(data["SE"]) & (data["SE"] > 0)
    n_invalid = (~valid).sum()
    if n_invalid > 0:
        print(
            f"Dataset {label}: {n_invalid}({n_invalid/nrows*100:.3f}%) variants with missing or non-finite BETA or SE, or SE <= 0, have been removed."
        )
    data = data.loc[valid].copy()

    if merge_keys == ["CHR", "POS"]:
        check_int_column(data, "CHR")
        check_int_column(data, "POS")
        missing_key = data["CHR"].isna() | data["POS"].isna()
        key_description = "CHR or POS"
    else:
        missing_key = missing_snp_ids(data["SNP"])
        key_description = "SNP ID (missing or placeholder such as '.')"
    n_missing_key = missing_key.sum()
    if n_missing_key > 0:
        print(f"Dataset {label}: {n_missing_key} variants without a valid {key_description} have been removed.")
    data = data.loc[~missing_key].copy()

    duplicate_keys = list(merge_keys)
    if "EA" in data.columns and "NEA" in data.columns:
        for col in ["EA", "NEA"]:
            alleles = data[col].astype(str).str.strip().str.upper()
            data[col] = alleles.where(alleles.str.fullmatch("[ACGT]+"))
        data["_allele_min"] = np.minimum(data["EA"].fillna(""), data["NEA"].fillna(""))
        data["_allele_max"] = np.maximum(data["EA"].fillna(""), data["NEA"].fillna(""))
        duplicate_keys += ["_allele_min", "_allele_max"]
    # Keep, for each duplicated variant, the row with the lowest P-value (largest |BETA/SE|; first row if tied)
    data = data.reset_index(drop=True)
    duplicated = (
        data.assign(_neg_abs_z=-(data["BETA"] / data["SE"]).abs())
        .sort_values("_neg_abs_z", kind="stable")
        .duplicated(subset=duplicate_keys, keep="first")
        .sort_index()
    )
    n_duplicated = duplicated.sum()
    if n_duplicated > 0:
        print(
            f"Warning: dataset {label} contains {n_duplicated} duplicated variants (same {'/'.join(merge_keys)}"
            f"{' and alleles' if len(duplicate_keys) > len(merge_keys) else ''}). For each variant, only the row with the lowest P-value (largest |BETA/SE|) is kept."
        )
    data = data.loc[~duplicated].drop(columns=["_allele_min", "_allele_max"], errors="ignore")
    return data


def _align_alleles(merged_data):
    """
    Keep the merged pairs whose alleles designate the same variant: {EA_1, NEA_1} = {EA_2, NEA_2} in either orientation,
    or, for non-palindromic SNVs, after taking the complementary strand of the dataset 2 alleles.
    BETA_2 is flipped when the effect alleles differ. Direct matches are placed first so that they are preferred over
    strand-complement matches when de-duplicating.

    Returns:
        tuple: (aligned DataFrame, number of removed pairs, number of strand-complement matches)
    """
    ea1, nea1, ea2, nea2 = (merged_data[col] for col in ["EA_1", "NEA_1", "EA_2", "NEA_2"])
    same = (ea1 == ea2) & (nea1 == nea2)
    swapped = (ea1 == nea2) & (nea1 == ea2)

    ea2_comp = ea2.str.translate(COMPLEMENT)
    nea2_comp = nea2.str.translate(COMPLEMENT)
    snv = (ea2.str.len() == 1) & (nea2.str.len() == 1)
    palindromic = ea2_comp == nea2
    complement_allowed = snv & ~palindromic & ~same & ~swapped
    same_comp = complement_allowed & (ea1 == ea2_comp) & (nea1 == nea2_comp)
    swapped_comp = complement_allowed & (ea1 == nea2_comp) & (nea1 == ea2_comp)

    keep = same | swapped | same_comp | swapped_comp
    merged_data.loc[swapped | swapped_comp, "BETA_2"] = -merged_data.loc[swapped | swapped_comp, "BETA_2"]
    merged_data["_complement_match"] = same_comp | swapped_comp
    aligned = merged_data.loc[keep].sort_values("_complement_match", kind="stable")
    n_complement = int(aligned["_complement_match"].sum())
    return aligned.drop(columns=["_complement_match"]), int((~keep).sum()), n_complement


def _resolve_sdY(data, sample_size, label):
    """
    Estimate sdY of a quantitative trait from the variance of the effect estimates, the EAF and the sample size
    of all the variants of the dataset (as coloc's sdY.est). EAF values outside (0, 1) are not used.
    If sdY cannot be estimated, sdY = 1 is assumed and a warning is printed.
    """
    assumption = (
        f"Assuming sdY{label} = 1, i.e. that the effects (BETA, SE) of trait {label} are on a standardized scale "
        f"(per standard deviation of the trait). If they are in other units, the prior on the effect sizes is misspecified "
        f"and the posterior probabilities can be wrong: provide sdY{label}, or n{label} with an EAF column so that sdY{label} can be estimated."
    )
    has_eaf = "EAF" in data.columns
    if sample_size is None or not has_eaf:
        if sample_size is None and not has_eaf:
            reason = f"n{label} is not provided and the trait {label} data has no EAF column"
        elif sample_size is None:
            reason = f"n{label} is not provided"
        else:
            reason = f"the trait {label} data has no EAF column"
        print(f"Warning: sdY{label} is not provided and cannot be estimated ({reason}). {assumption}")
        return 1

    eaf = pd.to_numeric(data["EAF"], errors="coerce")
    valid = (eaf > 0) & (eaf < 1)
    if not valid.any():
        print(f"Warning: sdY{label} is not provided and cannot be estimated (no EAF value strictly between 0 and 1 in the trait {label} data). {assumption}")
        return 1
    n_excluded = (~valid).sum()
    if n_excluded > 0:
        print(f"{n_excluded} variants of dataset {label} with a missing EAF or an EAF outside (0, 1) are not used to estimate sdY{label}.")

    sdY = sdY_est(data.loc[valid, "SE"] ** 2, eaf[valid], sample_size)
    if not np.isfinite(sdY):
        print(f"Warning: the estimated sdY{label} is not finite. {assumption}")
        return 1

    print(f"Using EAF and n{label} to estimate sdY{label}: {sdY:.2f}")
    return sdY


def _adjust_prior(p, nsnps, name, description):
    """Adjust a prior probability for very large regions, as coloc's adjust_prior."""
    if nsnps * p >= 1:
        adjusted = 1 / (nsnps + 1)
        print(f"{name} * number of {description} ({nsnps}) >= 1: setting {name} = 1/(nsnps + 1) = {adjusted:.3g}.")
        return adjusted
    return p


def coloc_abf_func(data1, data2, trait1_type="quant", trait2_type="quant",
                   sdY1=None, sdY2=None, n1=None, n2=None,
                   p1=1e-4, p2=1e-4, p12=1e-5, merge_on_snp=False):
    """
    Perform colocalization analysis between two GWAS datasets using approximate Bayes factors.
    Corresponds to the :meth:`Geno.colocalize` method.

    Args:
        data1: DataFrame containing GWAS results for trait 1
        data2: DataFrame containing GWAS results for trait 2
        trait1_type: Type of trait 1 ("quant" for quantitative traits or "cc" for case-control traits), default is "quant"
        trait2_type: Type of trait 2 ("quant" for quantitative traits or "cc" for case-control traits), default is "quant"
        sdY1: Standard deviation of trait 1 (required for quantitative traits, estimated from EAF and n1 or assumed to be 1 otherwise)
        sdY2: Standard deviation of trait 2 (required for quantitative traits, estimated from EAF and n2 or assumed to be 1 otherwise)
        n1: Sample size for trait 1 (used to estimate sdY if not provided)
        n2: Sample size for trait 2 (used to estimate sdY if not provided)
        p1: Prior probability SNP associated with trait 1, default is 1e-4
        p2: Prior probability SNP associated with trait 2, default is 1e-4
        p12: Prior probability SNP associated with both traits, default is 1e-5
        merge_on_snp: If True, merge the datasets on SNP column. If False, first attempt to merge on CHR and POS columns.

    Notes:
        As in coloc.abf, each dataset is checked before merging, sdY is estimated from all the variants of its dataset,
        p1 and p2 are adjusted with the number of variants of their dataset, and p12 with the number of overlapping variants.
    """
    for label, trait_type in (("1", trait1_type), ("2", trait2_type)):
        if trait_type not in ("quant", "cc"):
            raise ValueError(
                f"Invalid trait{label}_type: {trait_type!r}. Must be 'quant' (quantitative trait) or 'cc' (case-control trait)."
            )

    data1 = data1.reset_index(drop=True)
    data2 = data2.reset_index(drop=True)

    # Ensure that the BETA columns are preprocessed
    check_beta_column(data1, 'BETA')
    check_beta_column(data2, 'BETA')

    # Determine if we can merge on position, otherwise try SNP
    if all(col in data1.columns for col in ['CHR', 'POS']) and \
       all(col in data2.columns for col in ['CHR', 'POS']) and not merge_on_snp:
        print("Merging datasets using genomic positions (CHR, POS)")
        merge_keys = ["CHR", "POS"]
    elif 'SNP' in data1.columns and 'SNP' in data2.columns:
        if merge_on_snp:
            print("Merging datasets using SNP IDs.")
        else:
            print("Position columns (CHR, POS) not present in both datasets. Merging datasets using SNP IDs.")
        merge_keys = ["SNP"]
    else:
        raise ValueError("At least CHR/POS or SNP columns must be present in both datasets for colocalization analysis")

    # Restrict each dataset to coloc-relevant columns and check it before merging.
    data1 = _prepare_coloc_dataset(_subset_coloc_columns(data1, merge_keys), merge_keys, "1")
    data2 = _prepare_coloc_dataset(_subset_coloc_columns(data2, merge_keys), merge_keys, "2")
    n_snps1, n_snps2 = len(data1), len(data2)

    # Estimate sdY if not provided for quantitative traits (from all the variants of each dataset)
    if trait1_type == "quant" and sdY1 is None:
        sdY1 = _resolve_sdY(data1, n1, "1")
    if trait2_type == "quant" and sdY2 is None:
        sdY2 = _resolve_sdY(data2, n2, "2")

    data1 = data1.rename(columns={'EAF': 'EAF_1'})
    data2 = data2.rename(columns={'EAF': 'EAF_2'})
    merged_data = pd.merge(
        data1,
        data2,
        on=merge_keys,
        how='inner',
        suffixes=('_1', '_2')
    )

    # After merging, check if we can align alleles
    if all(col in merged_data.columns for col in ['EA_1', 'NEA_1', 'EA_2', 'NEA_2']):
        print("Aligning effect alleles between datasets")
        merged_data, n_unmatched, n_complement = _align_alleles(merged_data)
        if n_unmatched > 0:
            print(
                f"{n_unmatched} merged pairs have been removed because their alleles (EA, NEA) are missing or do not designate the same variant in both datasets."
            )
        if n_complement > 0:
            print(f"{n_complement} variants have been matched on the complementary strand (non-palindromic SNVs).")
    else:
        print("Allele columns (EA, NEA) not present in both datasets. "
              "This might lead to incorrect results if the effect estimates (BETA) were not obtained with the same reference allele in both datasets.")

    # One pair per variant: a variant (merge key and, when alleles are available, unordered allele pair) of either
    # dataset is used in a single pair. Multi-allelic sites keep one pair per matched variant. Direct allele matches
    # are placed first by _align_alleles and are preferred over strand-complement matches.
    if all(col in merged_data.columns for col in ['EA_1', 'NEA_1', 'EA_2', 'NEA_2']):
        n_before = len(merged_data)
        for k in ["1", "2"]:
            variant_keys = merge_keys + [f"_allele_min_{k}", f"_allele_max_{k}"]
            merged_data[f"_allele_min_{k}"] = np.minimum(merged_data[f"EA_{k}"], merged_data[f"NEA_{k}"])
            merged_data[f"_allele_max_{k}"] = np.maximum(merged_data[f"EA_{k}"], merged_data[f"NEA_{k}"])
            merged_data = merged_data.loc[~merged_data.duplicated(subset=variant_keys, keep="first")].copy()
        merged_data = merged_data.drop(columns=["_allele_min_1", "_allele_max_1", "_allele_min_2", "_allele_max_2"])
        n_duplicated = n_before - len(merged_data)
        description = f"{'/'.join(merge_keys)} and alleles"
    else:
        duplicated = merged_data.duplicated(subset=merge_keys, keep="first")
        n_duplicated = duplicated.sum()
        merged_data = merged_data.loc[~duplicated]
        description = '/'.join(merge_keys)
    if n_duplicated > 0:
        print(
            f"{n_duplicated} pairs with the same variant ({description}) as a previous pair have been removed (only the first matching pair is kept)."
        )

    # Clean up columns
    merged_data.drop(columns=["EA_2", "NEA_2", "SNP_2", "CHR_2", "POS_2"], inplace=True, errors='ignore')
    merged_data.rename(columns={"SNP_1": "SNP", "CHR_1": "CHR", "POS_1": "POS"}, inplace=True, errors='ignore')

    if merged_data.empty:
        raise ValueError("No overlapping variants found between the datasets")

    print(f"Using {len(merged_data)} overlapping variants for colocalization analysis")

    # Calculate Bayes factors for each dataset
    lABF_1 = approx_bf_estimates(merged_data['BETA_1'], merged_data['SE_1']**2,
                                trait_type=trait1_type, sdY=sdY1)
    lABF_2 = approx_bf_estimates(merged_data['BETA_2'], merged_data['SE_2']**2,
                                trait_type=trait2_type, sdY=sdY2)

    # Adjust priors based on the number of SNPs (as coloc.abf: dataset sizes for p1 and p2, overlap for p12)
    n_snps = len(merged_data)
    p1 = _adjust_prior(p1, n_snps1, "p1", "variants in dataset 1")
    p2 = _adjust_prior(p2, n_snps2, "p2", "variants in dataset 2")
    p12 = _adjust_prior(p12, n_snps, "p12", "overlapping variants")

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

def logsum_distinct_pairs(l1, l2):
    """
    Calculate log(sum over i != j of exp(l1_i + l2_j)) = log(sum_i exp(l1_i) * sum_{j != i} exp(l2_j)) without cancellation.

    The leave-one-out sums of exp(l2) are computed relative to the largest term of l2 (which they contain, so no
    cancellation occurs), and the leave-one-out sum of the largest term is computed directly. The result is -inf
    when there is a single variant.
    """
    l1 = np.asarray(l1, dtype=float)
    l2 = np.asarray(l2, dtype=float)
    if len(l2) < 2:
        return -np.inf
    k = int(np.argmax(l2))
    scaled = np.exp(l2 - l2[k])
    with np.errstate(divide="ignore"):  # the k-th value (possibly log(0)) is replaced below
        leave_one_out = l2[k] + np.log(np.sum(scaled) - scaled)
    leave_one_out[k] = logsum(np.delete(l2, k))
    return logsum(l1 + leave_one_out)

def combine_abf(l1, l2, p1, p2, p12):
    """Calculate posterior probabilities for different hypotheses"""
    lsum = l1 + l2

    lH0_abf = 0
    lH1_abf = np.log(p1) + logsum(l1)
    lH2_abf = np.log(p2) + logsum(l2)
    # H3: distinct causal variants i != j. Computed directly: as the difference logsum(l1) + logsum(l2) - logsum(l1 + l2),
    # it is lost by cancellation when one variant dominates both traits.
    lH3_abf = np.log(p1) + np.log(p2) + logsum_distinct_pairs(l1, l2)
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
