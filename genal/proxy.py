import pandas as pd
import numpy as np
import os
import subprocess
import re
import uuid

from .tools import get_reference_panel_path, get_plink_path, run_plink_command

## TO DO: accept lists of CHR/POS instead of SNP names for these functions


def query_outcome_proxy(df, ld, snps_to_extract, snps_df=[]):
    """
    Extract the best proxies from a dataframe, as well as specific SNPs.

    Given a dataframe `df` (originating from Geno.data) and a dataframe of potential proxies
    (output from `find_proxies`), this function extracts the best proxies from `df` as well as
    the SNPs specified in `snps_to_extract`.
    This is suited for querying outcome data.

    Args:
        df (pd.DataFrame): Dataframe of SNP information with the usual Geno columns
                           (SNP, BETA, SE, EAF, EA, NEA). EAF is not necessary.
        ld (pd.DataFrame): Dataframe of proxies (output from `find_proxies`).
        snps_to_extract (list): List of SNPs to extract in addition to the proxies.
        snps_df (list, optional): List of SNPs to choose the proxy from. Should be the list of
                                  SNPs in df. Can be provided to avoid recomputing it. Defaults to an empty list.

    Returns:
        pd.DataFrame: Dataframe with queried SNPs and their proxies.
    """
    # If ld is None
    if not isinstance(ld, pd.DataFrame):
        raise ValueError("ld is None (The SNPs to be proxied were not found in the reference panel)")

    # If snps_df is empty, populate it with SNPs from df
    if not snps_df:
        snps_df = df.SNP.values

    # Filter proxies that are present in df
    ld = ld[ld.SNP_B.isin(snps_df)]

    # Remove original SNPs
    ld = ld[ld["SNP_A"] != ld["SNP_B"]]

    # Sort by r and select the best proxy for each SNP
    ld = ld.reindex(ld["R"].abs().sort_values(ascending=False).index)
    ld = ld.groupby("SNP_A").first().reset_index(drop=False)

    # Determine SNPs to query
    snps_to_query = set(snps_to_extract) | set(ld.SNP_B.values)
    df_queried = df[df.SNP.isin(snps_to_query)]

    # Merge dataframes and identify proxies
    output = df_queried.merge(ld, how="left", left_on="SNP", right_on="SNP_B")
    output["proxy"] = output["SNP_B"].notnull()

    # In the plink output, the alleles taken as reference for proxying are "MAJ_A" and "MAJ_B" (major alleles in the reference panel)
    # We want to use as effect allele for the original SNP its minor allele in the reference panel
    # So, we flip BETA if the proxied SNP's effect allele is the major allele in the reference panel
    conditions = [
        output["EA"] == output["MAJ_B"],
        output["EA"] == output["NONMAJ_B"],
        ~output["proxy"],
        True,
    ]
    choices = [
        -output["BETA"],  # if EA == MAJ_B, flip the sign of BETA 
        output["BETA"],  # if EA == NONMAJ_B, BETA does not change
        output["BETA"],  # if SNP_B is NaN (The original SNP was not proxied), BETA does not change
        np.nan,  # if the original SNP was proxied but "EA" is neither "MAJ_A" nor "NONMAJ_A", BETA is NaN
    ]
    output["BETA"] = np.select(conditions, choices)

    # Flip BETA if the sign of R is negative: indicates that the positive correlation corresponds to MAJ_A with NONMAJ_B
    sign_r = np.sign(output["R"]) # Sign of R
    output["BETA"] = np.where(sign_r == -1, -output["BETA"], output["BETA"])

    # Delete SNPs with mismatched alleles
    nrow = output.shape[0]
    output = output.dropna(subset=["BETA"])
    if output.shape[0] < nrow:
        print(
            f"Deleted {nrow-output.shape[0]} base SNPs that did not have matching alleles in reference data."
        )
    print(f"Found proxies for {output['proxy'].sum()} SNPs.")

    # Replace the proxied SNPs with the position and alleles of the original SNPs
    output["SNP"] = np.where(output["proxy"], output["SNP_A"], output["SNP"])
    output["POS"] = np.where(output["proxy"], output["BP_A"], output["POS"])
    output["CHR"] = np.where(output["proxy"], output["CHR_A"], output["CHR"])
    output["EA"] = np.where(output["proxy"], output["NONMAJ_A"], output["EA"])
    output["NEA"] = np.where(output["proxy"], output["MAJ_A"], output["NEA"])
    if "EAF" in output.columns:
        output["EAF"] = np.where(output["proxy"], output["NONMAJ_FREQ_A"], output["EAF"])

    # Drop columns related to ld
    output = output.drop(columns=ld.columns)

    return output


def apply_proxies(df, ld, searchspace=None):
    """
    Given a dataframe (coming from GENO.data attribute) and a dataframe of proxies
    (output from find_proxies), replace the SNPs in df with their best proxies, if they exist.
    This function is suited for exposure data (before running a PRS for instance).

    Args:
        df (DataFrame): Dataframe of SNP information with the usual GENO columns (SNP, BETA, SE, EAF, EA, NEA). EAF is not necessary.
        ld (DataFrame): Dataframe of proxies (output from find_proxies).
        searchspace (list, optional): List of SNPs to restrict the list of potential proxies. By default, includes all the proxies found. Using a searchspace can be done either at the find_proxies step or at this step, but it is much faster to use it at this step.

    Returns:
        DataFrame: A DataFrame with SNPs replaced by their best proxies, if they exist.
    """
    # If ld is None
    if not isinstance(ld, pd.DataFrame):
        raise ValueError("ld is None (The SNPs to be proxied were not found in the reference panel)")
    
    # Check mandatory columns
    mandatory_cols = ["EA", "SNP", "BETA"]
    for col in mandatory_cols:
        if col not in df.columns:
            raise ValueError(f"The column {col} is not found in the data!")
    
    # Filter by searchspace if provided
    if searchspace:
        print("Filtering the potential proxies with the searchspace provided.")
        ld = ld[ld.SNP_B.isin(searchspace)]

    # Remove original SNPs and sort by r
    ld = ld[ld["SNP_A"] != ld["SNP_B"]]
    ld = ld.reindex(ld["R"].abs().sort_values(ascending=False).index)

    # Select the best proxy for each SNP
    ld = ld.groupby("SNP_A").first().reset_index(drop=False)

    # Merge the dataframes
    output = df.merge(ld, how="left", left_on="SNP", right_on="SNP_A")
    output["proxy"] = pd.notnull(output["SNP_B"])

    # In the plink output, the alleles taken as reference for proxying are "MAJ_A" and "MAJ_B" (major alleles in the reference panel)
    # We want to use as effect allele for the proxy SNP its minor allele in the reference panel
    # So, we flip BETA if the original SNP's effect allele is the major allele in the reference panel
    conditions = [
        output["EA"] == output["MAJ_A"],
        output["EA"] == output["NONMAJ_A"],
        ~output["proxy"],
        True,
    ]
    choices = [
        -output["BETA"],  # if EA == MAJ_A, flip the sign of BETA 
        output["BETA"],  # if EA == NONMAJ_A, BETA does not change
        output["BETA"],  # if SNP_B is NaN (The original SNP was not proxied), BETA does not change
        np.nan,  # if the original SNP was proxied but "EA" is neither "MAJ_A" nor "NONMAJ_A", BETA is NaN
    ]
    output["BETA"] = np.select(conditions, choices)

    # Flip BETA if the sign of R is negative: indicates that the positive correlation corresponds to MAJ_A with NONMAJ_B
    sign_r = np.sign(output["R"]) # Sign of R
    output["BETA"] = np.where(sign_r == -1, -output["BETA"], output["BETA"])

    # Delete SNPs with mismatched alleles
    nrow = output.shape[0]
    output = output.dropna(subset=["BETA"])
    if output.shape[0] < nrow:
        print(
            f"Deleted {nrow-output.shape[0]} base SNPs that did not have matching alleles in reference data."
        )
    print(f"Found proxies for {output['proxy'].sum()} SNPs.")

    # Replace the original SNPs with their proxy (if proxied)
    # As said above, we use as effect allele the minor allele in the reference panel
    output["SNP"] = np.where(output["proxy"], output["SNP_B"], output["SNP"])
    output["EA"] = np.where(output["proxy"], output["NONMAJ_B"], output["EA"])
    if "POS" in output.columns:
        output["POS"] = np.where(output["proxy"], output["BP_B"], output["POS"])
    if "CHR" in output.columns:
        output["CHR"] = np.where(output["proxy"], output["CHR_B"], output["CHR"])
    if "NEA" in output.columns:
        output["NEA"] = np.where(output["proxy"], output["MAJ_B"], output["NEA"])
    if "EAF" in output.columns:
        output["EAF"] = np.where(output["proxy"], output["NONMAJ_FREQ_B"], output["EAF"])

    # Drop ld columns
    output.drop(columns=ld.columns, inplace=True)

    return output


def find_proxies(
    snp_list,
    searchspace=None,
    reference_panel="EUR_37",
    kb=5000,
    r2=0.8,
    window_snps=1000000,
    threads=1,
    name=None
):
    """
    Given a list of SNPs, return a table of proxies using PLINK 2.0.

    Args:
        snp_list (list): List of rsids.
        searchspace (list, optional): List of SNPs to include in the search. By default, includes the whole reference panel.
        reference_panel (str, optional): The reference population to get linkage disequilibrium values and find proxies.
            Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
            Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
            Default is "EUR_37".
        kb (int, optional): Width of the genomic window to look for proxies. Defaults to 5000.
        r2 (float, optional): Minimum linkage disequilibrium value with the main SNP for a proxy to be included. Defaults to 0.8.
        window_snps (int, optional): Compute the LD value for SNPs that are not more than x SNPs apart from the main SNP. Defaults to 1000000 (equivalent to infinity).
        threads (int, optional): Number of threads to use. Defaults to 1.

    Returns:
        DataFrame: A DataFrame containing the proxies. Only biallelic SNPs are returned.
    """
    # Ensure tmp_GENAL directory exists
    os.makedirs(f"tmp_GENAL/", exist_ok=True)

    # Generate a default name if none is provided
    if name is None:
        name = str(uuid.uuid4())[:8]

    # Convert snp_list to numpy array
    snp_list = np.array(list(snp_list))

    # Check if searchspace is provided
    if searchspace is None:
        extract_arg = ""
    else:
        print("Searching proxies in the provided searchspace.")
        with open(f"tmp_GENAL/{name}_searchspace.txt", "w") as file:
            for s in searchspace + snp_list:
                file.write(str(s) + "\n")
        extract_arg = "--extract tmp_GENAL/{name}_searchspace.txt"

    # Save snp_list to a file
    np.savetxt(f"tmp_GENAL/{name}_snps_to_proxy.txt", snp_list, fmt="%s", delimiter=" ")

    # Get reference panel path and type
    ref_path, filetype = get_reference_panel_path(reference_panel)

    # Construct base command based on filetype
    base_cmd = f"{get_plink_path()}"
    if filetype == "bed":
        base_cmd += f" --bfile {ref_path}"
    else:  # pgen
        base_cmd += f" --pfile {ref_path}"

    # Construct base command based on filetype
    base_cmd = f"{get_plink_path()}"
    if filetype == "bed":
        base_cmd += f" --bfile {ref_path}"
    else:  # pgen
        base_cmd += f" --pfile {ref_path}"

    # Construct and execute the plink2 command
    command = (
        f"{base_cmd} {extract_arg} "
        f"--r-unphased 'cols=chrom,pos,id,maj,nonmaj,freq' "
        f"--ld-snp-list tmp_GENAL/{name}_snps_to_proxy.txt "
        f"--ld-window-kb {kb} "
        f"--ld-window-r2 {r2} "
        f"--ld-window {window_snps} "
        f"--threads {threads} "
        f"--out tmp_GENAL/{name}_proxy.targets"
    )
    
    run_plink_command(command)

    # Read log file to return amount of SNPs to be proxied present in the ref panel
    log_path = os.path.join("tmp_GENAL", f"{name}_proxy.targets.log")
    log_content = open(log_path).read()
    match = re.search(r'(\d+) variant[s] remaining', log_content)
    if match:
        n_present = int(match.group(1))
        if n_present == 0:
            print("None of the SNPs to be proxied are present in the reference panel.")
            return None
        else:
            print(f"{n_present} SNPs to be proxied are present in the reference panel.")

    # Read and process the output
    try:
        ld = pd.read_csv(f"tmp_GENAL/{name}_proxy.targets.vcor", sep="\s+")
    except FileNotFoundError:
        print("No proxies found that meet the specified criteria.")
        return None
      
    # Rename columns to match the expected format
    ld.rename(columns={
        'ID_A': 'SNP_A',
        'ID_B': 'SNP_B',
        '#CHROM_A': 'CHR_A',
        'CHROM_B': 'CHR_B',
        'POS_A': 'BP_A',
        'POS_B': 'BP_B',
        'UNPHASED_R': 'R',
    }, inplace=True)

    # Create PHASE column for compatibility
    #ld['PHASE'] = ld['A1'] + ld['B1'] + ld['A2'] + ld['B2']

    # Filter out multiallelic SNPs
    #ld = ld[ld["PHASE"].str.len() == 4]
    #ld = ld.reset_index(drop=True)

    # Convert integer columns to Int64 type
    for int_col in ["CHR_A", "CHR_B", "BP_A", "BP_B"]:
        ld[int_col] = ld[int_col].astype("Int64")

    return ld
