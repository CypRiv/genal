import pandas as pd
import numpy as np
import warnings
import time
import datetime
import pysam
import os
import subprocess
from tqdm import tqdm
import scipy.stats as st
from pandas.api.types import is_numeric_dtype
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
import glob
import copy

from io import StringIO

from .config import *


def query_outcome_proxy(df, ld, snps_to_extract, snps_df=[]):
    """
    Given a dataframe df (coming from GENO.data) and a dataframe of potential proxies (output from find_proxies), extract the best proxies from df as well as the SNPs in snps_to_extract. (Suited for outcome data.)
    df: dataframe of SNP information with the usual GENO columns (SNP,BETA,SE,EAF,EA,NEA) - EAF is not necessary
    ld: dataframe of proxies (output from find_proxies)
    snps_present: list of SNPs to extract in addition to the proxies
    snps_df: optional, list of SNPs to choose the proxy from. Should be the list of SNPs in df. Can be provided to avoid recomputing it.
    """
    if len(snps_df) == 0: snps_df = df.SNP.values
    ld = ld[ld.SNP_B.isin(snps_df)] ##The proxies should be present in df
    ld = ld.reindex(ld['R'].abs().sort_values(ascending=False).index) #Sort by r2
    ld = ld.groupby("SNP_A").first().reset_index(drop=False) #Select the best proxy for each SNP
    snps_to_query = set(snps_to_extract) | set(ld.SNP_B.values)
    df_queried = df[df.SNP.isin(snps_to_query)]
    output = df_queried.merge(ld, how="left", left_on="SNP", right_on="SNP_B")
    output["proxy"]=np.where(pd.notnull(output["SNP_B"]),True,False)
    ld_columns=ld.columns
    #Flip BETA if the proxied SNP alleles are switched in the reference panel
    conditions = [
        (output['EA'] == output['B2']),
        (output['EA'] == output['B1']),
        (~output['proxy']),
        ((output['EA'].isin([output['B1'], output['B2']]) == False) & (output['proxy']))
    ]
    choices = [
        -output['BETA'], # if EA == B2, flip the sign of BETA
        output['BETA'],  # if EA == B1, BETA does not change
        output['BETA'],  # if the original SNP was not proxied, BETA does not change
        np.nan       # if the original SNP was proxied but EA is neither 'B1' nor 'B2', BETA is NaN
    ]
    output['BETA'] = np.select(conditions, choices) # Applying the conditions and choices
    nrow=output.shape[0]
    output=output.dropna(subset=["BETA"]) 
    if output.shape[0]<nrow:
        print(f"Deleted {nrow-output.shape[0]} base SNPs that did not have matching alleles in reference data.")
    print(f"Found proxies for {output['proxy'].sum()} SNPs.")
    
    #Replace the original SNPs with their proxy (if proxied)
    output["SNP"]=np.where(output["proxy"],output["SNP_A"],output["SNP"])
    output["POS"]=np.where(output["proxy"],output["BP_A"],output["POS"])
    output["CHR"]=np.where(output["proxy"],output["CHR_A"],output["CHR"])
    output["EA"]=np.where(output["proxy"],output["A1"],output["EA"])
    output["NEA"]=np.where(output["proxy"],output["A2"],output["NEA"])
    if "EAF" in output.columns: output["EAF"]=np.where(output["proxy"],output["MAF_A"],output["EAF"])
    
    output=output.drop(columns=ld_columns) #Drop ld columns
    
    return output

def apply_proxies(df, ld, searchspace=np.empty(0)):
    """
    Given a dataframe (coming from GENO.data or GENO.data_clumped attributes) and a dataframe of proxies (output from find_proxies), replace the SNPs in df with their best proxies, if it exists. (Suited for exposure data.)
    searchspace=[]: list of SNPs to restrict the list of potential proxies. By default, include all the proxies found. Using a searchspace can be done either at the find_proxies step or at this step, but it is much faster to use it here.
    df: dataframe of SNP information with the usual GENO columns (SNP,BETA,SE,EAF,EA,NEA) - EAF is not necessary
    ld: dataframe of proxies (output from find_proxies)
    """
    if len(searchspace) != 0: ## Filter by searchspace
        print ("Filtering the potential proxies with the searchspace provided.")
        ld = ld[ld.SNP_B.isin(searchspace)]
    ld = ld.reindex(ld['R'].abs().sort_values(ascending=False).index) #Sort by r2
    ld = ld.groupby("SNP_A").first().reset_index(drop=False) #Select the best proxy for each SNP
    output = df.merge(ld,how="left",left_on="SNP",right_on="SNP_A") #Merge
    output["proxy"]=np.where(pd.notnull(output["SNP_B"]),True,False)
    ld_columns=ld.columns
    #Flip BETA if the original SNP alleles are switched in the reference panel
    conditions = [
        (output['EA'] == output['A2']),
        (output['EA'] == output['A1']),
        (~output['proxy']),
        ((output['EA'].isin([output['A1'], output['A2']]) == False) & (output['proxy']))
    ]
    choices = [
        -output['BETA'], # if EA == A2, flip the sign of BETA
        output['BETA'],  # if EA == A1, BETA does not change
        output['BETA'],  # if SNP_A is NaN (The original SNP was not proxied), BETA does not change
        np.nan       # if the original SNP was proxied but EA is neither 'A1' nor 'A2', BETA is NaN
    ]
    output['BETA'] = np.select(conditions, choices) # Applying the conditions and choices
    nrow=output.shape[0]
    output=output.dropna(subset=["BETA"]) 
    if output.shape[0]<nrow:
        print(f"Deleted {nrow-output.shape[0]} base SNPs that did not have matching alleles in reference data.")
    print(f"Found proxies for {output['proxy'].sum()} SNPs.")
    
    #Replace the original SNPs with their proxy (if proxied)
    output["SNP"]=np.where(output["proxy"],output["SNP_B"],output["SNP"])
    output["POS"]=np.where(output["proxy"],output["BP_B"],output["POS"])
    output["CHR"]=np.where(output["proxy"],output["CHR_B"],output["CHR"])
    output["EA"]=np.where(output["proxy"],output["B1"],output["EA"])
    output["NEA"]=np.where(output["proxy"],output["B2"],output["NEA"])
    if "EAF" in output.columns: output["EAF"]=np.where(output["proxy"],output["MAF_B"],output["EAF"])
    
    output=output.drop(columns=ld_columns) #Drop ld columns
    
    return output

def find_proxies(snp_list, searchspace=np.empty(0), ancestry="EUR", kb=5000,r2=0.6, window_snps=5000, threads=1):
    """
    Given a list of SNPs, return a table of proxies.
    snp_list: list of rsids
    searchspace=[]: list of SNPs to include in the search. By default, include the whole reference panel.
    ancestry="EUR"
    kb=5000: width of the genomic window to look for proxies
    r2=0.6: minimum linkage disequilibrium value with the main SNP for a proxy to be included 
    window_snps=5000: compute the LD value for SNPs that are not more than x SNPs apart from the main SNP
    threads=1: number of threads to use
    
    Return only biallelic SNPs
    """
    os.makedirs(f"tmp_GENAL/", exist_ok=True) ## create tmp folder 
    snp_list=np.array(list(snp_list)) #To array
    #searchspace=np.array(list(searchspace)) #To array
    if len(searchspace) == 0: ## Write the searchspace if not empty
        print("Searching proxies in the whole reference data.")
        extract_arg=""
    else:
        print("Searching proxies in the searchspace provided.")
        with open("tmp_GENAL/searchspace.txt", "w") as file:
            for s in searchspace:
                file.write(str(s) + '\n')
            for s in snp_list:
                file.write(str(s) + '\n')
        #np.savetxt("tmp_GENAL/searchspace.txt",searchspace, fmt="%s", delimiter='\n')
        extract_arg="--extract tmp_GENAL/searchspace.txt"
    np.savetxt("tmp_GENAL/snps_to_proxy.txt",snp_list, fmt="%s", delimiter=' ') ## Write the snp_list
    ## declare plink command
    command=f"{plink19_path} --bfile {ref_3kg_path}{ancestry} {extract_arg} --keep-allele-order --r in-phase with-freqs gz --ld-snp-list tmp_GENAL/snps_to_proxy.txt --ld-window-kb {kb} --ld-window-r2 {r2} --ld-window {window_snps} --out tmp_GENAL/proxy.targets --threads {threads}"
    output=subprocess.run(command,shell=True,capture_output=True,text=True,check=True) ## execute
    
    ## read output
    cmd = f'gunzip -c tmp_GENAL/proxy.targets.ld.gz'
    unzipped_content = subprocess.check_output(cmd, shell=True).decode('utf-8')
    ld = pd.read_csv(StringIO(unzipped_content), sep="\s+")
    
    # Cleaning 
    ld['PHASE'] = ld['PHASE'].str.replace("/", "")
    ld = ld[ld['PHASE'].apply(len) == 4] #Delete multiallelic SNPs
    ld = ld[ld['SNP_A'] != ld['SNP_B']] #Delete original SNPs
    ld = ld.reset_index(drop=True)
    temp = pd.DataFrame(ld['PHASE'].apply(list).to_list(), columns=['A1', 'B1', 'A2', 'B2']) #Split the 'PHASE' column into separate characters
    ld = pd.concat([ld, temp], axis=1)
    
    # Put the integer columns as Int64
    for int_col in ["CHR_A","CHR_B","BP_A","BP_B",]:
        ld[int_col]=ld[int_col].astype('Int64')

    return ld





























