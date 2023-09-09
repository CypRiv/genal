import pandas as pd
import numpy as np
import datetime 
import os, subprocess
import scipy.stats as st
from pandas.api.types import is_numeric_dtype

from .proxy import *
from .MR import *
from .MRpresso import *


def mrpresso_func(data, action, eaf_threshold, n_iterations, outlier_test, distortion_test, significance_p, cpus):
    """
    Function for the .MRpresso GENO method.
    """
    ## Check that action argument is a correct input
    if (action not in [1,2,3]):
        raise ValueError("The action argument only takes 1,2 or 3 as value")
        
    df_exposure = data[0]
    df_outcome = data[1]
    name_outcome = data[2] 
    ## Check EAF columns if action = 2
    if action == 2:
        if "EAF" not in df_exposure.columns: 
            print("Warning: action = 2 but EAF column is missing from exposure data: palindromic SNPs will be deleted (action set to 3).")
            action = 3
        elif "EAF" not in df_outcome.columns: 
            print("Warning: action = 2 but EAF column is missing from outcome data: palindromic SNPs will be deleted (action set to 3).")
            action = 3
            
    #Harmonize exposure and outcome data
    df_mr = harmonize_MR(df_exposure, df_outcome, action = action, eaf_threshold = eaf_threshold)

    # Re-check that there are no NAs
    na_rows = df_mr[df_mr[["BETA_e","SE_e","BETA_o","SE_o"]].isna().any(axis=1)]
    if len(na_rows) > 0:
        print(f"Deleting {len(na_rows)} SNPs with NA values in exposure or outcome BETA/SE columns: {na_rows['SNP']}")
        df_mr = df_mr[["BETA_e","SE_e","BETA_o","SE_o"]].dropna(inplace=True)
    else:
        df_mr = df_mr[["BETA_e","SE_e","BETA_o","SE_o"]]

    return mr_presso(df_mr, ["BETA_e"], n_iterations, outlier_test, distortion_test, significance_p, cpus)

def MR_func(data, methods, action , heterogeneity, eaf_threshold, nboot, penk, name_exposure, cpus):
    """
    Function for the .MR GENO method.
    """
    ## Check that action argument is a correct input
    if (action not in [1,2,3]):
        raise ValueError("The action argument only takes 1,2 or 3 as value")

    ## Check the methods argument
    valid_methods = ["IVW","IVW-FE","UWR","WM","WM-pen","Simple-median","Sign","Egger","Egger-boot"]
    if not all(m in valid_methods for m in methods):
        raise ValueError(f"The list of methods can only contain strings in {valid_methods}")
        
    df_exposure = data[0]
    df_outcome = data[1]
    name_outcome = data[2]
    ## Check EAF columns if action = 2
    if action == 2:
        if "EAF" not in df_exposure.columns: 
            print("Warning: action = 2 but EAF column is missing from exposure data: palindromic SNPs will be deleted (action set to 3).")
            action = 3
        elif "EAF" not in df_outcome.columns: 
            print("Warning: action = 2 but EAF column is missing from outcome data: palindromic SNPs will be deleted (action set to 3).")
            action = 3

    print(f"Running Mendelian Randomization with {name_exposure} as exposure and {name_outcome} as outcome.")
    #Harmonize exposure and outcome data
    df_mr = harmonize_MR(df_exposure, df_outcome, action = action, eaf_threshold = eaf_threshold)

    # Re-check that there are no NAs
    na_rows = df_mr[df_mr[["BETA_e","SE_e","BETA_o","SE_o"]].isna().any(axis=1)]
    if len(na_rows) > 0:
        print(f"Deleting {len(na_rows)} SNPs with NA values in exposure or outcome BETA/SE columns: {na_rows['SNP']}")
        df_mr = df_mr[["BETA_e","SE_e","BETA_o","SE_o"]].dropna(inplace=True)
    else:
        df_mr = df_mr[["BETA_e","SE_e","BETA_o","SE_o"]]

    # Prepare values for MR methods
    BETA_e = df_mr["BETA_e"]
    BETA_o = df_mr["BETA_o"]
    SE_e = df_mr["SE_e"]
    SE_o = df_mr["SE_o"]

    ## Mapping the methods passed as argument to the corresponding functions and freeze arguments
    function_map = {"Egger": partial(mr_egger_regression, BETA_e, SE_e, BETA_o, SE_o),
                    "Egger-boot": partial(mr_egger_regression_bootstrap, BETA_e, SE_e, BETA_o, SE_o, nboot, cpus),
                    "WM": partial(mr_weighted_median, BETA_e, SE_e, BETA_o, SE_o, nboot), 
                    "WM-pen": partial(mr_pen_wm, BETA_e, SE_e, BETA_o, SE_o, nboot, penk),
                    "Simple-median": partial(mr_simple_median, BETA_e, SE_e, BETA_o, SE_o, nboot),
                    "IVW": partial(mr_ivw, BETA_e, SE_e, BETA_o, SE_o), 
                    "IVW-RE": partial(mr_ivw_re, BETA_e, SE_e, BETA_o, SE_o), 
                    "IVW-FE": partial(mr_ivw_fe, BETA_e, SE_e, BETA_o, SE_o), 
                    "UWR": partial(mr_uwr, BETA_e, SE_e, BETA_o, SE_o), 
                    "Sign": partial(mr_sign, BETA_e, BETA_o)
                   }

    results = []
    for method in methods:
        func = function_map.get(method, None)
        result = func()
        results.extend(result)
    res = pd.DataFrame(results)
    res["exposure"] = name_exposure
    res["outcome"] = name_outcome
    if not heterogeneity:
        res = res[["exposure", "outcome", "method", "nSNP", "b", "se", "pval"]]
    else:
        res = res[["exposure", "outcome", "method", "nSNP", "b", "se", "pval", "Q", "Q_df", "Q_pval"]]
        res["Q_df"] = res["Q_df"].astype("Int64")

    return res

def query_outcome_func(data, outcome, name, proxy, reference_panel, kb, r2, window_snps, cpus):
    """
    Function for the query_outcome GENO method.
    """
    ## Check required columns in the exposure data
    for column in ["SNP","BETA","SE","EA","NEA"]:
        if not(column in data.columns):
            raise ValueError(f"The column {column} is not found in the data and is necessary.")

    ## Load the outcome dataframe (to be queried)
    if str(type(outcome))=="<class 'genal.GENO.GENO'>": ## Check the type
        df_outcome=outcome.data
        if name == "": name = outcome.name
        print(f"Outcome data successfully loaded from '{outcome.name}' geno object.")
    elif type(outcome)!=str: ## Check the path provided
        raise ValueError("You need to provide either a GENO object or filepath  string to the outcome variable.")
    elif not os.path.isfile(outcome):
        raise ValueError("The path you provided doesn't lead to a file.")
    elif not (outcome.endswith(".h5") or outcome.endswith(".hdf5")):
        raise ValueError("The file provided needs to be in .h5 or .hdf5 format. You can load GWAS summary stats in a GENO object and call the .save() attribute to create one.")
    else:
        df_outcome=pd.read_hdf(outcome,key="data")
        if name == "": name = os.path.splitext(os.path.basename(outcome))[0]
        print(f"Outcome data successfully loaded from path provided.")

    ## Check necessary columns from outcome
    for column in ["SNP","BETA","SE","EA","NEA"]:
        if not(column in df_outcome.columns):
            raise ValueError(f"The column {column} is not found in the outcome data and is necessary.")

    ## Get the list of SNPs from the outcome data
    print("Identifying the exposure SNPs present in the outcome data...")
    outcome_snps=set(df_outcome.SNP.values)
    n_snps=len(outcome_snps)
    ## identify the SNPs in the exposure data that are present and find proxies for the others
    exposure_snps = set(data.SNP.values)
    snps_present = exposure_snps & outcome_snps
    print(f"{len(snps_present)} SNPs out of {len(exposure_snps)} are present in the outcome data.")
    if proxy and (len(exposure_snps) - len(snps_present) > 0):
        snps_absent = exposure_snps - snps_present
        print(f"Searching proxies for {len(snps_absent)} SNPs...")
        ld = find_proxies(snps_absent, reference_panel=reference_panel, kb=kb,r2=r2, window_snps=window_snps, threads=cpus) #Find proxies for absent SNPs
        outcome = query_outcome_proxy(df_outcome, ld, snps_present, outcome_snps) #Query GWAS with proxying
        exposure = data[data.SNP.isin(outcome.SNP)] #Build final exposure dataframe
    else:
        exposure = data[data.SNP.isin(snps_present)]
        outcome = df_outcome[df_outcome.SNP.isin(snps_present)]
    exposure.reset_index(drop=True, inplace=True)
    outcome.reset_index(drop=True, inplace=True)
    
    print(f"(Exposure data, Outcome data, Outcome name) stored in the .outcome attribute.")
    
    return exposure, outcome, name

def harmonize_MR(df_exposure, df_outcome, action=2, eaf_threshold=0.42):
    """
    Harmonize exposure and outcome for MR analyses.
    Expects data as dataframes with GENO names: "SNP","BETA","SE","EA","NEA" and "EAF" if action=2
    action=2: Determines how to treat palindromes in the harmonizing step between exposure and outcome data. 
        =1: Doesn't attempt to flip them (= Assume all alleles are coded on the forward strand)
        =2: Use allele frequencies (EAF) to attempt to flip them (conservative, default)
        =3: Remove all palindromic SNPs (very conservative).
        =4: Keep all palindromes, determine if outcome GWAS is coded forward or backward and switch all of them accordingly.
    eaf_threshold=0.42: Maximal effect allele frequency accepted when attempting to flip palindromic SNPs (only relevant if action=2)
    """
    def flip_alleles(x):
        x = x.str.upper()
        x = x.replace("C", "g").replace("G", "c").replace("A", "t").replace("T", "a")
        x = x.str.upper()        
        return x
    
    #Rename columns
    df_exposure=df_exposure.rename(columns = {"EA":"EA_e","NEA":"NEA_e","EAF":"EAF_e","BETA":"BETA_e","SE":"SE_e"}, errors="ignore")
    df_outcome=df_outcome.rename(columns={"EA":"EA_o","NEA":"NEA_o","EAF":"EAF_o","BETA":"BETA_o","SE":"SE_o"},errors="ignore")
    df_outcome=df_outcome[df_outcome.columns.intersection(["SNP","EA_o","NEA_o","EAF_o","BETA_o","SE_o"])]
    df=df_exposure.merge(df_outcome,how="left",on="SNP") #Merge
    #Create EAF columns if they do not exist. Does not change the results (they need to be present for action=2 and they are ignored in the other cases).
    df['EAF_e'] = df.get('EAF_e', 0.5)
    df['EAF_o'] = df.get('EAF_o', 0.5) 
    
    # Identify palindromes
    condition1 = ((df['EA_e'] == 'A') & (df['NEA_e'] == 'T')) | ((df['EA_e'] == 'T') & (df['NEA_e'] == 'A'))
    condition2 = ((df['EA_e'] == 'C') & (df['NEA_e'] == 'G')) | ((df['EA_e'] == 'G') & (df['NEA_e'] == 'C'))
    df['palindrome'] = condition1 | condition2
    
    # Align effect alleles between exposure and outcome
    df["aligned"]=np.where((df.EA_e == df.EA_o) & (df.NEA_e == df.NEA_o), True, False) #Already aligned
    df["inverted"]=np.where((df.EA_e == df.NEA_o) & (df.NEA_e == df.EA_o), True, False) #Inverted
    df["to_flip"]=np.where((df.palindrome == False) & (df.aligned == False) & (df.inverted == False), True, False) #Neither aligned nor inverted nor palindromic
    # Flip the SNPs to be flipped
    if df.to_flip.sum() > 0:
        to_flip_idx = df[df["to_flip"]].index #Get indices of SNPs to be flipped
        df.loc[to_flip_idx, 'EA_o'] = flip_alleles(df.loc[to_flip_idx, 'EA_o'])
        df.loc[to_flip_idx, 'NEA_o'] = flip_alleles(df.loc[to_flip_idx, 'NEA_o'])
    df["inverted"]=np.where((df.EA_e == df.NEA_o) & (df.NEA_e == df.EA_o), True, False) #Recheck Inverted
    # Align the inverted SNPs
    if df.inverted.sum() > 0: 
        inverted_idx = df[df["inverted"]].index #Get indices of inverted SNPs
        df.loc[inverted_idx, ['EA_o', 'NEA_o']] = df.loc[inverted_idx, ['NEA_o', 'EA_o']].values #Swap outcome EA and NEA values
        df.loc[inverted_idx, 'BETA_o'] *= -1 #Invert outcome BETA
        df.loc[inverted_idx, 'EAF_o'] = 1 - df.loc[inverted_idx, 'EAF_o'] #Invert outcome EAF
    df["aligned"]=np.where((df.EA_e == df.EA_o) & (df.NEA_e == df.NEA_o), True, False) #Recheck aligned
    df["allele_mismatch"] = np.where(df.aligned == False, True, False) #If still not aligned: requires exclusion due to allele mismatch
    nrow = df.shape[0]
    df = df[~df["allele_mismatch"]]
    df.reset_index(inplace=True,drop=True)
    diff =  nrow - df.shape[0]
    if diff > 0: print(f"{diff} SNPs have been excluded due to a mismatch between the exposure and outcome alleles data.")
    
    #Treat palindromes if action == 2 or 3
    if action == 3: #Simply delete them
        snps_deleted = df[df.palindrome].SNP.values
        df = df[~df.palindrome]
        df.reset_index(drop=True, inplace=True)
        print(f"Action = 3: excluding {len(snps_deleted)} palindromic SNPs: {snps_deleted}")

    elif action==2: #Use EAF_e and EAF_o to align them if both EAF are below the given eaf_threshold
        # Identify ambiguous palindromes
        df["EAF_e"] = np.where(df.EAF_e.isna(), 0.5, df.EAF_e) #If EAF is nan for a SNP, it will be removed
        df["EAF_o"] = np.where(df.EAF_o.isna(), 0.5, df.EAF_o)
        minf = np.minimum(eaf_threshold,1-eaf_threshold) #Set the boundaries for intermediate frequencies
        maxf = 1 - minf
        df["ambiguous"] = (df['palindrome'] & (((minf <= df['EAF_e']) & (df['EAF_e'] <= maxf)) | ((minf <= df['EAF_o']) & (df['EAF_o'] <= maxf))))
        snps_deleted = df[df.ambiguous].SNP.values
        df = df[~df.ambiguous]
        diff = len(snps_deleted)
        if diff > 0: 
            print(f"Action = 2: {diff} SNPs excluded for being palindromic with intermediate allele frequencies: {snps_deleted}.")
        else:
            print(f"Action = 2: None of the SNPs are palindromic with intermediate allele frequency, keeping all of them.")
        # Identify palindromes that need to be flipped
        df["to_flip"] = (df['palindrome'] & ((df.EAF_e - 0.5) * (df.EAF_o - 0.5) < 0))
        if df.to_flip.sum() > 0:
            to_flip_idx = df[df["to_flip"]].index #Get indices of SNPs to be flipped
            df.loc[to_flip_idx, 'BETA_o'] *= -1 #Invert outcome BETA
            df.loc[to_flip_idx, 'EAF_o'] = 1 - df.loc[to_flip_idx, 'EAF_o'] #Invert outcome EAF
            print(f"Action = 2: {df.to_flip.sum()} palindromic SNPs have been flipped.")
        df.reset_index(drop=True, inplace=True)
    elif action==1:
        print("Action = 1: Keeping all palindromic SNPs without attempting to flip them.")

    return df











































