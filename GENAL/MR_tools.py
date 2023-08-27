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
import glob
import copy

from .config import *






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











































