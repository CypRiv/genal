import pandas as pd
import numpy as np
from pandas.api.types import is_numeric_dtype
import scipy.stats as st
import os

from .extract_prs import check_bfiles
from .tools import *


def association_test_func(data, covar_list, standardize, name, data_pheno, pheno_type):
    
    
    ## Check that extract_snps has been called.
    genetic_path = os.path.join("tmp_GENAL", f"{name}_allchr")
    if not check_bfiles(genetic_path):
        raise FileNotFoundError("You first need to run the extract_snps() method before performing association tests.")
        
    ## Set the phenotype in the FAM file and adjusts depending whether it's binary or continuous
    fam=pd.read_csv(genetic_path + ".fam",header=None,delimiter=" ")
    data_pheno_trait=data_pheno[["IID","PHENO"]].rename(columns={"IID":0}).copy()
    fam=fam.merge(data_pheno_trait,how="left",on=0)
    fam[5]=fam.PHENO
    fam=fam.drop(axis=1,columns="PHENO")
    if pheno_type=="binary":
        fam[5]=fam[5]+1
        fam[5]=fam[5].astype("Int64")
    if (pheno_type=="quant") & (standardize==True):
        print("Standardizing the phenotype to approximate a normal distribution. Use standardize = False if you do not want to standardize.")
        fam[5]=(fam[5]-fam[5].mean())/fam[5].std()
    fam[5]=fam[5].fillna(-9)
    fam.to_csv(genetic_path + ".fam",header=None,index=False,sep=" ")
    
    ## Set the covariate file
    if len(covar_list) > 0:
        for col in covar_list:
            if col not in data_pheno.columns:
                raise TypeError(f"The {col} column is not found in the .phenotype dataframe.")
        data_pheno = data_pheno[["IID", "IID"]+covar_list]
        data_pheno = data_pheno.rename(columns={data_pheno.columns[1]: "FID"})
        covar_filename = os.path.join("tmp_GENAL", f"{name}_covar.cov")
        data_pheno.to_csv(covar_filename, sep = " ", header = True, index = False)
        covar = True
    else:
        covar = False
        
    ## Run plink association test with the correct arguments corresponding to the type of the phenotype (binary vs continuous) and adding covar arguments if adjusting for covariates
    method= "logistic" if pheno_type=="binary" else "linear"
    covar_argument=f"--covar {covar_filename} --hide-covar" if covar==True else ""
    print(f"Running single-SNP {method} regression tests {f'adjusting for: {covar_list}' if covar==True else 'without covariates'}.")
    output = os.path.join("tmp_GENAL", name)
    command=f"{get_plink19_path()} --bfile {genetic_path} --{method} {covar_argument} --out {output}"
    subprocess.run(command,shell=True,capture_output=True,text=True,check=True)

    ## Read the results; log-transform OR if logistic; merge with clumped data to update BETA and P columns; flip the betas if the A1 allele from plink match the NEA and not the EA allele from the base data ; update SE column ; drop and print the number of SNPs removed due to A1 from plink matching neither EA nor NEA from the base data
    results_path = output+f".assoc."+method
    assoc=pd.read_csv(results_path,delimiter="\s+")
    assoc["BETA"]=np.log(assoc.OR) if pheno_type=="binary" else assoc.BETA
    data=data.drop(axis=1,columns=["BETA",'CHR','P'],errors="ignore").merge(assoc,how="inner",on="SNP")
    data["BETA"]=np.where(data.EA==data.A1,data.BETA,np.where(data.NEA==data.A1,-data.BETA,np.nan))
    data["SE"]=np.abs(data.BETA/st.norm.ppf(data.P/2))
    data=data.drop(axis=1,columns=["A1","TEST","NMISS","OR","STAT","BP"],errors="ignore")
    nrow_previous=data.shape[0]
    data=data.dropna(subset="BETA")
    delta_nrow=nrow_previous-data.shape[0]
    if delta_nrow>0 : print (f"{delta_nrow} rows were deleted due to a mismatch between the effect allele and the allele columns of the genetic data.") 
    return data
    


def set_phenotype_func(data, PHENO, PHENO_type, IID, alternate_control):
    """ Set a phenotype dataframe containing individual IDs and phenotype columns formatted for single-SNP association testing.
        data: pandas dataframe containing at least an individual IDs column and one phenotype column 
        IID: name of the individual IDs column in data
        PHENO: name of the phenotype column in data
        PHENO_type="" The function will try to guess if the phenotype is binary or quantitative. To avoid the guessing, specify "quant" or "binary". If binary, the function will code it as 0 for controls, 1 for cases.
        alternate_control=False. The function assumes that for a binary trait, the controls are coded with the most frequent value. If that is not the case, use = True.
    """
    ## Make sure the columns passed as arguments exist in the dataframe. If there exists columns whose name matches our standard name for PHENO and IID and which are different from the ones passed as arguments, drop them. Rename IID and PHENO to our standard names.
    for column in [PHENO,IID]:
        if column is None:
            raise ValueError(f"You need to provide a name for the {column} variable.")
        if not(column in data.columns):
            raise ValueError(f"The column {column} is not found in the data and is mandatory!")
    if PHENO!="PHENO":    
        data=data.drop(axis=1,columns=["PHENO"],errors="ignore")
    if IID!="IID":
        data=data.drop(axis=1,columns=["IID"],errors="ignore")
    data=data.rename(columns={IID:"IID",PHENO:"PHENO"})

    ## Guess the type of the phenotype
    if (PHENO_type==""):
        if len(np.unique(data.PHENO.dropna()))==2:
            PHENO_type="binary"
            print("The PHENO column looks like a binary trait. If that is not the case, use PHENO_type='quant'")
        else:
            PHENO_type="quant"
            print("The PHENO column looks like a quantitative trait. If that is not the case, use PHENO_type='binary'")

    ## Verify that the binary trait is indeed binary and code it in 0/1 by guessing the control/cases. If the guess is wrong, can be changed using alternate_control=True
    ## Verify that the quantitative trait is numeric
    if PHENO_type=="binary":
        if len(np.unique(data.PHENO.dropna())) != 2:
            raise ValueError(f"The {PHENO} column is not binary, it contains more than 2 distinct values!")
        code_control=data.PHENO.value_counts().index[0]
        code_case=data.PHENO.value_counts().index[1]
        print(f"We are assuming the value coding for controls in the PHENO column is {code_control}. If that is not the case, use alternate_control=True")
        if alternate_control==False:
            data=data.replace({"PHENO":{code_control:0,code_case:1}}) 
        else:
            data=data.replace({"PHENO":{code_control:1,code_case:0}})
    elif PHENO_type=="quant":
        if not is_numeric_dtype(data.PHENO):
            raise ValueError(f"The {PHENO} column is not numeric!")
    else:
        raise ValueError(f"The only passible values for the argument PHENO_type are 'binary' or 'quant'")

    ## Report the number of NAs in both columns.
    nrows = data.shape[0]
    n_nan_id = data.IID.isna().sum()
    n_nan_pheno = data.PHENO.isna().sum()
    if n_nan_id > 0:
        print(f"The phenotype dataframe contains {n_nan_id}({n_nan_id/nrows*100:.3f}%) NA values in the ID column. These rows will be dropped during association testing.")
    if n_nan_pheno > 0:
        print(f"The phenotype dataframe contains {n_nan_pheno}({n_nan_pheno/nrows*100:.3f}%) NA values in the PHENO column. These rows will be dropped during association testing.")
        
    print("The phenotype data is stored in the .phenotype attribute.")
    return data, PHENO_type
    
    
    
    
    
    
    
    

    