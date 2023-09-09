import pandas as pd
import numpy as np
import scipy.stats as st
import os, subprocess
from collections import Counter

from .tools import *

def remove_na(data):
    """
    Identify the columns containing NA values. Delete corresponding rows.
    """
    nrows = data.shape[0]
    columns_na = data.columns[data.isna().any()].tolist()
    data.dropna(inplace=True)
    n_del = nrows - data.shape[0]
    if n_del>0:
        print(f"Deleted {n_del}({n_del/nrows*100:.3f}%) rows containing NA values in columns {columns_na}. Use preprocessing = 1 to keep the rows containing NaN values.")
    return data

def check_snp_column(data):
    """
    Remove the duplicates in the SNP column.
    """
    duplicates = data.duplicated(subset=["SNP"], keep="first")
    data = data[~duplicates]
    n_del = duplicates.sum()
    if n_del > 0:
        print(f"{n_del}({n_del/data.shape[0]*100:.3f}%) duplicated SNPs have been removed. Use keep_dups=True to keep them.")  
    return data

def check_allele_column(data, allele_col, keep_multi):
    """
    Verify that the corresponding allele column is upper case strings. Set to nan if not formed with A, T, C, G letters. 
    Set to nan if multiallelic unless keep_multi=True.
    """
    nrows = data.shape[0]
    data[allele_col] = data[allele_col].astype(str).str.upper()
    atcg_condition = data[allele_col].str.contains('^[ATCG]+$', na=False)
    atcg_count = nrows - atcg_condition.sum()
    if atcg_count>0:
        data.loc[~atcg_condition, allele_col] = np.nan
        print(f"{atcg_count}({atcg_count/nrows*100:.3f}%) rows contain non A, T, C, G values in the {allele_col} column and are set to nan.")
    if not keep_multi:
        nrows = data.shape[0]
        multi_condition = (data[allele_col].str.len()>1)
        multi_count = multi_condition.sum()
        if multi_count > 0:
            data.loc[multi_condition, allele_col] = np.nan
            print(f"{multi_count}({multi_count/nrows*100:.3f}%) rows containing multiallelic values in the {allele_col} column and are set to nan. Use keep_multi=True to keep them.")
    return data

def fill_se_p(data):
    """
    If one column among P, SE is missing but the other and BETA are present, fill it.
    """
    # If SE is missing
    if (("P" in data.columns) & ("BETA" in data.columns) & ("SE" not in data.columns)):
        data["SE"] = np.where(data["P"]<1, 
                                   np.abs(data.BETA/st.norm.ppf(data.P/2)), 
                                   0)
        print("The SE (Standard Error) column has been created.")
    # If P is missing
    if (("SE" in data.columns) & ("BETA" in data.columns) & ("P" not in data.columns)):
        data["P"] = np.where(data["SE"]>0, 2*(1-st.norm.cdf(np.abs(data.BETA)/data.SE)), 1)
        print("The P (P-value) column has been created.")
    return data

def check_p_column(data):
    """
    Verify that the P column contains numeric values in the range [0,1]. Delete rows with non-appropriate values (depending on preprocessing value).
    Based on the maximal p-value, guess if the data is clumped or not.
    """
    nrows = data.shape[0]
    data["P"] = pd.to_numeric(data["P"], errors="coerce")
    data.loc[(data['P'] < 0) & (data['P'] > 1)] = np.nan
    n_missing = data["P"].isna().sum()
    if n_missing > 0:
        print(f"{n_missing}({n_missing/nrows*100:.3f}%) values in the P column have been set to nan for being missing, non numeric or out of range [0,1].")
    return data

def check_beta_column(data, effect_column, preprocessing):
    """
    If the BETA column is a column of odds ratios, log-transform it.
    If no effect_column argument is specified, guess if the BETA column are beta estimates or odds ratios.
    """
    if effect_column is None:
        if preprocessing == 0:
            return data
        median=np.median(data.BETA)
        if 0.5 < median < 1.5:
            effect_column="OR"
            print("The BETA column looks like Odds Ratios. Use effect_column='BETA' if it is a column of Beta estimates.")
        else:
            effect_column="BETA"
            print("The BETA column looks like Betas. Use effect_column='OR' if it is a column of Odds Ratios.")

    ## Log transform the effect column if appropriate
    if effect_column not in ["BETA","OR"]:
        raise ValueError("The argument effect_column accepts only 'BETA' or 'OR' as values.")
    if effect_column=="OR":
        data["BETA"]=np.log(data["BETA"])
        data.drop(columns="SE", errors="ignore", inplace=True)
        print("The BETA column has been log-transformed to obtain Beta estimates.")
    return data

def fill_ea_nea(data, reference_panel_df):
    """
    Fill in the EA and NEA columns based on reference data.
    """
    if "BETA" in data.columns:
        print("Warning: You have specified an effect (BETA) column but no effect allele (EA) column. An effect estimate is only meaningful if paired with its corresponding allele.")
    data = data.merge(reference_panel_df[["CHR","POS","A1","A2"]],on=["CHR","POS"],how="left")
    n_missing = data["A1"].isna().sum()
    data.rename(columns={"A1":"EA","A2":"NEA"}, inplace=True)
    print(f"Alleles columns created: effect (EA) and non-effect allele (NEA). {n_missing}({n_missing/data.shape[0]*100:.3f}%) values are set to nan because SNPs were not found in the reference data.")
    return data

def fill_nea(data, reference_panel_df):
    """
    Fill in the NEA column based on reference data.
    """
    data = data.merge(reference_panel_df[["CHR","POS","A1","A2"]], on=["CHR","POS"], how="left")
    conditions = [
        data["EA"] == data["A1"],
        data["EA"] == data["A2"]]
    choices = [data["A2"], data["A1"]]
    data["NEA"] = np.select(conditions, choices, default=np.nan)
    n_missing = data["NEA"].isna().sum()
    data.drop(columns=["A1","A2"], inplace=True)
    print(f"The NEA (Non Effect Allele) column has been created. {n_missing}({n_missing/data.shape[0]*100:.3f}%) values are set to nan because SNPs were not found in the reference data.")
    return data

def fill_coordinates_func(data, reference_panel_df):
    """
    Fill in the CHR/POS columns based on reference data.
    """
    if not "SNP" in data.columns:
        raise ValueError(f"The SNP column is not found in the data and is mandatory to fill coordinates!")
    data.drop(columns=["CHR", "POS"], inplace=True, errors="ignore")
    data = data.merge(reference_panel_df[["CHR","POS","SNP"]],on="SNP",how="left")
    n_missing = data["CHR"].isna().sum()
    data["CHR"] = data["CHR"].astype('Int32')
    data["POS"] = data["POS"].astype('Int32')
    print(f"The coordinates columns (CHR for chromosome and POS for position) have been created. {n_missing}({n_missing/data.shape[0]*100:.3f}%) values are set to nan because SNPs were not found in the reference data.") 
    return data

def fill_snpids_func(data, reference_panel_df):
    """
    Fill in the SNP column based on reference data.
    If not all SNPids are present in the reference panel, the original ones will be used to fill the missing values.
    If some SNPids are still missing, they will be replaced by a standard name: CHR:POS:EA 
    """
    for column in ["CHR","POS"]:
        if not(column in data.columns):
            raise ValueError(f"The column {column} is not found in the data and is mandatory to fill snpID!")
    data.rename(columns={"SNP":"SNP_original"}, inplace=True, errors="ignore")
    data = data.merge(reference_panel_df[["CHR","POS","SNP"]], on=["CHR","POS"], how="left")
    n_missing = data["SNP"].isna().sum()
    print(f"The SNP column (rsID) has been created. {n_missing}({n_missing/data.shape[0]*100:.3f}%) SNPs were not found in the reference data.")
    snp_original = False
    if "SNP_original" in data.columns:
        snp_original = True
        data["SNP"]=data["SNP"].combine_first(data["SNP_original"])
        data.drop(columns="SNP_original", inplace=True)
    if ("EA" in data.columns):
        missing_snp_condition = data["SNP"].isna()
        data.loc[missing_snp_condition, "SNP"] = (
            data.loc[missing_snp_condition, "CHR"].astype(str) + ":" + 
            data.loc[missing_snp_condition, "POS"].astype(str) + ":" + 
            data.loc[missing_snp_condition, "EA"].astype(str)
        )
    n_missing_2 = data["SNP"].isna().sum()
    missing_diff = n_missing - n_missing_2
    if missing_diff > 0:
        print (f"Of them, {missing_diff}({missing_diff/n_missing*100:.3f}%) SNP values were set to{' the original SNP value or' if snp_original else ''} CHR:POS:EA.")
    return data

def check_int_column(data, int_col):
    """
    Set the type of the int_col column to Int64.
    """
    nrows = data.shape[0]
    data[int_col] = pd.to_numeric(data[int_col], errors="coerce")
    data[int_col] = data[int_col].round(0).astype('Int32')
    n_nan = data[int_col].isna().sum()
    if n_nan > 0:
        print(f"The {int_col} column contains {n_nan}({n_nan/nrows*100:.3f}%) values set to NaN (due to being missing or non-integer).")
    return data

def adjust_column_names(data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns):
    """
    Rename columns to the standard names making sure that there are no duplicated names.
    Delete other columns if keep_columns=False, keep them if True.
    """
    rename_dict = {CHR:"CHR", POS:"POS", SNP:"SNP", EA:"EA", NEA:"NEA", BETA:"BETA", SE:"SE", P:"P", EAF:"EAF"}
    if not keep_columns:
        data = data.loc[:,data.columns.isin([CHR,POS,SNP,EA,NEA,BETA,SE,P,EAF])]
    data.rename(columns = rename_dict, inplace=True)
    #Check duplicated column names
    column_counts = Counter(data.columns)
    duplicated_columns = [col for col, count in column_counts.items() if (count > 1) and (col in rename_dict.values())]
    if duplicated_columns:
        raise ValueError(f"After adjusting the column names, the resulting dataframe has duplicated columns. Make sure your dataframe does not have a different column named as {duplicated_columns}.")
    return data

def check_arguments(df, preprocessing, reference_panel, clumped, effect_column, keep_columns, fill_snpids, fill_coordinates, keep_multi, keep_dups):
    """
    Verify that the arguments passed to the GENO initialization are valid.
    Apply the logic between the preprocessing argument and customization arguments. 
    If the keep_columns, keep_multi, keep_dups are None: set to True if preprocessing = 0,1 and False if preprocessing = 2.
    If the fill_snpids and fill_coordinates columns are None: set to False if preprocessing = 0.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df needs to be a pandas dataframe.")
    if preprocessing not in [0,1,2]:
        raise TypeError("The preprocessing argument takes value in 0, 1 or 2. 0: the dataframe is not modified; 1: missing columns are added based on reference data but no rows are deleted; 2: missing columns are added and rows with missing, duplicated or inappropriate values are deleted. Other arguments allow for more customization (fill_nipids, fill_corrdinates, keep_multi, keep_dups, keep_columns).")
    if not ((effect_column is None) or (effect_column in ["OR", "BETA"])):
        raise TyepError("The effect_column argument only takes values in [None, 'OR', 'BETA'].")
    variables = {
    "clumped": clumped,
    "keep_columns": keep_columns,
    "fill_snpids": fill_snpids,
    "fill_coordinates": fill_coordinates,
    "keep_multi": keep_multi,
    "keep_dups": keep_dups}
    for name, value in variables.items():
        if not (value is None or isinstance(value, bool)):
            raise TypeError(f"{name} only takes values in None, True, or False.")
    # Apply the preprocessing logic.
    def keeptype_column(arg):
        if arg is None:
            if preprocessing < 2:
                return True
            else:
                return False
        else:
            return arg
    def filltype_column(arg):
        if arg is None and preprocessing == 0:
            return False
        else:
            return arg
    keep_columns = keeptype_column(keep_columns)
    keep_multi = keeptype_column(keep_multi)
    keep_dups = keeptype_column(keep_dups)
    fill_snpids = filltype_column(fill_snpids)
    fill_coordinates = filltype_column(fill_coordinates)
    return keep_columns, keep_multi, keep_dups, fill_snpids, fill_coordinates


def save_data(data, name, path="", fmt="h5", sep="\t", header=True):
    if path!="":
        path_name=f"{path}/{name}.{fmt}"
    else:
        path_name=f"{name}.{fmt}"
        
    if fmt=="h5":
        df = data.copy()
        for col in df.select_dtypes(include='integer').columns:
            df[col] = df[col].astype('float64')
        df.to_hdf(path_name,mode="w",key="data") 
        
    elif fmt in ["csv", "txt"]:
        data.to_csv(path_name,sep=sep,header=header,index=False)
            
    elif fmt in ["vcf", "vcf.gz"]:
        ## to do
        return
    else:
        raise ValueError("The fmt argument takes value in (h5 (default), csv, txt, vcf, vcf.gz).")
    print (f"Data saved to {path_name}")
    return


def Combine_GENO(Gs,name="noname",clumped=False, skip_checks=False):
    #Combine a list of different GWAS objects into one
    C=pd.DataFrame()
    if clumped==True:
        for G in Gs:
            C=pd.concat([C,G.data_clumped])
    else:
        for G in Gs:
            C=pd.concat([C,G.data])
    C=C.reset_index(drop=True)
    return(GENO(C,name=name,clumped=clumped, skip_checks=skip_checks))

def delete_tmp():
    """
    Delete the tmp folder.
    """
    if os.path.isdir("tmp_GENAL"):
        subprocess.run("rm -r tmp_GENAL",shell=True,check=False)
        print ("The tmp_GENAL folder has been successfully deleted.")
    else:
        print("There is no tmp_GENAL folder to delete in the current directory.")
    return