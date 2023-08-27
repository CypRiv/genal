import pandas as pd
import numpy as np
import warnings
import time
import datetime
import pysam
import os
import subprocess
from collections import Counter
from tqdm import tqdm
import scipy.stats as st
from pandas.api.types import is_numeric_dtype
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from GENAL import PRS
import glob
import copy
#import dask.array as da
#import dask.dataframe as dd
#import h5py
from functools import partial

from .config import *
from .proxy import *
from .MR import *
from .MR_tools import *
from .clump import *
from .lift import *

## Create folders for temporary files and GENAL_MR
if not os.path.exists("tmp_GENAL"):
    os.makedirs("tmp_GENAL")

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

class COV:
    def __init__(self,data,name="covar_file",FID="FID",IID="IID"):
        ## Check the mandatory columns
        for column in [FID,IID]:
            if not(column in data.columns):
                raise ValueError("The column {column} is not found in the data!".format(column=column))
        self.name=name
        data=data.rename(columns={FID:"FID",IID:"IID"})
        self.data=data
        
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")
        os.chdir("tmp_GENAL")
        data.to_csv(f"{name}.cov",index=False,sep=" ")
        os.chdir("..")


class GENO:
    REF_DATASETS = ["multi","eur","sas","eas","amr"]
    REF_PANELS = ["EUR","SAS","EAS","AMR"]
    
    def __init__(self, df, name="noname", CHR="CHR",POS="POS",SNP="SNP",EA="EA",NEA="NEA",BETA="BETA",SE="SE",P="P",EAF="EAF", fill_snpids=None, fill_coordinates=None, keep_na=False, keep_multi=False, keep_dups=False, reference_panel="multi",  effect_column=None, clumped=None, reference_df=None, skip_checks=False, keep_columns=False):
        """Declare a GENO object used to store and transform data from a GWAS or derived from a GWAS. It does not include information related to individuals, only related to SNPs.
        To declare a GENO object, one needs a dataframe where each line is a SNP. Column names are specified by passing the corresponding arguments: CHR for chromosome, POS for genomic position, SNP for rsID, EA for effect allele, NEA for non-effect allele, BETA for the effect column, SE for standard error, P for p-value, EAF for effect allele frequency.
        The presence of all columns is not required to declare a GENO object, but certain methods require some of them.
        If you wish to update the rsID from CHR and POS based on a 1k genome reference, use fill_snpid=True.
        If you wish to update the CHR and/or POS from rsID based on a 1k genome reference, use fill_coordinates=True.
        If you wish to keep the NA values in the data, use keep_na=True, otherwise they will be deleted.
        GENO will try to determine whether the effect column are Betas or Odds Ratios and log-transform appropriately. If you want to override this guessing, specify effect_column with "BETA" or "OR" values.
        GENO will try to determine whether the data is clumped or not. If you want to override this guessing, specify clumped argument (True or False).
        keep_multi=False to remove multiallelic SNPs. 
        keep_dups=False to delete rows with duplicated SNP ids and keep the first line.
        reference_panel="multi" The reference panel to use for SNP adjustments. "eur", "amr", "sas", "eas" according to the 1k genome classification. "multi" combines the reference panels and contains the largest number of SNPs.
        reference_df=None A dataframe to use as reference panel. Must use the same column names as GENO objects. If provided, ref_panel argument is not considered.
        skip_checks=False To skip all formatting checks if value is True.
        keep_columns=False Will delete all columns which are not included in (CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF). Can avoid inconsistencies in some methods.
        """
        #self.check_arguments() todo
        
        self.data=df.copy()

        ## Skip checks
        if not skip_checks:
            
            ## Keep only the main columns and rename them to the standard names
            self.adjust_column_names(CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns)

            ## Check that the CHR and POS columns are integers (pandas int) and not float or else
            for int_col in ["CHR", "POS"]:
                if int_col in self.data.columns:
                    self.check_int_column(int_col)

            ## If SNP missing but CHR and POS present (or fill_snpid=True): fill the SNP column based on reference data 
            if ((("CHR" in self.data.columns) & ("POS" in self.data.columns) & ("SNP" not in self.data.columns)) or fill_snpids==True) and (not fill_snpids==False):
                self.fill_snpids(self.get_reference_panel(reference_panel, reference_df))
                
            ## If CHR and/or POS columns missing and SNP present (or fill_coordinates=True): fill CHR/POS based on reference data
            if (((("CHR" not in self.data.columns) | ("POS" not in self.data.columns)) & ("SNP" in self.data.columns)) or fill_coordinates==True) and (not fill_coordinates==False):
                self.fill_coordinates(self.get_reference_panel(reference_panel, reference_df))
            
            ## If NEA column is missing but EA and CHR/POS are present: fill it based on reference data
            if (("CHR" in self.data.columns) & ("POS" in self.data.columns) & ("NEA" not in self.data.columns) & ("EA" in self.data.columns)):
                self.fill_nea(self.get_reference_panel(reference_panel, reference_df))

            ## If NEA and EA columns are missing but CHR/POS are present: fill them based on reference data
            if (("CHR" in self.data.columns) & ("POS" in self.data.columns) & ("NEA" not in self.data.columns) & ("EA" not in self.data.columns)):
                self.fill_ea_nea(self.get_reference_panel(reference_panel, reference_df))

            ## Transform the effect column to Beta estimates if necessary
            if "BETA" in self.data.columns:
                self.check_beta_column(effect_column)

            ## Make sure the P column has appropriate values.
            if "P" in self.data.columns:
                self.check_p_column()
                
            ## If one of SE or P column missing but the other and BETA are present: fill it
            self.fill_se_p()

            ## Make sure the alleles columns are upper case strings and delete non-letters rows. Also delete multiallelic SNPs unless keep_multi=True.
            for allele_col in ["EA","NEA"]:
                if allele_col in self.data.columns:
                    self.check_allele_column(allele_col, keep_multi)    
                    
            ## Check whether some SNPs are present more than once based on the SNP column. Delete them if keep_dups=False.
            if "SNP" in self.data.columns and not keep_dups:
                self.check_snp_column()
         
            ## Check the presence of the main columns and throw warnings if they are not present
            for column in ["CHR","POS","SNP","EA","NEA","BETA","SE","P"]:
                if not(column in self.data.columns):
                    print(f"Warning: the data doesn't include a {column} column. This may become an issue later on.")       
                    
            ## NA handling: if keep_na==False: return the number of rows with missing values and delete them.
            if keep_na==False:
                self.remove_na()
                
            ## Reset index
            self.data.reset_index(drop=True, inplace=True)
            
        ##Initialize attributes
        self.init_attributes(name, clumped)
        
        return

    def init_attributes(self, name, clumped):
        """
        Initialize the name and outcome attributes.
        Guess if the data is clumped or not and set the data_clumped attribute accordingly.
        Initialize the ram and cpus attributes.
        """
        self.name = name
        if name == "noname":
            print("You haven't passed a specific name to GENO. It is recommended if you are working with several GENO objects.")
        self.outcome = []
        # Guess if clumped or not
        if clumped is None:
            if self.data["P"].max()<0.1:
                print("This data looks clumped. Use clumped=False if it is not.")
                clumped=True
            else:
                print("This data looks not clumped. Use clumped=True if it is.")
                clumped=False
        ## If the data is already clumped (clumped=True): also assign the data to self.data_clumped
        if clumped==True:
            self.data_clumped=self.data
        ## Set the maximal amount of ram/cpu to be used by the methods and dask chunksize
        ram = subprocess.run(["grep MemTotal /proc/meminfo"], shell=True, capture_output=True, text=True, check=True)
        ram = [int(s) for s in ram.stdout.split() if s.isdigit()][0]
        self.ram = round((ram/1000),-3)-1000
        self.cpus = os.cpu_count()
        self.chunksize = round((self.ram*0.8*1024**2) / (1000*self.cpus) / 100000)*100000
        return

    def check_snp_column(self):
        """
        Remove the duplicates in the SNP column.
        """
        duplicates = self.data.duplicated(subset=["SNP"], keep="first")
        self.data = self.data[~duplicates]
        n_del = duplicates.sum()
        if n_del > 0:
            print(f"{n_del}({n_del/self.data.shape[0]*100:.3f}%) duplicated SNPs have been removed. Use keep_dups=True to keep them.")  
        return
        
    def check_allele_column(self, allele_col, keep_multi):
        """
        Verify that the corresponding allele column is upper case strings. Delete strings which are not formed with A, T, C, G letters. 
        Delete multiallelic SNPs unless keep_multi=True.
        """
        self.data[allele_col] = self.data[allele_col].astype(str).str.upper()
        nrows = self.data.shape[0]
        self.data = self.data[self.data[allele_col].str.contains('^[ATCG]+$')]
        n_del = nrows-self.data.shape[0]
        if n_del>0:
            print(f"Deleting {n_del}({n_del/nrows*100:.3f}%) rows containing values different than strings formed with A, T, C, G letters in column {allele_col}.")
        if not keep_multi:
            nrows = self.data.shape[0]
            self.data = self.data[self.data[allele_col].str.len()==1]
            n_del = nrows-self.data.shape[0]
            if n_del>0:
                print(f"Deleting {n_del}({n_del/nrows*100:.3f}%) rows containing multiallelic SNPs in {allele_col} column. Use keep_multi=True to keep them.")
        return
        
    def remove_na(self):
        """
        Identify the columns containing NA values. Delete corresponding rows.
        """
        nrows = self.data.shape[0]
        columns_na = self.data.columns[self.data.isna().any()].tolist()
        self.data.dropna(inplace=True)
        n_del = nrows - self.data.shape[0]
        if n_del>0:
            print(f"Deleted {n_del}({n_del/nrows*100:.3f}%) rows containing NaN values in columns {columns_na}. Use keep_na=True to keep the rows containing NaN values.")
        return
        
    def check_beta_column(self, effect_column=None):
        """
        If the BETA column is a column of odds ratios, log-transform it.
        If no effect_column argument is specified, guess if the BETA column are beta estimates or odds ratios.
        """
        if effect_column is None:
            median=np.median(self.data.BETA)
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
            self.data["BETA"]=np.log(self.data["BETA"])
            self.data.drop(columns="SE", errors="ignore")
            print("The BETA column has been log-transformed to obtain Beta estimates.")
        return
        
    def fill_se_p(self):
        """
        If one column among P, SE is missing but the other and BETA are present, fill it.
        """
        # If SE is missing
        if (("P" in self.data.columns) & ("BETA" in self.data.columns) & ("SE" not in self.data.columns)):
            self.data["SE"] = np.where(self.data["P"]<1, 
                                       np.abs(self.data.BETA/st.norm.ppf(self.data.P/2)), 
                                       0)
            print("The SE (Standard Error) column has been created.")
        # If P is missing
        if (("SE" in self.data.columns) & ("BETA" in self.data.columns) & ("P" not in self.data.columns)):
            self.data["P"] = np.where(self.data["SE"]>0, 2*(1-st.norm.cdf(np.abs(self.data.BETA)/self.data.SE)), 1)
            print("The P (P-value) column has been created.")
        return
        
    def check_p_column(self):
        """
        Verify that the P column contains numeric values in the range [0,1]. Delete rows with non-appropriate values. 
        Based on the maximal p-value, guess if the data is clumped or not.
        """
        nrows = self.data.shape[0]
        self.data["P"] = pd.to_numeric(self.data["P"], errors="coerce")
        self.data.dropna(subset=["P"], inplace=True)
        self.data = self.data[(self.data['P'] >= 0) & (self.data['P'] <= 1)]
        n_del = nrows - self.data.shape[0]
        if n_del > 0:
            print(f"Deleting {n_del}({n_del/nrows*100:.3f}%) rows due to 'P' column values being non numeric or out of range [0,1].")
        return
        
    def fill_ea_nea(self, reference_panel_df):
        """
        Fill in the EA and NEA columns based on reference data.
        """
        if "BETA" in self.data.columns:
            print("Warning: You have specified an effect (BETA) column but no effect allele (EA) column. An effect estimate is only meaningful if paired with its corresponding allele.")
        self.data = self.data.merge(reference_panel_df[["CHR","POS","A1","A2"]],on=["CHR","POS"],how="left")
        self.data.rename(columns={"A1":"EA","A2":"NEA"}, inplace=True)
        print("Added alleles columns: effect (EA) and non-effect allele (NEA).")
        return
        
    def fill_nea(self, reference_panel_df):
        """
        Fill in the NEA column based on reference data.
        """
        self.data = self.data.merge(reference_panel_df[["CHR","POS","A1","A2"]], on=["CHR","POS"], how="left")
        conditions = [
            self.data["EA"] == self.data["A1"],
            self.data["EA"] == self.data["A2"]]
        choices = [self.data["A2"], self.data["A1"]]
        self.data["NEA"] = np.select(conditions, choices, default=np.nan)
        self.data.drop(columns=["A1","A2"], inplace=True)
        print("The NEA (Non Effect Allele) column has been created.")
        return
        
    def fill_coordinates(self, reference_panel_df):
        """
        Fill in the CHR/POS columns based on reference data.
        """
        if not "SNP" in self.data.columns:
            raise ValueError(f"The SNP column is not found in the data and is mandatory to fill coordinates!")
        self.data.drop(columns=["CHR", "POS"], inplace=True, errors="ignore")
        self.data = self.data.merge(reference_panel_df[["CHR","POS","SNP"]],on="SNP",how="left")
        self.data["CHR"] = self.data["CHR"].astype('Int32')
        self.data["POS"] = self.data["POS"].astype('Int32')
        print("The coordinates columns (CHR for chromosome and POS for position) have been created.") 
        return
 
    def fill_snpids(self, reference_panel_df):
        """
        Fill in the SNP column based on reference data.
        If not all SNPids are present in the reference panel, the original ones will be used to fill the missing values.
        If some SNPids are still missing, they will be replaced by a standard name: CHR:POS:EA 
        """
        for column in ["CHR","POS"]:
            if not(column in self.data.columns):
                raise ValueError(f"The column {column} is not found in the data and is mandatory to fill snpID!")
        self.data.rename(columns={"SNP":"SNP_original"}, inplace=True, errors="ignore")
        self.data = self.data.merge(reference_panel_df[["CHR","POS","SNP"]], on=["CHR","POS"], how="left")
        if "SNP_original" in self.data.columns:
            self.data["SNP"]=self.data["SNP"].combine_first(self.data["SNP_original"])
            self.data.drop(columns="SNP_original", inplace=True)
        if ("EA" in self.data.columns):
            missing_snp_condition = self.data["SNP"].isna()
            self.data.loc[missing_snp_condition, "SNP"] = (
                self.data.loc[missing_snp_condition, "CHR"].astype(str) + ":" + 
                self.data.loc[missing_snp_condition, "POS"].astype(str) + ":" + 
                self.data.loc[missing_snp_condition, "EA"].astype(str)
            )
        print("The SNP column has been created.")
        return

    def get_reference_panel(self, reference_panel="multi", reference_df=None):
        """
        Return the reference panel dataframe. 
        First, try to use the reference_df argument if passed. Otherwise, load it with the panel specified in reference_panel.
        """
        if not hasattr(self, "reference_panel"):
            if reference_df is not None:
                print("Using reference_df passed as argument as the reference panel.")
                self.reference_panel = reference_df
            else:
                reference_panel = reference_panel.lower()
                if reference_panel not in self.REF_DATASETS:
                    raise ValueError(f"The reference_panel argument can only take values in {self.REF_DATASETS} depending on the reference panel to be used.")
                self.reference_panel = pd.read_hdf(f"{ref_3kg_path}{reference_panel}.h5", key="data")
                self.reference_panel.columns = ["CHR","SNP","POS","A1","A2"]
                print(f"Loading the {reference_panel} reference panel.")
        return self.reference_panel
    
    def check_int_column(self, int_col):
        """
        Set the type of the int_col column to Int64 and delete rows with non-integer values.
        """
        nrows = self.data.shape[0]
        self.data[int_col] = pd.to_numeric(self.data[int_col], errors="coerce")
        self.data.dropna(subset=[int_col], inplace=True)
        self.data = self.data[self.data[int_col].apply(lambda x: x == int(x))]
        self.data[int_col] = self.data[int_col].astype('Int32')
        n_del = nrows - self.data.shape[0]
        if n_del > 0:
            print(f"Deleting {n_del}({n_del/nrows*100:.3f}%) rows with non-integer values in the {int_col} column.")
        return

    def adjust_column_names(self, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns):
        """
        Rename columns to the standard names making sure that there are no duplicated names.
        Delete other columns if keep_columns=False
        """
        rename_dict = {CHR:"CHR", POS:"POS", SNP:"SNP", EA:"EA", NEA:"NEA", BETA:"BETA", SE:"SE", P:"P", EAF:"EAF"}
        if not keep_columns:
            self.data = self.data.loc[:,self.data.columns.isin([CHR,POS,SNP,EA,NEA,BETA,SE,P,EAF])]
        self.data.rename(columns = rename_dict, inplace=True)
        #Check duplicated column names
        column_counts = Counter(self.data.columns)
        duplicated_columns = [col for col, count in column_counts.items() if (count > 1) and (col in rename_dict.values())]
        if duplicated_columns:
            raise ValueError(f"The resulting dataframe has duplicated columns. Make sure your dataframe does not have a different column named as {duplicated_columns}.")
        return

    def copy(self):
        """
        Return a deep copy of the GENO instance.
        """
        Copy=copy.deepcopy(self)
        ## Nested attributes don't seem to be copied by the copy.deepcopy() function, so we have to copy them manually (There is probably an automatic way I am not aware of.)
        if hasattr(self,"phenotype") and hasattr(self.phenotype,"type"): Copy.phenotype.type=self.phenotype.type
        return Copy
        
    
    def vcf(self,name="",clumped=False,gz=False):
        """
        Save the current data to vcf and vcf.gz/vcf.gz.tbi in the current directory to be used as MR outcome. 
        Delete rows with non-ATCG values in EA or NEA columns.
        clumped=True will save the clumped_data, otherwise the original data.
        name="" to define the name of vcf file. Otherwise, the name of the GENO object is used.
        gz=False whether the vcf file is to be converted to .vcf.gz/.vcf.gz.tbi pair or not
        """
        ## Select the name for the vcf file
        name=name.replace(".vcf","")
        name=name if name!="" else self.name
        path=f"{os.getcwd()}/{name}.vcf"
        
        ## Check that the EA and NEA columns contain only appropriate values and delete otherwise
        for allele_col in ["EA","NEA"]:
            if allele_col in self.data.columns:
                if clumped:
                    nrows=self.data_clumped.shape[0]
                    self.data_clumped = self.data_clumped[self.data_clumped[allele_col].isin(["A","T","C","G"])]
                    n_deleted=nrows-self.data_clumped.shape[0]
                    if n_deleted>0:
                        print(f"Deleting {n_deleted}({n_deleted/nrows*100:.3f}%) rows due to {allele_col} containing non-ATCG values")
                else:
                    nrows=self.data.shape[0]
                    self.data=self.data[self.data[allele_col].isin(["A","T","C","G"])]
                    n_deleted=nrows-self.data.shape[0]
                    if n_deleted>0:
                        print(f"Deleting {n_deleted}({n_deleted/nrows*100:.3f}%) rows due to {allele_col} containing non-ATCG values")

        
        ## Setup pandas2ri and call the tovcf R function
        ro.pandas2ri.activate()
        r=ro.r
        r['source'](f'{genal_path}GENAL_tovcf.R')
        tovcf_r=ro.globalenv['tovcf']
        if clumped: 
            tovcf_r(self.data_clumped,path,name)
            print("Saving the clumped data to vcf format.")
        else: 
            tovcf_r(self.data,path,name)
            print("Saving the unclumped data to vcf format.")
        
        ## Check the existence of the .vcf file and convert it to .vcf.gz/.vcf.gz.tbi if specified
        if os.path.isfile(path):
            print(f"{name}.vcf file has been successfully written.")
            if gz:
                pysam.tabix_compress(f"{name}.vcf",f"{name}.vcf.gz",force=True)
                pysam.tabix_index(f"{name}.vcf.gz",preset="vcf",force=True)
                print(f"{name}.vcf.gz and {name}.vcf.gz.tbi files have been successfully written.")
        else:
            print("There was a problem in the creation of the .vcf file.")
        
     
    
    def association_test(self,covar=True,standardize=True):
        """
        Perform single-SNP association testing for the clumped_data against the phenotype set with the set_phenotype() method.
        Requires extract_ukb method to be called before.
        Replace the BETA, SE and P columns but does not modify the phenotype attribute.
        covar=True to adjust the association tests with the standard covariates: sex, age, PC1-4
        standardize=True to standardize a quantitative phenotype before association testing (this is recommended to make the results more understandable).
        """
        ##Check that a phenotype has been set with the set_phenotype function.
        if not(hasattr(self,"phenotype")):
            raise ValueError("You first need to set a phenotype with .set_phenotype(data,PHENO,PHENO_type,IID)!") 
            
        ## Check that extract_ukb has been called.
        if not os.path.isfile(f"tmp_GENAL/{self.name}_merge_allchr.bed"):
            raise ValueError("You first need to run the extract_ukb() method before performing association tests.")
            
        ## Set the phenotype in the FAM file and adjusts depending whether it's binary or continuous
        fam=pd.read_csv(f"tmp_GENAL/{self.name}_merge_allchr.fam",header=None,delimiter=" ")
        data=self.phenotype[["IID","PHENO"]].rename(columns={"IID":0}).copy()
        fam=fam.merge(data,how="left",on=0)
        fam[5]=fam.PHENO
        fam=fam.drop(axis=1,columns="PHENO")
        if self.phenotype.type=="binary":
            fam[5]=fam[5]+1
            fam[5]=fam[5].astype("Int64")
        if (self.phenotype.type=="quant") & (standardize==True):
            fam[5]=(fam[5]-fam[5].mean())/fam[5].std()
        fam[5]=fam[5].fillna(-9)
        fam.to_csv(f"tmp_GENAL/{self.name}_merge_allchr.fam",header=None,index=False,sep=" ")
        
        ## Run plink association test with the correct arguments corresponding to the type of the phenotype (binary vs continuous) and adding covar arguments if adjusting for covariates
        method= "logistic" if self.phenotype.type=="binary" else "linear"
        covar_argument=f"--covar {genal_files_path}UKB_covar_file.cov --hide-covar" if covar==True else ""
        print(f"Running single-SNP {method} regression tests {'with' if covar==True else 'without'} covariates.")
        command=f"{plink19_path} --bfile tmp_GENAL/{self.name}_merge_allchr --{method} {covar_argument} --out tmp_GENAL/{self.name}"
        subprocess.run(command,shell=True,capture_output=True,text=True,check=True)
        
        ## Read the results; log-transform OR if logistic; merge with clumped data to update BETA and P columns; flip the betas if the A1 allele from plink match the NEA and not the EA allele from the base data ; update SE column ; drop and print the number of SNPs removed due to A1 from plink matching neither EA nor NEA from the base data
        assoc=pd.read_csv(f"tmp_GENAL/{self.name}.assoc.{method}",delimiter="\s+")
        assoc["BETA"]=np.log(assoc.OR) if self.phenotype.type=="binary" else assoc.BETA
        data=self.data_clumped
        data=data.drop(axis=1,columns=["BETA",'CHR','P'],errors="ignore").merge(assoc,how="inner",on="SNP")
        data["BETA"]=np.where(data.EA==data.A1,data.BETA,np.where(data.NEA==data.A1,-data.BETA,np.nan))
        data["SE"]=np.abs(data.BETA/st.norm.ppf(data.P/2))
        data=data.drop(axis=1,columns=["A1","TEST","NMISS","OR","STAT","BP"],errors="ignore")
        nrow_previous=data.shape[0]
        data=data.dropna(subset="BETA")
        delta_nrow=nrow_previous-data.shape[0]
        if delta_nrow>0 : print (f"{delta_nrow} rows were deleted due to the plink effective allele not found in the EA or NEA columns.") 
        self.data_clumped=data
        return 
       

            
    def set_phenotype(self, data, PHENO='', PHENO_type='', IID="person_id",alternate_control=False):
        """ Set an attribute .phenotype which is a dataframe with at least IID and PHENO columns
        data is a pandas dataframe containing an the IID and PHENO column specified in the respective arguments
        PHENO_type="" The function will try to guess if the phenotype is binary or quantitative. To avoid the guessing, specify "quant" or "binary". If binary, code it as 0 for control, 1 for cases.
        The function will determine if the IID column corresponds to the genomic IDs or to the new phenotype IDs. If it's the latter, replace it with the old genomic IDs. (The aim of the phenotype is to be used with the genetic data.)
        The function assumes that FID==IID, which is the case in UKB.
        alternate_control=False assumes that for a binary trait, the controls are coded with the most frequent value. If that is not the case, use =True.
        """
      
        ## Make sure the columns passed as arguments exist in the dataframe. If there exists columns whose name matches our standard name for PHENO and IID and which are different from the ones passed as arguments, drop them. Rename IID and PHENO to our standard names.
        for column in [PHENO,IID]:
            if not(column in data.columns):
                raise ValueError("The column {column} is not found in the data and is mandatory!".format(column=column))
        if PHENO!="PHENO":    
            data=data.drop(axis=1,columns=["PHENO"],errors="ignore")
        if IID!="IID":
            data=data.drop(axis=1,columns=["IID"],errors="ignore")
        data=data.rename(columns={IID:"IID",PHENO:"PHENO"})
        
        nrow_initial=data.shape[0]
        
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
                raise ValueError(f"The {PHENO} column is not binary!")
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
                
        ## Determine if the IID column corresponds to genomic IDs or to the new phenotype IDs. If necessary, replace it with the genomic IDs.
        bridge=pd.read_csv(f"{genal_files_path}UKB_PROJECTS_BRIDGE.txt",delimiter=" ")
        Pheno_ID=set(data.IID)
        if len(Pheno_ID.intersection(bridge.IID_old))<len(Pheno_ID.intersection(bridge.IID_new)):
            bridge["IID"]=bridge.IID_new
            bridge=bridge[["IID_old","IID"]]
            data=data.merge(bridge,how="inner",on="IID")
            data=data.drop(axis=1,columns=["IID"])
            data=data.rename(columns={"IID_old":"IID"})
        
        ## Create FID column (necessary for PRSice)
        data["FID"]=data.IID
            
        ## Select only useful columns and order them correctly (to satisfy PRSice requirements)
        data=data[["FID","IID","PHENO"]]
        
        ## Compare the number of rows and print the number of deleted ones if any.
        nrow_delta=nrow_initial-data.shape[0]
        if nrow_delta>0:
            print(f"{nrow_delta} rows ({nrow_delta/nrow_initial:.3f}%) have been deleted because the IDs provided were not the genomic ones and some of them were not present in the bridge file.")
        
        ## Set the attributes
        self.phenotype=data
        self.phenotype.type=PHENO_type
        
   
    
    def MR_all_outcomes(self,action=2,cpus=8,mem_per_cpu=20000, path=genal_mr_outcomes_path, pattern=".vcf.gz", partition="day"):
        """
        Perform MR with the data_clumped as exposure against all the standard outcomes located in Cyprien/Resources/GENAL_Outcomes, or a user-defined path to a folder containing outcome files (.vcf or .vcf.gz). The process is batch scripted and the results written to a .csv file in the local directory. To check the process and analyze the results, use MR_all_outcomes_result.
        action=2: how to treat palindromes in the harmonizing step between exposure and outcome data. =1 doesn't attempt to flip them, =2 use EAF to attempt to flip them (conservative, default), =3 remove all palindromic SNPs (very conservative).
        cpus=8: number of cpu cores used for the batschscript. 
        mem_per_cpu: quantity of ram (in mega) used PER cpu core.
        path=genal_mr_outcomes_path: the path to the folder containing all the outcomes
        pattern=".vcf.gz": the extension of the outcome files (either .vcf or .vcf.gz)
        """
        ## Check the existence of the required columns
        for column in ["SNP","BETA","SE","EAF","EA","NEA","P"]:
            if not(column in self.data.columns):
                if column=="EAF":
                    print("The EAF column is not present in the data. It is not necessary for MR unless you use action=2, which is, trying to adjust palindromic SNPs.")
                else:
                    raise ValueError(f"The column {column} is not found in the data!")
                
        ## Check that action argument is a correct input
        if (action not in (1,2,3)):
            raise ValueError("The action argument only takes 1,2 or 3 as value)")
            
        ## Check that the data_clumped attribute exists
        if not(hasattr(self,"data_clumped")):
            raise ValueError("The data needs to be clumped before using the MR_all_outcomes method! \
                \n Use .clump() first, or clumped=True when declaring the GENO.")
        
        ## Write the data to be used as exposure
        self.data_clumped.to_csv(f"tmp_GENAL/{self.name}_clumped.txt",sep="\t",index=False)
        
        ## Declare the name of the job, the name of the .sh file, cancel an existing job if it has the same name and launch the script.
        bash_name=f"MR_all_{self.name}"
        bash_script_name=f"tmp_GENAL/MR_all_outcomes_{self.name}.sh"
        if hasattr(self,"MR_all_jobid") and self.MR_all_jobid in subprocess.run(["squeue","--me"], capture_output=True, text=True).stdout:
            raise ValueError("The batch script was not launched because there is another MR_all_outcomes job for this GENO object currently running. If you wish to run several MR_all_outcomes jobs in parallel, do so from different GENO object and with different names. Otherwise, call the MRS_all_outcomes_results() attribute to inspect the progression of the job.")
            return
        
        command=f"Rscript {genal_path}GENAL_multiMR_script.R --path={os.getcwd()} --name={self.name} --genal_mr_outcomes_path={path} --plink19_path={plink19_path} --ref_3kg_path={ref_3kg_path} --cpus={cpus} --pattern={pattern} --action={action} > {os.getcwd()}/tmp_GENAL/{self.name}_R.output"
        
        ## Create the bashscript and run it
        create_bashscript(partition=partition, job_name=bash_name, bash_output_name=f"tmp_GENAL/{bash_name}", ntasks=1, cpus=cpus, mem_per_cpu=mem_per_cpu, time="10:00:00", bash_filename=bash_script_name, command=command)
        output=subprocess.run(["sbatch", f"{bash_script_name}"],check=True,text=True,capture_output=True)
        self.MR_all_jobid=output.stdout.strip().replace("Submitted batch job ","")
        self.MR_all_startime=datetime.datetime.now()
        print(f"Batch script launched: jobid {self.MR_all_jobid}.")
        return
    
    
    
    def MR_all_outcomes_results(self,p=0.05,results="p"):
        """Function to check the status of the batchscript launched by MR_all_outcomes. 
        If it has ended properly, load the results and return the significant ones in a table.
        result="p" to return only results with p_value<p as specificed by the p argument, "all" for all of them
        """
        ##Check if the jobid corresponding to the MR_all_outcomes job is still running. If still running: print the time it has been running.
        if hasattr(self,"MR_all_jobid") and self.MR_all_jobid in subprocess.run(["squeue","--me"], capture_output=True, text=True).stdout:
            print(f"The job is still running. It has been running for: {datetime.datetime.now()-self.MR_all_startime}")
        ##If not: if the results file is empty --> the job was not successful. Otherwise check the number of results, return them and print the significant ones.
        else:
            result=pd.read_csv(f"tmp_GENAL/Results_MR_all_outcomes_{self.name}.csv")
            result_len=result.shape[0]
            n_outcomes=len(glob.glob(f'{genal_mr_outcomes_path}*.vcf.gz'))
            if result_len==0:
                print("The MR_all_outcomes job has ended but was not successful. Refer to the .log file to see the error.")
            else:
                if n_outcomes>result_len/5:
                    print (f"The job has ended and was partially successfull: {int(result_len/5)} results for {n_outcomes} outcomes in total.")
                elif n_outcomes<result_len/5:
                    print("Result file is longer than expected, check what happened.")
                else:
                    print (f"The job has ended and was successfull! All {n_outcomes} outcomes yielded a result.")
                if results=="p":
                    result_s=result[(result.method=="Inverse variance weighted")&(result.pval<p)]
                    result_s=result_s.drop(columns=["id.exposure","id.outcome"])
                    result_s=result_s.sort_values("pval").reset_index(drop=True)
                    return result_s 
                elif results=="all":
                    result=result.drop(columns=["id.exposure","id.outcome"]).reset_index(drop=True)
                    return result
                else:
                    print("The results argument only takes values 'p' or 'all'")
            
    def save (self,path="",fmt="h5",sep="\t",header=True,clumped=False):
        """
        Save to a .h5 file (default) at the location specified by path or current folder if not.
        path: path of the folder to save the file
        fmt: format to use. Can be either .h5 (default), .csv, .txt, .vcf, .vcf.gz
        sep: delimiter to use (only for .csv and .txt), default is \t
        header=True: to save the column names or not (only for .csv and .txt)
        clumped=False to save the full data, or the clumped data (=True)
        """
        if path!="":
            path_name=f"{path}/{self.name}.{fmt}"
        else:
            path_name=f"{self.name}.{fmt}"
            
        if clumped and not(hasattr(self,"data_clumped")):
            raise ValueError("clumped=True was used but the data is not clumped. \
                \n Use .clump() first, or clumped=True when declaring the GENO.")
        
        if fmt=="h5":
            if clumped:
                for col in self.data_clumped.select_dtypes(include='Int64').columns:
                    self.data_clumped[col] = self.data_clumped[col].astype('float64')
                self.data_clumped.to_hdf(path_name,mode="w",key="data") 
            else:
                for col in self.data.select_dtypes(include='Int64').columns:
                    self.data[col] = self.data[col].astype('float64')      
                self.data.to_hdf(path_name,mode="w",key="data")

        elif fmt in ["csv", "txt"]:
            self.data_clumped.to_csv(path_name,sep=sep,header=header,index=False) if clumped else self.data.to_csv(path_name,sep=sep,header=header,index=False)
            
        elif fmt=="vcf":
            #to do
            return
        elif fmt=="vcf.gz":
            #to do
            return
        else:
            raise ValueError("fmt takes value in (.h5 (default), .csv, .txt, .vcf, .vcf.gz)")

        return
            
    def query_outcome (self, outcome, name="", proxy=True, ancestry="EUR", kb=5000, r2=0.6, window_snps=5000):
        """ 
        Query the outcome data, with or without proxying, to setup the dataframes required for running Mendelian Randomization (MR) with the data_clumped as exposure.
        Store a triple to the outcome attribute: (exposure_data, outcome_data, name) that is ready for MR.
        name="": the outcome data name
        #reset_outcome = True: Delete the outcomes stored in the outcome attribute. Otherwise, the processed data will be added to the list of already existing outcomes.
        Below are the parameters used when proxies are searched:
        ancestry="EUR"
        outcome: can be a GENO object (from a GWAS) or a filepath of the following types: .h5 or .hdf5 (output of the .save() attribute);
        kb=5000: width of the genomic window to look for proxies
        r2=0.6: minimum linkage disequilibrium value with the main SNP for a proxy to be included 
        window_snps=5000: compute the LD value for SNPs that are not more than x SNPs apart from the main SNP        
        """   
        ## Check that data_clumped attribute exists
        if not(hasattr(self,"data_clumped")):
            raise ValueError("The data needs to be clumped before using the quer_outcome method! \
                \n Use .clump() first, or clumped=True when declaring the GENO.")
            
        ## Check required columns in the exposure data
        for column in ["SNP","BETA","SE","EA","NEA"]:
            if not(column in self.data_clumped.columns):
                raise ValueError(f"The column {column} is not found in the data!")

        ## Load the outcome dataframe (to be queried)
        if str(type(outcome))=="<class 'GENAL.GENO.GENO'>": ## Check the type
            print(f"Loading outcome data from {outcome.name}.")
            df_outcome=outcome.data
            if name == "": name = outcome.name
        elif type(outcome)!=str: ## Check the path provided
            raise ValueError("You need to provide either a GENO object or filepath  string to the outcome variable.")
        elif not os.path.isfile(outcome):
            raise ValueError("The path you provided doesn't lead to a file.")
        elif not (outcome.endswith(".h5") or outcome.endswith(".hdf5")):
            raise ValueError("The file provided needs to be in .h5 or .hdf5 format. You can load GWAS summary stats in a GENO object and call the .save() attribute to create one.")
        else:
            df_outcome=pd.read_hdf(outcome,key="data")
            if name == "": name = outcome.split("/", expand=True)[-1].split(".h",expand=True)[1]
            print(f"Outcome data successfully loaded from path provided.")

        ## Check necessary columns from outcome
        for column in ["SNP","BETA","SE","EA","NEA"]:
            if not(column in df_outcome.columns):
                raise ValueError(f"The column {column} is not found in the outcome data.")

        ## Get the list of SNPs from the outcome data
        print("Identifying the exposure SNPs present in the outcome data...")
        outcome_snps=set(df_outcome.SNP.values)
        n_snps=len(outcome_snps)
        ## identify the SNPs in the exposure data that are present and find proxies for the others
        exposure_snps = set(self.data_clumped.SNP.values)
        snps_present = exposure_snps & outcome_snps
        print(f"{len(snps_present)} SNPs out of {len(exposure_snps)} are present in the outcome data.")
        if proxy and (len(exposure_snps) - len(snps_present) > 0):
            snps_absent = exposure_snps - snps_present
            print(f"Searching proxies for {len(snps_absent)} SNPs...")
            ld = find_proxies(snps_absent, ancestry=ancestry, kb=kb,r2=r2, window_snps=window_snps, threads=self.cpus) #Find proxies for absent SNPs
            outcome = query_outcome_proxy(df_outcome, ld, snps_present, outcome_snps) #Query GWAS with proxying
            exposure = self.data_clumped[self.data_clumped.SNP.isin(outcome.SNP)] #Build final exposure dataframe
        else:
            exposure = self.data_clumped[self.data_clumped.SNP.isin(snps_present)]
            outcome = df_outcome[df_outcome.SNP.isin(snps_present)]
        exposure.reset_index(drop=True, inplace=True)
        outcome.reset_index(drop=True, inplace=True)
        
        self.outcome = (exposure, outcome, name)
        
        return
    
    def MR_python (self, methods = ["IVW","WM","Egger"], action = 2, eaf_threshold=0.42, nboot = 10000):
        """
        Run a Mendelian Randomization with the clumped_data as exposure and the outcome(s) queried with the query_outcome() method.
        methods = ["IVW"]: a list of MR methods to run. Possible choices are "IVW": inverse variance-weighted; "WM": weighted median; "Egger": egger regression
        action=2: Determines how to treat palindromes in the harmonizing step between exposure and outcome data. 
            =1: Doesn't attempt to flip them (= Assume all alleles are coded on the forward strand)
            =2: Use allele frequencies (EAF) to attempt to flip them (conservative, default)
            =3: Remove all palindromic SNPs (very conservative).
        eaf_threshold=0.42: Maximal effect allele frequency accepted when attempting to flip palindromic SNPs (only relevant if action=2)
        nboot = 10000: number of bootstrap replications (for methods WM)

        return: a table with the MR results
        %%To do: implement cases when NEA is absent from exposure and/or outcome
        """  
        ## Check that action argument is a correct input
        if (action not in [1,2,3]):
            raise ValueError("The action argument only takes 1,2 or 3 as value")
            
        ## Check the methods argument
        valid_methods = ["IVW", "WM", "Egger"]
        if not all(m in valid_methods for m in methods):
            raise ValueError(f"The list of methods can only contain strings in {valid_methods}")
            
        ## Check that query_outcome has been called
        if not hasattr(self, "outcome"):
            raise ValueError("You must first call query_outcome() before running MR.")
            
        df_exposure = self.outcome[0]
        df_outcome = self.outcome[1]
        name = self.outcome[2]
        ## Check EAF columns if action = 2
        if action == 2:
            if "EAF" not in df_exposure.columns: 
                print("Warning: action = 2 but EAF column is missing from exposure data: palindromic SNPs will be deleted (action set to 3).")
                action = 3
            elif "EAF" not in df_outcome.columns: 
                print("Warning: action = 2 but EAF column is missing from outcome data: palindromic SNPs will be deleted (action set to 3).")
                action = 3
        
        df_harmonized = harmonize_MR(df_exposure, df_outcome, action = action, eaf_threshold = eaf_threshold)
        df_mr = df_harmonized[["BETA_e","SE_e","BETA_o","SE_o"]].dropna()
                
        ## Mapping the methods passed as argument to the corresponding functions and freeze arguments
        function_map = {"IVW": partial(mr_ivw, df_mr), "WM": partial(mr_weighted_median, df_mr, nboot), "Egger": partial(mr_egger_regression, df_mr)}
        results = []
        for method in methods:
            func = function_map.get(method, None)
            result = func()
            results.extend(result)
        res = pd.DataFrame(results)
        res["exposure"] = self.name
        res["outcome"] = self.outcome[2]
        res = res[["exposure", "outcome", "method", "nSNP", "b", "se", "pval"]]
 
        return res
            
            
    
    def MR (self, path,action=2,sensitivity=False,n=10000):
        """ 
        Perform an MR with the data_clumped as exposure and the .vcf (or .vcf.gz) file or list of files provided by the path argument as the outcome.
        Save the results (Tables and plots) in a folder MR_{self.name}.
        action=2: how to treat palindromes in the harmonizing step between exposure and outcome data. =1 doesn't attempt to flip them, =2 use EAF to attempt to flip them (conservative, default), =3 remove all palindromic SNPs (very conservative).
        sensitivity=False to determine if MRPresso is run in case of a significant inverse variance weighted result. If sensitivity=True, the function returns 2 dataframes with the second one being the MRPresso results.
        n=10000 is the number of MRPresso permutations performed. This only applies if sensitivity=True
        """
        
        ## Check the existence of the required columns
        for column in ["SNP","BETA","SE","EAF","EA","NEA","P"]:
            if not(column in self.data.columns):
                if column=="EAF":
                    print("The EAF column is not present in the data. It is not necessary for MR unless you use action=2, which is, trying to adjust palindromic SNPs.")
                else:
                    raise ValueError(f"The column {column} is not found in the data!")
                    
        ## Check that action argument is a correct input
        if (action not in (1,2,3)):
            raise ValueError("The action argument only takes 1,2 or 3 as value)")
            
        ## Check that the data_clumped attribute exists
        if not(hasattr(self,"data_clumped")):
            raise ValueError("The data needs to be clumped before using the MR method! \
                \n Use .clump() first, or clumped=True when declaring the GENO.")
            
        ## Write the data to be used as exposure
        #self.data_clumped.to_csv(f"tmp_GENAL/{self.name}_clumped.txt",sep="\t",index=False)
        
        ## Transform the path into a list if it's a single path
        path=[path] if type(path)==str else path

        ## Setup pandas2ri and extract the GENAL_MR.R function
        ro.pandas2ri.activate()
        r=ro.r
        r['source'](f'{genal_path}GENAL_MR.R')
        MR_r=ro.globalenv['MR']
        
        results=pd.DataFrame()
        if sensitivity: results_sensi=pd.DataFrame()
        for path in path:
            ## Check the path provided
            if type(path)!=str:
                raise ValueError("You need to provide a string to the path variable.")
            if not os.path.isfile(path):
                raise ValueError("The path you provided doesn't lead to a file.")
            if not path.endswith(".vcf") and not path.endswith(".vcf.gz"):
                raise ValueError("The file provided needs to be in .vcf format.")

            ## call the R function
            res=MR_r(self.data_clumped,path,action,ref_3kg_path,sensitivity,n,plink19_path)
            
            ## Convert the output of the GENAL_MR.R script back to python and adjust it depending on sensitivity parameter
            if sensitivity:
                res_p=ro.conversion.rpy2py(res[0])
                ## If the MR_Presso was run:
                if res_p.loc["3","pval"]<0.05:
                    res_sensi_p=ro.conversion.rpy2py(res[1][0][0])
                    MRpresso_P=res[1][0][1][0][1][0]
                    res_sensi_p["MRpresso_P"]=MRpresso_P
                    res_sensi_p["exposure"]=self.name
                    res_sensi_p["outcome"]=res_p.loc["3","outcome"]
                    ## If the MR_Presso was significant:
                    if not np.isnan(res_sensi_p.iloc[1,4]):
                        res_sensi_p["Distortion_P"]=res[1][0][1][2][2][0]
                        res_sensi_p["N_outliers"]=len(res[1][0][1][2][0])
                        
                        ## Determine the SNP names of the outliers --> Needs to be rechecked to make sure the SNPs are correct!
                        dat=ro.conversion.rpy2py(res[2])
                        ids=res[1][0][1][2][0]
                        Outliers=dat.iloc[ids-1].SNP.values
                        
                    else:
                        res_sensi_p["Distortion_P"]=np.nan
                        res_sensi_p["N_outliers"]=np.nan
                        Outliers=None
                    results_sensi=pd.concat([results_sensi,res_sensi_p])
            else:
                res_p=ro.conversion.rpy2py(res)
                
            ## Concatenate the dataframe
            res_p=res_p.drop(columns=["id.exposure","id.outcome"])
            res_p["exposure"]=self.name
            results=pd.concat([results,res_p])
            

        if sensitivity: return (results,results_sensi,Outliers)
        else: return results
       

            
    def clump(self, kb=250, r2=0.1, p1=5e-8, p2=0.01, reference_panel="EUR"):
        """ Clump the data in .data and assign it to the .data_clumped attribute. The clumping is done with plink.
        kb sets in thousands the window for the clumping
        r2 sets the linkage disequilibrium threshold 
        p1 sets the p-value used during the clumping (the SNPs above this threshold are not considered)
        p2 sets the p-value used after the clumping to further filter the clumped SNPs (set p2<p1 to not use this
        reference_panel="EUR" The reference population to get linkage disequilibrium values. Takes values in "EUR", "SAS", "AFR", "EAS", "AMR".
        """
        #Check input 
        reference_panel = reference_panel.upper()
        if reference_panel not in self.REF_PANELS:
            raise ValueError(f"The reference_panel argument can only take values in {self.REF_PANELS} depending on the reference panel to be used.")
            
        print(f"Clumping using {reference_panel} population as reference.")
        clumped_data = clump_data(self.data, plink19_path, ref_3kg_path + reference_panel, kb, r2, p1, p2, self.name, self.ram)
        
        if clumped_data is not None:
            self.data_clumped = clumped_data
            print("The clumped data is stored in the data_clumped attribute.")
        return 

    
    def standardize(self):
        """
        Standardize the Betas and adjust the SE column accordingly.
        """
        for column in ["BETA","SE"]:
            if not(column in self.data.columns):
                raise ValueError("The column {column} is not found in the data!".format(column=column))
        self.data["BETA"]=(self.data.BETA-np.mean(self.data.BETA))/np.std(self.data.BETA)
        self.data["SE"]=np.abs(self.data.BETA/st.norm.ppf(self.data.P/2))
        print("The Beta column has been standardized and the SE column has been adjusted.")

        
    def sort_group(self,method="lowest_p"):
        """
        Ways to handle duplicate SNPs if the instance is a combination of different GENOs.
        method="lowest_p" to keep the lowest P for each SNP. 
        """
        if method=="lowest_p":
            self.data=self.data.sort_values(by=["P"])
            self.data=self.data.groupby(by=["SNP"]).first().reset_index(drop=False)
        return

            
    def extract_ukb(self,clumped=True):
        """
        Extract the list of SNPs present in .data (or in .data_clumped if clumped==True) from the UKB plink files. 
        The output is a bed/bim/fam triple called {name}_extract_allchr including the SNPs from the UKB.
        """
        
        ## Create a tmp folder if it doesn't exist. Declare the file and bash names. Write the list of SNPs.
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")
        snp_list_name=f"{self.name}_list.txt"
        bash_script_name=f"tmp_GENAL/{self.name}_extract.sh"
        bash_name=f"{self.name}_extract"
        bedlist_name=f"{self.name}_bedlist.txt"
        if clumped==False:
            self.data["SNP"].to_csv(f"tmp_GENAL/{snp_list_name}",sep=" ",index=False,header=None)
            nrow=self.data.shape[0]
        else:
            self.data_clumped["SNP"].to_csv(f"tmp_GENAL/{snp_list_name}",sep=" ",index=False,header=None)
            nrow=self.data_clumped.shape[0]
            
        ## Declare the plink command. Call the create_bashscript function and run it in parallel (1 job per chromosome).
        command="{} --bfile {}plinkfiltered_${{SLURM_ARRAY_TASK_ID}} --extract tmp_GENAL/{}_list.txt --max-alleles 2 " \
        "--make-bed --out tmp_GENAL/{}_chr${{SLURM_ARRAY_TASK_ID}}".format(plink2_path,ukb_geno_path,self.name, bash_name)
        create_bashscript(job_name=bash_name, bash_output_name=f"tmp_GENAL/{bash_name}", ntasks=1, cpus=3, mem_per_cpu=50000, 
                          time="06:00:00",bash_filename=bash_script_name, command=command)
        output=subprocess.run(["sbatch", "--array", "1-22", f"{bash_script_name}"],check=True, text=True, capture_output=True)
        jobid=output.stdout.strip().replace("Submitted batch job ","")
        print(f"Submitting the batch job extraction for each chromosome. Jobid = {jobid}")
        
        ## Check when the job is done (could be done more efficiently). Call the create_bedlist function to determine on which chromosomes lie the extracted SNPs. Merge these chromosomes files with plink. 
        while jobid in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
            time.sleep(1)
        create_bedlist(f"tmp_GENAL/{bedlist_name}", f"tmp_GENAL/{bash_name}")    
        bash_script_name=f"tmp_GENAL/{self.name}_merge.sh"
        bash_name=f"{self.name}_merge"
        command_merge = "{} --merge-list tmp_GENAL/{} --make-bed --out tmp_GENAL/{}_allchr".format(plink19_path,bedlist_name, bash_name)
        create_bashscript(job_name=bash_name, bash_output_name=f"tmp_GENAL/{bash_name}", ntasks=1, cpus=5, mem_per_cpu=50000, 
                          time="06:00:00",bash_filename=bash_script_name, command=command_merge)
        
        output=subprocess.run(["sbatch",  f"{bash_script_name}"], check=True,text=True,capture_output=True)
        jobid=output.stdout.strip().replace("Submitted batch job ","")
        print(f"Submitting the batch job to merge SNPs extracted from each chromosome. Jobid = {jobid}")
        while jobid in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
            time.sleep(1)
            
        ## Check for the most common error that can occur in merging: multiallelic variants. If the error is present, rerun the merge excluding them.
        if "with 3+ alleles present" in open(f"tmp_GENAL/{bash_name}_allchr.log").read():
            print("Multiallelic variants detected: removing them before merging.")
            snps_to_exclude=pd.read_csv(f"tmp_GENAL/{self.name}_merge_allchr-merge.missnp",header=None)
            for i in range(22):
                if os.path.isfile(f"tmp_GENAL/{self.name}_extract_chr{i+1}.bim"):
                    bim=pd.read_csv(f"tmp_GENAL/{self.name}_extract_chr{i+1}.bim",sep="\t",header=None)
                    if len(set(bim[1]).intersection(set(snps_to_exclude[0]))) > 0:
                        command_extract=f"{plink19_path} --bfile tmp_GENAL/{self.name}_extract_chr{i+1} --exclude tmp_GENAL/{self.name}_merge_allchr-merge.missnp --make-bed --out tmp_GENAL/{self.name}_extract_chr{i+1}"
                        subprocess.run(command_extract, shell=True,capture_output=True,text=True,check=False)
            output=subprocess.run(["sbatch", f"{bash_script_name}"], check=True,text=True, capture_output=True)
            jobid=output.stdout.strip().replace("Submitted batch job ","")
            print(f"Reattempting the merge after deletion of multiallelic variants. Jobid = {jobid}")
            while jobid in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
                time.sleep(1)
                  
        ## In rare cases, plink fails to identify all the multiallelic SNPs with the first pass and we have to exclude again.
            if "with 3+ alleles present" in open(f"tmp_GENAL/{bash_name}_allchr.log").read():
                print("Multiallelic variants still detected: removing them before merging.")
                snps_to_exclude=pd.read_csv(f"tmp_GENAL/{self.name}_merge_allchr-merge.missnp",header=None)
                for i in range(22):
                    if os.path.isfile(f"tmp_GENAL/{self.name}_extract_chr{i}.bim"):
                        bim=pd.read_csv(f"tmp_GENAL/{self.name}_extract_chr{i}.bim",sep="\t",header=None)
                        if len(set(bim[1]).intersection(set(snps_to_exclude[0]))) > 0:
                            command_extract=f"{plink19_path} --bfile tmp_GENAL/{self.name}_extract_chr{i} --exclude tmp_GENAL/{self.name}_merge_allchr-merge.missnp --make-bed --out tmp_GENAL/{self.name}_extract_chr{i}"
                            subprocess.run(command_extract, shell=True,capture_output=True,text=True,check=False)
                output=subprocess.run(["sbatch", f"{bash_script_name}"], check=True, text=True, capture_output=True)
                jobid=output.stdout.strip().replace("Submitted batch job ","")
                print(f"Reattempting the merge after deletion of multiallelic variants. Jobid = {jobid}")
                while jobid in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
                    time.sleep(1)
        
        ## Report the number of SNPs not found in UKB data
        delta_nrow=nrow-int(subprocess.check_output(['wc', '-l', f"tmp_GENAL/{bash_name}_allchr.bim"]).split()[0])
        if delta_nrow > 0: print(f"{delta_nrow} SNPs were not extracted because they are not present in the data.")
        return

        
    def prs(self,weighted=True,maf=None,phenotypic_ids=True,software="plink",study="UKB"):
        """Compute a PRS with PRSice in UKB data on already clumped data
        If the P column is not present in the clumped data, creates a P=1 columns
        weighted=False will put all betas to 1 to create an unweighted PRS 
        maf will threshold by minor allele frequency (only available with software="prsice"
        phenotypic_ids=True will add a column with phenotypic ID to facilitate the merge with penotype datasets
        software="plink" to choose which software to use as the results can vary, either "plink" or "prsice"
        
        study="UKB" to choose the study in which the PRS is computed. Possible choices are UKB and HRS (which does not work for prsice currently).
        """
        
        ## Verify that the data has been clumped (or at least assigned to the data_clumped attribute)
        if not(hasattr(self,"data_clumped")):
            raise ValueError("The data needs to be clumped before using the prs method! \
                \n Use .clump() first, or clumped=True when declaring the GENO.")
        
        ## Check the mandatory columns
        data_to_prs=self.data_clumped.copy()
        for column in ["SNP","EA","BETA"]:
            if not(column in data_to_prs.columns):
                raise ValueError("The column {column} is not found in the data!".format(column=column))
        
        ## Create a P column =1 if no P-value found in the data
        if "P" not in self.data_clumped.columns:
            data_to_prs["P"]=1
            print("No P-value column found in data: creating a column P = 1.")
            
        ## Set the BETAs to 1 if weighted==False
        if weighted==False:
            data_to_prs["BETA"]=1
            print("Computing an unweighted prs.")
            
           
        if software=="prsice":
            ## Write the data to_csv in the tmp folder, adjust the maf threshold, and call PRSice on it 
            command=f'Rscript {prsice_path}PRSice.R \
            --dir {prsice_path} \
            --prsice {prsice_path}PRSice_linux \
            --base tmp_GENAL/To_prs.txt --A1 EA --snp SNP --stat BETA \
            --target {ukb_geno_path}plinkfiltered_# \
            --type bed --out tmp_GENAL/prs_{self.name} --beta --fastscore \
            --missing SET_ZERO --thread 1 --memory {int(self.ram/1000)}Gb --seed 45 --no-clump --no-regress \
            --bar-levels 1 --score avg'
            if maf!=None:
                data_to_prs=data_to_prs.loc[:,data_to_prs.columns.isin(["SNP","P","EA","NEA","BETA","EAF"])]
                if "EAF" not in data_to_prs.columns:
                    print("A minor allele frequency column must be present in the data in order to use the maf threshold.")
                else:
                    command=command+f" --base-maf EAF:{maf}"
            else:
                data_to_prs=data_to_prs.loc[:,data_to_prs.columns.isin(["SNP","P","EA","NEA","BETA"])]

            data_to_prs.to_csv("tmp_GENAL/To_prs.txt",sep="\t",index=False,header=True)
            output=subprocess.run(command, shell=True,capture_output=True,text=True,check=False)

            ## Handles a common PRSice error: duplicated SNP ID which requires to extract the SNP list from the first run
            if "duplicated SNP ID detected out of" in output.stderr:
                command=command+" --extract prs_output.valid"
                print("Rerunning PRSice after excluding the duplicated SNPs")
                output=subprocess.run(command, shell=True,capture_output=True,text=True,check=True)

            ## Read the results file, change columns names so it's in a format ready for declaring an instance of the PRS class
            if os.path.isfile(f"tmp_GENAL/prs_{self.name}.all_score"):
                print("The PRS computation was successfull!")
                df_score=pd.read_csv(f"tmp_GENAL/prs_{self.name}.all_score",sep=" ")
                df_score.columns=["FID","IID","SCORE"]
                return df_score
            else:
                print("The PRS computation was not successfull.")
                return output.stdout
            
        elif software=="plink":
                    
            ## Check that ukb_extract has been called if using the UKB.
            if not os.path.isfile(f"tmp_GENAL/{self.name}_merge_allchr.bed") and study=="UKB":
                raise ValueError("You first need to run the extract_ukb() method before computing a prs in the UKB with plink.")
                
            ## Set the path to the genetic data depending on the study variable
            if study=="UKB":
                genetic_path=f"tmp_GENAL/{self.name}_merge_allchr"
            elif study=="HRS":
                genetic_path=f"{hrs_geno_path}HRS2"
            
            ## Write the data to_csv in the tmp folder and call plink on it
            data_to_prs=data_to_prs[["SNP","EA","BETA"]]
            data_to_prs.to_csv("tmp_GENAL/To_prs.txt",sep="\t",index=False,header=True)
            output=subprocess.run(f"{plink19_path} --memory {self.ram} --bfile {genetic_path} \
            --score tmp_GENAL/To_prs.txt 1 2 3 header --out tmp_GENAL/prs_{self.name}", shell=True, capture_output=True, text=True, check=True)
            
            ## Read the results file, change columns names so it's in a format ready for declaring an instance of the PRS class
            if os.path.isfile(f"tmp_GENAL/prs_{self.name}.profile"):
                print("The PRS computation was successfull!")
                df_score=pd.read_csv(f"tmp_GENAL/prs_{self.name}.profile",sep="\s+")
                df_score=df_score[["FID","IID","SCORE"]]
                return df_score
            else:
                print("The PRS computation was not successfull.")
                return output.stdout

                    
    def lift(self, clumped=True, start="hg19",end="hg38", replace=False, extraction_file=False, chain_file="", name="", liftover=False, liftover_path=""):
        """Perform a liftover from a genetic build to another.
        If the chain file required to do the liftover is not present, download it first. It is also possible to manually provide the path to the chain file.
        start="hg19": current build of the data
        end="hg38": build to be lifted to
        If clumped==True, lift only the clumped data, otherwise the main data
        replace=False: whether to change the GENO object and update the .data or .data_clumped attributes
        extraction_file==True, also print a CHR POS SNP space delimited file for extraction in All of Us (WES data)
        chain_file="": path to a local chain file to be used for the lift. If provided, the start and end arguments are not considered.
        name="": can be used to specify a filename or filepath (without extension) to save the lifted dataframe. If not provided, will be saved in the current folder as [name]_lifted.txt
        """
        #Selecting data to lift
        if clumped==False or not(hasattr(self,"data_clumped")):
            data = self.data if replace else self.data.copy()
        elif clumped==True:
            data = self.data_clumped if replace else self.data_clumped.copy()
            
        print(f"Lifting the {'' if clumped else 'un'}clumped data{' inplace' if replace else ''}. The .data{'_clumped' if clumped else ''} attribute will {'' if replace else 'not'} be modified. Use replace={'False' if replace else 'True'} to {'leave it as is' if replace else 'lift inplace'}.")
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = lift_data(data, start="hg19",end="hg38", replace=False, extraction_file=False, chain_file=chain_file, name=f"{self.name if name == '' else name}" , genal_tools_path = genal_tools_path)
        return data
                
        
def create_bashscript(job_name,bash_output_name,bash_filename,command,
                      partition="day",time="20:00:00",cpus=1,mem_per_cpu=50000,ntasks=1):
    cmd='''#!/bin/bash
#SBATCH --partition={}
#SBATCH --job-name={}
#SBATCH --output={}
#SBATCH --ntasks={}
#SBATCH --cpus-per-task={}
#SBATCH --mem-per-cpu={}
#SBATCH --time={}

module load miniconda
conda activate {}

{}'''
    
    with open(bash_filename, "w+") as script_file:
        script_file.write(cmd.format(partition,job_name,bash_output_name,ntasks,cpus,mem_per_cpu,time,environment_name,command))

def check_bfiles(filepath):
    if os.path.exists("{}.bed".format(filepath)) and os.path.exists("{}.bim".format(filepath)) and os.path.exists("{}.fam".format(filepath)):
        return True
    return False

def create_bedlist(bedlist, output_name):
    with open(bedlist, "w+") as bedlist_file:
        for i in range(1,23):
            if check_bfiles("{}_chr{}".format(output_name,i)):
                bedlist_file.write("{}_chr{}\n".format(output_name,i))
                print("Bfiles for {}_chr{} added.".format(output_name,i))
            else:
                print("Bfiles for {}_chr{} do not exist.".format(output_name,i))