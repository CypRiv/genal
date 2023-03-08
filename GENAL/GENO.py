import pandas as pd
import numpy as np
import warnings
import time
import pysam
import os
import subprocess
from tqdm import tqdm
import scipy.stats as st
from pandas.api.types import is_numeric_dtype
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from GENAL import PRS
import glob
import copy

## Global paths
genal_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Softwares/Genal/"
genal_mr_outcomes_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Resources/GENAL_Outcomes/"
plink2_path="/gpfs/ysm/project/gf272/plink2"
plink19_path="/gpfs/ysm/project/gf272/software/plink2"
prsice_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Softwares/PRSice/"
liftover_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Softwares/LiftOver/"
ukb_geno_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Resources/UKB_geno_files/"
ref_3kg_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Resources/Ref/"

## Create folders for temporary files and GENAL_MR
if not os.path.exists("tmp_GENAL"):
    os.makedirs("tmp_GENAL")

def Combine_GENO(Gs,name="noname",clumped=False):
    #Combine a list of different GWAS objects into one
    C=pd.DataFrame()
    if clumped==True:
        for G in Gs:
            C=pd.concat([C,G.data_clumped])
    else:
        for G in Gs:
            C=pd.concat([C,G.data])
    C=C.reset_index(drop=True)
    return(GENO(C,name=name,clumped=clumped))

def delete():
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
    def __init__(self, data, name="noname", CHR="CHR",POS="POS",SNP="SNP",EA="EA",NEA="NEA",BETA="BETA",SE="SE",P="P",EAF="EAF", adjust_snpids=False, adjust_coordinates=False, keep_na=False, effect_column="", clumped=""):
        """Declare a GENO object used to store and transform data from a GWAS or derived from a GWAS. It does not include information related to individuals, only related to SNPs.
        To declare a GENO object, one needs a dataframe where each line is a SNP. Column names are specified by passing the corresponding arguments: CHR for chromosome, POS for genomic position, SNP for rsID, EA for effect allele, NEA for non-effect allele, BETA for the effect column, SE for standard error, P for p-value, EAF for effect allele frequency.
        The presence of all columns is not required to declare a GENO object, but certain methods require some of them.
        If you wish to update the rsID from CHR and POS based on a 1k genome reference, use adjust_snpid=True.
        If you wish to update the CHR and/or POS from rsID based on a 1k genome reference, use adjust_coordinates=True.
        If you wish to keep the NA values in the data, use keep_na=True, otherwise they will be deleted.
        GENO will try to determine whether the effect column are Betas or Odds Ratios and log-transform appropriately. If you want to override this guessing, specify effect_column with "BETA" or "OR" values.
        GENO will try to determine whether the data is clumped or not. If you want to override this guessing, specify clumped argument (True or False).
        """
            
        ## Keep only the main columns and rename them to our standard names
        data=data.loc[:,data.columns.isin([CHR,POS,SNP,EA,NEA,BETA,SE,P,EAF])]
        data=data.rename(columns={CHR:"CHR",POS:"POS",SNP:"SNP",EA:"EA",NEA:"NEA",BETA:"BETA",SE:"SE",P:"P",EAF:"EAF"})
        
        ## If SNP missing but CHR and POS present: use them to fill in the SNP based on reference data
        ## Perform the same if adjust_snpid==True
        ## If not all snpids are present in the ref data, it will use the original one to fill the missing values
        ## And if some snpids are still missing, they will be replaced by a standard name: CHR:POS:EA 
        if (("CHR" in data.columns) & ("POS" in data.columns) & ("SNP" not in data.columns)) or adjust_snpids==True:
            for column in ["CHR","POS"]:
                if not(column in data.columns):
                    raise ValueError(f"The column {column} is not found in the data and is mandatory to adjust rsID!")
            ref=pd.read_csv(f"{ref_3kg_path}Bim_ref_UKB_3k.txt",sep="\t",header=None, names=["CHR","SNP","F","POS","A1","A2"])
            data=data.rename(columns={"SNP":"SNP_original"})
            data=data.merge(ref[["CHR","POS","SNP"]],on=["CHR","POS"],how="left")
            if "SNP_original" in data.columns:
                data["SNP"]=data["SNP"].combine_first(data["SNP_original"])
                data=data.drop(columns="SNP_original")
            if ("EA" in data.columns):
                data["SNP"]=np.where(data.SNP.isna(),data.CHR.astype(str)+":"+data.POS.astype(str)+":"+data.EA.astype(str),
                                     data.SNP)
            print("The SNP column has been created.")
            
        ## If CHR and/or POS columns missing but SNP present: use SNP to fill them based on reference data
        ## Perform the same if adjust_coordinates==True
        if ((("CHR" not in data.columns) | ("POS" not in data.columns)) & ("SNP" in data.columns)) or adjust_coordinates==True:
            if not "SNP" in data.columns:
                raise ValueError(f"The column SNP is not found in the data and is mandatory to adjust coordinates!")
            ref=pd.read_csv(f"{ref_3kg_path}Bim_ref_UKB_3k.txt",sep="\t",header=None, names=["CHR","SNP","F","POS","A1","A2"])
            for column in ["CHR","POS"]:
                if column in data.columns:
                    data=data.drop(columns=column)
            data=data.merge(ref[["CHR","POS","SNP"]],on="SNP",how="left")
            print("The coordinates columns (CHR for chromosome and POS for position) have been created.")
            
            
        ## If NEA column is missing but EA and CHR/POS are present: fill it based on reference data
        if (("CHR" in data.columns) & ("POS" in data.columns) & ("NEA" not in data.columns) & ("EA" in data.columns)):
            if "ref" not in locals():
                ref=pd.read_csv(f"{ref_3kg_path}Bim_ref_UKB_3k.txt",
                                    sep="\t",header=None, names=["CHR","SNP","F","POS","A1","A2"])
            data=data.merge(ref[["CHR","POS","A1","A2"]],on=["CHR","POS"],how="left")
            data["NEA"]=np.where(data["EA"]==data["A1"],data["A2"],np.where(data["EA"]==data["A2"],data["A1"],np.nan))     
            data=data.drop(columns=["A1","A2"])
            print("The NEA (Non Effect Allele) column has been created.")

        ## Make sure the CHR and POS columns are integer (pandas int) and not float or else
        for int_col in ["CHR","POS"]:
            if int_col in data.columns:
                data[int_col]=data[int_col].astype('Int64')
            
        ## If SE column missing but BETA and P present: use them to fill in SE
        if (("P" in data.columns) & ("BETA" in data.columns)& ("SE" not in data.columns)):
            data["SE"]=np.abs(data.BETA/st.norm.ppf(data.P/2))
            print("The SE (Standard Error) column has been created.")
            
        ## If P column missing but BETA and SE present: use them to fill in P
        if (("SE" in data.columns) & ("BETA" in data.columns)& ("P" not in data.columns)):
            data["P"]=2*(1-st.norm.cdf(np.abs(data.BETA)/data.SE))
            print("The P (P-value) column has been created.")
                
        ## Make sure the alleles columns are upper case strings and delete non-letters rows
        for allele_col in ["EA","NEA"]:
            if allele_col in data.columns:
                data[allele_col]=data[allele_col].astype(str).str.upper()
                nrows_before=data.shape[0]
                data=data[data[allele_col].str.isalpha()].reset_index(drop=True)
                n_nonletters=nrows_before-data.shape[0]
                if n_nonletters>0:
                    print(f"Deleting {n_nonletters}({n_nonletters/nrows_before*100:.3f}%) rows containing non letters values in allele columns (EA, NEA).")
                
        ## Guess if effect column is OR or BETA if no effect_column argument is passed
        if BETA in data.columns:
            if effect_column=="":
                median=np.median(data.BETA)
                if np.abs(median-1)<np.abs(median):
                    effect_column="OR"
                    print("The BETA column looks like Odds Ratios. Use effect_column='BETA' if it is a column of Betas.")
                else:
                    effect_column="BETA"
                    print("The BETA column looks like Betas. Use effect_column='OR' if it is a column of Odds Ratios.")
        
            ## Log transform the effect column if appropriate
            if effect_column not in ["BETA","OR"]:
                raise ValueError("The argument effect_column accepts only 'BETA' or 'OR' as values!")
            if effect_column=="OR":
                data["BETA"]=np.log(data["BETA"])
                print("The BETA column has been log-transformed to obtain Betas.")
                
        ## NA handling: if keep_na==False: return the number of rows with missing values and delete them.
        if keep_na==False:
            nrows=data.shape[0]
            columns_na=data.columns[data.isna().any()].tolist()
            data=data.dropna().reset_index(drop=True)
            n_na=nrows-data.dropna().shape[0]
            if n_na>0:
                print(f"Deleting {n_na}({n_na/nrows*100:.3f}%) rows containing NaN values in columns {columns_na}. If you wish to keep the rows containing NaN values, use keep_na=True.")
        
        ## Check the presence of the main columns and throw warnings if they are not present
        for column in ["CHR","POS","SNP","EA","NEA","BETA","SE","P"]:
            if not(column in data.columns):
                print(f"Warning: the data doesn't include a {column} column. This may become an issue later on.")
        self.name=name
        if name=="noname":
            print("You haven't passed a specific name to this GENO instance. Be careful with methods creating tmp files.")
            
        ## Check whether some SNPs are present more than once, based on SNP, and throw warning if that's the case.
        if "SNP" in data.columns:
            dup_snp=data.SNP.duplicated().sum()
            if dup_snp>0:
                print(f"Warning: {dup_snp} SNPs are duplicated based on the SNP column. This can lead to errors in analyses.")
        
        
        ## Guess if the data is clumped or not if no clumped argument is passed
        if (P in data.columns) and (clumped==""):
            if data["P"].max()<0.8:
                print("This data looks clumped. Use clumped=False if it is not.")
                clumped=True
            else:
                print("This data looks not clumped. Use clumped=True if it is.")
                clumped=False
                
        ## Assign the main attribute: .data
        self.data=data
        
        ## If the data is already clumped (clumped=True): also assign the data to self.data_clumped
        if clumped==True:
            self.data_clumped=data
                
        ## Set the maximal amount of ram to be used by the methods
        ram=subprocess.run(["grep MemTotal /proc/meminfo"],shell=True,capture_output=True,text=True,check=True)
        ram=[int(s) for s in ram.stdout.split() if s.isdigit()][0]
        self.ram=round((ram/1000),-3)-1000
        

    def copy(self):
        """
        Return a deep copy of the GENO instance.
        """
        Copy=copy.deepcopy(self)
        ## Nested attributes don't seem to be copied by the copy.deepcopy() function, so we have to copy them manually (There is probably an automatic way I am not aware of.)
        if hasattr(self,"phenotype") and hasattr(self.phenotype,"type"): Copy.phenotype.type=self.phenotype.type
        return Copy
        
    
    def vcf(self,clumped=True,name="",gz=False):
        """
        Save the current data to vcf and vcf.gz/vcf.gz.tbi in the current directory to be used as MR outcome. 
        clumped=True will save the clumped_data, otherwise the original data.
        name="" to define the name of vcf file. Otherwise, the name of the GENO object is used.
        gz=False whether the vcf file is to be converted to .vcf.gz/.vcf.gz.tbi pair or not
        """
        ## Select the name for the vcf file
        name=name.replace(".vcf","")
        name=name if name!="" else self.name
        path=f"{os.getcwd()}/{name}.vcf"
        
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
        covar_argument=f"--covar {genal_path}UKB_covar_file.cov --hide-covar" if covar==True else ""
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
        bridge=pd.read_csv(f"{genal_path}UKB_PROJECTS_BRIDGE.txt",delimiter=" ")
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
        
   
    
    def MR_all_outcomes(self,action=2,cpus=8,mem_per_cpu=30000,path_outcomes=genal_mr_outcomes_path,pattern=".vcf.gz"):
        """
        Perform MR with the data_clumped as exposure against all the standard outcomes located in Cyprien/Resources/GENAL_Outcomes, or a user-defined path to a folder containing outcome files (.vcf or .vcf.gz). The process is batch scripted and the results written to a .csv file in the local directory. To check the process and analyze the results, use MR_all_outcomes_result.
        action=2: how to treat palindromes in the harmonizing step between exposure and outcome data. =1 doesn't attempt to flip them, =2 use EAF to attempt to flip them (conservative, default), =3 remove all palindromic SNPs (very conservative).
        cpus=8: number of cpu cores used for the batschscript. 
        mem_per_cpu: quantity of ram (in mega) used PER cpu core.
        path_outcomes=genal_mr_outcomes_path: the path to the folder containing all the outcomes
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
                
        ## Write the local path as well as the parameters passed as arguments to a txt file that will be used by the R script.
        param_file=f"tmp_GENAL/MR_all_outcomes.param"
        df_param=pd.DataFrame(data={"path":[os.getcwd()],"name":self.name,"action":[action],"cpus":[cpus],
                                    "genal_mr_outcomes_path":path_outcomes,"plink19_path":plink19_path,
                                   "ref_3kg_path":ref_3kg_path,"pattern":pattern})
        df_param.to_csv(param_file,index=False,sep="\t")
                                    
        ## Change the first line of the R script so it reads the correct param file.
        R_script_path=f"{genal_path}GENAL_multiMR_outcomes.R"
        with open(R_script_path) as f:
            lines = f.readlines()       
        lines[0]=f'param=read.table("{os.getcwd()}/tmp_GENAL/MR_all_outcomes.param",header=T) \n' 
        with open(R_script_path, "w") as f:
            f.writelines(lines)
        
        ## Write the data to be used as exposure
        self.data_clumped.to_csv(f"tmp_GENAL/{self.name}_clumped.txt",sep="\t",index=False)
        
        ## Declare the name of the job, the name of the .sh file, cancel an existing job if it has the same name and launch the script.
        bash_name=f"MR_all_{self.name}"
        bash_script_name=f"tmp_GENAL/MR_all_outcomes_{self.name}.sh"
        if bash_name[:18] in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
            print("Canceling a job with the same name. If it was not intentional, be careful not to call this function several times in a row. Call MRS_all_outcomes_result to inspect the progression of the job.")
            subprocess.run(f"scancel --name={bash_name}",shell=True,check=True)
        
        command=f"module load R \nRscript --vanilla {genal_path}GENAL_multiMR_outcomes.R > {os.getcwd()}/tmp_GENAL/{self.name}_R.output"
        
        ## Create the bashscript and run it
        create_bashscript(partition="pi_falcone", job_name=bash_name, bash_output_name=f"tmp_GENAL/{bash_name}", ntasks=1, cpus=cpus, mem_per_cpu=mem_per_cpu, bash_filename=bash_script_name, command=command)
        subprocess.run(["sbatch", f"{bash_script_name}"],check=True,text=True)
        return
    
    
    
    def MR_all_outcomes_results(self,p_value=0.05,results="p"):
        """Function to check the status of the batchscript launched by MR_all_outcomes. 
        If it has ended properly, load the results and return the significant ones in a table.
        result="p" to return only results with p<p_value as specificed by the p_value argument, "all" for all of them
        """
        ##Get our jobs in queue and check if the name corresponds to the prs regress job is listed.
        bash_name=f"MR_all_{self.name}"
        status=subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout
        ##If still running: print the time it has been running as well as the last line of the job status file.
        if bash_name[:18] in status:
            m=status.split(bash_name[:18]+" ")[1].split("R    ")[1].split(":")
            print(f"The job is still running. It has been running for: {m[0][-3:]+':'+m[1][:3]}")
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
                    result_s=result[(result.method=="Inverse variance weighted")&(result.pval<p_value)]
                    result_s=result_s.drop(columns=["id.exposure","id.outcome"])
                    result_s=result_s.sort_values("pval").reset_index(drop=True)
                    return result_s 
                elif results=="all":
                    result=result.drop(columns=["id.exposure","id.outcome"]).reset_index(drop=True)
                    return result
                else:
                    print("The results argument only takes values 'p' or 'all'")
            
    
    
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
            res=MR_r(self.data_clumped,path,action,ref_3kg_path,sensitivity,n)
            
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
                        
                        ## Determine the SNP names of the outliers
                        dat=ro.conversion.rpy2py(res[2])
                        ids=res[1][0][1][2][0]
                        Outliers=dat.iloc[ids-1].SNP.values
                        
                    else:
                        res_sensi_p["Distortion_P"]=np.nan
                        res_sensi_p["N_outliers"]=np.nan
                    results_sensi=pd.concat([results_sensi,res_sensi_p])
            else:
                res_p=ro.conversion.rpy2py(res)
                
            ## Concatenate the dataframe
            res_p=res_p.drop(columns=["id.exposure","id.outcome"])
            res_p["exposure"]=self.name
            results=pd.concat([results,res_p])
            

        if sensitivity: return (results,results_sensi,Outliers)
        else: return results
        
    
    def csv (self,path="",sep="\t",header=True,index=False,clumped=False):
        """
        Save to a csv file at the location specified by path or current folder if not
        Save either the original data, or the clumped data (clumped=True)
        """
        if path!="":
            path_name=f"{path}/{self.name}.txt"
        else:
            path_name=f"{self.name}.txt"
        if clumped:
            self.data_clumped.to_csv(path_name,sep=sep,header=header,index=index)
        else:
            self.data.to_csv(path_name,sep=sep,header=header,index=index)

            
    def clump(self,kb=250, r2=0.1,p1=5e-8,p2=0.01):
        """ Clump the data in .data and assign it to the .data_clumped attribute. The clumping is done with plink.
        kb sets in thousands the window for the clumping
        r2 sets the linkage disequilibrium threshold 
        p1 sets the p-value used during the clumping (the SNPs above this threshold are not considered)
        p2 sets the p-value used after the clumping to further filter the clumped SNPs (set p2<p1 to not use this
        """
        ## Check the existence of the required columns. Create the tmp file if necessary.
        for column in ["SNP","P"]:
            if not(column in self.data.columns):
                raise ValueError("The column {column} is not found in the data!".format(column=column))
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")
        
        ## Write only the necessary column to file. Call plink with the correct arguments.
        self.data[["SNP","P"]].to_csv(f"tmp_GENAL/{self.name}_to_clump.txt",index=False,sep="\t")
        output=subprocess.run([f"{plink19_path} \
        --memory {self.ram} \
        --bfile {ref_3kg_path}EUR \
        --clump tmp_GENAL/{self.name}_to_clump.txt --clump-kb {kb} --clump-r2 {r2} --clump-p1 {p1} \
        --clump-p2 {p2} --out tmp_GENAL/{self.name}"], shell=True,capture_output=True,text=True,check=True)
        
        ## Print common outcomes. If no SNP were found: exit the function. 
        if "more top variant IDs missing" in output.stderr:
            n_missing=output.stderr.split('more top variant IDs missing')[0].split('\n')[-1]
            print(f"Warning: {n_missing}top variant IDs missing")
        if "No significant --clump results." in output.stderr:
            print("No SNPs remaining after clumping.")
            return
        clumped_message=output.stdout.split("--clump: ")[1].split("\n")[0]
        print(clumped_message)
        
        ## Get the list of clumped SNPs. Take the corresponding data subset from .data and attribute it to .data_clumped and return it too.
        with open(f"tmp_GENAL/{self.name}.list", 'wb') as f:
            subprocess.run([f"awk '{{print $3}}' tmp_GENAL/{self.name}.clumped"],shell=True,text=True,check=True,stdout=f)
        plink_clumped=pd.read_csv(f"tmp_GENAL/{self.name}.list",sep=" ")
        clumped_data=self.data[self.data["SNP"].isin(plink_clumped["SNP"])]
        self.data_clumped=clumped_data.reset_index(drop=True)
        return clumped_data

    
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
        To use if the instance is a combination of different GENOs: keep the lowest P for each SNP in .data.
        """
        if method=="lowest_p":
            self.data=self.data.sort_values(by=["P"])
            self.data=self.data.groupby(by=["SNP"]).first().reset_index(drop=True)

            
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
        create_bashscript(job_name=bash_name, bash_output_name=f"tmp_GENAL/{bash_name}", ntasks=1, cpus=5, mem_per_cpu=50000, 
                          time="06:00:00",bash_filename=bash_script_name, command=command)
        subprocess.run(["sbatch", "--array", "1-22", f"{bash_script_name}"],check=True, text=True)
        
        ## Check when the job is done (could be done more efficiently). Call the create_bedlist function to determine on which chromosomes lie the extracted SNPs. Merge these chromosomes files with plink. 
        while bash_name in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
            time.sleep(1)
        create_bedlist(f"tmp_GENAL/{bedlist_name}", f"tmp_GENAL/{bash_name}")    
        bash_script_name=f"tmp_GENAL/{self.name}_merge.sh"
        bash_name=f"{self.name}_merge"
        command_merge = "{} --merge-list tmp_GENAL/{} --make-bed --out tmp_GENAL/{}_allchr".format(plink19_path,bedlist_name, bash_name)
        create_bashscript(job_name=bash_name, bash_output_name=f"tmp_GENAL/{bash_name}", ntasks=1, cpus=5, mem_per_cpu=50000, 
                          time="06:00:00",bash_filename=bash_script_name, command=command_merge)
        subprocess.run(["sbatch", "--array", "1-22", f"{bash_script_name}"],check=True,text=True)
        while bash_name in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
            time.sleep(1)
        ## Check for the most common error that can occur in merging: multiallelic variants. If the error is present, rerun the merge excluding them.
        if "variants with 3+ alleles present" in open(f"tmp_GENAL/{bash_name}_allchr.log").read():
            snps_to_exclude=pd.read_csv(f"tmp_GENAL/{self.name}_merge_allchr-merge.missnp",header=None)
            for i in range(22):
                if os.path.isfile(f"tmp_GENAL/{self.name}_extract_chr{i}.bim"):
                    bim=pd.read_csv(f"tmp_GENAL/{self.name}_extract_chr{i}.bim",sep="\t",header=None)
                    if len(set(bim[1]).intersection(set(snps_to_exclude[0]))) > 0:
                        command_extract=f"{plink19_path} --bfile tmp_GENAL/{self.name}_extract_chr{i} --exclude tmp_GENAL/{self.name}_merge_allchr-merge.missnp --make-bed --out tmp_GENAL/{self.name}_extract_chr{i}"
                        subprocess.run(command_extract, shell=True,capture_output=True,text=True,check=False)
            subprocess.run(["sbatch", "--array", "1-22", f"{bash_script_name}"],check=True,text=True)
            while bash_name in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
                time.sleep(1)
        
        ## Report the number of SNPs not found in UKB data
        delta_nrow=nrow-int(subprocess.check_output(['wc', '-l', f"tmp_GENAL/{bash_name}_allchr.bim"]).split()[0])
        if delta_nrow > 0: print(f"{delta_nrow} SNPs were not extracted because they are not present in the data.")
        return

        
    def prs(self,weighted=True,maf=None,phenotypic_ids=True,software="plink"):
        """Compute a PRS with PRSice in UKB data on already clumped data
        If the P column is not present in the clumped data, creates a P=1 columns
        weighted=False will put all betas to 1 to create an unweighted PRS 
        maf will threshold by minor allele frequency (only available with software="prsice"
        phenotypic_ids=True will add a column with phenotypic ID to facilitate the merge with penotype datasets
        software="plink" to choose which software to use as the results can vary, either "plink" or "prsice"
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
                    
            ## Check that ukb_extract has been called.
            if not os.path.isfile(f"tmp_GENAL/{self.name}_merge_allchr.bed"):
                raise ValueError("You first need to run the extract_ukb() method before computing a prs with plink.")
                
            ## Write the data to_csv in the tmp folder and call plink on it
            data_to_prs=data_to_prs[["SNP","EA","BETA"]]
            data_to_prs.to_csv("tmp_GENAL/To_prs.txt",sep="\t",index=False,header=True)
            output=subprocess.run(f"{plink19_path} --memory {self.ram} --bfile tmp_GENAL/{self.name}_merge_allchr \
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
            
        

    
    def prs_regress(self, weighted=True,clumped=False,model="add",alternate_control=False,fastscore=False,error_1_code=2,
                   step=5e-5,lower=5e-8,upper=0.5,maf=None):
        """Run PRSice in UKB data with regression against the phenotype defined with the set_phenotype function and search for the ideal p-value threshold.
        The function launches a batchscript. The result can be retrieved from prs_regress_result function once the computation is done. prs_regress_result can also be used to monitor the progress of the batchscript.
        weighted=False will put all betas to 1 to create an unweighted PRS 
        By default, perform the search starting with the unclumped data, but it is possible to do it on the clumped data with clumped==True.
        model="add" is the genetic model used for regression: add for additive, dom for dominant, rec for recessive, het for heterozygous
        alternate_control=True if the control group is smaller than the case group (unlikely)
        fastscore=True to compute only at the main thresholds (and not in between)
        maf=None will threshold by minor allele frequency before performing the search
        error_1_code=2 allows to control the number of times PRSice should be rerun excluding duplicates. In many cases that will be 2 (default), but sometimes, only 1 time is necessary (then set error_1_code=1), or not at all (=0)
        step=5e-5 is the length of the step the software takes when iterating over the P theresholds, the lower the more precise but the more computationally intensive
        lower=5e-8 and upper=0.5 indicate the bounds of the lookup over the P thresholds. A common practice is to start by taking big steps on a big interval and then rerun the function with smaller steps on a smaller interval around the most significant threshold to obtain a precise result.
        
        """
        
        ##Check that a phenotype has been set with the set_phenotype function.
        if not(hasattr(self,"phenotype")):
            raise ValueError("You first need to set a phenotype with .set_phenotype(data,PHENO,PHENO_type,IID)!") 
        
        ## Set the BETAs to 1 if weighted==False
        data_to_prs=self.data.copy()
        if weighted==False:
            data_to_prs["BETA"]=1
            print("Computing an unweighted prs.")
            
        ## Check the mandatory columns
        for column in ["SNP","P","EA","NEA","BETA"]:
            if column not in data_to_prs.columns:
                raise ValueError("The column {column} is not found in the data!".format(column=column))
                
        ## For a binary trait: code it in 1/2 to comply with PRSice requirements. 
        phenotype_to_prs=self.phenotype.copy()
        if self.phenotype.type=="binary":
            phenotype_to_prs["PHENO"]+=1
            
        ## Go the tmp folder, write the data and phenotype to_csv
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")
        data_to_prs=data_to_prs.loc[:,data_to_prs.columns.isin(["SNP","P","EA","NEA","BETA","EAF"])]
        data_to_prs.to_csv(f"tmp_GENAL/To_prs_{self.name}.txt",sep="\t",index=False,header=True)
        phenotype_to_prs.to_csv(f"tmp_GENAL/PHENO_{self.name}.txt",sep="\t",index=False,header=True)
        
        ## Call PRSice. Add arguments if binary trait, if fastscore, if maf threshold. 
        current_path=os.getcwd()+"/tmp_GENAL"
        command=f'Rscript {prsice_path}PRSice.R \
        --dir {prsice_path} \
        --prsice {prsice_path}PRSice_linux \
        --base {current_path}/To_prs_{self.name}.txt --A1 EA --A2 NEA --pvalue P --snp SNP --stat BETA \
        --target {ukb_geno_path}plinkfiltered_# \
        --type bed --out {current_path}/prs_regress_{self.name} --beta \
        --missing SET_ZERO --score avg \
        --seed 33 --no-clump --thread 2 --memory 90Gb \
        --bar-levels 5e-08,1e-05,5e-05,0.001,0.05,1 --lower {lower} --upper {upper} --interval {step} \
        --model {model} \
        --pheno {current_path}/PHENO_{self.name}.txt --pheno-col PHENO'
        
        if self.phenotype.type=="binary":
            command=command+" --binary-target T"
        else:
            command=command+" --binary-target F"
        if fastscore:
            command=command+" --fastscore"
        if maf!=None:
            if "EAF" not in data_to_prs.columns:
                print("A minor allele frequency column must be present in the data in order to use the maf threshold.")
            else:
                command=command+f" --base-maf EAF:{maf}"
                
        ## Handles a common PRSice error: duplicated SNP ID which requires to extract the SNP list from the first run
        if error_1_code>0:
            output=subprocess.run(command, shell=True,capture_output=True,text=True,check=False)
            if "duplicated SNP ID detected out of" in output.stderr:
                print("Rerunning PRSice after excluding the duplicated SNPs")
                command=command+f" --extract {current_path}/prs_regress_{self.name}.valid"
                if error_1_code>1:
                    output2=subprocess.run(command, shell=True,capture_output=True,text=True,check=False)
                    if "duplicated SNP ID detected out of" in output2.stderr:
                        print("Rerunning PRSice after excluding the duplicated SNPs (2)")
            
        ## Declare the name of the job, the name of the .sh file, cancel an existing job if it has the same name and launch the script.
        bash_name=f"prs_regress_{self.name}"
        bash_script_name=f"prs_regress_{self.name}.sh"
        
        if bash_name[:18] in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
            print("Canceling a job with the same name. If it was not intentional, be careful not to call this function several times in a row. Call prs_regress_result to inspect the progression of the job.")
            subprocess.run(f"scancel --name={bash_name}",shell=True,check=True)
            
        #If we change the specs, don't forget to change the ram and threads in the PRSice command
        os.chdir("tmp_GENAL")
        create_bashscript(partition="pi_falcone", job_name=bash_name, bash_output_name=bash_name, ntasks=1, cpus=2, mem_per_cpu=50000, bash_filename=bash_script_name, command=command)
        subprocess.run(["sbatch", f"{bash_script_name}"],check=True,text=True)
        os.chdir("..")

        return 
    
    def prs_regress_result(self):
        """Function to check the status of the batchscript launched by prs_regress. 
        If it has ended properly, load the results in a PRS format.
        """
        ##Get our jobs in queue and check if the name corresponding to the prs regress job is listed.
        bash_name=f"prs_regress_{self.name}"
        summary_name=f"tmp_GENAL/prs_regress_{self.name}.summary"
        status=subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout
        ##If the job is still in queue: print the time it has been running as well as the last line of the job status file.
        if bash_name[:18] in status:
            m=status.split(bash_name+" ")[1].split("R    ")[1].split(":")
            print(f"The job is still running. It has been running for: {m[0][-3:]+':'+m[1][:3]}")
            out=subprocess.run(f"tail -n 1 tmp_GENAL/{bash_name}",shell=True,capture_output=True,text=True,check=True)
            update_string=out.stdout.split('\n')[-1]
            print(f"Progression update: {update_string}")
        ##If not: if the summary file doesn't exist --> the job was not successful. If it exists --> job was successful, summary file is loaded and printed, result file is loaded and returned.
        else:
            if not os.path.isfile(summary_name):
                print("The prs_regress job has ended but was not successful. Refer to the .log file to see the error.")
            elif os.path.isfile(summary_name):
                print ("The job has ended and was successfull!")
                summary=pd.read_csv(summary_name,sep="\t")
                print (f"The p-value threshold selected is {summary.Threshold}, it corresponds to {summary.Num_SNP} SNPs, and it reached a regression P-value of {summary.P} against the phenotype.")
                df_score=pd.read_csv(f"prs_regress_{self.name}.best",sep="\t")
                df_score=df_score[["FID","IID","PRS"]].rename(columns={PRS:"SCORE"})
                return df_score
                    

            
    def lift(self,clumped=True,extraction_file=True):
        """Perform a liftover from build 37 to 38 (possible to add other lifts later)
        Doesn't change the attributes. 
        The full data in build 38 is saved to the file name_38
        If extraction_file==True, also print a CHR POS SNP space delimited file for extraction in All of Us (WES data)
        If clumped==True, lift only the clumped data, otherwise the main data
        """
        ##Check that the mandatory columns are present, make sure the type is right, and create the tmp folder
        for column in ["CHR","POS"]:
            if not(column in self.data.columns):
                raise ValueError("The column {column} is not found in the data!".format(column=column))
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")
        
        ## Decide whether to lift the clumped data or the base data based on argument clumped and the presence of attribute data_clumped. Write the correct data in the format needed by liftOver.
        if clumped==False or not(hasattr(self,"data_clumped")):
            print("Lifting the unclumped data.")
            data=self.data.copy()
            data["CHR"]="chr"+data.CHR.astype(str)
            data[["CHR","POS","POS"]].to_csv(f"tmp_GENAL/{self.name}.prelift",sep=" ",index=False,header=False)
        elif clumped==True:
            print("Lifting the clumped data.")
            data=self.data_clumped.copy()
            data["CHR"]="chr"+data.CHR.astype(str)
            data[["CHR","POS","POS"]].to_csv(f"tmp_GENAL/{self.name}.prelift",sep=" ",index=False,header=False)
            
        ## Call the liftOver software.
        command=f'{liftover_path}liftOver tmp_GENAL/{self.name}.prelift \
        {liftover_path}hg19ToHg38.over.chain \
        tmp_GENAL/{self.name}.postlift tmp_GENAL/unMapped'
        output=subprocess.run(command, shell=True,capture_output=True,text=True,check=True)
        
        ## Read the output, print the number of unlifted SNPs and remove them from the prelift data. 
        df_post=pd.read_csv(f"tmp_GENAL/{self.name}.postlift",sep="\t",header=None)
        unMapped = open('tmp_GENAL/unMapped', 'r')
        Lines = unMapped.readlines()
        if len(Lines)>0:
            print(f"{int(len(Lines)/2)} SNPs could not been lifted.")
        indices=list()
        for i in range(1,len(Lines),2):
            c=Lines[i].strip()
            (chrom,pos,pos)=c.split("\t")
            indices.append(str(chrom)+":"+str(pos))
        data["SNP_IDS"]=data.CHR.astype(str)+":"+data.POS.astype(str)
        data=data[~(data.SNP_IDS.isin(indices))].drop(columns="SNP_IDS").reset_index(drop=True)
        
        ## Merge prelift and postlift data. Unknown chr from the output of liftOver are assigned the value 99. SNPs mapped to unknown chr are deleted from the final data and their number printed.
        data["POS"]=df_post[1].astype(int)
        data["CHR"]=df_post[0].str.split("chr",expand=True)[1].str.split("_",expand=True)[0].replace({"X":99,"Un":99}).astype(int)
        nrow_before=data.shape[0]
        data=data[data.CHR!=99]
        nrow_diff=nrow_before-data.shape[0]
        if nrow_diff>0:
            print(f"{nrow_diff} SNPs were lifted to an unknown chromosome and deleted from the final files.")
        
        ## Save files: whole data and also the extraction file to be used in All of Us if extraction_file=True
        data.to_csv(f"{self.name}_38.txt",sep="\t",header=True,index=False)
        if extraction_file:
            if not("SNP" in data.columns):
                data["SNP"]=data.CHR.astype(str)+":"+data.POS.astype(str)
            data[["CHR","POS","SNP"]].to_csv(f"{self.name}_38_extraction.txt",sep=" ",header=False,index=False)
        return data


        
def create_bashscript(job_name,bash_output_name,bash_filename,command,
                      partition="pi_falcone,general",time="100:00:00",cpus=1,mem_per_cpu=50000,ntasks=1):
    cmd='''#!/bin/bash
#SBATCH --partition={}
#SBATCH --job-name={}
#SBATCH --output={}
#SBATCH --ntasks={}
#SBATCH --cpus-per-task={}
#SBATCH --mem-per-cpu={}
#SBATCH --time={}

{}'''
    
    with open(bash_filename, "w+") as script_file:
        script_file.write(cmd.format(partition,job_name,bash_output_name,ntasks,cpus,mem_per_cpu,time,command))

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