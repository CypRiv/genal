import pandas as pd
import numpy as np
import warnings
import time
import os
import subprocess
import scipy.stats as st
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter

## Global paths
genal_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Softwares/Genal/"
genal_mr_outcomes_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Resources/GENAL_Outcomes/"
plink2_path="/gpfs/ysm/project/gf272/plink2"
plink19_path="/gpfs/ysm/project/gf272/software/plink2"
prsice_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Softwares/PRSice/"
liftover_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Softwares/LiftOver/"
ukb_geno_path="/SAY/standard2/falconelab-CC0841-MEDNEU/UKB_Filtered/"
ref_3kg_path="/gpfs/gibbs/pi/falcone/LabMembers/Cyprien/Resources/Ref/"

class PRS:
    def __init__(self,data,name="noname",IID="IID",SCORE="SCORE",standardize=True):
        """
        Declare a PRS object used to store, transform and analyze data of a Polygenic Risk Score (PRS) in a population.
        It does not include information related to individual SNPs, only about the aggregate PRS measure.
        data is a dataframe with each line representing an individual a containing at least the IID, and SCORE columns. It can contains more columns to be used as covariates for certain methods.
        name="noname" is the name used in the different methods when creating files or printing results and plots.
        standardize=True to standardize the PRS (recommended)
        The initialization will create the most common ntiles of the PRS: tertiles, quintiles ... and store them in appropriate variables.
        This version is meant to be used in the UKB so it only handles person ids (IID) and not family ids (FID)
        """
        
        ## Check the presence of the main columns and throw errors if they are not present. 
        ## Throw warning if name is not well defined.
        for column in [IID,SCORE]:
            if not(column in data.columns):
                raise ValueError(f"The column {column} is not found in the data!")
        self.name=name
        if name=="noname":
            print("You haven't passed a specific name to this PRS instance. Be careful with methods creating tmp files.")
            
        ## Rename the main columns to our standard names
        data=data.rename(columns={IID:"IID",SCORE:"SCORE"})
        
        ## Standardize the PRS
        if standardize: 
            data["SCORE"]=(data.SCORE-np.mean(data.SCORE))/np.std(data.SCORE)
            print("The PRS has been standardized to have a mean of 0 and standard deviation of 1.")
        
        ## Create the most common ntiles of the PRS and add their names to the repartition list attribute
        data["SCORE_tert"]=pd.qcut(data.SCORE,q=3,labels=range(1,4)).astype("Int64")
        data["SCORE_quart"]=pd.qcut(data.SCORE,q=4,labels=range(1,5)).astype("Int64")
        data["SCORE_quin"]=pd.qcut(data.SCORE,q=5,labels=range(1,6)).astype("Int64")
        data["SCORE_dec"]=pd.qcut(data.SCORE,q=10,labels=range(1,11)).astype("Int64")
        data["SCORE_twent"]=pd.qcut(data.SCORE,q=20,labels=range(1,21)).astype("Int64")
        data["SCORE_2080"]=pd.qcut(data.SCORE,q=[0,0.20,0.80,1],labels=range(1,4)).astype("Int64")
        self.repartition_list=["SCORE_tert","SCORE_quart","SCORE_quin","SCORE_dec","SCORE_twent","SCORE_2080"]
        print("The following PRS repartitions have been added to the data: tertiles, quartiles, quintiles, deciles, twentiles, 0-20-80-1.")
        
        ## Assign the main attribute: .data
        self.data=data
        
    def ntile(self,array,name,add_to_list=True):
        """
        Create a new variable with a custom n-tile repartition of the SCORE.
        array is an array representing the repartition. For instance, [0,0.33,0.66,1] for tertiles. Can also be an integer, and in that case it's going to be the number of categories, evenly distributed.
        name is the name to be given to the variable. Be careful that if the name already exists in the data, the corresponding variable values will be replaced.
        add_to_list=True will add the variable name to the attribute repartition_list representing different repartitions of the SCORE. This list is used in some methods to compare the repartitions.
        """
        ## Check if the name argument already exists in the column names.
        if name in self.data.columns:
            print("Be careful: this name is already the name of a column and it will be replaced!")
        ## Add the column in the data
        self.data[name]=pd.qcut(self.data.SCORE,q=array,labels=range(1,len(array))).astype("Int64")
        ## Add the name to the list if add_to_list==True
        if add_to_list: self.repartition_list.append(name)
            
          
    def add_variables(self, data, IID="person_id"):
        """
        Add phenotype variables to the PRS data. The data passed as argument is left joined on the PRS data.
        data: dataframe with an IID column that should correspond to the genomic or the phenotypic IDs. 
        IID="IID" indicates the name of the column containing the individual IDs.
        """
         ## Make sure the IID column passed as argument exists in the dataframe. Drop an eventual column named "IID" if IID!="IID". Rename IID column.
        if not(IID in data.columns):
            raise ValueError("The column {column} is not found in the data and is mandatory!".format(column=column))
        if IID!="IID":
            data=data.drop(axis=1,columns=["IID"],errors="ignore")
        data=data.rename(columns={IID:"IID"})
        
        ## Make sure that none of the columns of the passed dataframe already exists in our data.
        i=0
        for column in data.columns:
            if column in self.data.columns and column != "IID":
                i=1
                print(f"A variable {column} already exists in the data, change its name to allow the merge.")
        if i!=0: raise ValueError("Some variable names were already present in the PRS object data. Please change their names and try again.")
        
        ## Determine if the IID column corresponds to genomic IDs or to the new phenotype IDs. If necessary, replace it with the genomic IDs. (We assume the PRS ID column corresponds to genomic IDs.)
        bridge=pd.read_csv(f"{genal_path}UKB_PROJECTS_BRIDGE.txt",delimiter=" ")
        Pheno_ID=set(data.IID)
        nrow_initial=data.shape[0]
        if len(Pheno_ID.intersection(bridge.IID_old))<len(Pheno_ID.intersection(bridge.IID_new)):
            bridge["IID"]=bridge.IID_new
            bridge=bridge[["IID_old","IID"]]
            data=data.merge(bridge,how="inner",on="IID")
            data=data.drop(axis=1,columns=["IID"])
            data=data.rename(columns={"IID_old":"IID"})
    
        ## Compare the number of rows and print the number of deleted ones if any.
        nrow_delta=nrow_initial-data.shape[0]
        if nrow_delta>0:
            print(f"{nrow_delta} rows ({nrow_delta/nrow_initial:.3f}%) have been deleted because the IDs provided were not the genomic ones and some of them were not present in the bridge file.")
            
        ## Merge with the existing data
        self.data=self.data.merge(data,on="IID",how="left")
        
            

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        