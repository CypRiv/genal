import pandas as pd
import numpy as np
import warnings
import os, subprocess
import subprocess
import copy

from .proxy import *
from .MR import *
from .MR_tools import *
from .clump import *
from .lift import *
from .tools import *
from .geno_tools import *
from .association import *
from .extract_prs import *

# Add proxying function (input is df + searchspace and returns proxied df)
# Multi-MR with python MR
# Don't forget to change the init to put .genal folder and rewrite it only if not present
# Warning that users might not have shell (for the .ram attribute)


class GENO:
    
    def __init__(self, df, name="noname", CHR="CHR",POS="POS",SNP="SNP",EA="EA",NEA="NEA",BETA="BETA",SE="SE",P="P",EAF="EAF", preprocessing=1, reference_panel="eur", clumped=None, effect_column=None, keep_columns=None, keep_multi=None, keep_dups=None, fill_snpids=None, fill_coordinates=None):
        """Declare a GENO object used to store and transform Single Nucleotide Polymorphisms (SNP) data. It is intended to handle data from a GWAS or derived from a GWAS with information such as SNP name, position on the genome, SNP-trait effects and corresponding effect allele, population allele frequency.
        df: a pandas dataframe where each line is a SNP. Column names are specified by passing the corresponding arguments. The presence of all columns is not required to declare a GENO object, but certain methods require some of them. 
        CHR = "CHR": name of the chromosome column
        POS = "POS": name of the genomic position column
        SNP = "SNP": name of the rsID column 
        EA = "EA": name of the effect allele column
        NEA = "NEA": name of the non-effect allele column
        BETA = "BETA": name of the effect estimate column
        SE = "SE": name of the effect standard error column
        P = "P": name of the p-value column
        EAF = "EAF": name of the effect allele frequency column. 
        preprocessing = 1. Standard level of preprocessing to apply to all processing steps. 0: the dataframe is not modified; 1: missing columns are added based on reference data and invalid values set to nan but no rows are deleted; 2: missing columns are added and rows with missing, duplicated or invalid values are deleted. Other arguments allow for more customization (fill_nipids, fill_coordinates, keep_multi, keep_dups, keep_columns).
        reference_panel = "eur". The reference panel to use for SNP adjustments. Takes value in "eur", "amr", "sas", "eas", "multi" following the 1k genome classification. "multi" combines the reference panels from the four ancestries and contains the largest number of SNPs. Can also be: a dataframe with columns ["CHR","SNP","POS","A1","A2"] to use as reference panel. Or a path to a .bim file to use as reference. 
        clumped = None. genal will try to determine whether the data is clumped or not. If you want to override this guess, specify if True (data is clumped) or False (data is not clumped).
        effect_column = None. genal will try to determine whether the effect column are Betas or Odds Ratios and log-transform appropriately. If you want to override this guessing, specify "BETA" or "OR".
        keep_columns = None. False: All columns which are not included in (CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF) will be deleted (Can avoid inconsistencies in some methods). True: These columns will be kept. None: defers to the preprocessing value (0,1: True; 2: False).
        keep_multi = None: False: multiallelic SNPs will be removed. True: These rows will be kept. None: defers to the preprocessing value (0,1: True; 2: False).
        keep_dups = None. False: rows with duplicated SNP ids will be deleted and only the first line conserved. True: These rows will be kept. None: defers to the preprocessing value (0,1: True; 2: False).
        fill_snpids = None. True: The SNP (rsID) column will be created or replaced based on CHR/POS columns and a reference genome. False: Do not attempt to update the rsID column. None: Create the SNP column if absent and CHR/POS columns are present and preprocessing = 1,2.
        fill_coordinates = None. True: The CHR and/or POS will be created or replaced based on SNP column and a reference genome. False: Do not attempt to update the CHR/POS columns. None: Create the CHR/POS columns if absent and SNP column is present and preprocessing = 1,2.
        """
        #Check arguments and apply preprocessing logic.
        keep_columns, keep_multi, keep_dups, fill_snpids, fill_coordinates = check_arguments(df, preprocessing, reference_panel, clumped, effect_column, keep_columns, fill_snpids, fill_coordinates, keep_multi, keep_dups)
        
        data=df.copy()
        
        self.checks=[] #Keep track of the checks performed, to avoid doing them again in the methods
            
        ## Keep only the main columns and rename them to the standard names
        data = adjust_column_names(data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns)

        ## Check that the CHR and POS columns are integers (pandas int) and not float or else
        for int_col in ["CHR", "POS"]:
            if int_col in data.columns and preprocessing > 0:
                data = check_int_column(data, int_col)
                self.checks.append(int_col)

        ## If SNP missing but CHR and POS present (or fill_snpid=True): fill the SNP column based on reference data 
        if ((("CHR" in data.columns) & ("POS" in data.columns) & ("SNP" not in data.columns)) or fill_snpids==True) and (not fill_snpids==False):
            data = fill_snpids_func(data, self.get_reference_panel(reference_panel))

        ## If CHR and/or POS columns missing and SNP present (or fill_coordinates=True): fill CHR/POS based on reference data
        if (((("CHR" not in data.columns) | ("POS" not in data.columns)) & ("SNP" in data.columns)) or fill_coordinates==True) and (not fill_coordinates==False):
            data = fill_coordinates_func(data, self.get_reference_panel(reference_panel))

        ## If NEA column is missing but EA and CHR/POS are present: fill it based on reference data
        if (("CHR" in data.columns) & ("POS" in data.columns) & ("NEA" not in data.columns) & ("EA" in data.columns)) and (preprocessing > 0):
            data = fill_nea(data, self.get_reference_panel(reference_panel))

        ## If NEA and EA columns are missing but CHR/POS are present: fill them based on reference data
        if (("CHR" in data.columns) & ("POS" in data.columns) & ("NEA" not in data.columns) & ("EA" not in data.columns)) and (preprocessing > 0):
            data = fill_ea_nea(data, self.get_reference_panel(reference_panel))

        ## Transform the effect column to Beta estimates if necessary
        if "BETA" in data.columns:
            data = check_beta_column(data, effect_column, preprocessing)

        ## Make sure the P column has appropriate values.
        if "P" in data.columns and preprocessing > 0:
            data = check_p_column(data)
            self.checks.append("P")

        ## If one of SE or P column missing but the other and BETA are present: fill it
        if preprocessing > 0:
            data = fill_se_p(data)

        ## Make sure the alleles columns are upper case strings and delete non-letters rows. Also set multiallelic SNPs to nan if not keep_multi=True.
        for allele_col in ["EA","NEA"]:
            if (allele_col in data.columns) and (preprocessing > 0):
                data = check_allele_column(data, allele_col, keep_multi)    
                self.checks.append(allele_col)

        ## Check whether some SNPs are present more than once based on the SNP column. Delete them if keep_dups=False.
        if "SNP" in data.columns and not keep_dups:
            data = check_snp_column(data)
            self.checks.append("SNP")

        ## Check the presence of the main columns and throw warnings if they are not present
        for column in ["CHR","POS","SNP","EA","NEA","BETA","SE","P"]:
            if not(column in data.columns):
                print(f"Warning: the data doesn't include a {column} column. This may become an issue later on.")       

        ## NA handling: if preprocessing = 2 return the number of rows with missing values and delete them.
        if preprocessing == 2:
            data = remove_na(data)
            self.checks.append("NA_removal")

        ## Reset index
        #self.data.reset_index(drop=True, inplace=True)
            
        ##Initialize attributes
        self.data = data
        self.init_attributes(name, clumped)
        
        return

    def init_attributes(self, name, clumped):
        """
        Initialize the name and outcome attributes.
        Guess if the data is clumped or not and set the data_clumped attribute accordingly.
        Initialize the ram and cpus attributes.
        Create the tmp_GENAL folder.
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
        self.cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=os.cpu_count()))
        #self.chunksize = round((self.ram*0.8*1024**2) / (1000*self.cpus) / 100000)*100000
        
        ## Create folders for temporary files
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")
        return    

    def get_reference_panel(self, reference_panel="eur"):
        """
        Return the reference panel dataframe. 
        """
        if not hasattr(self, "reference_panel"):
            if isinstance(reference_panel, pd.DataFrame):
                for col in ["CHR","SNP","POS","A1","A2"]:
                    if col not in reference_panel.columns:
                        raise ValueError(f"The {col} column is not present in the reference_panel provided and is necessary.")
                print("Using the provided reference_panel dataframe as the reference panel.")
                self.reference_panel = reference_df
            else:
                self.reference_panel = load_reference_panel(reference_panel)
        return self.reference_panel

    def copy(self):
        """
        Return a deep copy of the GENO instance.
        """
        Copy=copy.deepcopy(self)
        return Copy
        
    def save (self,path="",fmt="h5",sep="\t", header=True, clumped=True):
        """
        Save to a .h5 file in the "data" field (default) at the location specified by path or current folder if not.
        path: path of the folder to save the file
        fmt: format to use. Can be either .h5 (default), .csv, .txt
        sep: delimiter to use (only for .csv and .txt), default is \t
        header=True: to save the column names or not (only for .csv and .txt)
        clumped=False to save the full data, or the clumped data (=True)
        %To Do: add .vcf and .vcf.gz
        """
        if clumped:
            if not(hasattr(self,"data_clumped")):
                raise ValueError("clumped=True was used but the data is not clumped. \
                \n Use .clump() first, or clumped=False.")
            else:
                save_data(self.data_clumped, name = self.name+"_clumped", path=path, fmt=fmt, sep=sep, header=header)
        else:
            save_data(self.data, name = self.name, path=path, fmt=fmt, sep=sep, header=header)        
        return
    
    def set_phenotype(self, data, IID=None, PHENO=None, PHENO_type='', alternate_control=False):
        """ Set an attribute .phenotype which is a dataframe containing individual IDs and phenotype columns. This is required to run single-SNP association tests with the association_test method.
        data: pandas dataframe containing at least an individual IDs column and one phenotype column 
        IID: name of the individual IDs column in data. The IDs must match the genetic IDs of the fam file that will be used for association testing.
        PHENO: name of the phenotype column in data (the trait that will be the dependent variable for the association tests)
        PHENO_type="" The function will try to guess if the phenotype is binary or quantitative. To avoid the guessing, specify "quant" or "binary". If binary, the function will code it as 0 for controls, 1 for cases.
        alternate_control=False. The function assumes that for a binary trait, the controls are coded with the most frequent value. If that is not the case, use True.
        """
        data, phenotype_type = set_phenotype_func(data, PHENO, PHENO_type, IID, alternate_control)
        ## Set the attributes
        phenotype = (data, phenotype_type)
        self.phenotype = phenotype
        return     
    
    def association_test(self, covar=[], standardize=True, clumped=True):
        """
        Perform single-SNP association testing against a phenotype.
        Requires to set the phenotype with set_phenotype() and to call the extract_snps method to extract the variants from genomic files.
        Update the BETA, SE and P columns with the association tests results. 
        covar=[]. List of columns in the phenotype dataframe to use as covariates in the association tests.
        standardize=True to standardize a quantitative phenotype before association testing (usually done to make the results more understandable).
        clumped=True to run the association tests for the clumped list of SNPs. False to use the unclumped list.
        """
        ##Check that a phenotype has been set with the set_phenotype function.
        if not(hasattr(self,"phenotype")):
            raise ValueError("You first need to set a phenotype with .set_phenotype(data,PHENO,PHENO_type,IID)!") 
            
        data = self.data_clumped if clumped else self.data
        data = association_test_func(data, covar, standardize, self.name, self.phenotype[0], self.phenotype[1])
        #Update attributes
        if clumped: self.data_clumped = data
        else: self.data = data
        print(f"The BETA, SE, P columns of the .data{'_clumped' if clumped else ''} attribute have been updated.")
        return
            
    def query_outcome (self, outcome, name="", proxy=True, reference_panel="eur", kb=5000, r2=0.6, window_snps=5000):
        """ 
        Query the outcome data, with or without proxying, to setup the dataframes required for running Mendelian Randomization (MR) with the data_clumped as exposure.
        Store a triple to the outcome attribute: (exposure_data, outcome_data, name) that is ready to run MR methods.
        name="": the outcome data name
        #reset_outcome = True: Delete the outcomes stored in the outcome attribute. Otherwise, the processed data will be added to the list of already existing outcomes.
        Below are the parameters used when proxies are searched:
        reference_panel="EUR" The reference population to get linkage disequilibrium values and find proxies (only used if proxy=True). Takes values in "EUR", "SAS", "AFR", "EAS", "AMR". It is also possible to provide a path leading to a specific bed/bim/fam reference panel.
        outcome: can be a GENO object (from a GWAS) or a filepath of the following types: .h5 or .hdf5 (output of the .save() attribute);
        kb=5000: width of the genomic window to look for proxies
        r2=0.6: minimum linkage disequilibrium value with the main SNP for a proxy to be included 
        window_snps=5000: compute the LD value for SNPs that are not more than x SNPs apart from the main SNP        
        """   
        ## Check that data_clumped attribute exists
        if not(hasattr(self,"data_clumped")):
            raise ValueError("The data needs to be clumped before using the query_outcome method! \
                \n Use .clump() first, or clumped=True when declaring the GENO.")
        
        exposure, outcome, name = query_outcome_func(self.data_clumped, outcome, name, proxy, reference_panel, kb, r2, window_snps, self.cpus)     
        self.outcome = (exposure, outcome, name) 
        return
    
    def MR (self, methods = ["IVW","IVW-FE","UWR", "WM","WM-pen","Simple-median","Sign","Egger","Egger-boot"], action = 2, heterogeneity = False, eaf_threshold=0.42, nboot = 10000, penk = 20):
        """
        Run a Mendelian Randomization with the .clumped_data as exposure and the outcome queried with the .query_outcome() method.
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
        ## Check that query_outcome has been called
        if not hasattr(self, "outcome"):
            raise ValueError("You must first call query_outcome() before running MR.")
                               
        return MR_func(self.outcome, methods, action , heterogeneity, eaf_threshold, nboot, penk, self.name, self.cpus)  
    
    def MRpresso(self, action = 2, eaf_threshold=0.42, n_iterations = 10000, outlier_test = True, distortion_test = True, significance_p = 0.05, cpus = 1):
        """
        Perform the MR-PRESSO Mendelian Randomization algorithm for detection of horizontal pleiotropy with the .clumped_data as exposure and the outcome queried with the .query_outcome() method. 
        action=2: Determines how to treat palindromes in the harmonizing step between exposure and outcome data. 
        =1: Doesn't attempt to flip them (= Assume all alleles are coded on the forward strand)
        =2: Use allele frequencies (EAF) to attempt to flip them (conservative, default)
        =3: Remove all palindromic SNPs (very conservative).
        eaf_threshold=0.42: Maximal effect allele frequency accepted when attempting to flip palindromic SNPs (only relevant if action=2)
        n_iterations: number of steps performed (random data generation). Increase this number for better stability of the results.
        outlier_test: identification of outlier SNPs responsible for horizontal pleiotropy (performed if global test p_value < significance_p)
        distortion_test: testing of significant distortion in the causal estimates before and after outlier removal (performed if global test p_value < significance_p)
        significance_p: statistical significance threshold for the detection of horizontal pleiotropy (both for the global test and outlier identification)

        return: [mod_table, GlobalTest, OutlierTest, DistortionTest]
            mod_table: table with the original (before outlier removal) and outlier-corrected (after outlier removal) inverse variance-weighted MR results 
            GlobalTest: p-value of the global MR-PRESSO test indicating the presence of horizontal pleiotropy
            OutlierTest: table assigning one p-value to each SNP representing the likelihood of this SNP being responsible for the global pleiotropy (set to nan if global test p_value > SignifThreshold)
            DistortionTest: p-value for the distortion test 
        """ 
        ## Check that query_outcome has been called
        if not hasattr(self, "outcome"):
            raise ValueError("You must first call query_outcome() before running MR.")

        return mrpresso_func(self.outcome, action, eaf_threshold, n_iterations, outlier_test, distortion_test, significance_p, self.cpus)
            
    def clump(self, kb=250, r2=0.1, p1=5e-8, p2=0.01, reference_panel="eur"):
        """ Clump the data in .data and assign it to the .data_clumped attribute. The clumping is done with plink.
        kb=250. clumping window in thousands SNPs
        r2=0.1. linkage disequilibrium threshold (between 0 and 1)
        p1=5e-8. p-value used during the clumping (the SNPs above this threshold are not considered)
        p2=0.01 p-value used after the clumping to further filter the clumped SNPs (if p2<p1, it won't be considered)
        reference_panel="eur" The reference population to get linkage disequilibrium values. Takes values in "eur", "sas", "afr", "eas", "amr". It is also possible to provide a path leading to a specific bed/bim/fam reference panel.
        """
        #Check input 
        clumped_data, checks = clump_data(self.data, reference_panel, get_plink19_path(), kb, r2, p1, p2, self.name, self.ram, self.checks) 
        self.checks = checks
        if clumped_data is not None:
            self.data_clumped = clumped_data
            print("The clumped data is stored in the .data_clumped attribute.")
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

            
    def extract_snps(self, clumped=True, path=""):
        """
        Extract the list of SNPs present in .data (or in .data_clumped if clumped=True) from the provided files. 
        path: path to the genomic files. If the files are split by chromosomes, replace the chromosome number by '$'. For instance path = ukb_chr$_file. The provided path is saved and if you call this function again, you don't need to specify the path (provided you want to use the same genomic files).
        #slurm=False: if True, the extraction and merge jobs will be launched as batch scripts
        The output is a bed/bim/fam triple called {name}_extract_allchr including the SNPs from the UKB.
        """
        snp_list = self.data_clumped["SNP"] if clumped else self.data["SNP"]
        extract_snps_func(snp_list, self.name, path)
        return

        
    def prs(self, name =None, clumped = True, weighted=True, path = None):
        """Compute a PRS and save it to a .csv file in the current directory.
        name = None. Name or path of the saved prs file. If not given, use the name of the GENO object.
        clumped = True. To use the data contained in the .data_clumped attribute (after using the .clump() method or clump=True when initializing the GENO object.) False to use the unclumped data (in the .data attribute)
        weighted=True to perform a PRS weighted by the BETA column estimates. False for an unweighted PRS (equivalent to BETAs == 1)
        path = None. Can be used to provide a path to a bed/bim/fam set of genetic files to use for PRS calculation. If not provided, will use the genetic data extracted with the .extract_snps method.
        """
        ## Verify that the data has been clumped (or at least assigned to the data_clumped attribute)
        if not(hasattr(self,"data_clumped")) and clumped:
            raise ValueError("clumped = True argument but the data is not clumped. \ \n Use .clump() first to clump the data or use clumped = False to perform PRS on the unclumped data.")
        if clumped:
            print("Using clumped data for PRS calculation.")
            prs_data = prs_func(self.data_clumped, weighted, path, checks=self.checks, ram=self.ram, name=self.name)       
        else:
            print("Using unclumped data for PRS calculation.")
            prs_data = prs_func(self.data, weighted, path, checks=self.checks, ram=self.ram, name=self.name)

        if name is None:
            prs_filename = f"{self.name}_prs.csv"
        else:
            prs_filename = os.path.splitext(name)[0] + ".csv"
        prs_data.to_csv(prs_filename, index=False, header=True)
        print (f"PRS data saved to {prs_filename}")
        return prs_data

                    
    def lift(self, clumped=True, start="hg19",end="hg38", replace=False, extraction_file=False, chain_file="", name=None, liftover=False, liftover_path=""):
        """Perform a liftover from a genetic build to another.
        If the chain file required to do the liftover is not present, download it first. It is also possible to manually provide the path to the chain file.
        start="hg19": current build of the data
        end="hg38": build to be lifted to
        If clumped==True, lift the clumped data, otherwise the unclumped data
        replace=False: whether to change the GENO object and update the .data or .data_clumped attributes
        extraction_file==True, also print a CHR POS SNP space delimited file for extraction in All of Us (WES data)
        chain_file="": path to a local chain file to be used for the lift. If provided, the start and end arguments are not considered.
        name=None: can be used to specify a filename or filepath (without extension) to save the lifted dataframe. If not provided, will be saved in the current folder as [name]_lifted.txt
        """
        #Selecting data to lift
        if clumped==False or not(hasattr(self,"data_clumped")):
            data = self.data if replace else self.data.copy()
        elif clumped==True:
            data = self.data_clumped if replace else self.data_clumped.copy()
            
        print(f"Lifting the {'' if clumped else 'un'}clumped data{' inplace' if replace else ''}. The .data{'_clumped' if clumped else ''} attribute will {'' if replace else 'not'} be modified. Use replace={'False' if replace else 'True'} to {'leave it as is' if replace else 'lift inplace'}.")
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = lift_data(data, start=start, end=end, replace=replace, extraction_file=extraction_file, chain_file=chain_file, name=f"{self.name if name is None else name}")
        return data
                

