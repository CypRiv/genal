import pandas as pd
import numpy as np
import warnings
import os, subprocess
import subprocess
import copy
import psutil

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
# Get proxies (simply return a list of proxies)
# Multi-MR with python MR
# Warning that users might not have shell (for the .ram attribute)
# Phenoscanner


class GENO:
    """
    A class to handle GWAS-derived data, including SNP rsID, genome position, 
    SNP-trait effects, and effect allele frequencies.

    Attributes:
        name (str): Name of the object.
        data (pd.DataFrame): Main DataFrame containing SNP data.
        data_clumped (pd.DataFrame, optional): DataFrame of clumped data. Initialized after 
            'clump' method execution or when using clumped=True during initialization.
        phenotype (pd.DataFrame, str): Tuple with a DataFrame of individual-level phenotype 
            data and a string representing the phenotype trait column. Initialized after 
            running the 'set_phenotype' method.
        MR_data (pd.DataFrame, pd.DataFrame, str): Tuple containing DataFrames for associations 
            with exposure and outcome, and a string for the outcome name. Initialized after 
            running the 'query_outcome' method.
        ram (int): Available memory.
        cpus (int): Number of available CPUs.
        checks (list): List of checks performed on the main DataFrame.
        reference_panel (pd.DataFrame): Reference population SNP data used for SNP info 
            adjustments. Initialized when first needed.

    Methods:
        clump():
            Clumps the main data and stores the result in data_clumped.

        prs():
            Computes Polygenic Risk Score on genomic data.

        set_phenotype():
            Assigns a DataFrame with individual-level data and a phenotype trait to 
            the phenotype attribute.

        association_test():
            Computes SNP-trait effect estimates, standard errors, and p-values.

        query_outcome():
            Extracts SNPs from outcome data with proxying and initializes MR_data.

        MR():
            Performs Mendelian Randomization between SNP-exposure and SNP-outcome data.

        MRpresso():
            Executes the MR-PRESSO algorithm for horizontal pleiotropy correction between 
            SNP-exposure and SNP-outcome data.

        lift():
            Lifts SNP data from one genomic build to another.
    """
    
    def __init__(self, df, name="noname", CHR="CHR",POS="POS",SNP="SNP",EA="EA",NEA="NEA",BETA="BETA",SE="SE",P="P",EAF="EAF", preprocessing=1, reference_panel="eur", clumped=None, effect_column=None, keep_columns=None, keep_multi=None, keep_dups=None, fill_snpids=None, fill_coordinates=None):
        """
        Initializes the GENO object used to store and transform Single Nucleotide Polymorphisms (SNP) data.
        
        Args:
            df (pd.DataFrame): DataFrame where each row represents a SNP.
            name (str, optional): Name for the object. Defaults to "noname".
            CHR (str, optional): Column name for chromosome. Defaults to "CHR".
            POS (str, optional): Column name for genomic position. Defaults to "POS".
            SNP (str, optional): Column name for SNP identifier. Defaults to "SNP".
            EA (str, optional): Column name for effect allele. Defaults to "EA".
            NEA (str, optional): Column name for non-effect allele. Defaults to "NEA".
            BETA (str, optional): Column name for effect estimate. Defaults to "BETA".
            SE (str, optional): Column name for effect standard error. Defaults to "SE".
            P (str, optional): Column name for p-value. Defaults to "P".
            EAF (str, optional): Column name for effect allele frequency. Defaults to "EAF".
            preprocessing (int, optional): Level of preprocessing to apply. Options include:
                - 0: The dataframe is not modified.
                - 1: Missing columns are added based on reference data and invalid values set to NaN, but no rows are deleted.
                - 2: Missing columns are added, and rows with missing, duplicated, or invalid values are deleted.
                Defaults to 1.
            reference_panel (str or pd.DataFrame, optional): Reference panel for SNP adjustments. Can be a string representing ancestry classification ("eur", "afr", "eas", "sas", "amr") or a DataFrame with ["CHR","SNP","POS","A1","A2"] columns or a path to a .bim file. Defaults to "eur".
            clumped (bool, optional): Specifies if the data is already clumped. If None, the method tries to determine this. Defaults to None.
            effect_column (str, optional): Specifies the type of effect column ("BETA" or "OR"). If None, the method tries to determine it. Odds Ratios will be log-transformed and the standard error adjusted. Defaults to None.
            keep_columns (bool, optional): Determines if non-main columns should be kept. If None, defers to preprocessing value. Defaults to None.
            keep_multi (bool, optional): Determines if multiallelic SNPs should be kept. If None, defers to preprocessing value. Defaults to None.
            keep_dups (bool, optional): Determines if rows with duplicate SNP IDs should be kept. If None, defers to preprocessing value. Defaults to None.
            fill_snpids (bool, optional): Decides if the SNP (rsID) column should be created or replaced based on CHR/POS columns and a reference genome. If None, defers to preprocessing value. Defaults to None.
            fill_coordinates (bool, optional): Decides if CHR and/or POS should be created or replaced based on SNP column and a reference genome. If None, defers to preprocessing value. Defaults to None.
        """
        
        #Check arguments and solve arguments logic.
        keep_columns, keep_multi, keep_dups, fill_snpids, fill_coordinates = check_arguments(df, preprocessing, reference_panel, clumped, effect_column, keep_columns, fill_snpids, fill_coordinates, keep_multi, keep_dups)
        
        # Create a copy of input dataframe
        data=df.copy()
        
        # List to keep track of performed checks
        self.checks=[]
            
        # Standardize column names based on provided parameters +/- delete other columns
        data = adjust_column_names(data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns)

        # Ensure CHR and POS columns are integers if preprocessing is enabled
        for int_col in ["CHR", "POS"]:
            if int_col in data.columns and preprocessing > 0:
                data = check_int_column(data, int_col)
                self.checks.append(int_col)

        # Fill missing SNP column from reference data if necessary
        should_fill_snpids = (("CHR" in data.columns) and ("POS" in data.columns) and 
                              ("SNP" not in data.columns)) or fill_snpids
        if should_fill_snpids and fill_snpids is not False:
            data = fill_snpids_func(data, self.get_reference_panel(reference_panel))

        # Fill missing CHR/POS columns from reference data if necessary
        should_fill_coordinates = ((("CHR" not in data.columns) or ("POS" not in data.columns)) and 
                                  ("SNP" in data.columns)) or fill_coordinates
        if should_fill_coordinates and fill_coordinates is not False:
            data = fill_coordinates_func(data, self.get_reference_panel(reference_panel))

        # Fill missing NEA column from reference data if necessary and preprocessing is enabled
        
        missing_nea_condition = ("CHR" in data.columns and "POS" in data.columns and 
                                    "NEA" not in data.columns and "EA" in data.columns)
        if missing_nea_condition and preprocessing > 0:
            data = fill_nea(data, self.get_reference_panel(reference_panel))

        # Fill missing EA and NEA columns from reference data if necessary and preprocessing is enabled
        missing_ea_nea_condition = ("CHR" in data.columns and "POS" in data.columns and 
                                    "NEA" not in data.columns and "EA" not in data.columns)
        if missing_ea_nea_condition and preprocessing > 0:
            data = fill_ea_nea(data, self.get_reference_panel(reference_panel))

        # Convert effect column to Beta estimates if present
        if "BETA" in data.columns:
            data = check_beta_column(data, effect_column, preprocessing)

        # Ensure P column contains valid values
        if "P" in data.columns and preprocessing > 0:
            data = check_p_column(data)
            self.checks.append("P")

        # Fill missing SE or P columns if necessary
        if preprocessing > 0:
            data = fill_se_p(data)

        # Process allele columns
        for allele_col in ["EA","NEA"]:
            if (allele_col in data.columns) and (preprocessing > 0):
                data = check_allele_column(data, allele_col, keep_multi)    
                self.checks.append(allele_col)

        # Check for and handle duplicate SNPs if necessary
        if "SNP" in data.columns and not keep_dups:
            data = check_snp_column(data)
            self.checks.append("SNP")

        # Warn if essential columns are missing
        for column in ["CHR","POS","SNP","EA","NEA","BETA","SE","P"]:
            if not(column in data.columns):
                print(f"Warning: the data doesn't include a {column} column. This may become an issue later on.")       

        # Remove missing values if preprocessing level is set to 2
        if preprocessing == 2:
            data = remove_na(data)
            self.checks.append("NA_removal")

        ## Reset index
        #self.data.reset_index(drop=True, inplace=True)
            
        # Set object attributes
        self.data = data
        self.init_attributes(name, clumped)
        
        return
  

    def init_attributes(self, name, clumped):
        """
        Initializes several attributes for the GENO object, including the `name`,
        `outcome`, and memory/CPU related attributes. Also guesses whether the provided
        data is clumped based on the `P` column.

        Args:
            name (str): Name for the GENO object.
            clumped (bool, optional): Specifies if the data is already clumped. If None,
                the method tries to determine this based on the `P` column.

        Attributes:
            name (str): Name for the GENO object.
            outcome (list): List of outcomes (initialized as empty).
            data_clumped (pd.DataFrame): Dataframe of clumped data if `clumped` is True.
            cpus (int): Number of CPUs to be used.
            ram (int): Amount of RAM to be used in MBs.
        """

        self.name = name
        if name == "noname":
            print("You haven't passed a specific name to GENO. It is recommended if you are working with several GENO objects.")
        self.outcome = []

        # Guess if clumped or not
        if clumped is None:
            if self.data["P"].max() < 0.1:
                print("This data looks clumped. Use clumped=False if it is not.")
                clumped = True
            else:
                print("This data looks not clumped. Use clumped=True if it is.")
                clumped = False

        # If the data is already clumped (clumped=True): also assign the data to self.data_clumped
        if clumped:
            self.data_clumped = self.data

        # Set the maximal amount of ram/cpu to be used by the methods and dask chunksize
        self.cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=os.cpu_count()))
        non_hpc_ram_per_cpu = psutil.virtual_memory().available / (1024 ** 2 * self.cpus)
        ram_per_cpu = int(os.environ.get('SLURM_MEM_PER_CPU', default=non_hpc_ram_per_cpu))
        self.ram = int(ram_per_cpu * self.cpus * 0.8)
        #self.chunksize = round((self.ram*0.8*1024**2) / (1000*self.cpus) / 100000)*100000

        # Create folders for temporary files
        if not os.path.exists("tmp_GENAL"):
            try:
                os.makedirs("tmp_GENAL")
            except OSError:
                raise OSError("Unable to create the 'tmp_GENAL' directory. Check permissions.")

            
    def get_reference_panel(self, reference_panel="eur"):
        """
        Retrieve or set the reference panel for the GENO object.

        If the GENO object does not have a reference panel attribute set,
        this method will try to set it based on the provided `reference_panel`
        argument. This can be either a string indicating a predefined reference panel 
        or a DataFrame with specific columns or a path to a .bim file.

        Args:
            reference_panel (str or pd.DataFrame, optional): Either a string indicating a predefined
                reference panel (default is "eur") or a DataFrame with necessary columns or a valid path 
                to a .bim file

        Returns:
            pd.DataFrame: The reference panel DataFrame for the GENO object.

        Raises:
            ValueError: If the provided DataFrame doesn't have the necessary columns.
        """

        # Check if the object already has a reference panel set
        if not hasattr(self, "reference_panel"):

            # If the provided reference_panel is a DataFrame, verify its structure and dtypes
            if isinstance(reference_panel, pd.DataFrame):
                required_cols = ["CHR", "SNP", "POS", "A1", "A2"]

                for col in required_cols:
                    if col not in reference_panel.columns:
                        raise ValueError(f"The {col} column is not present in the reference_panel provided and is necessary.")   

                print("Using the provided reference_panel dataframe as the reference panel.")
                self.reference_panel = reference_panel.copy()
            else:
                # Load the reference panel based on the provided string identifier
                self.reference_panel = load_reference_panel(reference_panel)

        return self.reference_panel

    
    def clump(self, kb=250, r2=0.1, p1=5e-8, p2=0.01, reference_panel="eur"):
        """
        Clump the data based on linkage disequilibrium and assign it to the .data_clumped attribute. 
        The clumping process is executed using plink.

        Args:
            kb (int, optional): Clumping window in terms of thousands of SNPs. Default is 250.
            r2 (float, optional): Linkage disequilibrium threshold, values between 0 and 1. Default is 0.1.
            p1 (float, optional): P-value threshold during clumping. SNPs above this value are not considered. Default is 5e-8.
            p2 (float, optional): P-value threshold post-clumping to further filter the clumped SNPs. If p2 < p1, it won't be considered. Default is 0.01.
            reference_panel (str, optional): The reference population for linkage disequilibrium values. Accepts values "eur", "sas", "afr", "eas", "amr". Alternatively, a path leading to a specific bed/bim/fam reference panel can be provided. Default is "eur".

        Returns:
            None: The result is stored in the .data_clumped attribute.
        """

        # Validate input and clump the data using the specified parameters
        clumped_data, checks = clump_data(self.data, reference_panel, get_plink19_path(), kb, r2, p1, p2, self.name, self.ram, self.checks) 

        # Update checks attribute with latest checks
        self.checks = checks

        # If clumped data is successfully generated, assign it to the object's attribute
        if clumped_data is not None:
            self.data_clumped = clumped_data
            print("The clumped data is stored in the .data_clumped attribute.")
  

    def extract_snps(self, clumped=True, path=None):
        """
        Extract the list of SNPs present in the data.

        Args:
            clumped (bool, optional): If True, SNPs will be extracted from .data_clumped. 
                                      If False, from .data. Default is True.
            path (str, optional): Path to a bed/bim/fam set of genetic files. 
                If files are split by chromosomes, replace the chromosome number with '$'. 
                For instance: path = "ukb_chr$_file". Default is None.

        Returns:
            None: The output is a bed/bim/fam triple in the tmp_GENAL folder 
            with the format "{name}_extract_allchr" which includes the SNPs from the UKB.

        Notes:
            The provided path is saved to the config file. If this function is called again, 
            you don't need to specify the path if you want to use the same genomic files.
        """

        # Extract the relevant SNP list based on whether clumped data is desired or not
        snp_list = self.data_clumped["SNP"] if clumped else self.data["SNP"]

        # Extract SNPs using the provided path and SNP list
        extract_snps_func(snp_list, self.name, path)

        
    def prs(self, name=None, clumped=True, weighted=True, path=None):
        """
        Compute a Polygenic Risk Score (PRS) and save it as a CSV file in the current directory.

        Args:
            name (str, optional): Name or path of the saved PRS file. If not given, 
                                  it defaults to the name of the GENO object.
            clumped (bool, optional): If True, uses data in .data_clumped. If False, 
                                      uses data in .data. Default is True.
            weighted (bool, optional): If True, performs a PRS weighted by the BETA column estimates. 
                                       If False, performs an unweighted PRS. Default is True.
            path (str, optional): Path to a bed/bim/fam set of genetic files for PRS calculation.
                                  If files are split by chromosomes, replace the chromosome number 
                                  with '$'. For instance: path = "ukb_chr$_file". 
                                  If not provided, it will use the genetic data extracted 
                                  with the .extract_snps method. Default is None.

        Returns:
            pd.DataFrame: The computed PRS data.

        Raises:
            ValueError: If the data hasn't been clumped and 'clumped' parameter is True.
        """

        # Check if clumped data exists when the clumped argument is True
        if not hasattr(self, "data_clumped") and clumped:
            raise ValueError("clumped=True but data isn't clumped. Use .clump() first or set clumped=False.")

        # Compute PRS based on whether clumped or unclumped data is desired
        if clumped:
            print("Using clumped data for PRS calculation.")
            prs_data = prs_func(self.data_clumped, weighted, path, checks=self.checks, ram=self.ram, name=self.name)       
        else:
            print("Using unclumped data for PRS calculation.")
            prs_data = prs_func(self.data, weighted, path, checks=self.checks, ram=self.ram, name=self.name)

        # Save the computed PRS data as a CSV file
        if name is None:
            prs_filename = f"{self.name}_prs.csv"
        else:
            prs_filename = os.path.splitext(name)[0] + ".csv"
        prs_data.to_csv(prs_filename, index=False, header=True)
        print(f"PRS data saved to {prs_filename}")

        return prs_data
  

    def set_phenotype(self, data, IID=None, PHENO=None, PHENO_type=None, alternate_control=False):
        """
        Assign a phenotype dataframe to the .phenotype attribute.

        This method sets the .phenotype attribute which is essential to perform 
        single-SNP association tests using the association_test method.

        Args:
            data (pd.DataFrame): DataFrame containing individual-level row data with at least an individual IDs column 
                                 and one phenotype column.
            IID (str, optional): Name of the individual IDs column in 'data'. These IDs should 
                                 correspond to the genetic IDs in the FAM file that will be used for association testing.
            PHENO (str, optional): Name of the phenotype column in 'data' which will be used 
                                   as the dependent variable for association tests.
            PHENO_type (str, optional): If not specified, the function will try to infer if 
                                        the phenotype is binary or quantitative. To bypass this, 
                                        use "quant" for quantitative or "binary" for binary phenotypes.
                                        Default is None.
            alternate_control (bool, optional): By default, the function assumes that for a binary 
                                                trait, the controls have the most frequent value. 
                                                Set to True if this is not the case. Default is False.

        Returns:
            None: Sets the .phenotype attribute for the instance.
        """
        processed_data, inferred_pheno_type = set_phenotype_func(data, PHENO, PHENO_type, IID, alternate_control)

        # Assign the processed data and inferred phenotype type to the .phenotype attribute
        self.phenotype = (processed_data, inferred_pheno_type)


    def association_test(self, covar=[], standardize=True, clumped=True):
        """
        Conduct single-SNP association tests against a phenotype.

        This method requires the phenotype to be set using the set_phenotype() function. 
        The method also expects the extract_snps method to have been called prior to this.

        Args:
            covar (list, optional): List of columns in the phenotype dataframe to be used 
                                    as covariates in the association tests. Default is an empty list.
            standardize (bool, optional): If True, it will standardize a quantitative phenotype 
                                          before performing association tests. This is typically done 
                                          to make results more interpretable. Default is True.
            clumped (bool, optional): If True, association tests will be run for the clumped 
                                      list of SNPs. If False, it will use the unclumped list. Default is True.

        Returns:
            None: Updates the BETA, SE, and P columns of the data attribute based on the results 
                  of the association tests.
        """

        # Ensure that the phenotype has been set using set_phenotype
        if not hasattr(self, "phenotype"):
            raise ValueError("You first need to set a phenotype using .set_phenotype(data, PHENO, PHENO_type, IID)!")

        # Determine which data to use (clumped or unclumped)
        target_data = self.data_clumped if clumped else self.data

        # Perform the association test
        updated_data = association_test_func(target_data, covar, standardize, self.name, self.phenotype[0], self.phenotype[1])

        # Update the instance data based on whether we're using clumped data or not
        if clumped:
            self.data_clumped = updated_data
        else:
            self.data = updated_data

        print(f"The BETA, SE, P columns of the .data{'_clumped' if clumped else ''} attribute have been updated.")
  

    def query_outcome(self, outcome, name=None, proxy=True, reference_panel="eur", 
                      kb=5000, r2=0.6, window_snps=5000):
        """
        Prepares dataframes required for Mendelian Randomization (MR) with `data_clumped` as exposure.

        Queries the outcome data, with or without proxying, and assigns a tuple to 
        the outcome attribute: (exposure_data, outcome_data, name) ready for MR methods.

        Args:
            outcome: Can be a GENO object (from a GWAS) or a filepath of types: .h5 or .hdf5 (created with the :meth:`GENO.save` method.
            name (str, optional): Name for the outcome data. Defaults to None.
            proxy (bool, optional): If true, proxies are searched. Default is True.
            reference_panel (str, optional): The reference population to get linkage 
                disequilibrium values and find proxies (only if proxy=True). Acceptable values 
                include "EUR", "SAS", "AFR", "EAS", "AMR" or a path to a specific bed/bim/fam panel. 
                Default is "EUR".
            kb (int, optional): Width of the genomic window to look for proxies. Default is 5000.
            r2 (float, optional): Minimum linkage disequilibrium value with the main SNP 
                for a proxy to be included. Default is 0.6.
            window_snps (int, optional): Compute the LD value for SNPs that are not 
                more than x SNPs away from the main SNP. Default is 5000.

        Returns:
            None: Sets the `MR_data` attribute for the instance.
        """

        # Ensure that data_clumped attribute exists
        if not hasattr(self, "data_clumped"):
            raise ValueError("The data needs to be clumped before using the query_outcome method."
                             " Use .clump() or set clumped=True when initializing the GENO.")

        exposure, outcome_data, outcome_name = query_outcome_func(
            self.data_clumped, outcome, name, proxy, reference_panel, kb, r2, window_snps, self.cpus
        )

        # Assign the processed data to the MR_data attribute
        self.MR_data = (exposure, outcome_data, outcome_name)


    def MR(self, methods=["IVW","IVW-FE","UWR","WM","WM-pen","Simple-median","Sign","Egger","Egger-boot"], 
           action=2, eaf_threshold=0.42, heterogeneity=False, nboot=10000, penk=20):
        """
        Executes Mendelian Randomization (MR) using the `data_clumped` attribute as exposure data and `MR_data` attribute as outcome data queried using the `query_outcome` method.

        Args:
            methods (list, optional): List of MR methods to run. Possible options include:
                "IVW": inverse variance-weighted with random effects and under-dispersion correction
                "IVW-FE": inverse variance-weighted with fixed effects
                "IVW-RE": inverse variance-weighted with random effects and without under-dispersion correction
                "UWR": unweighted regression
                "WM": weighted median (bootstrapped standard errors)
                "WM-pen": penalised weighted median (bootstrapped standard errors)
                "Simple-median": simple median (bootstrapped standard errors)
                "Sign": sign concordance test
                "Egger": egger regression
                "Egger-boot": egger regression with bootstrapped standard errors
                Default is ["IVW","IVW-FE","UWR","WM","WM-pen","Simple-median","Sign","Egger","Egger-boot"].
            action (int, optional): How to treat palindromes during harmonizing between 
                exposure and outcome data. Accepts:
                1: Doesn't flip them (Assumes all alleles are on the forward strand)
                2: Uses allele frequencies to attempt to flip (conservative, default)
                3: Removes all palindromic SNPs (very conservative)
            eaf_threshold (float, optional): Max effect allele frequency accepted when 
                flipping palindromic SNPs (relevant if action=2). Default is 0.42.
            heterogeneity (bool, optional): If True, includes heterogeneity tests in the results (Cochran's Q test).Default is False.
            nboot (int, optional): Number of bootstrap replications for methods with bootstrapping. Default is 10000.
            penk (int, optional): Penalty value for the WM-pen method. Default is 20.

        Returns:
            pd.DataFrame: A table with MR results.
        """

        # Ensure that query_outcome has been previously called
        if not hasattr(self, "MR_data"):
            raise ValueError("You must first call query_outcome() before running MR.")

        return MR_func(
            self.MR_data, methods, action, heterogeneity, eaf_threshold, nboot, penk, self.name, self.cpus
        )
 

    def MRpresso(self, action=2, eaf_threshold=0.42, n_iterations=10000, outlier_test=True, 
                 distortion_test=True, significance_p=0.05, cpus=-1):
        """
        Executes the MR-PRESSO Mendelian Randomization algorithm for detection and correction of horizontal pleiotropy.

        Args:
            action (int, optional): Treatment for palindromes during harmonizing between 
                exposure and outcome data. Options:
                - 1: Don't flip (assume all alleles are on the forward strand)
                - 2: Use allele frequencies to flip (default)
                - 3: Remove all palindromic SNPs
            eaf_threshold (float, optional): Max effect allele frequency when flipping 
                palindromic SNPs (relevant if action=2). Default is 0.42.
            n_iterations (int, optional): Number of random data generation steps for 
                improved result stability. Default is 10000.
            outlier_test (bool, optional): Identify outlier SNPs responsible for horizontal 
                pleiotropy if global test p_value < significance_p. Default is True.
            distortion_test (bool, optional): Test significant distortion in causal estimates 
                before and after outlier removal if global test p_value < significance_p. 
                Default is True.
            significance_p (float, optional): Statistical significance threshold for 
                horizontal pleiotropy detection (both global test and outlier identification).
                Default is 0.05.
            cpus (int, optional): number of cpu cores to be used for the parallel random data generation.

        Returns: 
            list: Contains the following elements:
                - mod_table: DataFrame containing the original (before outlier removal) 
                             and outlier-corrected (after outlier removal) inverse variance-weighted MR results.
                - GlobalTest: p-value of the global MR-PRESSO test indicating the presence of horizontal pleiotropy.
                - OutlierTest: DataFrame assigning a p-value to each SNP representing the likelihood of this 
                               SNP being responsible for the global pleiotropy. Set to NaN if global test p_value > significance_p.
                - DistortionTest: p-value for the distortion test.
        """

        if not hasattr(self, "outcome"):
            raise ValueError("You must first call query_outcome() before running MR.")
        cpus = self.cpus if cpus == -1 else cpus
        
        return mrpresso_func(
            self.MR_data, action, eaf_threshold, n_iterations, outlier_test, distortion_test, significance_p, cpus
        )


    def lift(self, clumped=True, start="hg19", end="hg38", replace=False, 
             extraction_file=False, chain_file=None, name=None, liftover=False, 
             liftover_path=None):
        """
        Perform a liftover from one genetic build to another.

        Args:
            clumped (bool, optional): If True, uses data in .data_clumped. If False, 
                                          uses data in .data. Default is True.
            start (str, optional): Current build of the data. Default is "hg19".
            end (str, optional): Target build for the liftover. Default is "hg38".
            replace (bool, optional): If True, updates the object data attributes in place. Default is False.
            extraction_file (bool, optional): If True, prints a CHR POS SNP space-delimited 
                file. Default is False.
            chain_file (str, optional): Path to a local chain file for the lift. 
                If provided, `start` and `end` arguments are not considered. Default is None.
            name (str, optional): Filename or filepath (without extension) to save the lifted dataframe. 
                Default saves as [name]_lifted.txt in the current directory.

        Returns:
            pd.DataFrame: Data after being lifted.
        """

        data_type = 'clumped' if clumped else 'unclumped'
        data = getattr(self, f"data_{data_type}" if clumped else "data")
        if not replace:
            data = data.copy()

        print(f"Lifting the {data_type} data{' inplace' if replace else ''}. "
              f"The .data{'_clumped' if clumped else ''} attribute will {'' if replace else 'not '}be modified. "
              f"Use replace={'False' if replace else 'True'} to {'leave it as is' if replace else 'lift inplace'}.")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = lift_data(data, start=start, end=end, extraction_file=extraction_file, 
                             chain_file=chain_file, name=f"{self.name if name is None else name}")

        return data 
   

    def standardize(self):
        """
        Standardize the Betas and adjust the SE column accordingly.

        Raises:
            ValueError: If the required columns are not found in the data.
        """
        required_columns = ["BETA", "SE"]
        for column in required_columns:
            if column not in self.data.columns:
                raise ValueError(f"The column {column} is not found in the data!")

        self.data["BETA"] = (self.data.BETA - np.mean(self.data.BETA)) / np.std(self.data.BETA)
        self.data["SE"] = np.abs(self.data.BETA / st.norm.ppf(self.data.P / 2))
        print("The Beta column has been standardized and the SE column has been adjusted.")


    def sort_group(self, method="lowest_p"):
        """
        Handle duplicate SNPs. Useful if the instance combines different GENOs.

        Args:
            method (str, optional): How to handle duplicates. Default is "lowest_p", 
                                    which retains the lowest P-value for each SNP.

        Returns:
            None
        """
        if method == "lowest_p":
            self.data = self.data.sort_values(by=["P"])
            self.data = self.data.groupby(by=["SNP"]).first().reset_index(drop=False)
        return


    def copy(self, name=None):
        """
        Create a deep copy of the GENO instance. 
        
        Args:
            name (str, optional): Name to give the copied GENO object

        Returns:
            GENO: A deep copy of the instance.
        """
        copied_geno = copy.deepcopy(self)
        if name is not None: copied_geno.name = name
        return copied_geno


    def save(self, path="", fmt="h5", sep="\t", header=True, clumped=True):
        """
        Save the GENO data to a file.

        Args:
            path (str, optional): Folder path to save the file. Defaults to the current directory.
            fmt (str, optional): File format. Options: .h5 (default), .csv, .txt. Future: .vcf, .vcf.gz.
            sep (str, optional): Delimiter for .csv and .txt formats. Default is tab.
            header (bool, optional): Save column names for .csv and .txt formats. Default is True.
            clumped (bool, optional): If True, save clumped data. Otherwise, save full data. Default is True.

        Raises:
            ValueError: If clumped data is requested but data is not clumped.
        """
        if clumped:
            if not hasattr(self, "data_clumped"):
                raise ValueError("clumped=True was used but the data is not clumped. "
                                 "Use .clump() first, or set clumped=False.")
            else:
                save_data(self.data_clumped, name=self.name+"_clumped", path=path, fmt=fmt, sep=sep, header=header)
        else:
            save_data(self.data, name=self.name, path=path, fmt=fmt, sep=sep, header=header)
        return
