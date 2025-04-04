import pandas as pd
import numpy as np
import warnings
import os
import copy
import psutil
import uuid
from functools import partial
from concurrent.futures import ProcessPoolExecutor
import scipy.stats as st
from plotnine import ggplot, aes, geom_point, geom_errorbarh, geom_errorbar, theme, element_text, geom_abline, labs, expand_limits, geom_vline


from .proxy import find_proxies, apply_proxies
from .MR_tools import query_outcome_func, MR_func, mrpresso_func
from .clump import clump_data_plink2
from .lift import lift_data
from .tools import create_tmp, load_reference_panel, setup_genetic_path, check_reference_panel
from .geno_tools import (
    save_data,
    check_arguments,
    adjust_column_names,
    check_int_column,
    fill_snpids_func,
    fill_coordinates_func,
    fill_nea,
    fill_ea_nea,
    check_beta_column,
    check_p_column,
    fill_se_p,
    check_allele_column,
    check_snp_column,
    remove_na,
    filter_by_gene_func
)
from .association import set_phenotype_func, association_test_func_plink2
from .extract_prs import extract_snps_func, prs_func
from .snp_query import async_query_gwas_catalog
from .constants import STANDARD_COLUMNS, REF_PANEL_COLUMNS, CHECKS_DICT, MR_METHODS_NAMES
from .colocalization import coloc_abf_func

# Do all the MR steps (query_outcome, harmonize etc) based on CHR/POS and not SNPs
# Add proxying function (input is df + searchspace (list of SNP or path to .bim, can be separated by chromosomes) and returns proxied df)
# Get proxies (simply return a list of proxies)
# Include proxying option to association_test
# Multivariable-MR
# Check stability with variants on sexual chromosomes
# Check the build of user data (potentially with a list of SNPs with different positions)
# update_snpids function: take alleles into account during the merge if they are present in the user data


class Geno:
    """
    A class to handle GWAS-derived data, including SNP rsID, genome position,
    SNP-trait effects, and effect allele frequencies.

    Attributes:
        data (pd.DataFrame): Main DataFrame containing SNP data.
        phenotype (pd.DataFrame, str): Tuple with a DataFrame of individual-level phenotype
            data and a string representing the phenotype trait column. Initialized after
            running the 'set_phenotype' method.
        MR_data (pd.DataFrame, pd.DataFrame, str): Tuple containing DataFrames for associations
            with exposure and outcome, and a string for the outcome name. Initialized after
            running the 'query_outcome' method.
        MR_results (pd.DataFrame, pd.DataFrame, str, str): Contains an MR results dataframe, a dataframe of harmonized SNPs, an exposure name, an outcome name. Assigned after calling the MR method and used for plotting with the MR_plot method.
        ram (int): Available memory.
        cpus (int): Number of available CPUs.
        checks (dict): Dictionary of checks performed on the main DataFrame.
        name (str): ID of the object (for internal reference and debugging purposes).
        reference_panel (pd.DataFrame): Reference population SNP data used for SNP info
            adjustments. Initialized when first needed.
        reference_panel_name (str): string to identify the reference_panel (path or population string)

    Methods:
        preprocess_data: Clean and preprocess the 'data' attribute (the main dataframe of SNP-level data).
        clump: Clump the main data based on reference panels and return a new Geno object with the clumped data.
        prs: Computes Polygenic Risk Score on genomic data.
        set_phenotype: Assigns a DataFrame with individual-level data and a phenotype trait to the 'phenotype' attribute.
        association_test: Computes SNP-trait effect estimates, standard errors, and p-values.
        query_outcome: Extracts SNPs from SNP-outcome association data and stores it in the 'MR_data' attribute.
        MR: Performs Mendelian Randomization between the SNP-exposure and SNP-outcome data stored in the 'MR_data' attribute. Stores the results in the 'MR_results' attribute.
        MR_plot: Plot the results of the MR analysis stored in the 'MR_results' attribute.
        MR_forest: Creates and returns a forest plot of MR results, with one row per method.
        MRpresso: Executes the MR-PRESSO algorithm for horizontal pleiotropy correction between the SNP-exposure and SNP-outcome data stored in the 'MR_data' attribute.
        lift: Lifts SNP data from one genomic build to another.
        query_gwas_catalog: Query the GWAS Catalog for SNP-trait associations.
    """

    def __init__(
        self,
        df,
        CHR="CHR",
        POS="POS",
        SNP="SNP",
        EA="EA",
        NEA="NEA",
        BETA="BETA",
        SE="SE",
        P="P",
        EAF="EAF",
        keep_columns=True,
    ):
        """
        Initializes the Geno object used to store and transform Single Nucleotide Polymorphisms (SNP) data.

        Args:
            df (pd.DataFrame): DataFrame where each row represents a SNP.
            CHR (str, optional): Column name for chromosome. Defaults to "CHR".
            POS (str, optional): Column name for genomic position. Defaults to "POS".
            SNP (str, optional): Column name for SNP identifier. Defaults to "SNP".
            EA (str, optional): Column name for effect allele. Defaults to "EA".
            NEA (str, optional): Column name for non-effect allele. Defaults to "NEA".
            BETA (str, optional): Column name for effect estimate. Defaults to "BETA".
            SE (str, optional): Column name for effect standard error. Defaults to "SE".
            P (str, optional): Column name for p-value. Defaults to "P".
            EAF (str, optional): Column name for effect allele frequency. Defaults to "EAF".
            keep_columns (bool, optional): Determines if non-main columns should be kept. Defaults to True.

        Attributes:
            name (str): Randomly generated ID for the Geno object.
            cpus (int): Number of CPUs to be used.
            ram (int): Amount of RAM to be used in MBs.
            checks (dict): Dictionary of checks performed on the main DataFrame.
            reference_panel (pd.DataFrame): Reference population SNP data used for SNP info
                adjustments. Initialized when first needed.
            reference_panel_name (str): string to identify the reference_panel (path or population string)
            phenotype (pd.DataFrame, str): Tuple with a DataFrame of individual-level phenotype
                data and a string representing the phenotype trait column. Initialized after
                running the 'set_phenotype' method.
            MR_data (pd.DataFrame, pd.DataFrame, str): Tuple containing DataFrames for associations
                with exposure and outcome, and a string for the outcome name. Initialized after
                running the 'query_outcome' method.
            MR_results (pd.DataFrame, pd.DataFrame, str, str): Contains an MR results dataframe, a dataframe of harmonized SNPs, an exposure name, an outcome name. 
                Assigned after calling the MR method and used for plotting with the MR_plot method.
            MRpresso_subset_data (pd.DataFrame, pd.DataFrame, str, str): Contains a dataframe of subsetted harmonized SNPs without outliers. 
                Assigned after calling the MRpresso method.
        """

        # Validate df type
        if not isinstance(df, pd.DataFrame):
            raise TypeError("df needs to be a pandas dataframe.")
        data = df.copy()

        # Standardize column names based on provided parameters +/- delete other columns
        data = adjust_column_names(
            data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns
        )

        # Set object attributes
        self.data = data
        self.name = str(uuid.uuid4())[:8]

        # List to keep track of checks performed
        self.checks = CHECKS_DICT.copy()

        # Set the maximal amount of ram/cpu to be used by the methods and dask chunksize
        self.cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", default=os.cpu_count()))
        non_hpc_ram_per_cpu = psutil.virtual_memory().total  / (
            1024**2 * self.cpus
        )
        ram_per_cpu = int(
            os.environ.get("SLURM_MEM_PER_CPU", default=non_hpc_ram_per_cpu)
        )
        self.ram = int(ram_per_cpu * self.cpus * 0.8)

        create_tmp()

        return

    def preprocess_data(
        self,
        preprocessing='Fill',
        reference_panel="37",
        effect_column=None,
        keep_indel=None,
        keep_dups=None,
        fill_snpids=None,
        fill_coordinates=None,
    ):
        """
        Clean and preprocess the main dataframe of Single Nucleotide Polymorphisms (SNP) data.

        Args:
            preprocessing (str, optional): Level of preprocessing to apply. Options include:
                - "None": The dataframe is not modified.
                - "Fill": Missing columns are added based on reference data and invalid values set to NaN, but no rows are deleted.
                - "Fill_delete": Missing columns are added, and rows with missing, duplicated, or invalid values are deleted.
                Defaults to 'Fill'.
            reference_panel (str or pd.DataFrame, optional): Reference panel for SNP adjustments. 
                Options include: "37", "38" which correspond to the 1000G phase 3 reference panels in build 37 and 38 respectively.
                Can also be a pd.DataFrame with the necessary columns (CHR, POS, SNP, A1, A2).
                Also accepts a path to a specific bed/bim/fam or pgen/pvar/psam panel.
            effect_column (str, optional): Specifies the type of effect column ("BETA" or "OR"). If None, the method tries to determine it. Odds Ratios will be log-transformed and the standard error adjusted. Defaults to None.
            keep_indel (bool, optional): Determines if insertions/deletions should be kept. If None, defers to preprocessing value. Defaults to None.
            keep_dups (bool, optional): Determines if rows with duplicate SNP IDs should be kept. If None, defers to preprocessing value. Defaults to None.
            fill_snpids (bool, optional): Decides if the SNP (rsID) column should be created or replaced based on CHR/POS columns and a reference genome. If None, defers to preprocessing value. Defaults to None.
            fill_coordinates (bool, optional): Decides if CHR and/or POS should be created or replaced based on SNP column and a reference genome. If None, defers to preprocessing value. Defaults to None.

        Note:
            If you pass a standard reference_panel name (e.g. "EUR_37"), it will be converted to "37".
        """

        data = self.data

        # Check arguments and solve arguments logic.
        keep_indel, keep_dups, fill_snpids, fill_coordinates = check_arguments(
            preprocessing,
            effect_column,
            fill_snpids,
            fill_coordinates,
            keep_indel,
            keep_dups,
        )

        # Ensure CHR and POS columns are integers if preprocessing is enabled
        for int_col in ["CHR", "POS"]:
            if int_col in data.columns and preprocessing in ['Fill', 'Fill_delete']:
                check_int_column(data, int_col)
                self.checks[int_col] = True

        # Fill missing SNP column from reference data if necessary
        should_fill_snpids = (
            ("CHR" in data.columns)
            and ("POS" in data.columns)
            and ("SNP" not in data.columns)
        ) or fill_snpids
        if should_fill_snpids and fill_snpids is not False:
            data = fill_snpids_func(data, self.get_reference_panel(reference_panel), keep_indel)
            # If EA and NEA were present, they have been preprocessed in fill_snpids_func
            if "EA" in data.columns and "NEA" in data.columns:
                self.checks["EA"] = True
                self.checks["NEA"] = True

        # Fill missing CHR/POS columns from reference data if necessary
        should_fill_coordinates = (
            (not ("CHR" in data.columns) or not ("POS" in data.columns))
            and ("SNP" in data.columns)
        ) or fill_coordinates
        if should_fill_coordinates and fill_coordinates is not False:
            data = fill_coordinates_func(
                data, self.get_reference_panel(reference_panel)
            )

        # Fill missing NEA column from reference data if necessary and preprocessing is enabled
        missing_nea_condition = (
            "CHR" in data.columns
            and "POS" in data.columns
            and "NEA" not in data.columns
            and "EA" in data.columns
        )
        if missing_nea_condition and preprocessing in ['Fill', 'Fill_delete']:
            check_allele_column(data, "EA", keep_indel)
            self.checks["EA"] = True
            data = fill_nea(data, self.get_reference_panel(reference_panel))

        # Fill missing EA and NEA columns from reference data if necessary and preprocessing is enabled
        missing_ea_nea_condition = (
            "CHR" in data.columns
            and "POS" in data.columns
            and "NEA" not in data.columns
            and "EA" not in data.columns
        )
        if missing_ea_nea_condition and preprocessing in ['Fill', 'Fill_delete']:
            data = fill_ea_nea(data, self.get_reference_panel(reference_panel))

        # Convert effect column to Beta estimates if present
        if "BETA" in data.columns:
            check_beta_column(data, effect_column, preprocessing)
            self.checks["BETA"] = True

        # Ensure P column contains valid values
        if "P" in data.columns and preprocessing in ['Fill', 'Fill_delete']:
            check_p_column(data)
            self.checks["P"] = True

        # Fill missing SE or P columns if necessary
        if preprocessing in ['Fill', 'Fill_delete']:
            fill_se_p(data)

        # Process allele columns
        for allele_col in ["EA", "NEA"]:
            check_allele_condition = (allele_col in data.columns) and (
                (preprocessing in ['Fill', 'Fill_delete']) or (not keep_indel)
            )
            if check_allele_condition and not self.checks[allele_col]:
                check_allele_column(data, allele_col, keep_indel)
                self.checks[allele_col] = True

        # Check for and handle duplicate SNPs if necessary
        if "SNP" in data.columns and not keep_dups:
            check_snp_column(data)
            self.checks["SNP"] = True

        # Warn if essential columns are missing
        missing_columns = [column for column in STANDARD_COLUMNS if column not in data.columns]
        if missing_columns:
            print(f"Warning: the data doesn't include the following columns: {', '.join(missing_columns)}. Some methods may need them.")

        # Remove missing values if preprocessing level is set to 'Fill_delete'
        if preprocessing == 'Fill_delete':
            remove_na(data)
            self.checks["NA_removal"] = True

        ## Reset index
        # self.data.reset_index(drop=True, inplace=True)

        self.data = data

    def get_reference_panel(self, reference_panel="37"):
        """
        Retrieve or set the reference panel for the Geno object.

        If the Geno object does not have a reference panel attribute set, this method will try to set it based on the provided `reference_panel` argument. 
        This can be either "37" or "38" for the provided reference panel (based on 1000G phase 3) 
        or a DataFrame with specific columns or a path to a .bim or .pvar file.

        Note:
            If you pass a standard reference_panel name (e.g. "EUR_37"), it will be converted to "37".

        Args:
            reference_panel (str or pd.DataFrame, optional): Either "37" or "38" for the provided reference panel (based on 1000G phase 3)
                or a DataFrame with necessary columns or a valid path to a .bim or .pvar file

        Returns:
            pd.DataFrame: The reference panel DataFrame for the Geno object.

        Raises:
            ValueError: If the provided DataFrame doesn't have the necessary columns.
        """
        # Check if the user provided a dataframe
        if isinstance(reference_panel, pd.DataFrame):
            # If the provided reference_panel is a DataFrame, verify its structure and dtypes
            for col in REF_PANEL_COLUMNS:
                if col not in reference_panel.columns:
                    raise ValueError(
                        f"The {col} column is not present in the reference_panel provided and is necessary."
                    )

            print(
                "Using the provided reference_panel dataframe as the reference panel."
            )
            # Check and standardize the reference panel
            reference_panel = check_reference_panel(reference_panel)
            self.reference_panel = reference_panel.copy()
            self.reference_panel_name = "USER_PROVIDED"
            
        # Else, check if there is already a reference_panel with the same ID. If not, load it based on provided string
        elif not (hasattr(self, "reference_panel") and 
                  hasattr(self, "reference_panel_name") and
                  self.reference_panel_name==reference_panel):
            self.reference_panel = load_reference_panel(reference_panel)
            self.reference_panel_name = reference_panel   

        return self.reference_panel

    def clump(self, kb=10000, r2=0.01, p1=5e-8, p2=0.01, reference_panel="EUR_37"):
        """
        Clump the data based on linkage disequilibrium and return another Geno object with the clumped data.
        The clumping process is executed using plink.

        Args:
            kb (int, optional): Clumping window in thousands of SNPs. Default is 10000 (very strict).
            r2 (float, optional): Linkage disequilibrium threshold, values between 0 and 1. Default is 0.01 (strict).
            p1 (float, optional): P-value threshold during clumping. SNPs with a P-value higher than this value are excluded. Default is 5e-8.
            p2 (float, optional): P-value threshold post-clumping to further filter the clumped SNPs. If p2 < p1, it won't be considered. Default is 0.01.
            reference_panel (str, optional): The reference population for linkage disequilibrium values. 
                Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
                Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
                Default is "EUR_37".

        Returns:
            genal.Geno: A new Geno object based on the clumped data.
        """
        
        # Ensure required columns exist in the data
        for column in ["SNP", "P"]:
            if column not in self.data.columns:
                raise ValueError(f"The column {column} is not found in the data")
                
        create_tmp() #Make sure temporary folder exists

        # Validate and process SNP and P columns, if not already done
        if not self.checks.get("SNP"):
            check_snp_column(self.data)
            self.checks["SNP"] = True

        if not self.checks.get("P"):
            check_p_column(self.data)
            self.checks["P"] = True

        initial_rows = self.data.shape[0]
        self.data.dropna(subset=["SNP", "P"], inplace=True)
        deleted_rows = initial_rows - self.data.shape[0]
        if deleted_rows > 0:
            print(
                f"{deleted_rows} ({deleted_rows/initial_rows*100:.3f}%) rows with NA values in columns SNP or P have been deleted."
            )

        # Create tmp directory if it doesn't exist
        create_tmp()

        # Clump the data using the specified parameters
        clumped_data = clump_data_plink2(
            self.data,
            reference_panel,
            kb,
            r2,
            p1,
            p2,
            self.name,
            self.ram,
        )

        # If clumped data is successfully generated, assign it to the object's attribute
        if clumped_data is not None:
            Clumped = self.copy(clumped_data)
            return Clumped
        return None

    def update_snpids(self, path=None, replace=False):
        """
        Update or create the column of SNP name based on genetic data and genomic position.
        
        Args:
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files.
                If files are split by chromosomes, replace the chromosome number with '$'.
                For instance: path = "ukb_chr$_file". Defaults to the path from the configuration.
            replace (bool, optional): To update the .data attribute with the updated SNP column or not.
                
        Returns:
            None: It updates the dataframe in the .data attribute.
            
        Note:
            This can be used before extracting SNPs from the genetic data if there is possibility of a mismatch between the SNP name contained in the Geno dataframe (SNP-level data) and the SNP name used in the genetic data (individual-level data). 
            Notably, this can avoid losing SNPs due to ID mismatch during polygenic risk scoring or single-SNP association testing.
        """
        
        # Check mandatory columns
        for column in ["CHR", "POS"]:
            if not (column in self.data.columns):
                raise ValueError(
                    f"The column {column} is not found in the data and is mandatory to update snpIDs!"
                )
        data = self.data.copy() # We don't want to modify the data attribute
        path, filetype = setup_genetic_path(path) #Verify the path and get filetype
        is_split = "$" in path #If data is split by chromosomes or not
        
        # Merge data with the variant info dataframe
        if not is_split:
            if filetype == "bed":
                variant_info = pd.read_csv(
                    path + ".bim", sep="\t", names=["CHR", "SNP", "F", "POS", "A1", "A2"]
                )
            else: # pgen
                variant_info = pd.read_csv(
                    path + ".pvar", sep="\t", comment="#", names=["CHR", "POS", "SNP", "A1", "A2"]
                )

            # Convert CHR to string and remove 'chr' prefix if present
            variant_info["CHR"] = variant_info["CHR"].astype(str).str.replace("^chr", "", regex=True)
            # Convert numeric values to int, leave X/Y as strings
            variant_info["CHR"] = variant_info["CHR"].apply(
                lambda x: int(float(x)) if x.replace(".", "").isdigit() else x
            )
            # Remove duplicates in the CHR/POS columns
            variant_info.drop_duplicates(subset=["CHR", "POS"], keep="first", inplace=True)
            variant_info.drop_duplicates(subset=["SNP"], keep="first", inplace=True)

            data = data.merge(
                variant_info[["CHR", "POS", "SNP"]], on=["CHR", "POS"], how="left", suffixes=('', '_new')
            )
        else:
            chr_dict = {k: v for k, v in data.groupby('CHR')} #Split the dataframe by chromosome
            partial_merge_command_parallel = partial(
                merge_command_parallel, path=path, filetype=filetype
            )  # Wrapper function
            with ProcessPoolExecutor() as executor: #Merge each dataframe subset
                results = list(executor.map(partial_merge_command_parallel, chr_dict.values()))
                data = pd.concat(results, ignore_index=True, axis=0) #And combine them again
                
        # Update the SNP column
        data['SNP'] = data['SNP_new'].fillna(data['SNP'])
        data.drop(columns = ['SNP_new'], inplace=True)
        
        if replace: self.data = data #Update attribute if replace argument
        return data
    
    def extract_snps(self, path=None):
        """
        Extract the list of SNPs of this Geno object from the genetic data provided.

        Args:
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files.
                If files are split by chromosomes, replace the chromosome number with '$'.
                For instance: path = "ukb_chr$_file". Default is None.

        Returns:
            None: The output is a pgen/pvar/psam triple in the tmp_GENAL folder
            with the format "{name}_extract_allchr" which includes the SNPs from the UKB.

        Notes:
            The provided path is saved to the config file. If this function is called again,
            you don't need to specify the path if you want to use the same genomic files.
        """

        create_tmp() #Make sure temporary folder exists
        
        # Extract the SNP list
        snp_list = self.data["SNP"]

        # Renaming to avoid conflicts with previous extraction
        self.name = str(uuid.uuid4())[:8]
        # Extract SNPs using the provided path and SNP list
        _ = extract_snps_func(snp_list, self.name, path)

        return

    def prs(self, name=None, 
            weighted=True, 
            path=None, 
            proxy=False,
            reference_panel="EUR_37",
            kb=5000,
            r2=0.8,
            window_snps=1000000,
            
           ):
        """
        Compute a Polygenic Risk Score (PRS) and save it as a CSV file in the current directory.

        Args:
            name (str, optional): Name or path of the saved PRS file.
            weighted (bool, optional): If True, performs a PRS weighted by the BETA column estimates.
                                       If False, performs an unweighted PRS. Default is True.
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files for PRS calculation.
                                  If files are split by chromosomes, replace the chromosome number
                                  with '$'. For instance: path = "ukb_chr$_file".
                                  If not provided, it will use the genetic path most recently used 
                                  (if any). Default is None.
            position (bool, optional): Use the genomic positions instead of the SNP names to find the
                                  SNPs in the genetic data (recommended).
            proxy (bool, optional): If true, proxies are searched. Default is True.
            reference_panel (str, optional): The reference population used to derive linkage disequilibrium values and find proxies (only if proxy=True). 
                Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
                Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
                Default is "EUR_37".
            kb (int, optional): Width of the genomic window to look for proxies. Default is 5000.
            r2 (float, optional): Minimum linkage disequilibrium value with the main SNP
                for a proxy to be included. Default is 0.8.
            window_snps (int, optional): Compute the LD value for SNPs that are not
                more than x SNPs away from the main SNP. Default is 1000000 (equivalent to infinity, plink default).

        Returns:
            pd.DataFrame: The computed PRS data.
        """
        
        path, filetype = setup_genetic_path(path) # Check path
        create_tmp() #Make sure temporary folder exists
        
        # Check for mandatory columns in data
        mandatory_cols = ["EA", "BETA"]
        for col in mandatory_cols:
            if col not in self.data.columns:
                raise ValueError(f"The column {col} is not found in the data!")
                
        # Based on column presents, run the PRS with SNP names or genomic positions (with preference for positions)
        if "CHR" in self.data.columns and "POS" in self.data.columns:
            print("CHR/POS columns present: SNPs searched based on genomic positions.")
            data_prs = self.update_snpids(path = path)
        elif "SNP" in self.data.columns:
            print("CHR/POS columns absent: SNPs searched based on SNP name.")
            data_prs = self.data.copy()
        else:
            raise ValueError("Either the SNP or the CHR/POS columns need to be present to run a PRS.")

        # Check SNP and EA columns
        if not self.checks.get("SNP"):
            check_snp_column(data_prs)
        if not self.checks.get("EA"):
            check_allele_column(data_prs, "EA", keep_indel=False)
        if not self.checks.get("BETA"):
            check_beta_column(data_prs, effect_column=None, preprocessing='Fill_delete')

        initial_rows = data_prs.shape[0]
        data_prs.dropna(subset=["SNP", "EA", "BETA"], inplace=True)
        deleted_rows = initial_rows - data_prs.shape[0]
        if deleted_rows > 0:
            print(
                f"{deleted_rows} ({deleted_rows/initial_rows*100:.3f}%) rows with NA values in columns SNP, EA, or BETA have been deleted."
            )
            
        # If proxy option
        if proxy:
            print("Identifying the SNPs present in the genetic data...")
            # Obtain the list of SNPs present in the genetic data
            if path.count("$") == 1: #If split: merge all SNP columns of the variant files in parallel
                genetic_snp_list = []
                for i in range(1,23):
                    path_i = path.replace("$", str(i))
                    if filetype == "bed":
                        variant_file = path_i + ".bim"
                        if os.path.exists(variant_file):
                            variants = pd.read_csv(
                                variant_file, sep="\t", 
                                names=["CHR", "SNP", "F", "POS", "A1", "A2"], 
                                usecols=["SNP"]
                            )
                            genetic_snp_list.extend(variants.SNP.tolist())
                    else: # filetype == "pgen"
                        variant_file = path_i + ".pvar"
                        if os.path.exists(variant_file):
                            variants = pd.read_csv(
                                variant_file, sep="\t", comment="#",
                                names=["CHR", "POS", "SNP", "A1", "A2"],
                                usecols=["SNP"]
                            )
                            genetic_snp_list.extend(variants.SNP.tolist())
            else: #If not split
                if filetype == "bed":
                    variant_file = path + ".bim"
                    if os.path.exists(variant_file):
                        variants = pd.read_csv(
                            variant_file, sep="\t",
                            names=["CHR", "SNP", "F", "POS", "A1", "A2"]
                        )
                    else:
                        raise ValueError(f"The file {variant_file} does not exist.")
                else: # filetype == "pgen"
                    variant_file = path + ".pvar"
                    if os.path.exists(variant_file):
                        variants = pd.read_csv(
                            variant_file, sep="\t", comment="#",
                            names=["CHR", "POS", "SNP", "A1", "A2"]
                        )
                    else:
                        raise ValueError(f"The file {variant_file} does not exist.")
                genetic_snp_list = variants.SNP.tolist()
                
            # Identify the SNPs already present in the genetic data
            genetic_snps = set(genetic_snp_list)
            exposure_snps = set(data_prs.SNP.values)
            snps_present = len(exposure_snps & genetic_snps)
            snps_absent = exposure_snps - genetic_snps

            print(
                f"{snps_present} SNPs out of {len(exposure_snps)} are present in the genetic data."
            )

            # Search proxies only if there are absent SNPs
            if snps_absent:
                print(f"Searching proxies for {len(snps_absent)} SNPs...")
                ld = find_proxies(
                    snps_absent,
                    reference_panel=reference_panel,
                    kb=kb,
                    r2=r2,
                    window_snps=window_snps,
                    threads=self.cpus,
                    name=self.name
                )
                # Apply proxies if found
                if isinstance(ld, pd.DataFrame) and not ld.empty:
                    data_prs = apply_proxies(data_prs, ld, searchspace=genetic_snps)
                    # Remove duplicates in the SNP column
                    duplicate_indices = data_prs[data_prs.duplicated(subset=["SNP"], keep="first")].index
                    n_del = len(duplicate_indices)
                    if n_del > 0:
                        data_prs.drop(index=duplicate_indices, inplace=True)
                        print(f"After proxying, {n_del} SNPs had duplicated IDs and were removed.")
        
        # Renaming to avoid conflicts with previous extraction
        self.name = str(uuid.uuid4())[:8]
        # Compute PRS
        prs_data = prs_func(data_prs, weighted, path, ram=self.ram, name=self.name)

        # Save the computed PRS data as a CSV file
        name = self.name if not name else name
        prs_filename = os.path.splitext(name)[0] + ".csv"
        prs_data.to_csv(prs_filename, index=False, header=True)
        print(f"PRS data saved to {prs_filename}")

        return

    def set_phenotype(
        self, data, IID=None, PHENO=None, PHENO_type=None, FID=None, alternate_control=False
    ):
        """
        Assign a phenotype dataframe to the .phenotype attribute.

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
            FID (str, optional): Name of the family ID column in 'data'. If not provided,
                                FID values will be set to the same as IID values.
            alternate_control (bool, optional): By default, the function assumes that for a binary
                                                trait, the controls have the most frequent value.
                                                Set to True if this is not the case. Default is False.
        """

        processed_data, inferred_pheno_type = set_phenotype_func(
            data, PHENO, PHENO_type, IID, FID, alternate_control
        )

        # Assign the processed data and inferred phenotype type to the .phenotype attribute
        self.phenotype = (processed_data, inferred_pheno_type, PHENO)

        return

    def association_test(self, path=None, covar=[], standardize=True):
        """
        Conduct single-SNP association tests against a phenotype.

        Args:
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files.
                If files are split by chromosomes, replace the chromosome number with '$'.
                For instance: path = "ukb_chr$_file". Default is None.
            covar (list, optional): List of columns in the phenotype dataframe to be used
                                    as covariates in the association tests. Default is an empty list.
            standardize (bool, optional): If True, it will standardize a quantitative phenotype
                                          before performing association tests. This is typically done
                                          to make results more interpretable. Default is True.

        Returns:
            None: Updates the BETA, SE, and P columns of the data attribute based on the results
                  of the association tests.

        Note:
            This method requires the phenotype to be set using the set_phenotype() function.
        """

        # Ensure that the phenotype has been set using set_phenotype
        if not hasattr(self, "phenotype"):
            raise ValueError(
                "You first need to set a phenotype using .set_phenotype(data, PHENO, PHENO_type, IID)!"
            )
        
        # Ensure that covar does not contain the PHENO column
        if self.phenotype[2] in covar:
            print(f"The phenotype column {self.phenotype[2]} is also present in the covariates list and will be removed from the association tests.")
            covar.remove(self.phenotype[2])

        create_tmp() #Make sure temporary folder exists
        
        n_original = self.data.shape[0]
        # Based on column presents, extract the SNPs based on names or genomic positions (with preference for positions)
        if "CHR" in self.data.columns and "POS" in self.data.columns:
            print("CHR/POS columns present: SNPs searched based on genomic positions.")
            data = self.update_snpids(path = path)
        elif "SNP" in self.data.columns:
            print("CHR/POS columns absent: SNPs searched based on SNP name.")
            data = self.data
        else:
            raise ValueError("Either the SNP or the CHR/POS columns need to be present to identify SNPs in genetic data.")
        
        # Extract the SNP list
        snp_list = data["SNP"]

        # Renaming to avoid conflicts with previous extraction
        self.name = str(uuid.uuid4())[:8]
        # Extract SNPs using the provided path and SNP list
        path = extract_snps_func(snp_list, self.name, path)
        if path == "FAILED":
            raise ValueError("No SNPs were extracted from the genetic data and the association test can't be run.")
        
        # Perform the association test
        updated_data = association_test_func_plink2(
            data,
            covar,
            standardize,
            self.name,
            self.phenotype[0],
            self.phenotype[1],
        )

        n_updated = updated_data.shape[0]
        print(f"The BETA, SE, P columns of the .data attribute have been updated for {n_updated} SNPs.")
        if n_updated < n_original:
            print(f"{n_original - n_updated}({(n_original - n_updated)/n_original*100:.3f}%) SNPs have been removed.")

        # Update the instance data
        self.data = updated_data
        return

    def query_outcome(
        self,
        outcome,
        name=None,
        proxy=True,
        reference_panel="EUR_37",
        kb=5000,
        r2=0.8,
        window_snps=1000000,
    ):
        """
        Prepares dataframes required for Mendelian Randomization (MR) with the SNP information in `data` as exposure.

        Queries the outcome data, with or without proxying, and assigns a tuple to
        the outcome attribute: (exposure_data, outcome_data, name) ready for MR methods.

        Args:
            outcome: Can be a Geno object (from a GWAS) or a filepath of types: .h5 or .hdf5 (created with the :meth:`Geno.save` method.
            name (str, optional): Name for the outcome data. Defaults to None.
            proxy (bool, optional): If true, proxies are searched. Default is True.
            reference_panel (str, optional): The reference population to get linkage disequilibrium values and find proxies (only if proxy=True). 
                Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
                Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
                Default is "EUR_37".
            kb (int, optional): Width of the genomic window to look for proxies. Default is 5000.
            r2 (float, optional): Minimum linkage disequilibrium value with the main SNP
                for a proxy to be included. Default is 0.8.
            window_snps (int, optional): Compute the LD value for SNPs that are not
                more than x SNPs away from the main SNP. Default is 1000000 (equivalent to infinity, plink default).

        Returns:
            None: Sets the `MR_data` attribute for the instance.
        """

        exposure, outcome_data, outcome_name = query_outcome_func(
            self.data,
            outcome,
            name,
            proxy,
            reference_panel,
            kb,
            r2,
            window_snps,
            self.cpus,
        )

        # Assign the processed data to the MR_data attribute
        self.MR_data = [exposure, outcome_data, outcome_name]
        return

    def MR(
        self,
        methods=["IVW","IVW-FE","WM","Simple-mode","Egger"],
        action=2,
        eaf_threshold=0.42,
        heterogeneity=False,
        nboot=1000,
        penk=20,
        phi=1,
        exposure_name=None,
        outcome_name=None,
        cpus=-1,
        odds=False,
        use_mrpresso_data=False
    ):
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
                "Simple-mode": simple mode method
                "Weighted-mode": weighted mode method
                Default is ["IVW","IVW-FE","WM","Simple-mode","Weighted-mode","Egger"].
            action (int, optional): How to treat palindromes during harmonizing between
                exposure and outcome data. Accepts:
                1: Doesn't flip them (Assumes all alleles are on the forward strand)
                2: Uses allele frequencies to attempt to flip (conservative, default)
                3: Removes all palindromic SNPs (very conservative)
            eaf_threshold (float, optional): Max effect allele frequency accepted when
                flipping palindromic SNPs (relevant if action=2). Default is 0.42.
            heterogeneity (bool, optional): If True, includes heterogeneity tests in the results (Cochran's Q test).Default is False.
            nboot (int, optional): Number of bootstrap replications for methods with bootstrapping. Default is 1000.
            penk (int, optional): Penalty value for the WM-pen method. Default is 20.
            phi (int, optional): Factor for the bandwidth parameter used in the kernel density estimation of the mode methods
            exposure_name (str, optional): Name of the exposure data (only for display purposes).
            outcome_name (str, optional): Name of the outcome data (only for display purposes).
            odds (bool, optional): If True, adds an odds ratio column with 95% confidence intervals. Default is False.
            use_mrpresso_data (bool, optional): If True and MRpresso has been run, uses the subset of instruments after outlier removal. Default is False.

        Returns:
            pd.DataFrame: A table with MR results.
        """

        # Ensure that query_outcome has been previously called
        if not hasattr(self, "MR_data"):
            raise ValueError("You must first call query_outcome() before running MR.")

        # Get MRpresso subset data if requested and available
        subset_data = None
        if use_mrpresso_data:
            if not hasattr(self, "MRpresso_subset_data"):
                raise ValueError("Use_mrpresso_data is set to True but MRpresso results not found. Please run MRpresso first.")
            elif self.MRpresso_subset_data is None:
                print("Use_mrpresso_data is set to True but MRpresso did not identify any outliers. Using all instruments.")
            else:
                subset_data = self.MRpresso_subset_data
                print(f"Using data for {len(subset_data)} instruments after MR-PRESSO outlier removal. The action is the one used in MR-PRESSO.")


        if outcome_name:
            self.MR_data[2] = outcome_name
        exp_name = exposure_name if exposure_name else self.name
        cpus = self.cpus if cpus == -1 else cpus
        res, df_mr = MR_func(
            self.MR_data,
            methods,
            action,
            eaf_threshold,
            nboot,
            penk,
            phi,
            exp_name,
            self.cpus,
            subset_data
        )

        if not heterogeneity:
            res = res.loc[:,["exposure", "outcome", "method", "nSNP", "b", "se", "pval"]]

        if odds and not res.empty:
            # Calculate odds ratios and confidence intervals using .loc
            res.loc[:,'OR_95CI'] = res.apply(lambda row: 
                f"{np.exp(row['b']):.3f} ({np.exp(row['b'] - 1.96*row['se']):.3f}-{np.exp(row['b'] + 1.96*row['se']):.3f})" 
                if not pd.isna(row['b']) and not pd.isna(row['se']) 
                else np.nan, axis=1)

        self.MR_results = (res, df_mr, exposure_name, outcome_name)

        return res
    
    def MR_plot(
        self,
        methods=[
            "IVW",
            "WM",
            "Simple-mode",
            "Egger",
        ],
        exposure_name=None,
        outcome_name=None,
        filename=None
    ):
        """
        Creates and returns a scatter plot of individual SNP effects with lines representing different Mendelian Randomization (MR) methods. Each MR method specified in the 'methods' argument is represented as a line in the plot.

        Args:
            methods (list of str, optional): A list of MR methods to be included in the plot. Default methods are "IVW", "WM", "Simple-median", and "Egger".
            exposure_name (str, optional): A custom label for the exposure effect axis. If None, uses the label provided in the MR function call or a default label.
            outcome_name (str, optional): A custom label for the outcome effect axis. If None, uses the label provided in the MR function call or a default label.
            filename (str, optional): The filename where the plot will be saved. If None, the plot is not saved.

        Returns:
            plotnine.ggplot.ggplot: A plotnine ggplot object representing the scatter plot of individual SNP effects with MR method lines.

        Raises:
            ValueError: If MR analysis has not been performed prior to calling this function.

        Note:
            This function requires prior execution of the `MR` method to compute MR results. Make sure the MR analysis is performed on the data before calling `MR_plot`.
        """
        if not hasattr(self, "MR_results"):
            raise ValueError("You need to run an MR analysis with the MR method before calling the MR_plot function.")
        
        ## Extract the previously computed MR results
        df_mr = self.MR_results[1]
        res = self.MR_results[0]
        exposure_name = self.MR_results[2] if not exposure_name else exposure_name
        exposure_name = "Effect on the exposure" if not exposure_name else f"Effect on {exposure_name}"
        outcome_name = self.MR_results[3] if not outcome_name else outcome_name
        outcome_name = "Effect on the outcome" if not outcome_name else f"Effect on {outcome_name}"
        
        ## Switch all exposure betas to >= 0
        df_mr['BETA_e'], df_mr['BETA_o'] = np.where(df_mr['BETA_e'] < 0, (-df_mr['BETA_e'], -df_mr['BETA_o']), (df_mr['BETA_e'], df_mr['BETA_o']))

        ## Create the scatter plot with error bars
        plot = (
            ggplot(df_mr, aes('BETA_e', 'BETA_o'))

            + geom_errorbarh(aes(xmin='BETA_e-SE_e', xmax='BETA_e+SE_e'), height=0, color="gray", size=0.1) 
            + geom_errorbar(aes(ymin='BETA_o-SE_o', ymax='BETA_o+SE_o'), width=0, color="gray", size=0.1)
            + geom_point(color='black', size=0.2) 
            + geom_abline(slope=0, intercept=0, color='black')
            + labs(x=exposure_name, y=outcome_name) 
            + theme(
                axis_title=element_text(size=12),
                axis_text=element_text(size=10),
                figure_size=(10,6)
            )
            + expand_limits(x=0)
        )
        
        ## Add the lines corresponding to the specified MR methods (if present in the computation)
        lines = []
        for method in methods:
            if method not in MR_METHODS_NAMES.keys():
                warnings.warn(f"{method} is not an appropriate MR method. MR methods can be IVW, WM, Egger... Please refer to the documentation for more.")
                continue
            ## If not an Egger method: only need the slope
            if not method.startswith("Egger"):
                method_name = MR_METHODS_NAMES[method]
                res_row = res[res.method == method_name]
                if res_row.shape[0] == 0:
                    warnings.warn(f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot.")
                elif res_row.shape[0] == 1:
                    lines.append({
                        'slope': res_row["b"].values[0], 
                        'intercept': 0, 
                        'MR Methods': method_name  # Use method_name as the color label
                    })
            ## For Egger methods: need to get the slope and the intercept
            else:
                method_name = MR_METHODS_NAMES[method][0]
                method_name_intercept = MR_METHODS_NAMES[method][1]
                res_row = res[res.method == method_name]
                res_row_intercept = res[res.method == method_name_intercept]
                if res_row.shape[0] == 0:
                    warnings.warn(f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot.")
                elif res_row.shape[0] == 1 and res_row_intercept.shape[0] == 1:
                    lines.append({
                        'slope': res_row["b"].values[0], 
                        'intercept': res_row_intercept["b"].values[0], 
                        'MR Methods': method_name  # Use method_name as the color label
                    })
        line_data = pd.DataFrame(lines)
        plot += geom_abline(aes(slope='slope', intercept='intercept', color='MR Methods'), data=line_data)
        
        ## Save plot if filename is specified
        if filename:
            plot.save(f"{filename}.png", dpi=500, width=10, height=6, verbose=False)
        
        return plot
    
    def MR_forest(
        self,
        methods=[
            "IVW",
            "WM", 
            "Simple-mode",
            "Egger",
        ],
        exposure_name=None,
        outcome_name=None,
        odds=False,
        filename=None
    ):
        """
        Creates and returns a forest plot of MR results, with one row per method.

        Args:
            methods (list of str, optional): A list of MR methods to be included in the plot. 
                Default methods are "IVW", "WM", "Simple-mode", and "Egger".
            exposure_name (str, optional): A custom label for the exposure. If None, uses 
                the label provided in the MR function call or a default label.
            outcome_name (str, optional): A custom label for the outcome. If None, uses 
                the label provided in the MR function call or a default label.
            odds (bool, optional): If True, plots odds ratios instead of betas. Default is False.
            filename (str, optional): The filename where the plot will be saved. If None, 
                the plot is not saved.

        Returns:
            plotnine.ggplot.ggplot: A plotnine ggplot object representing the forest plot 
                of MR results.

        Raises:
            ValueError: If MR analysis has not been performed prior to calling this function.
        """
        if not hasattr(self, "MR_results"):
            raise ValueError("You need to run an MR analysis with the MR method before calling the MR_forest function.")
        
        # Extract the previously computed MR results
        res = self.MR_results[0]
        exposure_name = self.MR_results[2] if not exposure_name else exposure_name
        exposure_name = "Exposure" if not exposure_name else exposure_name
        outcome_name = self.MR_results[3] if not outcome_name else outcome_name
        outcome_name = "Outcome" if not outcome_name else outcome_name
        
        # Create plotting dataframe
        plot_data = []
        for method in methods:
            if method not in MR_METHODS_NAMES.keys():
                warnings.warn(f"{method} is not an appropriate MR method. MR methods can be IVW, WM, Egger... Please refer to the documentation for more.")
                continue
                
            # Get the method name, beta estimate, standard error, and number of SNPs
            method_name = MR_METHODS_NAMES[method]
            
            # Handle Egger method which has tuple of names
            if isinstance(method_name, tuple):
                method_name = method_name[0]  # Use first element for Egger regression
                
            res_row = res[res.method == method_name]
            if res_row.shape[0] == 0:
                warnings.warn(f"The {method_name} ({method}) method was not included in the MR method call and will be excluded from the plot.")
            elif res_row.shape[0] == 1:
                plot_data.append({
                    'Method': method_name,
                    'Estimate': res_row["b"].values[0],
                    'SE': res_row["se"].values[0], 
                    'nSNP': res_row["nSNP"].values[0]
                })
        
        plot_df = pd.DataFrame(plot_data)
        
        # Calculate confidence intervals
        if odds:
            plot_df['CI_lower'] = np.exp(plot_df['Estimate'] - 1.96 * plot_df['SE'])
            plot_df['CI_upper'] = np.exp(plot_df['Estimate'] + 1.96 * plot_df['SE'])
            plot_df['Estimate'] = np.exp(plot_df['Estimate'])
            x_label = f"Odds Ratio for {outcome_name} per unit increase in {exposure_name}"
        else:
            plot_df['CI_lower'] = plot_df['Estimate'] - 1.96 * plot_df['SE']
            plot_df['CI_upper'] = plot_df['Estimate'] + 1.96 * plot_df['SE']
            x_label = f"Effect on {outcome_name} per unit increase in {exposure_name}"
        
        # Convert Method to categorical with reversed order to preserve row order in plot
        plot_df['Method'] = pd.Categorical(plot_df['Method'], categories=plot_df['Method'].values[::-1], ordered=True)
        
        # Create the forest plot
        plot = (
            ggplot(plot_df, aes(x='Estimate', y='Method'))
            + geom_point(size=3)
            + geom_errorbarh(aes(xmin='CI_lower', xmax='CI_upper'), height=0.2)
            + theme(
                axis_title=element_text(size=12),
                axis_text=element_text(size=10),
                figure_size=(10,6)
            )
            + labs(x=x_label, y="")
        )
        
        # Add vertical line at null effect
        if odds:
            plot += geom_vline(xintercept=1, linetype='dashed', color='gray')
        else:
            plot += geom_vline(xintercept=0, linetype='dashed', color='gray')
        
        # Save plot if filename is specified
        if filename:
            plot.save(f"{filename}.png", dpi=500, width=10, height=6, verbose=False)
        
        return plot

    def MRpresso(
        self,
        action=2,
        eaf_threshold=0.42,
        n_iterations=10000,
        outlier_test=True,
        distortion_test=True,
        significance_p=0.05,
        cpus=-1,
    ):
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
            tuple: Contains the following elements:
                - mod_table: DataFrame containing the original (before outlier removal)
                             and outlier-corrected (after outlier removal) inverse variance-weighted MR results.
                - GlobalTest: p-value of the global MR-PRESSO test indicating the presence of horizontal pleiotropy.
                - OutlierTest: DataFrame assigning a p-value to each SNP representing the likelihood of this
                               SNP being responsible for the global pleiotropy. Set to NaN if global test p_value > significance_p.
                - DistortionTest: p-value for the distortion test.
        """

        if not hasattr(self, "MR_data"):
            raise ValueError("You must first call query_outcome() before running MR.")
        cpus = self.cpus if cpus == -1 else cpus

        mod_table, GlobalTest, OutlierTest, BiasTest, subset_data = mrpresso_func(
            self.MR_data,
            action,
            eaf_threshold,
            n_iterations,
            outlier_test,
            distortion_test,
            significance_p,
            cpus,
        )
        if subset_data is not None:
            subset_data.drop(columns=["Weights"], inplace=True, errors='ignore')
        self.MRpresso_subset_data = subset_data
        self.MRpresso_results = mod_table, GlobalTest, OutlierTest, BiasTest

        return mod_table, GlobalTest, OutlierTest, BiasTest

    def filter_by_gene(self, gene_id, id_type="symbol", window_size=1000000, build="37", replace=False):
        """
        Filter the data to include only variants that are within a specified distance of a specific gene.

        Args:
            gene_id (str): Identifier for the gene/protein to filter variants around.
            id_type (str, optional): Type of identifier provided. Options are:
                - "symbol": Gene symbol (e.g., "APOE")
                - "HGNC": HGNC ID (e.g., "HGNC:613")
                - "name": Full gene name (e.g., "apolipoprotein E")
                - "Ensembl": Ensembl gene ID (e.g., "ENSG00000130203")
                - "NCBI": NCBI gene ID (e.g., "348")
                - "UCSC": UCSC gene ID (e.g., "uc001hbu.2")
                - "Vega": Vega gene ID (e.g., "OTTHUMG00000019505")
                Default is "symbol".
            window_size (int, optional): Size of the window around the gene in base pairs. Default is 1,000,000 (1Mb).
            build (str, optional): Genome build of the data. Default is "37".
            replace (bool, optional): If True, replace the existing data attribute with the filtered data. Default is True.
        Returns:
            if replace is True:
                pd.DataFrame: Filtered DataFrame containing only variants within the specified window 
                    around the gene, with additional column 'Distance'.
            if replace is False:
                genal.Geno: A new Geno object with the filtered data.

        Raises:
            ValueError: If required columns are missing, gene information cannot be found, or invalid id_type is provided.
            
        Notes:
            - Distance is calculated from the nearest gene boundary (start or end position)
            - Null distances indicate the variant is within the gene
        """
        # Check required columns
        for col in ["CHR", "POS"]:
            if col not in self.data.columns:
                raise ValueError(f"Column {col} must be present in the input data!")
            
        # Do the appropriate preprocessing on CHR and POS columns if not already done
        if not self.checks.get("CHR"):
            check_int_column(self.data, "CHR")
            self.checks["CHR"] = True
        if not self.checks.get("POS"):
            check_int_column(self.data, "POS")
            self.checks["POS"] = True

        filtered = filter_by_gene_func(self.data, gene_id, id_type, window_size, build)
        
        if replace:
            self.data = filtered
        else:
            Geno_filtered = self.copy(filtered)
            return Geno_filtered

    def colocalize(self, outcome, method="abf", trait1_type=None, trait2_type=None,
                   sdY1=None, sdY2=None, n1=None, n2=None, p1=1e-4, p2=1e-4, p12=1e-5, merge_on_snp=False):
        """
        Perform colocalization analysis between two GWAS datasets.
        
        Args:
            outcome: Another Geno object containing the outcome dataset
            method: Method to use for colocalization (default: "abf")
            trait1_type: Type of exposure trait ("quant" for quantitative traits or "cc" for case-control traits)
            trait2_type: Type of outcome trait ("quant" for quantitative traits or "cc" for case-control traits)
            sdY1: Standard deviation of exposure trait (required for quantitative traits, but can be estimated from EAF and sample size)
            sdY2: Standard deviation of outcome trait (required for quantitative traits, but can be estimated from EAF and sample size)
            n1: Sample size for exposure (used to estimate sdY1 if sdY1 is not provided)
            n2: Sample size for outcome (used to estimate sdY2 if sdY2 is not provided)
            p1: Prior probability SNP associated with exposure
            p2: Prior probability SNP associated with outcome
            p12: Prior probability SNP associated with both traits
            merge_on_snp: If True, merge the datasets on SNP column. If False, first attempt to merge on CHR and POS columns.
        """
        # Ensure required columns exist in both datasets
        required_cols = ['BETA', 'SE']
        for col in required_cols:
            if col not in self.data.columns:
                raise ValueError(f"Column {col} must be present in exposure dataset")
            if col not in outcome.data.columns:
                raise ValueError(f"Column {col} must be present in outcome dataset")
            
        if trait1_type is None:
            print("trait1_type not specified. Assuming trait 1 is a quantitative trait.")
            trait1_type = "quant"
        if trait2_type is None:
            print("trait2_type not specified. Assuming trait 2 is a quantitative trait.")
            trait2_type = "quant"
                
        # Make copies of the data to avoid modifying the original data
        data1 = self.data.copy()
        data2 = outcome.data.copy()
            
        # Call the implementation function
        return coloc_abf_func(data1, 
                             data2,
                             trait1_type=trait1_type,
                             trait2_type=trait2_type,
                             sdY1=sdY1, 
                             sdY2=sdY2,
                             n1=n1,
                             n2=n2,
                             p1=p1,
                             p2=p2,
                             p12=p12,
                             merge_on_snp=merge_on_snp)


    def lift(
        self,
        start="hg19",
        end="hg38",
        replace=False,
        extraction_file=False,
        chain_file=None,
        name=None,
        liftover_path=None,
    ):
        """
        Perform a liftover from one genetic build to another.

        Args:
            start (str, optional): Current build of the data. Default is "hg19".
            end (str, optional): Target build for the liftover. Default is "hg38".
            replace (bool, optional): If True, updates the data attribute in place. Default is False.
            extraction_file (bool, optional): If True, prints a CHR POS SNP space-delimited
                file. Default is False.
            chain_file (str, optional): Path to a local chain file for the lift.
                If provided, `start` and `end` arguments are not considered. Default is None.
            name (str, optional): Filename or filepath (without extension) to save the lifted dataframe.
                If not provided, the data is not saved.
            liftover_path (str, optional): Specify the path to the USCS liftover executable. If not provided, the lift will be done in python (slower for large amount of SNPs).

        Returns:
            pd.DataFrame: Data after being lifted.
        """
        # Ensure mandatory columns are present in the input data
        for column in ["CHR", "POS"]:
            if column not in self.data.columns:
                raise ValueError(f"The column {column} is not found in the data!")

        create_tmp()  # Create tmp folder if does not exist

        # Select appropriate data or copy of data depending on replace argument
        if not replace:
            data = self.data.copy()
        else:
            data = self.data

        # Do the appropriate preprocessing on CHR and POS columns if not already done
        if not self.checks.get("CHR"):
            check_int_column(data, "CHR")
        if not self.checks.get("POS"):
            check_int_column(data, "POS")

        # Update the checks if replace = True
        if replace:
            self.checks["CHR"] = True
            self.checks["POS"] = True

        print(
            f"Lifting the data{' inplace' if replace else ''}. "
            f"The .data attribute will {'' if replace else 'not '}be modified. "
            f"Use replace={'False' if replace else 'True'} to {'leave it as is' if replace else 'lift inplace'}."
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = lift_data(
                data,
                start=start,
                end=end,
                extraction_file=extraction_file,
                chain_file=chain_file,
                name=name,
                liftover_path=liftover_path,
                object_id=self.name,
            )

        return data
    
    def query_gwas_catalog(
        self,
        p_threshold=5e-8,
        return_p=False,
        return_study=False,
        replace=True,
        max_associations=None,
        timeout=-1,
    ):
        """
        Queries the GWAS Catalog Rest API and add an "ASSOC" column containing associated traits for each SNP.
        
        Args:
            p_threshold (float, optional): Only associations that are at least as significant are reported. Default is 5e-8.
            return_p (bool, optional): If True, include the p-value in the results. Default is False.
            return_study (bool, optional): If True, include the ID of the study from which the association is derived in the results. Default is False.
            replace (bool, optional): If True, updates the data attribute in place. Default is True.
            max_associations (int, optional): If not None, only the first `max_associations` associations are reported for each SNP. Default is None.
            timeout (int, optional): Timeout for each query in seconds. Default is -1 (custom timeout based on number of SNPs to query). Choose None for no timeout.

        Returns:
            pd.DataFrame: Data attribute with an additional column "ASSOC".
                The elements of this column are lists of strings or tuples depending on the `return_p` and `return_study` flags. If the SNP could not be queried, the value is set to "FAILED_QUERY".
        """
        # Ensure mandatory column is present in the input data
        if "SNP" not in self.data.columns:
            raise ValueError(f"The SNP column is necessary for the GWAS query!")
            
        # Select appropriate data or copy of data depending on replace argument
        if not replace:
            data = self.data.copy()
        else:
            data = self.data
            
        print(
            f"Querying the GWAS Catalog and creating the ASSOC column. \n"
            f"Only associations with a p-value <= {p_threshold} are reported. Use the p_threshold argument to change the threshold."
        )
        if max_associations:
            print(f"Reporting the first {max_associations} associations for each SNP.")
        if not return_p:
            print(f"To report the p-value for each association, use return_p=True.")
        if not return_study :
            print(f"To report the study ID for each association, use return_study=True.")
        print(
            f"The .data attribute will {'be' if replace else 'not be'} modified. "
            f"{'Use replace=False to leave it as is.' if replace else ''}"
            )
        
        # Estimate a reasonable timeout given the number of SNPs to query (45 SNPs per second)
        timeout = max(len(data) / 45, 30) if timeout == -1 else timeout

        # Call the async function to query all SNPs
        results_snps, errors, timeouts = async_query_gwas_catalog(
            data.SNP.to_list(), 
            p_threshold=p_threshold, 
            return_p=return_p, 
            return_study=return_study,
            max_associations=max_associations,
            timeout=timeout,
        )
        
        # Create the column
        data["ASSOC"] = data['SNP'].map(results_snps).fillna("FAILED_QUERY")
        data.loc[data["SNP"].isin(timeouts), "ASSOC"] = "TIMEOUT"
        
        print("The ASSOC column has been successfully created.")
        print(f"{len(errors)} ({len(errors)/len(data)*100:.2f}%) SNPs failed to query (not found in GWAS Catalog) and {len(timeouts)} ({len(timeouts)/len(data)*100:.1f}%) SNPs timed out after {timeout:.2f} seconds." 
              f" You can increase the timeout value with the timeout argument.")
            
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

        self.data["BETA"] = (self.data.BETA - np.mean(self.data.BETA)) / np.std(
            self.data.BETA
        )
        self.data["SE"] = np.abs(self.data.BETA / st.norm.ppf(self.data.P / 2))
        print(
            "The Beta column has been standardized and the SE column has been adjusted."
        )

    def sort_group(self, method="lowest_p"):
        """
        Handle duplicate SNPs. Useful if the instance combines different Genos.

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

    def copy(self, data):
        """
        Create another Geno instance with the updated data attribute.
        The relevant attributes are copied as well (checks, phenotype, reference_panel, reference_panel_name).
        Attributes that are not copied are MR_data, MR_results, MRpresso_subset_data, MRpresso_results.

        Returns:
            Geno: A deep copy of the instance.
        """
        Geno_copy = Geno(data, keep_columns=True)
        Geno_copy.checks = self.checks.copy()
        if hasattr(self, "phenotype"):
            Geno_copy.phenotype = self.phenotype
            if hasattr(self, "reference_panel"):
                Geno_copy.reference_panel = self.reference_panel
        if hasattr(self, "reference_panel_name"):
            Geno_copy.reference_panel_name = self.reference_panel_name
        return Geno_copy

    def save(self, path="", fmt="h5", sep="\t", header=True):
        """
        Save the Geno data to a file.

        Args:
            path (str, optional): Folder path to save the file. Defaults to the current directory.
            fmt (str, optional): File format. Options: .h5 (default), .csv, .txt. Future: .vcf, .vcf.gz.
            sep (str, optional): Delimiter for .csv and .txt formats. Default is tab.
            header (bool, optional): Save column names for .csv and .txt formats. Default is True.
        """
        save_data(self.data, name=self.name, path=path, fmt=fmt, sep=sep, header=header)
        return

def merge_command_parallel(df_subset, path, filetype):
    """Helper function of the update_snpids method to update SNP in parallel when genetic data is split by chromosome."""
    chr_number = df_subset.iloc[0]['CHR']
    
    if filetype == "bed":
        variant_file = path.replace("$", str(chr_number)) + ".bim"
        if not os.path.exists(variant_file):
            return
        variants = pd.read_csv(
            variant_file, sep="\t", names=["CHR", "SNP", "F", "POS", "A1", "A2"]
        )
    else:  # filetype == "pgen"
        variant_file = path.replace("$", str(chr_number)) + ".pvar"
        if not os.path.exists(variant_file):
            return
        variants = pd.read_csv(
            variant_file, sep="\t", comment="#",
            names=["CHR", "POS", "SNP", "A1", "A2"]
        )

    # Convert CHR to string and remove 'chr' prefix if present
    variants["CHR"] = variants["CHR"].astype(str).str.replace("^chr", "", regex=True)
    # Convert numeric values to int, leave X/Y as strings
    variants["CHR"] = variants["CHR"].apply(
        lambda x: int(float(x)) if x.replace(".", "").isdigit() else x
    )
    
    variants.drop_duplicates(subset=["CHR", "POS"], keep='first', inplace=True)
    variants.drop_duplicates(subset=["SNP"], keep='first', inplace=True)

    return df_subset.merge(
        variants[["CHR", "POS", "SNP"]], on=["CHR", "POS"], how="left", suffixes=('', '_new')
    )