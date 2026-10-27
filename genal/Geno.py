import pandas as pd
import numpy as np
import warnings
import os
import copy
import psutil
import uuid
import scipy.stats as st

from .tools import (
    create_tmp,
    load_reference_panel,
    setup_genetic_path,
    check_reference_panel
)
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
    fill_fstatistic,
    check_allele_column,
    check_snp_column,
    check_snp_ids,
    check_eaf_column,
    remove_na,
    update_eaf_func,
)
from .constants import STANDARD_COLUMNS, REF_PANEL_COLUMNS, CHECKS_DICT, MR_METHODS_NAMES

# Do all the MR steps (query_outcome, harmonize etc) based on CHR/POS and not SNPs
# Add proxying function (input is df + searchspace (list of SNP or path to .bim, can be separated by chromosomes) and returns proxied df)
# Get proxies (simply return a list of proxies)
# Include proxying option to association_test
# Multivariable-MR
# Check stability with variants on sexual chromosomes
# Check the build of user data (potentially with a list of SNPs with different positions)
# Consider how update_snpids could be replaced by the extract range function in plink2 (to gain speed)


class Geno:
    """
    A class to handle GWAS-derived data, including SNP rsID, genome position,
    SNP-trait effects, and effect allele frequencies.

    Attributes:
        data (pd.DataFrame): Main DataFrame containing SNP data.
        phenotype (tuple[pd.DataFrame, str, str]): Tuple set by :meth:`Geno.set_phenotype`:
            (phenotype_dataframe, phenotype_type, original_PHENO_column_name).
        MR_data (list): [exposure DataFrame, outcome DataFrame, outcome name]. Set by the
            :meth:`Geno.query_outcome` method.
        MR_results (pd.DataFrame, pd.DataFrame, str, str): Contains an MR results dataframe, a dataframe of harmonized SNPs, an exposure name, an outcome name. Assigned after calling the MR method and used for plotting with the MR_plot method.
        MR_loo_results (tuple): (loo_table, all_instruments_row, exposure_name, outcome_name, method, method_name,
            odds, heterogeneity, mrpresso_subset_used). Set by the :meth:`Geno.MR_loo` method and used by the MR_loo_plot method.
        ram (int): RAM in MB used by PLINK commands. Defaults to 80% of the system RAM.
        cpus (int): Number of CPUs used by parallel steps. Defaults to the number of CPUs available to the process
            (its CPU affinity, e.g. a SLURM allocation) minus 1, and at least 1.
        checks (dict): Dictionary of checks performed on the main DataFrame.
        name (str): ID of the object (for internal reference and debugging purposes).
        reference_panel (pd.DataFrame): Reference population SNP data used for SNP info
            adjustments. Initialized when first needed.
        reference_panel_name (str): string to identify the reference_panel (path or population string)
        MRpresso_subset_data (pd.DataFrame): Harmonized SNPs without the MR-PRESSO outliers (None if no outliers were
            removed). Assigned after calling the MRpresso method.
        MRpresso_results (tuple): (mod_table, GlobalTest, OutlierTest, BiasTest) returned by the MRpresso method.
        MRpresso_outliers (list): SNP IDs flagged by the MR-PRESSO outlier test. Assigned after calling the MRpresso method.
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
        FSTAT="FSTAT",
    ):
        """
        Initializes the Geno object used to store and transform Single Nucleotide Polymorphisms (SNP) data.

        Args:
            df (pd.DataFrame): DataFrame where each row represents a SNP. The index of df is not kept: the `data`
                attribute is a copy of df with a RangeIndex (0 to n-1), so that repeated index labels cannot affect row operations.
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
            FSTAT (str, optional): Column name for F-statistic. Defaults to "FSTAT".
        """

        # Validate df type
        if not isinstance(df, pd.DataFrame):
            raise TypeError("df needs to be a pandas dataframe.")
        data = df.reset_index(drop=True)

        # Standardize column names based on provided parameters +/- delete other columns
        data = adjust_column_names(
            data, CHR, POS, SNP, EA, NEA, BETA, SE, P, EAF, keep_columns, FSTAT
        )

        # Set object attributes
        self.data = data
        self.name = str(uuid.uuid4())[:8]

        # List to keep track of checks performed
        self.checks = CHECKS_DICT.copy()

        # Set the maximal amount of ram/cpu to be used by the methods
        self.cpus = max(1, available_cpus() - 1)
        self.ram = int(psutil.virtual_memory().total / 1024**2 * 0.8)

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
        fill_f=False,
    ):
        """
        Clean and preprocess the main dataframe of Single Nucleotide Polymorphisms (SNP) data.

        Args:
            preprocessing (str, optional): Level of preprocessing to apply. Defaults to 'Fill'. Options include:

                - "None": The dataframe is not modified.
                - "Fill": Missing columns are added based on reference data and invalid values set to NaN, but no rows are deleted.
                - "Fill_delete": Missing columns are added, and rows with missing, duplicated, or invalid values are deleted.
            reference_panel (str or pd.DataFrame, optional): Reference panel for SNP adjustments. 
                Options include: "37", "38" which correspond to the 1000G phase 3 reference panels in build 37 and 38 respectively.
                Can also be a pd.DataFrame with the necessary columns (CHR, POS, SNP, A1, A2).
                Also accepts a path to a specific bed/bim/fam or pgen/pvar/psam panel.
            effect_column (str, optional): Specifies the type of effect column ("BETA" or "OR"). Defaults to None.

                - "BETA": the BETA column is used as is.
                - "OR": the BETA column contains odds ratios and is log-transformed (odds ratios <= 0 give NaN).
                  If SE and P are available, SE is kept when abs(log(OR))/SE reproduces the z-statistics implied by P
                  (SE of log(OR)); otherwise (e.g. SE of the odds ratio) it is removed and recomputed from BETA and P.
                  Without P, SE is kept and assumed to be the SE of log(OR); without SE, SE is computed from BETA and P.
                - None: the scale is determined from the data. BETA with negative or zero values are Beta estimates.
                  If all BETA are positive and SE and P are available, the scale (Beta estimates or odds ratios) is the
                  one for which abs(BETA)/SE or abs(log(BETA))/SE is closer to the z-statistics implied by P. Otherwise, a heuristic
                  is used: odds ratios if the median of BETA is between 0.5 and 1.5. The decision and its evidence are printed,
                  and odds ratios are then handled as with "OR".
            keep_indel (bool, optional): Determines if insertions/deletions should be kept. If None, defers to preprocessing value. Defaults to None.
            keep_dups (bool, optional): Determines if rows with duplicate SNP IDs should be kept. If None, defers to preprocessing value. Defaults to None.
                When duplicates are removed, the first row of each SNP is kept (after the removal of rows with missing values under "Fill_delete"),
                and rows with a missing SNP ID are never considered duplicates. Placeholder SNP IDs ('', '.', 'NA', 'nan', 'None')
                are set to NaN (missing IDs, deleted under "Fill_delete").
            fill_snpids (bool, optional): Decides if the SNP (rsID) column should be created or replaced based on CHR/POS columns and a reference genome. If None, defers to preprocessing value. Defaults to None.
                When EA and NEA are present, the reference variant must have the same alleles; SNPs not found get the ID
                CHR:POS:allele1:allele2 (alleles in alphabetical order). Without EA/NEA, positions with several reference variants are left as NaN.
            fill_coordinates (bool, optional): Decides if CHR and/or POS should be created or replaced based on SNP column and a reference genome. If None, defers to preprocessing value. Defaults to None.
                SNP IDs found at several positions in the reference are left as NaN.
            fill_f (bool, optional): If True, force recomputation/overwrite of the FSTAT column even if it
                already exists. If False (default), FSTAT is created if missing, or only missing values
                are filled if the column already exists.

        Note:
            If you pass a standard reference_panel name (e.g. "EUR_37"), it will be converted to "37".
            CHR values "X", "Y", "XY", "MT"/"M" (with or without "chr") are converted to 23, 24, 25, 26. POS values equal
            to 0 are invalid (positions are 1-based) and set to NaN.
            EAF values outside [0, 1] (or non-numeric) are set to NaN.
            Missing EA and NEA columns are filled from the reference only for data without a BETA column: with a BETA
            column, a ValueError is raised since the allele to which BETA refers cannot be inferred (provide EA; NEA can then be filled).
            A missing P column is computed from BETA/SE as 2*sf(abs(BETA/SE)) (NaN when SE is missing or zero), and a missing
            SE column from BETA/P (NaN when BETA = 0, P = 0 or P = 1). Odds ratios <= 0 give a NaN BETA.
            The BETA column is converted only once: calling preprocess_data again does not transform it again.
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

        # Placeholder SNP IDs ('', '.', 'NA', ...) are missing IDs
        if "SNP" in data.columns and preprocessing in ['Fill', 'Fill_delete']:
            check_snp_ids(data)

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

        # Convert effect column to Beta estimates if present (only once: the BETA column of a preprocessed object is already on the Beta scale)
        if "BETA" in data.columns and preprocessing in ['Fill', 'Fill_delete']:
            if self.checks["BETA"]:
                print("The BETA column has already been processed by a previous preprocess_data call and is not transformed again.")
            else:
                check_beta_column(data, effect_column)
                self.checks["BETA"] = True

        # Ensure P column contains valid values
        if "P" in data.columns and preprocessing in ['Fill', 'Fill_delete']:
            check_p_column(data)
            self.checks["P"] = True

        # Ensure EAF column contains valid values
        if "EAF" in data.columns and preprocessing in ['Fill', 'Fill_delete']:
            check_eaf_column(data)

        # Fill missing SE or P columns if necessary
        if preprocessing in ['Fill', 'Fill_delete']:
            fill_se_p(data)

        # Compute or fill the FSTAT column
        # Under normal preprocessing: compute/create/fill FSTAT
        # If preprocessing == "None": only compute FSTAT when fill_f=True
        if preprocessing in ['Fill', 'Fill_delete'] or fill_f:
            fill_fstatistic(data, overwrite=fill_f)

        # Process allele columns
        for allele_col in ["EA", "NEA"]:
            check_allele_condition = (allele_col in data.columns) and (
                (preprocessing in ['Fill', 'Fill_delete']) or (not keep_indel)
            )
            if check_allele_condition and not self.checks[allele_col]:
                check_allele_column(data, allele_col, keep_indel)
                self.checks[allele_col] = True

        # Warn if essential columns are missing
        missing_columns = [column for column in STANDARD_COLUMNS if column not in data.columns]
        if missing_columns:
            print(f"Warning: the data doesn't include the following columns: {', '.join(missing_columns)}. Some methods may need them.")

        # Remove missing values if preprocessing level is set to 'Fill_delete'
        # (before removing duplicates, so that a complete copy of a duplicated SNP is kept rather than an incomplete one)
        if preprocessing == 'Fill_delete':
            remove_na(data)
            self.checks["NA_removal"] = True

        # Check for and handle duplicate SNPs if necessary (the first row of each SNP is kept)
        if "SNP" in data.columns and not keep_dups:
            check_snp_column(data)
            self.checks["SNP"] = True

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
            kb (int, optional): Clumping window in kilobases (kb): SNPs less than kb kilobases away from an index SNP can be clumped with it. Default is 10000 (very strict).
            r2 (float, optional): Linkage disequilibrium threshold, values between 0 and 1. SNPs with an r2 above this value with an index SNP are clumped with it. Default is 0.01 (strict).
            p1 (float, optional): P-value threshold for index SNPs. SNPs with a P-value higher than this value are excluded. Default is 5e-8.
            p2 (float, optional): Secondary P-value threshold of PLINK (--clump-p2): the clumped SNPs with a P-value below it are listed in the SP2 column of the PLINK report.
                It does not change the index SNPs, which are the SNPs kept in the clumped data. Default is 0.01.
            reference_panel (str, optional): The reference population for linkage disequilibrium values. 
                Acceptable populations are "EUR", "SAS", "AFR", "EAS", "AMR" and available builds are 37 and 38 ("EUR_38" or "AFR_37" etc...)
                Also accepts or a path to a specific bed/bim/fam or pgen/pvar/psam panel.
                Default is "EUR_37".

        Returns:
            genal.Geno or None: A new Geno object with the clumped data (the index SNPs), or None if no SNP remains after clumping.

        Note:
            Rows of the data attribute with a missing SNP or P value are deleted, and duplicated SNP IDs are removed
            (the first row is kept) if preprocess_data has not already removed them.
            Only the SNPs with a P-value below p1 or p2 are given to PLINK (the others cannot change the clumped data):
            the TOTAL, NONSIG and S* counts of the PLINK report (.clumps file in the tmp_GENAL folder) are counts of these SNPs.
        """
        # Validate parameter types
        if not isinstance(kb, int):
            raise TypeError("The 'kb' parameter must be an integer.")
        if not isinstance(r2, (int, float)):
            raise TypeError("The 'r2' parameter must be a float or an integer.")
        if not isinstance(p1, (int, float)):
            raise TypeError("The 'p1' parameter must be a float or an integer.")
        if not isinstance(p2, (int, float)):
            raise TypeError("The 'p2' parameter must be a float or an integer.")

        # Validate parameter values
        if not 0 <= r2 <= 1:
            raise ValueError("The 'r2' parameter must be between 0 and 1.")
        
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
        from .clump import clump_data_plink2

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
        Update or create the column of SNP name based on genetic data, genomic position and alleles.

        Each row is matched to the variant of the genetic data at the same CHR/POS whose alleles are EA and NEA
        (in any order). For SNVs that are not palindromic (not A/T or C/G), the strand-complemented allele pair is also
        accepted, and the EA and NEA of these rows are complemented to the strand of the genetic data. If NEA is missing,
        EA must be an allele of the variant and the variant must be the only one at this position.
        Rows whose SNP ID is present in the genetic data (for data split by chromosome: in the file of their chromosome)
        keep it, even if their CHR/POS match a variant with another ID.
        Rows without a single matching variant ID keep their SNP value (or get a missing SNP value if there is no SNP column).
        A warning is printed when less than 90% of the rows with a variant at their position match its alleles
        (likely genome build mismatch between the data and the genetic files).

        Args:
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files.
                If files are split by chromosomes, replace the chromosome number with '$'.
                For instance: path = "ukb_chr$_file". Defaults to the path from the configuration.
                Files of chromosomes 23-26 can be named with their number or with X, Y, XY, MT.
            replace (bool, optional): To update the .data attribute with the updated SNP column or not.

        Returns:
            pd.DataFrame: A copy of the data with an updated `SNP` column matching the genotype dataset.
            If `replace=True`, the `.data` attribute is also updated in place.
            
        Note:
            This can be used before extracting SNPs from the genetic data if there is possibility of a mismatch between the SNP name contained in the Geno dataframe (SNP-level data) and the SNP name used in the genetic data (individual-level data). 
            Notably, this can avoid losing SNPs due to ID mismatch during polygenic risk scoring or single-SNP association testing.
        """
        from .tools import read_variant_file, COMPLEMENT
        from .geno_tools import missing_snp_ids
        
        # Check mandatory columns
        for column in ["CHR", "POS"]:
            if not (column in self.data.columns):
                raise ValueError(
                    f"The column {column} is not found in the data and is mandatory to update snpIDs!"
                )
        data = self.data.copy() # We don't want to modify the data attribute
        path, filetype = setup_genetic_path(path) #Verify the path and get filetype
        is_split = "$" in path #If data is split by chromosomes or not
        match_columns = [col for col in ["CHR", "POS", "SNP", "EA", "NEA"] if col in data.columns]
        
        # Match the rows of data with the variants of the genetic data
        if not is_split:
            variants = read_variant_file(path + (".bim" if filetype == "bed" else ".pvar"))
            snp_new, flipped, counts = match_variant_ids(data[match_columns], variants)
        else:
            chr_rows = list(data.groupby("CHR", dropna=False, sort=False).indices.values()) #Row positions of each chromosome
            results = [match_chromosome_variant_ids(data[match_columns].iloc[rows], path, filetype) for rows in chr_rows]
            snp_new = np.full(data.shape[0], np.nan, dtype=object)
            flipped = np.zeros(data.shape[0], dtype=bool)
            counts = dict.fromkeys(MATCH_COUNTS, 0)
            for rows, (snp_chr, flipped_chr, counts_chr) in zip(chr_rows, results):
                snp_new[rows] = snp_chr
                flipped[rows] = flipped_chr
                counts = {key: counts[key] + counts_chr[key] for key in counts}

        # Complement the alleles of the rows matched on the opposite strand
        if flipped.any():
            for allele in ["EA", "NEA"]:
                data.loc[flipped, allele] = data.loc[flipped, allele].astype(str).str.upper().map(COMPLEMENT).to_numpy()

        # Update the SNP column
        has_snp = 'SNP' in data.columns
        if has_snp:
            original = data['SNP']
            data['SNP'] = pd.Series(snp_new, index=data.index).fillna(original)
            n_updated = int(((data['SNP'] != original) & data['SNP'].notna()).sum())
        else:
            data['SNP'] = snp_new
            n_updated = int(pd.notna(snp_new).sum())

        # Report the matching outcome
        n_rows = max(1, data.shape[0])
        n_opposite = counts["opposite_strand"]
        print(
            f"{counts['matched']}({counts['matched']/n_rows*100:.3f}%) rows match a variant of the genetic data on CHR/POS and alleles"
            + (f", including {n_opposite} on the opposite strand." if n_opposite > 0 else ".")
        )
        n_unmatched = data.shape[0] - counts["matched"]
        if n_unmatched > 0:
            print(
                f"{n_unmatched}({n_unmatched/n_rows*100:.3f}%) rows are not matched: {counts['no_position_match']} without variant at their position, "
                f"{counts['no_allele_match']} without variant with matching alleles at their position, {counts['ambiguous']} ambiguous (several variants at their position could match)."
            )
        if counts["own_id"] > 0:
            print(
                f"{counts['own_id']} rows keep their SNP ID, which is present in the genetic data"
                + (f" ({counts['own_id_other_match']} of them match a variant with a different ID on CHR/POS)." if counts["own_id_other_match"] > 0 else ".")
            )
        if flipped.any():
            print(f"The EA and NEA alleles of {flipped.sum()} rows matched on the opposite strand have been complemented to the strand of the genetic data.")
        n_missing = int(missing_snp_ids(data['SNP']).sum())
        print(
            f"The SNP column has been updated for {n_updated} rows. "
            + ("The other rows keep their SNP value" if has_snp else "The other rows have a missing SNP value")
            + (f" ({n_missing} rows have a missing SNP ID and can't be found in the genetic data)." if has_snp and n_missing > 0 else ".")
        )
        # Allele concordance at the positions present in the genetic data (low with a genome build mismatch)
        n_concordant = counts["matched"] + counts["ambiguous"]
        n_at_position = n_concordant + counts["no_allele_match"]
        if n_at_position > 0 and n_concordant / n_at_position < 0.9:
            print(
                f"Warning: only {n_concordant/n_at_position*100:.3f}% of the {n_at_position} rows with a variant of the genetic data at their position match its alleles. "
                "The data and the genetic files are likely on different genome builds."
            )
        if replace: self.data = data #Update attribute if replace argument
        return data
    
    def extract_snps(self, path=None):
        """
        Extract the list of SNPs of this Geno object from the genetic data provided.

        Args:
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files.
                If files are split by chromosomes, replace the chromosome number with '$'.
                For instance: path = "ukb_chr$_file". Default is None.
                Files of chromosomes 23-26 can be named with their number or with X, Y, XY, MT.

        Returns:
            None: The output is a pgen/pvar/psam triple in the tmp_GENAL folder
            with the format "{name}_allchr" which includes the extracted SNPs.

        Notes:
            The provided path is saved to the config file. If this function is called again,
            you don't need to specify the path if you want to use the same genomic files.
            When several variants of the genetic data share an ID (e.g. split multiallelic variants), the variant
            matching the EA and NEA alleles of the data is extracted.
        """

        create_tmp() #Make sure temporary folder exists
        
        # Extract the SNP list
        snp_list = self.data["SNP"]

        # Renaming to avoid conflicts with previous extraction
        self.name = str(uuid.uuid4())[:8]
        # Extract SNPs using the provided path and SNP list
        from .extract_prs import extract_snps_func

        _ = extract_snps_func(snp_list, self.name, path, self.ram, self.cpus, data=self.data)

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
                                       BETA is used as given; for odds ratios, run preprocess_data(effect_column="OR") first.
                                       If False, performs an unweighted PRS: each SNP contributes the dosage of its
                                       trait-increasing allele (EA if BETA >= 0, NEA if BETA < 0) with a weight of 1
                                       (SNPs with a negative BETA and no NEA contribute their EA with a weight of -1). Default is True.
            path (str, optional): Path to a bed/bim/fam or pgen/pvar/psam set of genetic files for PRS calculation.
                                  If files are split by chromosomes, replace the chromosome number
                                  with '$'. For instance: path = "ukb_chr$_file".
                                  Files of chromosomes 23-26 can be named with their number or with X, Y, XY, MT.
                                  If not provided, it will use the genetic path most recently used 
                                  (if any). Default is None.
            proxy (bool, optional): If True, search proxy SNPs when variants are missing from the
                target genotype dataset. Default is False.
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
            None: Writes `<name>.csv` (or `<self.name>.csv` if `name=None`) and stores intermediate
            PLINK outputs under `tmp_GENAL/`. The CSV file has the columns FID, IID, ALLELE_CT (number of
            non-missing alleles), NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG (SCORE1_SUM divided by the number of alleles
            of the SNPs used, i.e. 2 per SNP) and SCORE1_SUM (sum over the SNPs used of the weight times the dosage
            of the scored allele, missing genotypes being mean-imputed).
        """
        from .tools import read_variant_file, split_chromosome_path, SPLIT_CHROMOSOMES
        from .geno_tools import missing_snp_ids
        
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
        # BETA is used as given (odds ratios must be log-transformed with preprocess_data(effect_column="OR"))
        data_prs["BETA"] = pd.to_numeric(data_prs["BETA"], errors="coerce")
        if not weighted and "NEA" in data_prs.columns and not self.checks.get("NEA"):
            check_allele_column(data_prs, "NEA", keep_indel=False)

        initial_rows = data_prs.shape[0]
        data_prs.dropna(subset=["SNP", "EA", "BETA"], inplace=True)
        data_prs = data_prs[~missing_snp_ids(data_prs["SNP"]).to_numpy()]
        deleted_rows = initial_rows - data_prs.shape[0]
        if deleted_rows > 0:
            print(
                f"{deleted_rows} ({deleted_rows/initial_rows*100:.3f}%) rows with missing values in columns SNP (NA or placeholder IDs such as '.'), EA, or BETA have been deleted."
            )
            
        from .extract_prs import prs_func

        # If proxy option
        if proxy:
            from .proxy import apply_proxies, find_proxies

            print("Identifying the SNPs present in the genetic data...")
            # Obtain the list of SNPs present in the genetic data
            extension = ".bim" if filetype == "bed" else ".pvar"
            if path.count("$") == 1: #If split: merge all SNP columns of the variant files
                genetic_snp_list = []
                for i in SPLIT_CHROMOSOMES:
                    prefix = split_chromosome_path(path, i, filetype)
                    if prefix is not None:
                        genetic_snp_list.extend(read_variant_file(prefix + extension).SNP.tolist())
            else: #If not split
                variant_file = path + extension
                if not os.path.exists(variant_file):
                    raise ValueError(f"The file {variant_file} does not exist.")
                genetic_snp_list = read_variant_file(variant_file).SNP.tolist()
                
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
                    # Remove duplicates in the SNP column (by position, keeping the first row)
                    duplicated = data_prs.duplicated(subset=["SNP"], keep="first").to_numpy()
                    n_del = int(duplicated.sum())
                    if n_del > 0:
                        data_prs = data_prs[~duplicated]
                        print(f"After proxying, {n_del} SNPs had duplicated IDs and were removed.")
        
        # Renaming to avoid conflicts with previous extraction
        self.name = str(uuid.uuid4())[:8]
        # Compute PRS
        prs_data = prs_func(data_prs, weighted, path, ram=self.ram, cpus=self.cpus, name=self.name)

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
            IID (str): Name of the individual IDs column in 'data' (required). These IDs should
                                correspond to the genetic IDs in the FAM file that will be used for association testing.
            PHENO (str): Name of the phenotype column in 'data' (required) which will be used
                                as the dependent variable for association tests.
            PHENO_type (str, optional): If not specified, the function will try to infer if
                                        the phenotype is binary or quantitative. To bypass this,
                                        use "quant" for quantitative or "binary" for binary phenotypes.
                                        Default is None.
            FID (str, optional): Name of the family ID column in 'data'. If not provided,
                                FID values will be set to the same as IID values.
            alternate_control (bool, optional): For a binary trait, values 0/1 are read as 0 = control, 1 = case,
                                                values 1/2 as 1 = control, 2 = case, and for any other pair of
                                                values the most frequent value is the control.
                                                Set to True to invert the control and case values. Default is False.
        """

        from .association import set_phenotype_func

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
                Files of chromosomes 23-26 can be named with their number or with X, Y, XY, MT.
            covar (list, optional): List of columns in the phenotype dataframe to be used
                                    as covariates in the association tests. Default is an empty list.
                                    Text or categorical columns are encoded as indicator columns (k - 1 for
                                    k levels, the first level in sorted order being the reference); numeric
                                    columns are continuous (convert numeric codes, e.g. of assessment centres,
                                    to str or category). Individuals with a missing covariate are excluded,
                                    and covariates with a single value among the analysed individuals are removed.
            standardize (bool, optional): If True, it will standardize a quantitative phenotype
                                          before performing association tests. This is typically done
                                          to make results more interpretable. Default is True.

        Returns:
            None: Replaces the data attribute with the tested SNPs, with BETA, SE, P and FSTAT columns updated from the
                  results of the association tests. SNPs with a missing SNP ID or not extracted from the genetic data, with alleles not matching
                  the tested variant, with an NA result, or at multiallelic variants where the contrast tested by PLINK
                  is not EA versus NEA (or without NEA) are removed.

        Note:
            This method requires the phenotype to be set using the set_phenotype() function.
            The data must contain the EA column. The PLINK results are matched to the data by variant ID and aligned on EA
            (using NEA when present, otherwise the alleles of the genetic variant).
            Phenotype columns of the genetic data (.psam/.fam file) are ignored, and covariates are
            variance-standardized by PLINK (this does not change the SNP effect estimates).
        """

        # Ensure that the phenotype has been set using set_phenotype
        if not hasattr(self, "phenotype"):
            raise ValueError(
                "You first need to set a phenotype using .set_phenotype(data, PHENO, PHENO_type, IID)!"
            )
        
        if "EA" not in self.data.columns:
            raise ValueError("The column EA is not found in the data and is mandatory for association tests.")

        # Ensure that covar does not contain the PHENO column
        if self.phenotype[2] in covar:
            print(f"The phenotype column {self.phenotype[2]} is also present in the covariates list and will be removed from the association tests.")
            covar = [col for col in covar if col != self.phenotype[2]]

        create_tmp() #Make sure temporary folder exists
        from .association import association_test_func_plink2
        from .extract_prs import extract_snps_func
        
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

        # Rows without SNP ID can't be found in the genetic data
        from .geno_tools import missing_snp_ids
        missing = missing_snp_ids(data["SNP"])
        if missing.any():
            print(f"{missing.sum()} rows with a missing SNP ID (NA or placeholder such as '.') can't be tested and are removed.")
            data = data[~missing.to_numpy()]
        
        # Extract the SNP list
        snp_list = data["SNP"]

        # Renaming to avoid conflicts with previous extraction
        self.name = str(uuid.uuid4())[:8]
        # Extract SNPs using the provided path and SNP list
        path = extract_snps_func(snp_list, self.name, path, self.ram, self.cpus, data=data)
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

        # Recompute FSTAT to be consistent with the updated BETA/SE/P values
        fill_fstatistic(updated_data, overwrite=True)

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

        Queries the outcome data, with or without proxying, and stores the list
        [exposure_data, outcome_data, outcome_name] in the `MR_data` attribute, ready for MR methods.
        Rows with a missing SNP ID (NaN or placeholders such as '', '.', 'NA') are excluded from both datasets.
        MR-PRESSO results stored from a previous outcome (`MRpresso_results`, `MRpresso_subset_data`,
        `MRpresso_outliers`) are cleared.

        Args:
            outcome: Can be a Geno object (from a GWAS) or a `.parquet` filepath created with the :meth:`Geno.save` method.
            name (str, optional): Name for the outcome data. Defaults to None: the name of the outcome Geno object,
                or the file name (without extension) for a filepath.
            proxy (bool, optional): If true, proxies are searched. Proxies must be in LD with the missing SNP
                (r2 threshold), non-palindromic, allele-compatible with the outcome data, with valid BETA/SE, and not
                themselves exposure SNPs present in the outcome data. They are assigned greedily by decreasing abs(r)
                (a proxy serves at most one missing SNP). The proxied association is reported for the missing SNP, with
                its reference panel minor allele as effect allele and its reference panel frequency as EAF (palindromic
                SNPs are then treated like any other SNP by the `action` of the MR methods). Default is True.
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

        from .MR_tools import query_outcome_func

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

        # MR-PRESSO results refer to the previous MR_data: clear them
        cleared = False
        for attr in ["MRpresso_subset_data", "MRpresso_results", "MRpresso_outliers"]:
            if hasattr(self, attr):
                delattr(self, attr)
                cleared = True
        if cleared:
            print("MR-PRESSO results of the previous outcome have been cleared (run MRpresso again for this outcome).")
        return

    def MR(
        self,
        methods=["IVW","IVW-FE","WM","Simple-mode","Egger"],
        action=2,
        eaf_threshold=0.42,
        heterogeneity=False,
        nboot=10000,
        penk=20,
        phi=1,
        exposure_name=None,
        outcome_name=None,
        cpus=None,
        odds=False,
        use_mrpresso_data=False,
        seed=None,
    ):
        """
        Executes Mendelian Randomization (MR) using exposure/outcome data stored in `MR_data`
        (created by calling :meth:`Geno.query_outcome`).

        If a single instrument remains after harmonization, the Wald ratio (BETA_o / BETA_e, with standard
        error SE_o / abs(BETA_e)) is returned instead of the requested methods.

        Args:
            methods (list, optional): List of MR methods to run. Default is ["IVW","IVW-FE","WM","Simple-mode","Egger"].
                Possible options include:

                - "IVW": inverse variance-weighted with random effects and under-dispersion correction
                - "IVW-FE": inverse variance-weighted with fixed effects
                - "IVW-RE": inverse variance-weighted with random effects and without under-dispersion correction
                - "UWR": unweighted regression (ordinary least squares standard error)
                - "WM": weighted median (bootstrapped standard errors)
                - "WM-pen": penalised weighted median (bootstrapped standard errors)
                - "Simple-median": simple median (bootstrapped standard errors)
                - "Sign": sign concordance test
                - "Egger": egger regression
                - "Egger-boot": egger regression with bootstrapped standard errors (regular MR Egger estimates with the
                  standard deviation of bootstrap replicates as standard error and a normal p-value; the bootstrap
                  SE ignores between-variant heterogeneity and can be smaller than the "Egger" SE)
                - "Simple-mode": simple mode method
                - "Weighted-mode": weighted mode method
                - "all": all the methods above
            action (int, optional): How to treat palindromes during harmonizing between
                exposure and outcome data. Accepts:
                1: Doesn't flip them (Assumes all alleles are on the forward strand)
                2: Uses allele frequencies to attempt to flip (conservative, default)
                3: Removes all palindromic SNPs (very conservative)
            eaf_threshold (float, optional): Max effect allele frequency accepted when
                flipping palindromic SNPs (relevant if action=2). Palindromic SNPs with an EAF strictly
                between eaf_threshold and 1 - eaf_threshold (or a missing EAF) are removed. Default is 0.42.
            heterogeneity (bool, optional): If True, includes heterogeneity tests in the results (Cochran's Q test). Default is False.
            nboot (int, optional): Number of bootstrap replications for methods with bootstrapping. Default is 10000.
            penk (int, optional): Penalty value for the WM-pen method. Default is 20.
            phi (int, optional): Factor for the bandwidth parameter used in the kernel density estimation of the mode methods
            exposure_name (str, optional): Name of the exposure data (only for display purposes).
                Defaults to the name of the Geno object.
            outcome_name (str, optional): Name of the outcome data (only for display purposes).
            cpus (int, optional): Not used: the MR methods run on a single core. A message is printed when it is given.
            odds (bool, optional): If True, adds an odds ratio column with 95% confidence intervals
                (not computed for the Egger intercept rows, which are not causal effects). The intervals use the
                distribution of each method's p-value: t with nSNP - 2 degrees of freedom for MR Egger, t with nSNP - 1
                for the mode methods, normal otherwise. Default is False.
            use_mrpresso_data (bool, optional): If True and MRpresso has been run, uses the subset of instruments after outlier removal
                (all instruments if MR-PRESSO identified no outliers; a ValueError is raised if all instruments are outliers). Default is False.
            seed (int, optional): Random seed for the bootstrap steps (WM, WM-pen, Simple-median, Egger-boot,
                Simple-mode, Weighted-mode). Default None (non-reproducible).

        Returns:
            pd.DataFrame: A table with MR results.
        """

        # Ensure that query_outcome has been previously called
        if not hasattr(self, "MR_data"):
            raise ValueError("You must first call query_outcome() before running MR.")

        from .MR_tools import MR_func, mrpresso_subset, odds_ratio_ci

        # Get MRpresso subset data if requested and available
        subset_data = None
        if use_mrpresso_data:
            if not hasattr(self, "MRpresso_subset_data"):
                raise ValueError("Use_mrpresso_data is set to True but MRpresso results not found. Please run MRpresso first.")
            subset_data = mrpresso_subset(self.MRpresso_subset_data, getattr(self, "MRpresso_outliers", []))

        if outcome_name:
            self.MR_data[2] = outcome_name
        exp_name = exposure_name if exposure_name else self.name
        if cpus is not None:
            print("MR is now single-core only thanks to its optimization: the cpus argument is ignored.")

        res, df_mr = MR_func(
            self.MR_data,
            methods,
            action,
            eaf_threshold,
            nboot,
            penk,
            phi,
            exp_name,
            1,
            subset_data,
            seed,
        )

        if not heterogeneity:
            res = res.loc[:,["exposure", "outcome", "method", "nSNP", "b", "se", "pval"]]

        if odds and not res.empty:
            # Odds ratios and 95% confidence intervals (same distribution as the p-value of each method;
            # intercepts are not causal effects: no odds ratio)
            intercept_methods = [MR_METHODS_NAMES["Egger"][1], MR_METHODS_NAMES["Egger-boot"][1]]
            res.loc[:,'OR_95CI'] = res.apply(lambda row:
                odds_ratio_ci(row['b'], row['se'], row['method'], row['nSNP'])
                if row['method'] not in intercept_methods
                else np.nan, axis=1)

        self.MR_results = (res, df_mr, exposure_name, outcome_name)

        return res

    def MR_loo(
        self,
        method="IVW",
        action=2,
        eaf_threshold=0.42,
        heterogeneity=False,
        nboot=1000,
        penk=20,
        phi=1,
        exposure_name=None,
        outcome_name=None,
        cpus=None,
        odds=False,
        use_mrpresso_data=False,
        seed=None,
    ):
        """
        Executes leave-one-out Mendelian Randomization (MR) using exposure/outcome data stored in `MR_data`
        (created by calling :meth:`Geno.query_outcome`). The MR method is fitted on all instruments and on
        each subset with one instrument removed. At least 3 instruments are required.

        Args:
            method (str, optional): A single MR method key (see :meth:`Geno.MR`; "all" is not accepted). Default is "IVW".
            action (int, optional): How to treat palindromes during harmonizing (1, 2 or 3, see :meth:`Geno.MR`). Default is 2.
            eaf_threshold (float, optional): Max effect allele frequency accepted when
                flipping palindromic SNPs (relevant if action=2). Default is 0.42.
            heterogeneity (bool, optional): If True, includes heterogeneity statistics (Q, Q_df, Q_pval). Default is False.
            nboot (int, optional): Number of bootstrap replications for methods with bootstrapping. Default is 1000.
            penk (int, optional): Penalty value for the WM-pen method. Default is 20.
            phi (int, optional): Factor for the bandwidth parameter used in the kernel density estimation of the mode methods. Default is 1.
            exposure_name (str, optional): Name of the exposure data (only for display purposes).
                Defaults to the name of the Geno object.
            outcome_name (str, optional): Name of the outcome data (only for display purposes).
            cpus (int, optional): Not used: the MR methods run on a single core. A message is printed when it is given.
            odds (bool, optional): If True, adds an odds ratio column with 95% confidence intervals (same distribution
                as the p-value of the method, see :meth:`Geno.MR`). Default is False.
            use_mrpresso_data (bool, optional): If True and MRpresso has been run, uses the subset of instruments after outlier removal
                (all instruments if MR-PRESSO identified no outliers; a ValueError is raised if all instruments are outliers). Default is False.
            seed (int, optional): Random seed for the bootstrap steps of bootstrap-based methods (each fit uses its own
                random generator derived from the seed). Default None (non-reproducible).

        Returns:
            pd.DataFrame: A tidy table with leave-one-out estimates.
        """
        if not hasattr(self, "MR_data"):
            raise ValueError("You must first call query_outcome() before running MR.")

        from .MR_tools import MR_loo_func, mrpresso_subset, odds_ratio_ci

        subset_data = None
        if use_mrpresso_data:
            if not hasattr(self, "MRpresso_subset_data"):
                raise ValueError(
                    "Use_mrpresso_data is set to True but MRpresso results not found. Please run MRpresso first."
                )
            subset_data = mrpresso_subset(self.MRpresso_subset_data, getattr(self, "MRpresso_outliers", []))

        if outcome_name:
            self.MR_data[2] = outcome_name

        exp_name = exposure_name if exposure_name else self.name
        outcome_name_final = self.MR_data[2]
        if cpus is not None:
            print("MR is now single-core only thanks to its optimization: the cpus argument is ignored.")

        loo_df, all_row, method_display_name = MR_loo_func(
            self.MR_data,
            method,
            action,
            eaf_threshold,
            heterogeneity,
            nboot,
            penk,
            phi,
            exp_name,
            1,
            subset_data,
            seed,
        )

        if odds and not loo_df.empty:
            loo_df.loc[:, "OR_95CI"] = loo_df.apply(
                lambda row: odds_ratio_ci(row["b"], row["se"], row["method"], row["nSNP"]),
                axis=1,
            )

        self.MR_loo_results = (
            loo_df,
            all_row,
            exp_name,
            outcome_name_final,
            method,
            method_display_name,
            odds,
            heterogeneity,
            subset_data is not None,  # True if MR-PRESSO outliers were removed
        )
        self.MR_loo_config = {
            "action": action,
            "eaf_threshold": eaf_threshold,
            "nboot": nboot,
            "penk": penk,
            "phi": phi,
            "seed": seed,
        }

        return loo_df

    def MR_loo_plot(
        self,
        snps_per_page=40,
        page=None,
        top_influential=True,
        filename=None,
        exposure_name=None,
        outcome_name=None,
        figure_size=None,
        *,
        methods=None,
        use_mrpresso_data=None,
    ):
        """
        Create a leave-one-out forest plot using results stored from :meth:`Geno.MR_loo`.

        Args:
            snps_per_page (int, optional): Maximum number of instruments per page (at least 5). Default is 40.
            page (int, optional): Page to render (starting at 1) when top_influential is False. Default is None (all pages).
            top_influential (bool, optional): If True, plots on a single page the `snps_per_page` instruments whose removal
                changes the estimate the most. If False, plots all the instruments over several pages. Default is True.
            filename (str, optional): Filename prefix where the plot is saved (`<filename>.png`, or `<filename>_pageNN.png`
                for each page of a multi-page plot). If None, the plot is not saved.
            exposure_name (str, optional): Exposure label. Defaults to the exposure name used in MR_loo.
            outcome_name (str, optional): Outcome label. Defaults to the outcome name used in MR_loo.
            figure_size (tuple[float, float] | None): Figure size in inches for both the returned
                plot and the saved PNG (when `filename` is provided). In small figures, text is reduced (to at least 6 points) and, below 4.5 inches wide, legends are placed below the plot.
            methods (str or list of str, optional): Keyword-only. Additional MR method keys (see :meth:`Geno.MR`) shown as
                overall estimate rows, computed with the settings of the MR_loo call. Default is None.
            use_mrpresso_data (bool, optional): Keyword-only. If True, highlights the SNPs flagged as outliers by
                :meth:`Geno.MRpresso` (and adds the MR-PRESSO corrected estimate when MR_loo used all instruments).
                Default is None: the use_mrpresso_data value of the MR_loo call.

        Returns:
            plotnine.ggplot.ggplot or list: The plot, or a list of plots (one per page) when several pages are rendered.

        Raises:
            ValueError: If MR_loo has not been run or returned no results.
        """
        from .plots import mr_loo_plot

        return mr_loo_plot(
            getattr(self, "MR_loo_results", None),
            snps_per_page=snps_per_page,
            page=page,
            top_influential=top_influential,
            methods=methods,
            use_mrpresso_data=use_mrpresso_data,
            mr_data=getattr(self, "MR_data", None),
            mrpresso_outliers=getattr(self, "MRpresso_outliers", None),
            mrpresso_subset_data=getattr(self, "MRpresso_subset_data", None),
            mr_loo_config=getattr(self, "MR_loo_config", None),
            filename=filename,
            exposure_name=exposure_name,
            outcome_name=outcome_name,
            figure_size=figure_size,
        )
    
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
        filename=None,
        figure_size=None,
        *,
        use_mrpresso_data=False,
    ):
        """
        Creates and returns a scatter plot of individual SNP effects with lines representing different Mendelian Randomization (MR) methods. Each MR method specified in the 'methods' argument is represented as a line in the plot.

        Args:
            methods (list of str, optional): A list of MR methods to be included in the plot.
                Default methods are "IVW", "WM", "Simple-mode", and "Egger".
            exposure_name (str, optional): A custom label for the exposure effect axis. If None, uses the label provided in the MR function call or a default label.
            outcome_name (str, optional): A custom label for the outcome effect axis. If None, uses the label provided in the MR function call or a default label.
            filename (str, optional): The filename prefix where the plot will be saved (a `.png` extension is added).
                If None, the plot is not saved.
            figure_size (tuple[float, float] | None): Figure size in inches for both the returned
                plot and the saved PNG (when `filename` is provided). Defaults to (10, 6). In small figures, text is reduced (to at least 6 points) and, below 4.5 inches wide, legends are placed below the plot.
            use_mrpresso_data (bool, optional): If True, highlights the SNPs flagged as outliers by
                :meth:`Geno.MRpresso` (requires MRpresso to be run first). Default is False.

        Returns:
            plotnine.ggplot.ggplot: A plotnine ggplot object representing the scatter plot of individual SNP effects with MR method lines.

        Raises:
            ValueError: If MR analysis has not been performed prior to calling this function.

        Note:
            This function requires prior execution of the `MR` method to compute MR results. Make sure the MR analysis is performed on the data before calling `MR_plot`.
        """
        from .plots import mr_plot

        return mr_plot(
            getattr(self, "MR_results", None),
            methods=methods,
            exposure_name=exposure_name,
            outcome_name=outcome_name,
            use_mrpresso_data=use_mrpresso_data,
            mrpresso_outliers=getattr(self, "MRpresso_outliers", None),
            filename=filename,
            figure_size=figure_size,
        )
    
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
        filename=None,
        figure_size=None,
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
            figure_size (tuple[float, float] | None): Figure size in inches for both the returned
                plot and the saved PNG (when `filename` is provided). Default: 7 inches wide, height depending on
                the number of methods. In small figures, text is reduced (to at least 6 points).

        Returns:
            plotnine.ggplot.ggplot: A plotnine ggplot object representing the forest plot 
                of MR results.

        Raises:
            ValueError: If MR analysis has not been performed prior to calling this function.
        """
        from .plots import mr_forest

        return mr_forest(
            getattr(self, "MR_results", None),
            methods=methods,
            exposure_name=exposure_name,
            outcome_name=outcome_name,
            odds=odds,
            filename=filename,
            figure_size=figure_size,
        )

    def MR_funnel(
        self,
        methods="IVW",
        exposure_name=None,
        outcome_name=None,
        se_method="nome",
        symmetric_x=True,
        add_null_line=False,
        figure_size=(10, 6),
        filename=None,
        *,
        use_mrpresso_data=False,
    ):
        """
        Creates and returns a funnel plot of single-SNP causal estimates (Wald ratios).

        Args:
            methods (None | str | list[str] | tuple[str], optional): MR method(s) used to draw
                vertical reference lines. Default is "IVW". Accepts:

                - None: no vertical lines
                - str: a single method key (e.g., "IVW")
                - list/tuple of str: multiple method keys (e.g., ["WM", "IVW", "Egger"])
            exposure_name (str, optional): Exposure label (not displayed on the funnel plot).
            outcome_name (str, optional): Outcome label (not displayed on the funnel plot).
            se_method (str, optional): Standard error method for Wald ratios. One of {"nome", "delta"}.
                Default is "nome".
            symmetric_x (bool, optional): If True, apply symmetric x-limits centered on the first
                resolved MR method line (or on the median ratio if no method line is available).
            add_null_line (bool, optional): If True, adds a subtle vertical line at x=0.
            figure_size (tuple[float, float], optional): Figure size in inches for both the returned
                plot and the saved PNG (when `filename` is provided). Defaults to (10, 6). In small figures, text is reduced (to at least 6 points) and, below 4.5 inches wide, legends are placed below the plot.
            filename (str, optional): Filename prefix where the plot will be saved (a `.png`
                extension is added). If None, the plot is not saved.
            use_mrpresso_data (bool, optional): If True, highlights the SNPs flagged as outliers by
                :meth:`Geno.MRpresso` (requires MRpresso to be run first). Default is False.

        Returns:
            plotnine.ggplot.ggplot: A plotnine ggplot object representing the funnel plot.

        Raises:
            ValueError: If MR analysis has not been performed prior to calling this function.
        """
        from .plots import mr_funnel

        return mr_funnel(
            getattr(self, "MR_results", None),
            methods=methods,
            exposure_name=exposure_name,
            outcome_name=outcome_name,
            se_method=se_method,
            symmetric_x=symmetric_x,
            add_null_line=add_null_line,
            figure_size=figure_size,
            use_mrpresso_data=use_mrpresso_data,
            mrpresso_outliers=getattr(self, "MRpresso_outliers", None),
            filename=filename,
        )

    def MRpresso(
        self,
        action=2,
        eaf_threshold=0.42,
        n_iterations=10000,
        outlier_test=True,
        distortion_test=True,
        significance_p=0.05,
        cpus=None,
        seed=None,
    ):
        """
        Executes the MR-PRESSO Mendelian Randomization algorithm for detection and correction of horizontal pleiotropy.
        It computes the same tests as the MR-PRESSO R package, over 99.7% faster (about 0.3 s for 812 instruments and
        10,000 iterations, on a single core). At least 4 instruments are required after harmonization: otherwise empty results are returned and no MR-PRESSO
        results are stored (MR methods with `use_mrpresso_data=True` then require a successful MRpresso call).

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
            cpus (int, optional): Not used: MR-PRESSO runs on a single core. A message is printed when it is given.
            seed (int, optional): Random seed for the simulation steps (random data generation and distortion test).
                Default None (non-reproducible).

        Returns:
            tuple: (mod_table, GlobalTest, OutlierTest, DistortionTest), with:

            - mod_table: DataFrame containing the original (before outlier removal) and outlier-corrected
              (after outlier removal) inverse variance-weighted MR results. The outlier-corrected row is computed
              whenever the outlier test flags at least one (but not all) SNPs, and is NaN otherwise.
            - GlobalTest: dict with the observed residual sum of squares ("RSSobs") and the p-value of the global
              MR-PRESSO test of horizontal pleiotropy ("global_test_p").
            - OutlierTest: DataFrame (indexed by SNP) with the Bonferroni-corrected p-value of each SNP ("Pvalue")
              for being responsible for the global pleiotropy. Empty if outlier_test is False or if the global test
              p-value is not below significance_p.
            - DistortionTest: dict with the outlier SNP IDs ("outliers_indices"), the distortion coefficient
              ("distortion_test_coefficient") and the p-value of the distortion test ("distortion_test_p").
              Empty if the distortion test is not run.

            Distortion test (as in the MR-PRESSO R package): the distortion coefficient is
            100 * (b_all - b_no_outliers) / abs(b_no_outliers), with b the inverse variance-weighted estimate. Its null
            distribution is obtained from n_iterations IVW fits on J - k instruments: the k outliers plus J - 2k
            non-outliers drawn at random with replacement (J being the total number of instruments), and the p-value
            is the proportion of null coefficients larger in absolute value than the observed one. The test is skipped
            (NaN coefficient and p-value) when no SNP or all SNPs are outliers, when at least half of the instruments
            are outliers, or when the outlier-corrected estimate is 0.
            Empirical p-values reported as 0 mean p < 1/n_iterations (global and distortion tests) or
            p < nSNP/n_iterations (Bonferroni-corrected outlier test).
            The SNP IDs flagged by the outlier test are stored in `MRpresso_outliers` (used by plots with
            `use_mrpresso_data=True`). Whenever at least one (but not all) SNPs are flagged, the outlier-removed harmonized
            data is stored in `MRpresso_subset_data` (used by MR methods with `use_mrpresso_data=True`), whether or not
            the distortion test is run.
        """

        if not hasattr(self, "MR_data"):
            raise ValueError("You must first call query_outcome() before running MR.")

        from .MR_tools import mrpresso_func

        mod_table, GlobalTest, OutlierTest, BiasTest, subset_data = mrpresso_func(
            self.MR_data,
            action,
            eaf_threshold,
            n_iterations,
            outlier_test,
            distortion_test,
            significance_p,
            cpus,
            seed,
        )

        # Results of a previous MRpresso call are replaced; nothing is stored when MR-PRESSO could not be run
        for attr in ["MRpresso_subset_data", "MRpresso_results", "MRpresso_outliers"]:
            if hasattr(self, attr):
                delattr(self, attr)
        if mod_table.empty:
            return mod_table, GlobalTest, OutlierTest, BiasTest

        if subset_data is not None:
            subset_data.drop(columns=["Weights"], inplace=True, errors='ignore')
        self.MRpresso_subset_data = subset_data
        self.MRpresso_results = mod_table, GlobalTest, OutlierTest, BiasTest
        self.MRpresso_outliers = (
            [] if OutlierTest.empty
            else OutlierTest.index[OutlierTest["Pvalue"] <= significance_p].astype(str).tolist()
        )

        return mod_table, GlobalTest, OutlierTest, BiasTest

    def filter_by_gene(self, gene, id_type="symbol", window_size=1000000, *, build, replace=False):
        """
        Filter the data to include only variants that are within a specified distance of a specific gene.

        Args:
            gene (str or int): Identifier for the gene/protein to filter variants around.
            id_type (str, optional): Type of identifier provided. Default is "symbol". Text identifiers are matched
                case-insensitively. Options are:

                - "symbol": Gene symbol (e.g., "APOE")
                - "HGNC": HGNC ID (e.g., "HGNC:613")
                - "name": Full gene name (e.g., "apolipoprotein E")
                - "Ensembl": Ensembl gene ID (e.g., "ENSG00000130203")
                - "NCBI": NCBI gene ID (e.g., "348" or 348)
                - "UCSC": UCSC gene ID (e.g., "uc001hbu.2")
                - "Vega": Vega gene ID (e.g., "OTTHUMG00000019505")
            window_size (int, optional): Total window size in base pairs (±window_size/2 around the gene).
                Default is 1,000,000 (±500kb).
            build (str): Genome build of the CHR/POS coordinates of the data, "37" or "38" (required, keyword only).
            replace (bool, optional): If True, filter the current object in place. Default is False.

        Returns:
            Geno or None: If replace is False (default), a new :class:`~genal.Geno` with the filtered data; if replace
            is True, None (the `.data` attribute is replaced). The filtered data has a Distance column (in base pairs,
            to the nearest gene boundary: negative before the gene start, 0 within the gene, positive after the gene end).

        Raises:
            ValueError: If required columns are missing, gene information cannot be found (gene absent, or without
                chromosome or coordinates in the requested build in the gene information file), or invalid id_type is provided.
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

        from .genes import filter_by_gene_func

        filtered = filter_by_gene_func(self.data, gene, id_type, window_size, build=build)
        
        if replace:
            self.data = filtered
        else:
            Geno_filtered = self.copy(filtered)
            return Geno_filtered

    def colocalize(self, outcome, method="abf", trait1_type=None, trait2_type=None,
                   sdY1=None, sdY2=None, n1=None, n2=None, p1=1e-4, p2=1e-4, p12=1e-5, merge_on_snp=False):
        """
        Perform colocalization analysis between two GWAS datasets.

        Each dataset is checked before merging: variants with missing or non-finite BETA/SE or SE <= 0, with a missing
        merge key (CHR/POS, or a missing or placeholder SNP ID such as '.'), and duplicated variants (same merge key and
        alleles: only the row with the lowest P-value is kept) are removed. When both datasets have EA/NEA columns, a pair
        of variants is considered the same variant only if their alleles match (in either orientation, or on the
        complementary strand for non-palindromic SNVs), and each variant is used in a single pair (multi-allelic sites
        keep one pair per matched variant).

        Args:
            outcome: Another Geno object containing the outcome dataset
            method: Method to use for colocalization. Only "abf" (approximate Bayes factors, as coloc.abf) is available. Default is "abf".
            trait1_type: Type of exposure trait ("quant" for quantitative traits or "cc" for case-control traits), default is "quant"
            trait2_type: Type of outcome trait ("quant" for quantitative traits or "cc" for case-control traits), default is "quant"
            sdY1: Standard deviation of exposure trait (required for quantitative traits, but can be estimated from EAF and sample size).
                If neither sdY1 nor n1 with an EAF column is available, sdY1 = 1 is assumed, i.e. the effects are assumed
                to be on a standardized scale (per standard deviation of the trait).
            sdY2: Standard deviation of outcome trait (required for quantitative traits, but can be estimated from EAF and sample size).
                Same assumption as sdY1 if it cannot be estimated.
            n1: Sample size for exposure (used to estimate sdY1 from the EAF and SE of all exposure variants if sdY1 is not provided)
            n2: Sample size for outcome (used to estimate sdY2 from the EAF and SE of all outcome variants if sdY2 is not provided)
            p1: Prior probability SNP associated with exposure, default is 1e-4
            p2: Prior probability SNP associated with outcome, default is 1e-4
            p12: Prior probability SNP associated with both traits, default is 1e-5.
                As in coloc.abf, p1 (p2) is set to 1/(n + 1) when n * p1 >= 1, n being the number of variants of the
                exposure (outcome) dataset, and p12 likewise with n the number of overlapping variants.
            merge_on_snp: If True, merge the datasets on SNP column. If False, first attempt to merge on CHR and POS columns.

        Returns:
            dict: Number of overlapping variants used ('nsnps') and posterior probabilities 'PP.H0.abf' to 'PP.H4.abf'.
        """
        if method != "abf":
            raise ValueError(f"Invalid method: {method!r}. Only method='abf' is available for colocalization.")

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
            
        from .colocalization import coloc_abf_func

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

        Variants lifted to an alternative, unplaced or patch contig (e.g. ``chrUn_...``, ``..._alt``, ``..._random``, ``..._fix``)
        are treated as unlifted and removed. Chromosomes 23, 24 and 26 are lifted as chrX, chrY and chrM.
        For variants lifted to the opposite strand (inverted regions between builds), the EA and NEA alleles
        are reverse-complemented.

        Args:
            start (str, optional): Current build of the data: "hg17", "hg18", "hg19" or "hg38"
                ("GRCh37" and "GRCh38" are accepted for "hg19" and "hg38"). Default is "hg19".
            end (str, optional): Target build for the liftover (same accepted values as start). Default is "hg38".
            replace (bool, optional): If True, updates the data attribute in place. Default is False.
            extraction_file (bool, optional): If True, prints a CHR POS SNP space-delimited
                file named "{name}_lifted_extraction.txt" ("{object name}_lifted_extraction.txt" if name is not provided). Default is False.
            chain_file (str, optional): Path to a local chain file for the lift.
                If provided, `start` and `end` arguments are not considered. Default is None.
            name (str, optional): Filename or filepath (without extension) to save the lifted dataframe.
                If not provided, the data is not saved.
            liftover_path (str, optional): Path to the UCSC liftOver executable, to lift with it instead of in python. If not provided (default), the lift is done in python, which is faster.

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

        from .lift import lift_data

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

        if replace:
            self.data = data
        return data
    
    def query_gwas_catalog(
        self,
        p_threshold=5e-8,
        replace=True,
        max_associations=None,
        timeout=-1,
    ):
        """
        Queries the GWAS Catalog REST API and adds structured association annotations for each SNP.
        
        Args:
            p_threshold (float, optional): Only associations that are at least as significant are reported. Default is 5e-8.
                P-values are compared exactly, using the mantissa and exponent reported by the GWAS Catalog
                (the float `p_value` returned in ASSOC is 0.0 for p-values below ~1e-308).
            replace (bool, optional): If True, updates the data attribute in place. Default is True.
            max_associations (int, optional): If not None, only the `max_associations` qualifying trait-study associations with the smallest exact p-values are reported for each SNP after filtering and deduplication. Default is None.
            timeout (int, optional): Timeout for each HTTP request in seconds. Default is -1, which uses 30 seconds per request. Choose None for no timeout.

        Returns:
            pd.DataFrame: Data attribute with additional columns:
                - ASSOC: list of dictionaries with keys `trait`, `p_value`, and `study_accession`
                - ASSOC_STATUS: query outcome status
                - ASSOC_ERROR: short error detail or None
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
            f"Querying the GWAS Catalog and creating the ASSOC, ASSOC_STATUS, and ASSOC_ERROR columns. \n"
            f"Only associations with a p-value <= {p_threshold} are reported. Use the p_threshold argument to change the threshold."
        )
        if max_associations:
            print(
                f"Reporting at most {max_associations} qualifying trait-study associations for each SNP after filtering and deduplication."
            )
        print("Normalizing SNP identifiers to canonical rsIDs and rejecting non-rs inputs.")
        print(
            f"The .data attribute will {'be' if replace else 'not be'} modified. "
            f"{'Use replace=False to leave it as is.' if replace else ''}"
            )
        
        timeout = 30 if timeout == -1 else timeout

        from .snp_query import async_query_gwas_catalog

        # Call the async function to query all SNPs
        outcomes, summary = async_query_gwas_catalog(
            data.SNP.to_list(), 
            p_threshold=p_threshold, 
            max_associations=max_associations,
            timeout=timeout,
        )
        
        data["ASSOC"] = [outcome["associations"] for outcome in outcomes]
        data["ASSOC_STATUS"] = [outcome["status"] for outcome in outcomes]
        data["ASSOC_ERROR"] = [outcome["error"] for outcome in outcomes]

        status_counts = summary["status_counts"]
        retry_counts = summary["retry_counts"]
        deduplicated_queries = summary["valid_count"] - summary["unique_query_count"]

        print(
            f"Normalized {summary['requested_count']} SNPs to {summary['unique_query_count']} unique GWAS Catalog rsID queries."
        )
        print(
            f"Rejected {summary['invalid_count']} SNPs with invalid rsID format and avoided {deduplicated_queries} duplicate GWAS Catalog queries."
        )
        if summary["retry_total"]:
            retry_text = ", ".join(
                f"{status}={retry_counts[status]}"
                for status in sorted(retry_counts)
            )
            print(f"Retried GWAS Catalog requests {summary['retry_total']} times ({retry_text}).")
        else:
            print("No GWAS Catalog retries were needed.")

        print("The ASSOC, ASSOC_STATUS, and ASSOC_ERROR columns have been successfully created.")
        print(
            "Final GWAS Catalog status counts: "
            + ", ".join(
                f"{status}={status_counts[status]}"
                for status in [
                    "ok",
                    "no_associations",
                    "invalid_format",
                    "not_found",
                    "timeout",
                    "rate_limited",
                    "server_error",
                    "client_error",
                ]
            )
            + "."
        )
        if timeout is not None:
            print(
                f"Per-request timeout was set to {timeout:.2f} seconds. Increase the timeout argument if many requests finish with ASSOC_STATUS='timeout'."
            )
            
        return data
        
    def sort_group(self, method="lowest_p"):
        """
        Handle duplicate SNPs. Useful if the instance combines different Genos.

        Args:
            method (str, optional): How to handle duplicates. Default is "lowest_p",
                                    which keeps, for each SNP, the complete row with the lowest P-value
                                    (rows with a missing SNP ID are all kept). The data is then sorted by P-value.

        Returns:
            None
        """
        if method == "lowest_p":
            data = self.data.sort_values(by=["P"], kind="stable")
            duplicated = data["SNP"].notna() & data.duplicated(subset=["SNP"], keep="first")
            self.data = data[~duplicated].reset_index(drop=True)
            n_del = duplicated.sum()
            print(
                f"{n_del} rows with duplicated SNP IDs have been removed, keeping the row with the lowest P-value for each SNP."
            )
        return

    def copy(self, data=None):
        """
        Create another Geno instance from `data` (by default, from the data attribute).

        The checks, phenotype, reference_panel and reference_panel_name attributes are carried over (the phenotype and
        reference panel are shared with this instance, not copied). The analysis results (MR_data, MR_results,
        MR_loo_results, MRpresso_results, MRpresso_subset_data, MRpresso_outliers) are not carried over, and the created
        instance gets its own name.

        Args:
            data (pd.DataFrame, optional): Data of the created instance. Defaults to the data attribute.

        Returns:
            Geno: The created instance.
        """
        if data is None:
            data = self.data
        Geno_copy = Geno(data, keep_columns=True)
        Geno_copy.checks = self.checks.copy()
        if hasattr(self, "phenotype"):
            Geno_copy.phenotype = self.phenotype
        if hasattr(self, "reference_panel"):
            Geno_copy.reference_panel = self.reference_panel
        if hasattr(self, "reference_panel_name"):
            Geno_copy.reference_panel_name = self.reference_panel_name
        return Geno_copy

    def save(self, path="", fmt="parquet", sep="\t", header=True):
        """
        Save the Geno data to the file `<path>/<name>.<fmt>`, `name` being the name attribute of the object.

        Args:
            path (str, optional): Folder path to save the file. Defaults to the current directory.
            fmt (str, optional): File format: "parquet" (default), "csv" or "txt".
                Parquet is saved with `index=True`, `compression="zstd"`, `engine="pyarrow"`.
            sep (str, optional): Delimiter for .csv and .txt formats. Default is tab.
            header (bool, optional): Save column names for .csv and .txt formats. Default is True.
        """
        save_data(self.data, name=self.name, path=path, fmt=fmt, sep=sep, header=header)
        return
    
    def standardize_betas(self, outcome_sd=None, n_samples=None, replace=True):
        """
        Standardizes the BETA and SE columns of the Geno object's data.

        This method standardizes the effect sizes (BETA) and their standard errors (SE)
        using either a provided standard deviation of the outcome trait or by estimating
        it from the summary statistics using the total sample size.

        Args:
            outcome_sd (float, optional): The standard deviation of the outcome trait.
                If provided, this value is used to standardize BETA and SE.
                Defaults to None.
            n_samples (int, optional): The total number of samples in the GWAS.
                This is used to estimate the outcome's standard deviation if `outcome_sd`
                is not provided. `EAF` and `SE` columns must be present in the data.
                Defaults to None.
            replace (bool, optional): If True (default), updates the data attribute in place. If False, operates on a copy.

        Returns:
            pd.DataFrame or None: Data after being standardized. None (data not modified) if the standard deviation
            cannot be estimated (no row with both EAF and SE) or is 0.

        Raises:
            ValueError: If neither `outcome_sd` nor `n_samples` is provided.
            ValueError: If `n_samples` is provided, but 'EAF' or 'SE' columns
                are missing from the data.

        Note:
            Valid only for quantitative traits (a warning is printed).
        """
        print("Warning: Standardizing estimates for binary traits will make the MR results not interpretable as odds ratios.")
        if outcome_sd is not None:
            sd = outcome_sd
        elif n_samples is not None:
            required_cols = ['EAF', 'SE']
            if not all(col in self.data.columns for col in required_cols):
                raise ValueError("When `n_samples` is provided, 'EAF' and 'SE' columns are required.")

            temp_df = self.data[['EAF', 'SE']].dropna()

            if temp_df.empty:
                warnings.warn("No valid EAF and SE values to estimate outcome_sd, standardization skipped.")
                return

            var_snp = 2 * temp_df['EAF'] * (1 - temp_df['EAF'])
            var_outcome_estimates = temp_df['SE']**2 * n_samples * var_snp
            
            sd_outcome_estimates = np.sqrt(var_outcome_estimates)
            sd = np.median(sd_outcome_estimates)
            print(f"Estimated the trait standard deviation from EAF and SE columns: {sd:.4f}.")

        else:
            raise ValueError("Either `outcome_sd` or `n_samples` must be provided to standardize betas.")

        if sd == 0:
            warnings.warn("Outcome standard deviation is zero, standardization cannot be performed.")
            return
            
        if replace:
            self.data['BETA'] = self.data['BETA'] / sd
            self.data['SE'] = self.data['SE'] / sd
            return self.data
        else:
            data = self.data.copy()
            data['BETA'] = data['BETA'] / sd
            data['SE'] = data['SE'] / sd
            return data

    def update_eaf(self, reference_panel="EUR_37", replace=True, fill=True):
        """
        Update or create the EAF (Effect Allele Frequency) column using a reference panel.

        This method calculates allele frequencies from a specified reference panel using PLINK
        and updates the 'EAF' column for SNPs in the dataset. It matches SNPs based on
        CHR/POS (or SNP ID if CHR/POS are absent) and on the alleles: a reference variant is used only
        if {EA, NEA} = {REF, ALT} (or, if there is no NEA column, if EA is one of its alleles).
        SNPs matching several reference variants with different frequencies are left unmatched.
        The number of SNPs absent from the reference panel (no reference variant at their position, or SNP ID not found)
        is reported. A warning is issued when less than 90% of the SNPs found in the reference panel match its alleles
        (likely genome build mismatch).

        Args:
            reference_panel (str, optional): The reference panel to use for deriving EAF.
                Can be a standard name (e.g., "EUR_37", "AFR_38") or a path to a
                custom PLINK fileset (bed/bim/fam or pgen/pvar/psam).
                Defaults to "EUR_37".
            replace (bool, optional): If True (default), modifies the instance's `data` attribute
                in place. If False, operates on a copy.
            fill (bool, optional): If `True` (default), existing `EAF` values for SNPs not
                found in the reference panel will be preserved. If `False`, `EAF` values
                for unmatched SNPs will be set to `NaN`.

        Returns:
            pd.DataFrame: The data with the updated 'EAF' column (the data attribute itself if replace=True).

        Raises:
            ValueError: If the required columns ('EA' and either 'SNP' or 'CHR'/'POS')
                are not present in the data.
        """

        if 'EA' not in self.data.columns:
            raise ValueError("The 'EA' column is required to update EAF.")
        if not ('SNP' in self.data.columns or ('CHR' in self.data.columns and 'POS' in self.data.columns)):
            raise ValueError("Either 'SNP' or both 'CHR' and 'POS' columns are required.")

        data = self.data if replace else self.data.copy()

        data_updated = update_eaf_func(
            data=data,
            reference_panel=reference_panel,
            object_name=self.name,
            ram=self.ram,
            fill=fill
        )

        if replace:
            self.data = data_updated
        return data_updated



MATCH_COUNTS = ["matched", "opposite_strand", "ambiguous", "no_allele_match", "no_position_match", "own_id", "own_id_other_match"]


def match_variant_ids(data, variants):
    """
    Find the ID of the variant of the genetic data corresponding to each row of data, based on CHR/POS and alleles.

    A row matches a variant at the same CHR/POS whose alleles are EA and NEA (in any order). For SNVs that are not
    palindromic (not A/T or C/G), the strand-complemented pair is also accepted. If NEA is missing, EA must be an
    allele of the variant and the variant must be the only one at this position (without EA, the variant must be the
    only one at this position). Rows matching several variants with different IDs are ambiguous and left unmatched.
    Rows whose own SNP ID is present in the genetic data keep it, whatever the CHR/POS match.

    Args:
        data (pd.DataFrame): Data with the CHR and POS columns (SNP, EA and NEA are used if present).
        variants (pd.DataFrame): Variants with the CHR, POS, SNP, A1, A2 columns (see :func:`genal.tools.read_variant_file`).

    Returns:
        tuple: IDs (np.ndarray, NaN for rows without an ID from the genetic data), rows whose alleles must be
            complemented to the strand of the genetic data (np.ndarray of booleans) and counts of rows per
            matching outcome (dict).
    """
    from .tools import normalize_chr, COMPLEMENT

    n = data.shape[0]
    def alleles(col):
        if col not in data.columns:
            return np.full(n, "", dtype=object)
        return data[col].fillna("").astype(str).str.upper().to_numpy()
    rows = pd.DataFrame({
        "row": np.arange(n),
        "CHR": normalize_chr(data["CHR"]).to_numpy(dtype=float),
        "POS": pd.to_numeric(data["POS"], errors="coerce").to_numpy(dtype=float),
        "EA": alleles("EA"),
        "NEA": alleles("NEA"),
    })
    genetic_ids = variants["SNP"].dropna()
    # Only the variants at the positions of the rows can match: the others are dropped first (all the variants at these
    # positions are kept, so the number of variants per position is unchanged)
    variants = variants.loc[variants["POS"].isin(rows["POS"].dropna().unique()), ["CHR", "POS", "SNP", "A1", "A2"]]
    variants = variants.astype({"CHR": float, "POS": float})
    variants["n_at_position"] = variants.groupby(["CHR", "POS"])["POS"].transform("size")
    variants = variants.dropna(subset=["SNP"])
    m = rows.merge(variants, on=["CHR", "POS"], how="inner")
    m["A1"] = m["A1"].astype(str).str.upper()
    m["A2"] = m["A2"].astype(str).str.upper()

    def in_variant(allele):
        """Whether each allele is one of the alleles of the variant (ALT can list several alleles)."""
        hit = (allele == m["A1"]) | (allele == m["A2"])
        multi = m["A2"].str.contains(",", regex=False)
        if multi.any():
            hit[multi] = [a in [ref, *alt.split(",")] for a, ref, alt in zip(allele[multi], m["A1"][multi], m["A2"][multi])]
        return hit & (allele != "")

    has_ea, has_nea = m["EA"] != "", m["NEA"] != ""
    ea_in, nea_in = in_variant(m["EA"]), in_variant(m["NEA"])
    single = m["n_at_position"] == 1
    direct = np.where(has_nea, ea_in & nea_in & (m["EA"] != m["NEA"]), np.where(has_ea, ea_in, True) & single)
    ambiguous_position = ~has_nea & np.where(has_ea, ea_in, True) & ~single
    comp_ea, comp_nea = m["EA"].map(COMPLEMENT).fillna(""), m["NEA"].map(COMPLEMENT).fillna("")
    flippable = has_nea & (comp_ea != "") & (comp_nea != "") & (comp_ea != m["NEA"])
    opposite = flippable & in_variant(comp_ea) & in_variant(comp_nea) & ~direct

    def unique_id(mask):
        """IDs of the rows matching variants with a single ID, and rows matching variants with several IDs."""
        ids = m.loc[mask].groupby("row")["SNP"]
        n_ids = ids.nunique()
        return ids.first()[n_ids == 1], n_ids.index[n_ids > 1]

    direct_ids, direct_ambiguous = unique_id(direct)
    undecided = ~m["row"].isin(direct_ids.index) & ~m["row"].isin(direct_ambiguous)
    opposite_ids, opposite_ambiguous = unique_id(opposite & undecided)

    snp_ids = np.full(n, np.nan, dtype=object)
    snp_ids[direct_ids.index] = direct_ids.to_numpy()
    snp_ids[opposite_ids.index] = opposite_ids.to_numpy()
    flipped = np.zeros(n, dtype=bool)
    flipped[opposite_ids.index] = True

    n_matched = len(direct_ids) + len(opposite_ids)
    n_ambiguous = len(set(direct_ambiguous) | set(opposite_ambiguous) | set(m.loc[ambiguous_position, "row"]))
    n_at_position = m["row"].nunique()
    counts = {
        "matched": n_matched,
        "opposite_strand": len(opposite_ids),
        "ambiguous": n_ambiguous,
        "no_allele_match": n_at_position - n_matched - n_ambiguous,
        "no_position_match": n - n_at_position,
        "own_id": 0,
        "own_id_other_match": 0,
    }

    # Rows whose own SNP ID is present in the genetic data keep it (and its strand); missing or placeholder IDs never match
    if "SNP" in data.columns:
        from .geno_tools import missing_snp_ids
        own_ids = data["SNP"].to_numpy(dtype=object)
        own_kept = (data["SNP"].isin(genetic_ids) & ~missing_snp_ids(data["SNP"])).to_numpy()
        same_id = pd.notna(snp_ids) & (snp_ids == own_ids)
        counts["own_id"] = int(own_kept.sum())
        counts["own_id_other_match"] = int((own_kept & pd.notna(snp_ids) & ~same_id).sum())
        flipped &= ~own_kept | same_id
        snp_ids = np.where(own_kept, own_ids, snp_ids)
    return snp_ids, flipped, counts


def available_cpus():
    """
    Number of CPUs the process may use, at least 1: the CPUs of its affinity mask (e.g. the CPUs allocated by SLURM)
    where Python exposes it (Python 3.13+, or Linux), otherwise all the CPUs of the machine.
    """
    if hasattr(os, "process_cpu_count"):
        n_cpus = os.process_cpu_count()
    elif hasattr(os, "sched_getaffinity"):
        n_cpus = len(os.sched_getaffinity(0))
    else:
        n_cpus = os.cpu_count()
    return max(1, n_cpus or 1)


def match_chromosome_variant_ids(df_subset, path, filetype):
    """Helper function of the update_snpids method to match the SNPs of one chromosome when genetic data is split by chromosome."""
    from .tools import read_variant_file, normalize_chr, split_chromosome_path

    chrom = normalize_chr(df_subset["CHR"].iloc[:1]).iloc[0]
    prefix = split_chromosome_path(path, int(chrom), filetype) if pd.notna(chrom) else None
    if prefix is None:
        n = df_subset.shape[0]
        counts = dict.fromkeys(MATCH_COUNTS, 0)
        counts["no_position_match"] = n
        return np.full(n, np.nan, dtype=object), np.zeros(n, dtype=bool), counts
    return match_variant_ids(df_subset, read_variant_file(prefix + (".bim" if filetype == "bed" else ".pvar")))
