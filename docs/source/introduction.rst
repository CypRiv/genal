============
Introduction
============

Genal is a python module designed to make it easy to run genetic risk scores and mendelian randomization analyses. It integrates a collection of tools that facilitate the cleaning of single nucleotide polymorphisms (SNP) data and enable the execution of clinical population genetic workflows. The functionalities provided by genal include clumping, lifting, association testing, Polygenic Risk Scoring, and Mendelian Randomization analyses (including a parallel python implementation of MR-PRESSO), all within a single Python module.

The module prioritizes user-friendliness and intuitive operation, aiming to reduce the complexity of genetic analysis for researchers. Despite its focus on simplicity, Genal does not sacrifice customization or the precision of analysis. Researchers can expect to maintain analytical rigour while benefiting from the streamlined python experience.

This page provides installation steps followed by a short tutorial in using genal.

============
Installation
============

The genal package can be easily installed with pip::

    pip install genal

The main genal functionalities require a working installation of PLINK v1.9 (and not 2.0 as certain functionalities have not been updated yet) that can be downloaded here: https://www.cog-genomics.org/plink/ 
Once downloaded, the path to the plink 1.9 executable should be set with::

    genal.set_plink(path="/path/to/plink/executable/file")


========
Tutorial
========

For the purpose of this tutorial, we are going to build a PRS of systolic blood pressure (SBP) and investigate the genetically-determined effect of SBP on the risk of stroke. We will use both summary statistics from Genome-Wide Association Studies (GWAS) and individual-level data from the UK Biobank as our test population. We are going to go through the following steps:

* Data loading 
* Data preprocessing
    * Perform a series of checks and cleaning operations on the SNP-level data to ensure it is in appropriate format
* Clump the data to obtain genome-wide significant and independent SNPs that will serve as our genetic instruments for SBP
* Build a genomic risk score of SBP in a test population
    * Build a risk score with proxies
* Perform Mendelian Randomization
    * Run a mendelian randomization analysis of SBP as exposure and acute ischemic stroke as outcome (using summary statistics for both SBP and stroke)
    * Run sensitivity analyses of the previous mendelian randomization using the weighted median, MR-Egger, and MR-PRESSO 
* Calibrate SNP-trait weights with individual level genetic data
    * Run single-SNP association tests to calibrate the SBP genetic instruments in a test population
* Show some utility tools
    * Lift the data to another genomic build
        * In pure python
        * Using LiftOver
    * Phenoscanner (todo)


Data loading
============

We start this tutorial with publicly available summary statistics data from a large GWAS of systolic blood pressure (https://www.nature.com/articles/s41588-018-0205-x). After downloading and unzipping the summary statistics, we load the data in a pandas dataframe::

    import pandas as pd
    sbp_gwas = pd.read_csv("Evangelou_30224653_SBP.txt",sep=" ")

The loaded dataframe::

    >>> sbp_gwas.head(5)
             MarkerName Allele1 Allele2   Freq1  Effect  StdErr         P  TotalSampleSize  N_effective
    0  10:100000625:SNP       a       g  0.5660  0.0523  0.0303  0.083940           738170       736847
    1  10:100000645:SNP       a       c  0.7936  0.0200  0.0372  0.591100           738168       735018
    2  10:100003242:SNP       t       g  0.8831  0.1417  0.0469  0.002526           738168       733070
    3  10:100003304:SNP       a       g  0.9609  0.0245  0.0838  0.769800           737054       663809
    4  10:100003785:SNP       t       c  0.6406 -0.0680  0.0313  0.029870           738169       735681

We can now load this data in a :class:`~genal.Geno` object. The :class:`~genal.Geno` class is the central piece of the package. It is made to store Single Nucleotide Polymorphisms (SNP) data and make it easy to preprocess and clean. 
The :class:`~genal.Geno` takes as input a pandas dataframe where each row corresponds to a SNP, and with columns describing the position and possibly the effect of the SNP for the given trait (SBP in our case). To indicate the name of the columns, the following arguments can be passed:
* CHR: Column name for chromosome. Defaults to "CHR".
* POS: Column name for genomic position. Defaults to "POS".
* SNP: Column name for SNP identifier (rsid). Defaults to "SNP".
* EA: Column name for effect allele. Defaults to "EA".
* NEA: Column name for non-effect allele. Defaults to "NEA".
* BETA: Column name for effect estimate. Defaults to "BETA".
* SE: Column name for effect standard error. Defaults to "SE".
* P: Column name for effect p-value. Defaults to "P".
* EAF: Column name for effect allele frequency. Defaults to "EAF".

In our case, and after inspecting the dataframe, we must first extract the chromosome and position information from the "MarkerName" column into two new columns "CHR" and "POS"::

    sbp_gwas[["CHR", "POS", "Filler"]] = sbp_gwas["MarkerName"].str.split(":", expand=True)
    
The resulting dataframe has now separate columns for the chromosome and genomic position::

    >>> sbp_gwas.head(5)
             MarkerName Allele1 Allele2   Freq1  Effect  StdErr         P  TotalSampleSize  N_effective CHR        POS Filler
    0  10:100000625:SNP       a       g  0.5660  0.0523  0.0303  0.083940           738170       736847  10  100000625    SNP
    1  10:100000645:SNP       a       c  0.7936  0.0200  0.0372  0.591100           738168       735018  10  100000645    SNP
    2  10:100003242:SNP       t       g  0.8831  0.1417  0.0469  0.002526           738168       733070  10  100003242    SNP
    3  10:100003304:SNP       a       g  0.9609  0.0245  0.0838  0.769800           737054       663809  10  100003304    SNP
    4  10:100003785:SNP       t       c  0.6406 -0.0680  0.0313  0.029870           738169       735681  10  100003785    SNP

and it can now be loaded in a :class:`~genal.Geno` object::
    
    import genal
    SBP_Geno = genal.Geno(sbp_gwas, CHR = "CHR", POS = "POS", EA = "Allele1", NEA = "Allele2", BETA = "Effect", SE = "StdErr", P = "P", EAF = "Freq1", keep_columns = False)
    
The last argument (keep_columns = False) indicates that we do not wish to keep the other (non-main) columns in the dataframe.

.. note::

   Make sure to read the readme file usually provided with the summary statistics to identify the correct columns. It is particularly important to correctly identify the allele that represents the effect allele. Also, you do not need all columns to move forward, as some can be inputed as we will see next.
   
   
Data preprocessing
==================

Now that we have loaded the data in a :class:`~genal.Geno` object, we can clean and format it. Methods such as Polygenic Risk Scoring or Mendelian Randomization require the SNP data to be in a specific format. Also, raw summary statistics can sometimes contain missing or invalid values that need to be handled. Finally, some columns may be missing from the data (such as the SNP rsid column, or the non-effect allele column...) and these columns can be created based on existing ones and a reference panel. 
Genal can run all the basic cleaning and preprocessing steps in one command::

    SBP_Geno.preprocess_data(preprocessing = 2)

The preprocessing argument specifies the global level of preprocessing applied to the data:
* preprocessing = 0 means that the data won't be modified.
* preprocessing = 1 means that missing columns will be added based on reference data and invalid values set to NaN, but no rows will be deleted.
* preprocessing = 2 means that missing columns will be added, and all rows containing missing, duplicated, or invalid values will be deleted. This option is recommended before running genetic methods.

By default, and depending on the global preprocessing level (0, 1, 2) chosen, the :meth:`~genal.Geno.preprocess_data` method will run the following checks:
* Make sure the CHR (chromosome) and POS (genomic position) columns are integers.
* Ensure the EA (effect allele) and NEA (non-effect allele) columns are upper case characters containing A, T, C, G letters. Set multiallelic values to NaN.
* Ensure the P (p-value) column contains valid values.
* Ensure there are no duplicated SNPs based on rsid.
* Determines if the BETA (effect) column contains beta estimates or odds ratios. Log-transform odds ratios if necessary.
* Create SNP column using a reference panel if CHR and POS columns are present.
* Create CHR and/or POS column using a reference panel if SNP column is present.
* Create NEA (non-effect allele) column using a reference panel if EA (effect allele) column is present.
* Create the SE (standard-error) column if the BETA and P (p-value) columns are present.
* Create the P column if the BETA and SE columns are present.

If you do not wish to run certain steps, or wish to run only certain steps, you can use additional arguments. Please refer to the API for more info: :meth:`~genal.Geno.preprocess_data`.

In our case, the SNP column (for SNP identifier - rsid) was missing from our dataframe and has been added based on a 1000 genome reference panel::

    Using the EUR reference panel.
    The SNP column (rsID) has been created. 197511(2.787%) SNPs were not found in the reference data and their ID set to CHR:POS:EA.
    The BETA column looks like Beta estimates. Use effect_column='OR' if it is a column of Odds Ratios.
    
You can always check the data of a :class:`~genal.Geno` object by accessing the 'data' attribute::

    >>> SBP_Geno.data
            EA NEA     EAF    BETA      SE         P  CHR        POS         SNP
    0        A   G  0.5660  0.0523  0.0303  0.083940   10  100000625   rs7899632
    1        A   C  0.7936  0.0200  0.0372  0.591100   10  100000645  rs61875309
    2        T   G  0.8831  0.1417  0.0469  0.002526   10  100003242  rs12258651
    3        A   G  0.9609  0.0245  0.0838  0.769800   10  100003304  rs72828461
    4        T   C  0.6406 -0.0680  0.0313  0.029870   10  100003785   rs1359508
    ...     ..  ..     ...     ...     ...       ...  ...        ...         ...
    7088116  A   G  0.4474  0.0114  0.0307  0.710600    9   99997049  rs10817273
    7088117  T   C  0.4529 -0.0151  0.0307  0.623900    9   99997707  rs11794422
    7088118  T   C  0.4511 -0.0224  0.0306  0.463900    9   99998403  rs10981296
    7088119  C   G  0.8094 -0.0204  0.0389  0.599800    9   99998646  rs10981297
    7088120  A   G  0.9028 -0.0184  0.0517  0.722300    9   99999468  rs10981301
    
And we see that the SNP column with the rsids has been added based on the reference data.
You do not need to obtain the 1000 genome reference panel yourself, Genal will download it the first time you use it. By default, the reference panel used is the european (eur) one. You can specify another valid reference panel (afr, eas, sas, amr) with the reference_panel argument::

    SBP_Geno.preprocess_data(preprocessing = 2, reference_panel = "afr")
    
You can also use a custom reference panel by specifying to the reference_panel argument a path to bed/bim/fam files (without the extension).


Clumping
========

Clumping is the step at which we select the SNPs that will be used as our genetic instruments in future Polygenic Risk Scores and Mendelian Randomization analyses. The process involves identifying the SNPs that are strongly associated with our trait of interest (systolic blood pressure in this tutorial) and are independent from each other. This second step ensures that selected SNPs are not highly correlated, (i.e. they are not in close linkage disequilibrium). For this step, we again need to use a reference panel. The SNP-data loaded in a :class:`~genal.Geno` object can be clumped using the :meth:`~genal.Geno.clump` method. It will return another :class:`~genal.Geno` object containing only the clumped data::

    SBP_clumped = SBP_Geno.clump(p1 = 5e-8, r2 = 0.1, kb = 250, reference_panel = "eur")
    
It will output the number of instruments obtained::

    Using the EUR reference panel.
    Warning: 760  top variant IDs missing
    1545 clumps formed from 73594 top variants.
    
You can specify the thresholds you want to use for the clumping with the following arguments:
* p1: P-value threshold during clumping. SNPs with a P-value higher than this value are excluded.
* r2: Linkage disequilibrium threshold for the independence check. Takes values between 0 and 1.
* kb: Genomic window used for the independence check (the unit is thousands of base-pair positions).
* reference_panel: The reference population used to derive linkage disequilibrium values and select independent SNPs

Polygenic Risk Scoring
======================

Computing a Polygenic Risk Score (PRS) can be done in one line with the :meth:`~genal.Geno.prs` method::

    SBP_clumped.prs(name = "SBP_prs", path = path/to/genetic/files)
    
The genetic files of the target population can be either one triple of bed/bim/fam files containing information for all SNPs, or they can be divided by chromosome (one bed/bim/fam triple for chr 1, another for chr 2 etc...). In the latter case, provide the path by replacing the chromosome number by '$' and Genal will extract the necessary SNPs from each chromosome and merge them before running the PRS. For instance, if the genetic files are named Pop_chr1.bed, Pop_chr1.bim, Pop_chr1.fam, Pop_chr2.bed ..., you can use::

    SBP_clumped.prs(name = "SBP_prs", path = "Pop_chr$")
    
The name argument specifies the name of the .csv file that will be saved with the individual risk scores. 
The output of the :meth:`~genal.Geno.prs` method will include how many SNPs were used to compute the risk score. It can happen that some of the SNPs are multiallelic in the genetic data (even if they are not multiallelic in our SNP data) and need to be excluded. It can also happen that some of the SNPs are missing from the genetic files of the target population (for instance if the data has not been imputed)::

    CHR/POS columns present: SNPs searched based on genomic positions.
    Extracting SNPs for each chromosome...
    SNPs extracted for chr1.
    SNPs extracted for chr2.
    SNPs extracted for chr3.
    SNPs extracted for chr4.
    SNPs extracted for chr5.
    SNPs extracted for chr6.
    SNPs extracted for chr7.
    SNPs extracted for chr8.
    SNPs extracted for chr9.
    SNPs extracted for chr10.
    SNPs extracted for chr11.
    SNPs extracted for chr12.
    SNPs extracted for chr13.
    SNPs extracted for chr14.
    SNPs extracted for chr15.
    SNPs extracted for chr16.
    SNPs extracted for chr17.
    SNPs extracted for chr18.
    SNPs extracted for chr19.
    SNPs extracted for chr20.
    SNPs extracted for chr21.
    SNPs extracted for chr22.
    Merging SNPs extracted from each chromosome...
    Created bed/bim/fam fileset with extracted SNPs: tmp_GENAL/4f4ce6a7_allchr
    Extraction completed. 786(50.874%) SNPs were not extracted from the genetic data.
    Computing a weighted PRS using tmp_GENAL/4f4ce6a7_allchr data.
    The PRS computation was successful and used 759/1545 (49.126%) SNPs.
    PRS data saved to SBP_prs.csv

Here, we see that about half of the SNPs were not extracted from the data. In such cases, we may want to try and salvage some of these SNPs by looking for proxies (SNPs in high linkage disequilibrium, i.e. highly correlated SNPs). This can be done by specifying the 'proxy = True' argument::

    SBP_clumped.prs(name = "SBP_prs" ,path = "Pop_chr$", proxy = True, reference_panel = "eur", r2=0.8, kb=5000, window_snps=5000)
    
and the output is::

    CHR/POS columns present: SNPs searched based on genomic positions.
    Identifying the SNPs present in the genetic data...
    759 SNPs out of 1545 are present in the genetic data.
    Searching proxies for 786 SNPs...
    Using the EUR reference panel.
    Filtering the potential proxies with the searchspace provided.
    Found proxies for 578 missing SNPs.
    7(0.455%) duplicated SNPs have been removed. Use keep_dups=True to keep them.
    Extracting SNPs for each chromosome...
    SNPs extracted for chr1.
    SNPs extracted for chr2.
    SNPs extracted for chr3.
    SNPs extracted for chr4.
    SNPs extracted for chr5.
    SNPs extracted for chr6.
    SNPs extracted for chr7.
    SNPs extracted for chr8.
    SNPs extracted for chr9.
    SNPs extracted for chr10.
    SNPs extracted for chr11.
    SNPs extracted for chr12.
    SNPs extracted for chr13.
    SNPs extracted for chr14.
    SNPs extracted for chr15.
    SNPs extracted for chr16.
    SNPs extracted for chr17.
    SNPs extracted for chr18.
    SNPs extracted for chr19.
    SNPs extracted for chr20.
    SNPs extracted for chr21.
    SNPs extracted for chr22.
    Merging SNPs extracted from each chromosome...
    Created bed/bim/fam fileset with extracted SNPs: tmp_GENAL/4f4ce6a7_allchr
    Extraction completed. 208(13.524%) SNPs were not extracted from the genetic data.
    Computing a weighted PRS using tmp_GENAL/4f4ce6a7_allchr data.
    The PRS computation was successful and used 1330/1538 (86.476%) SNPs.
    PRS data saved to SBP_prs.csv
    
In our case, we have been able to find proxies for 571 of the 786 SNPs that were missing in the population genetic data (7 potential proxies have been removed because they were identical to SNPs already present in our data).
You can customize how the proxies are chosen with the following arguments:
* reference_panel: The reference population used to derive linkage disequilibrium values and find proxies
* kb: Width of the genomic window to look for proxies (in thousands of base-pair positions)
* r2: Minimum linkage disequilibrium value with the original SNP for a proxy to be included
* window_snps: Width of the window to look for proxies (in number of SNPs)

.. note::

You can call the :meth:`~genal.Geno.prs` method on any Geno object (containing at least the EA, BETA and either SNP or CHR/POS columns). The data does not need to be clumped and there is no limit to the number of instruments used to compute the scores.


Mendelian Randomization
=======================

To run MR, we need to load both our exposure and outcome SNP-level data in :class:`~genal.Geno` objects. In our case, the genetic instruments of the MR are the SNPs associated with blood pressure at genome-wide significant levels resulting from the clumping of the blood pressure GWAS. They are stored in our SBP_clumped :class:`~genal.Geno` object which also include their association with the exposure trait (instrument-SBP estimates in the BETA column).

To get their association with the outcome trait (instrument-stroke estimates), we are going to use SNP-level data from a large GWAS of stroke performed by the GIGASTROKE consortium (https://www.nature.com/articles/s41586-022-05165-3)::

    stroke_gwas = pd.read_csv("GCST90104539_buildGRCh37.tsv",sep="\t")
    
We inspect it to determine the column names::

    >>> stroke_gwas.head(5)
       chromosome  base_pair_location  effect_allele_frequency    beta  standard_error  p_value  odds_ratio  ci_lower  ci_upper effect_allele other_allele
    0           5            29439275                   0.3569  0.0030          0.0070   0.6658    1.003005  0.989337  1.016861             T            C
    1           5            85928892                   0.0639 -0.0152          0.0137   0.2686    0.984915  0.958820  1.011720             T            C
    2          10           128341232                   0.4613  0.0025          0.0065   0.6998    1.002503  0.989812  1.015357             T            C
    3           3            62707519                   0.0536  0.0152          0.0152   0.3177    1.015316  0.985514  1.046019             T            C
    4           2            80464120                   0.9789  0.0057          0.0254   0.8223    1.005716  0.956874  1.057052             T            G
    
We load it in a :class:`~genal.Geno` object::

    Stroke_Geno = genal.Geno(stroke_gwas, CHR = "chromosome", POS = "base_pair_location", EA = "effect_allele", NEA = "other_allele", BETA = "beta", SE = "standard_error", P = "p_value", EAF = "effect_allele_frequency", keep_columns = False)
    
We preprocess it as well to put it in the correct format and make sure there is no invalid values::

    Stroke_Geno.preprocess_data(preprocessing = 2)
    
Now, we need to extract our instruments (SNPs of the SBP_clumped data) from the outcome data to obtain their association with the outcome trait (stroke). It can be done by calling the :meth:`~genal.Geno.query_outcome` method::

    SBP_clumped.query_outcome(Stroke_Geno, proxy = False)
    
Genal will print how many SNPs were successfully found and extracted from the outcome data::

    Outcome data successfully loaded from 'b352e412' geno object.
    Identifying the exposure SNPs present in the outcome data...
    1541 SNPs out of 1545 are present in the outcome data.
    (Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute.
    
Here as well you have the option to use proxies for the instruments that are not present in the outcome data::

    SBP_clumped.query_outcome(Stroke_geno, proxy = True, reference_panel = "eur", kb = 5000, r2 = 0.6, window_snps = 5000)

And Genal will print the number of missing instruments which have been proxied::

    Outcome data successfully loaded from 'b352e412' geno object.
    Identifying the exposure SNPs present in the outcome data...
    1541 SNPs out of 1545 are present in the outcome data.
    Searching proxies for 4 SNPs...
    Using the EUR reference panel.
    Found proxies for 4 SNPs.
    (Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute.
    
After extracting the instruments from the outcome data, the SBP_clumped :class:`~genal.Geno` object contains an 'MR_data' attribute containing the instruments-exposure and instruments-outcome associations necessary to run MR. Running MR is now as simple as calling the :meth:`~genal.Geno.MR` method of the SBP_clumped :class:`~genal.Geno` object::

    SBP_clumped.MR(action = 3, exposure_name = "SBP", outcome_name = "Stroke_eur")

The :meth:`~genal.Geno.MR` method returns a dataframe containing the estimates and p-values for different MR methods::

       exposure     outcome                                     method  nSNP         b        se      pval
    0       SBP  Stroke_eur                  Inverse-Variance Weighted  1312  0.023394  0.001132    <e-100
    1       SBP  Stroke_eur  Inverse Variance weighted (Fixed effects)  1312  0.023394  0.000807    <e-100
    2       SBP  Stroke_eur                      Unweighted regression  1312  0.021764  0.078648  0.781986
    3       SBP  Stroke_eur                            Weighted Median  1312  0.022891  0.001423    <e-100
    4       SBP  Stroke_eur                  Penalised weighted median  1312  0.021525  0.001432    <e-100
    5       SBP  Stroke_eur                              Simple median  1312  0.021480  0.001364    <e-100
    6       SBP  Stroke_eur                      Sign concordance test  1312  0.373476       NaN       0.0
    7       SBP  Stroke_eur                                   MR Egger  1312  0.029312  0.003063    <e-100
    8       SBP  Stroke_eur                            Egger Intercept  1312 -0.001799  0.000865  0.037777
    9       SBP  Stroke_eur                         MR Egger bootstrap  1312  0.030342  0.002093    <e-100
    10      SBP  Stroke_eur                  Egger Intercept bootstrap  1312 -0.002758  0.000740    0.0015
    
You can specify several arguments. We refer to the API for a full list, but the most important one is the 'action' argument. It determines how palindromic SNPs are treated during the exposure-outcome harmonization step. Palindromic SNPs are SNPs where the nucleotide change reads the same forward and backward on complementary strands of DNA (for instance EA = 'A' and NEA = 'T').
* action = 1: Palindromic SNPs are not treated (assumes all alleles are on the forward strand)
* action = 2: Uses effect allele frequencies to attempt to flip them (conservative, default)
* action = 3: Removes all palindromic SNPs (very conservative)
If you choose the option 2 or 3 (recommended), Genal will print the list of palindromic SNPs that have been removed from the analysis.
By default, all MR methods (inverse-variance weighted, weighted median, MR-Egger etc) are going to be run. But if you do not wish to run all of them, you can specify a 'methods' argument. More details in the :meth:`~genal.Geno.MR` API. 
For more fine-tuning, such as settings for the number of boostrapping iterations, please refer to the API.
    
If you wish to include the heterogeneity values (Cochran's Q), you can use the heterogeneity argument. Here, the heterogeneity for the inverse-variance weighted method::

    SBP_clumped.MR(action = 3, methods = ["IVW"], exposure_name = "SBP", outcome_name = "Stroke_eur", heterogeneity = True)
    
And that will give::

      exposure     outcome                     method  nSNP         b        se    pval            Q  Q_df  Q_pval
    0      SBP  Stroke_eur  Inverse-Variance Weighted  1312  0.023394  0.001132  <e-100  2583.740268  1311  <e-100
    
As expected, many MR methods indicate that SBP is strongly associated with stroke, but there are some signs of horizontal pleiotropy (instruments influencing the outcome through a different pathway than the one used as exposure) given the significant MR-Egger intercept p-value.
To investigate horizontal pleiotropy in more details, a very useful method is Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO). MR-PRESSO is a method designed to detect and correct for horizontal pleiotropy. It will identify which instruments are likely to be pleiotropic on their effect on the outcome, and it will rerun an inverse-variance weighted MR after excluding them. It can be run using the :meth:`~genal.Geno.MRpresso` method::

    SBP_clumped.MRpresso(action = 3)
    
As with the :meth:`~genal.Geno.MR` method, the action argument determines how the pleiotropic SNPs will be treated. The output is a list containing:
* A table containing the original and outlier-corrected inverse variance-weighted results
* The global test p-value indicating the presence of horizontal pleiotropy
* A dataframe of p-values, one for each instrument, representing the likelihood that this instrument is pleiotropic (only relevant if the global test is significant).
* A dictionary containing the outputs of the distortion test. This test assesses whether the removal of the pleiotropic instruments has significantly altered the original MR estimate.
    * An array containing the indices of the pleiotropic SNPs
    * The coefficient of the distortion test
    * The p-value of the distortion test


SNP-association testing
=======================

We may want to calibrate instrument-trait estimates in a specific population for which we have individual-level data (genetic files as well as phenotypic data). For instance, if the GWAS of SBP was done in a european population, we may want to adjust the estimates based on data coming from a population of a different ancestry. This can be done in 2 steps:
* Loading the phenotypic data in a dataframe and calling the :meth:`~genal.Geno.set_phenotype` method
* Calling the :meth:`~genal.Geno.association_test` method to run the association tests and update the estimates
Let's start by loading phenotypic data::

    df_pheno = pd.read_csv("path/to/trait/data")

.. note::

One important detail is to make sure that the individual IDs are identical between the phenotypic data and the genetic data for the target population.
Then, it is advised to make a copy of the :class:`~genal.Geno` object containing our instruments as we are going to update their coefficients and to avoid any confusion::

    SBP_adjusted = SBP_clumped.copy()

We can then call the :meth:`~genal.Geno.set_phenotype` method, specifying which column contains our trait of interest (for the association testing) and which column contains the IDs::

    SBP_adjusted.set_phenotype(df_pheno, PHENO = "htn", IID = "IID")

At this point, Genal will identify if the phenotype is binary or quantitative (to determine the regression model. If the phenotype is binary, it will assume that the most frequent value is coding for control (and the other value for case), this can be changed with 'alternate_control=True'::

    Detected a binary phenotype in the 'PHENO' column. Specify 'PHENO_type="quant"' if this is incorrect.
    Identified 0 as the control code in 'PHENO'. Set 'alternate_control=True' to inverse this interpretation.
    The phenotype data is stored in the .phenotype attribute.
    
We can then run the association tests, specifying the path to the genetic files in plink format, and any columns we may want to include as covariates in the regression tests::

    SBP_adjusted.association_test(covar=["age"], path = "path/to/genetic/files")

Genal will print information regarding the number of individuals used in the tests and the kind of tests performed. It is advised to make sure that these information are consistent with your data::

    CHR/POS columns present: SNPs searched based on genomic positions.
    Extracting SNPs for each chromosome...
    SNPs extracted for chr1.
    SNPs extracted for chr2.
    SNPs extracted for chr3.
    SNPs extracted for chr4.
    SNPs extracted for chr5.
    SNPs extracted for chr6.
    SNPs extracted for chr7.
    SNPs extracted for chr8.
    SNPs extracted for chr9.
    SNPs extracted for chr10.
    SNPs extracted for chr11.
    SNPs extracted for chr12.
    SNPs extracted for chr13.
    SNPs extracted for chr14.
    SNPs extracted for chr15.
    SNPs extracted for chr16.
    SNPs extracted for chr17.
    SNPs extracted for chr18.
    SNPs extracted for chr19.
    SNPs extracted for chr20.
    SNPs extracted for chr21.
    SNPs extracted for chr22.
    Merging SNPs extracted from each chromosome...
    Created bed/bim/fam fileset with extracted SNPs: tmp_GENAL/e415aab3_allchr
    39131 individuals are present in the genetic data and have a valid phenotype trait.
    Running single-SNP logistic regression tests on tmp_GENAL/e415aab3_allchr data with adjustment for: ['age'].
    The BETA, SE, P columns of the .data attribute have been updated.
    
The SBP_adjusted.data attribute has been updated in the BETA, SE, and P columns with the results of the association tests. 

Lifting
=======

It is sometimes necessary to lift the SNP data to a different build. For instance, if the genetic data of our target population is in build 38 (hg38), but the GWAS summary statistics are in build 37 (hg19).
This can easily be done in Genal using the :meth:`~genal.Geno.lift` method::

    SBP_clumped.lift(start = "hg19", end = "hg38", replace = False)
    
This output a table with the lifted SBP instruments (stored in the SBP_clumped object) from build 37 (hg19) to build 38 (hg38). We specified 'replace = False' to not modify the SBP_clumped.data attribute, but we may want to modify it (before running a PRS in a population stored in build 38 for instance). Genal will download the appropriate chain files required for the lift, and it will be done in pure python by default. However, if you plan to lift large datasets of SNPs (the whole summary statistics for instance), it may be useful to install the LiftOver executable that will run faster than the pure python version. It can be downloaded here: https://genome-store.ucsc.edu/ You will need to create an account, scroll down to "LiftOver program", add it to your cart, and declare that you are a non-profit user. 
You can need specify the path of the LiftOver executable to the 'liftover_path' argument::

    SBP_Geno.lift(start = "hg19", end = "hg38", replace = False, liftover_path = "path/to/liftover/exec")












































