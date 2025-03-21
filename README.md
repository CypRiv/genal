[![Python 3.8](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue)](https://www.python.org/downloads/release/python-3100/)

<img src="/genal_logo.png" data-canonical-src="/genal_logo.png" height="80" />  

<center><h1> genal: A Python Toolkit for Genetic Risk Scoring and Mendelian Randomization </h1></center>


# Table of contents
1. [Introduction](#introduction)
2. [Citation](#citation)
3. [Requirements for the genal module](#paragraph1)
4. [Installation and how to use genal](#paragraph2)
    1. [Installation](#paragraph2.1)
    2. [Documentation](#paragraph2.2)
5. [Tutorial and presentation of the main tools](#paragraph3)
    1. [Data loading](#paragraph3.1)
    2. [Data preprocessing](#paragraph3.2)
    3. [Clumping](#paragraph3.3)
    4. [Polygenic Risk Scoring](#paragraph3.4)
    5. [Mendelian Randomization](#paragraph3.5)
    6. [SNP-association testing](#paragraph3.6)
    7. [Lifting](#paragraph3.7)
    8. [GWAS Catalog](#paragraph3.8)


## Introduction <a name="introduction"></a>
Genal is a python module designed to make it easy and intuitive to run genetic risk scores and Mendelian Randomization analyses. The functionalities provided by genal include:

- Data preprocessing and cleaning of variant data (usually from GWAS summary statistics)
- Selection of independent genetic instruments through clumping
- Polygenic risk score calculations (currently C+T only, soon LDPRED2 and PRScs)
- More than 10 Mendelian Randomization methods, including heterogeneity and pleiotropy tests, with parallel processing support:
  - Inverse Variance Weighted (IVW) methods
  - MR-Egger methods
  - Weighted Median methods
  - Mode methods
  - MR-PRESSO in parallel
  - More to come...
- SNP-trait association testing
- Lifting of genetic data to another genomic build
- Variant-phenotype querying with the GWAS Catalog

### Key Features

- **Efficient Parallel Processing**: Parallel computation for bootstrapping-based MR methods and MR-PRESSO significantly reduces computation time compared to the original R packages (up to 85% faster for MR-PRESSO)
- **Flexible Data Handling**: Automatic formatting of variant data and summary statistics
- **Comprehensive MR Pipeline**: From data preprocessing to sensitivity analyses and plotting in a single package
- **Reference Panel Support**: Automatically download and use the latest 1000 Genomes reference panels in builds 37 and 38 with the option to use custom reference panels
- **Customizable**: Ability to choose all the parameters, but defaults are set to the most common values
- **Proxy SNP Support**: Includes functionality for finding and using proxy SNPs when instruments are missing (for polygenic risk scores, Mendelian Randomization)

The objective of genal is to bring the functionalities of well-established R packages such as TwoSampleMR, MR-Presso, MendelianRandomization, and gwasvcf, in a more user-friendly Python environment. This approach ensures that users have access to tried and tested techniques with the versatility of Python's data science tools. 

This package is still under development, feel free to report any issues or suggest improvements!

<img src="/Genal_flowchart.png" data-canonical-src="/Genal_flowchart.png" style="max-width:100%;" />

Genal flowchart. Created in https://www.BioRender.com
## Citation <a name="citation"></a> 
This project was developed by Cyprien A. Rivier.
If you're using genal, please cite the following paper:  
**Genal: A Python Toolkit for Genetic Risk Scoring and Mendelian Randomization.**  
Cyprien A. Rivier, Santiago Clocchiatti-Tuozzo, Shufan Huo, Victor Torres-Lopez, Daniela Renedo, Kevin N. Sheth, Guido J. Falcone, Julian N. Acosta.  
Bioinformatics Advances 2024.  
doi: https://doi.org/10.1093/bioadv/vbae207

If you're using methods derived from R packages, such as MR-PRESSO, please also cite the original papers.

## Requirements for the genal module <a name="paragraph1"></a> 
***Python 3.8 or later***. https://www.python.org/ <br> 


## Installation and How to use the genal module <a name="paragraph2"></a>

### Installation <a name="paragraph2.1"></a>

> **Note:**
> 
> **Optional**: It is recommended to create a new environment to avoid dependencies conflicts. Here, we create a new conda environment called 'genal_env'.
> ```
> conda create --name genal_env python=3.8
> conda activate genal_env
> ```


Download and install the package with pip:    
```
pip install genal-python
```
And import it in a python environment with:
```python
import genal
```

The main genal functionalities require a working installation of PLINK v2.0. 
If you have already installed plink v2.0, you can set the path to its executable with:

```
genal.set_plink(path="/path/to/plink/executable/file")
```

If plink is not installed, genal can install the correct version for your system with the following line:

```
genal.install_plink()
```

### Documentation <a name="paragraph2.2"></a>

For detailed information on how to use the functionalities of Genal, please refer to the documentation: https://genal.rtfd.io

The documentation covers:
- Installation 
- This tutorial
- The list of the main functions with complete description of their arguments
- An exhaustive API reference


## Tutorial <a name="paragraph3"></a>
For this tutorial, we will obtain genetic instruments for systolic blood pressure (SBP), compute a Polygenic Risk Score (PRS), and run a Mendelian Randomization analysis to investigate the genetically-determined effect of SBP on the risk of stroke. We will utilize summary statistics from Genome-Wide Association Studies (GWAS) and individual-level data from the UK Biobank. The steps include:

- Data loading
- Data preprocessing
  - Perform checks and cleaning operations on SNP-level data for appropriate formatting
- Clump data to identify significant, independent SNPs as genetic instruments for SBP
- Build a genomic risk score for SBP in a test population
  - Include risk score calculations with proxies
- Perform Mendelian Randomization
  - Analyze SBP as an exposure and acute stroke as an outcome
  - Plot the results
  - Conduct sensitivity analyses using the weighted median, MR-Egger, and MR-PRESSO methods
- Calibrate SNP-trait weights with individual-level genetic data
  - Execute single-SNP association tests for calibrating SBP genetic instruments
- Utility tools demonstration
  - Data lifting to another genomic build
    - In pure Python
    - Using LiftOver
  - Querying the GWAS Catalog

### Data loading <a name="paragraph3.1"></a>

We start this tutorial with publicly available summary statistics from a large GWAS study of systolic blood pressure. [Link to study](https://www.nature.com/articles/s41588-018-0205-x). [Download link](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz). After downloading and unzipping the summary statistics, we load them into a pandas DataFrame:

```python
import pandas as pd
sbp_gwas = pd.read_csv("Evangelou_30224653_SBP.txt", sep=" ")
sbp_gwas.head(5)
```
We can now load this data into a `genal.Geno` instance. The `genal.Geno` class is the central piece of the package. It is designed to store Single Nucleotide Polymorphisms (SNP) data and make it easy to preprocess and clean it.

The `genal.Geno` takes as input a pandas dataframe where each row corresponds to a SNP, with columns describing the position and possibly the effect of the SNP for the given trait (SBP in our case). To indicate the names of the columns, the following arguments can be passed:
- **CHR**: Column name for chromosome. Defaults to `'CHR'`.
- **POS**: Column name for genomic position. Defaults to `'POS'`.
- **SNP**: Column name for SNP identifier (rsid). Defaults to `'SNP'`.
- **EA**: Column name for effect allele. Defaults to `'EA'`.
- **NEA**: Column name for non-effect allele. Defaults to `'NEA'`.
- **BETA**: Column name for effect estimate. Defaults to `'BETA'`.
- **SE**: Column name for effect standard error. Defaults to `'SE'`.
- **P**: Column name for effect p-value. Defaults to `'P'`.
- **EAF**: Column name for effect allele frequency. Defaults to `'EAF'`.

> **Note:**
> 
> You do not need all columns to move forward, as not all columns are required by every function. Additionally, some columns can be imputed as we will see in the next paragraph.

After inspecting the dataframe, we first need to extract the chromosome and position information from the `MarkerName` column into two new columns `CHR` and `POS`:

```python
sbp_gwas[["CHR", "POS", "Filler"]] = sbp_gwas["MarkerName"].str.split(":", expand=True)
sbp_gwas.head(2)
```
| MarkerName        | Allele1 | Allele2 | Freq1 | Effect | StdErr | P        | TotalSampleSize | N_effective | CHR | POS       | Filler |
|-------------------|---------|---------|-------|--------|--------|----------|-----------------|-------------|-----|-----------|--------|
| 10:100000625:SNP  | a       | g       | 0.5660| 0.0523 | 0.0303 | 0.083940 | 738170          | 736847      | 10  | 100000625 | SNP    |
| 10:100000645:SNP  | a       | c       | 0.7936| 0.0200 | 0.0372 | 0.591100 | 738168          | 735018      | 10  | 100000645 | SNP    |

And it can now be loaded into a `genal.Geno` instance:

```python
import genal
SBP_Geno = genal.Geno(sbp_gwas, CHR="CHR", POS="POS", EA="Allele1", NEA="Allele2", BETA="Effect", SE="StdErr", P="P", EAF="Freq1", keep_columns=False)
```

The last argument (`keep_columns = False`) indicates that we do not wish to keep the other (non-main) columns in the dataframe. Defaults to `True`.

> **Note:**
> 
> Make sure to read the readme file usually provided with the summary statistics to identify the correct columns. It is particularly important to correctly identify the allele that represents the effect allele.

### Data preprocessing <a name="paragraph3.2"></a>

Now that we have loaded the data into a `genal.Geno` instance, we can begin cleaning and formatting it. Methods such as Polygenic Risk Scoring or Mendelian Randomization require the SNP data to be in a specific format. Also, raw summary statistics can sometimes contain missing or invalid values that need to be handled. Additionally, some columns may be missing from the data (such as the SNP rsid column, or the non-effect allele column) and these columns can be created based on existing ones and a reference panel.

Genal can run all the basic cleaning and preprocessing steps in one command:

```python
SBP_Geno.preprocess_data(preprocessing = 'Fill_delete')
```

The `preprocessing` argument specifies the global level of preprocessing applied to the data:
- `preprocessing = 'None'`: The data won't be modified.
- `preprocessing = 'Fill'`: Missing columns will be added based on reference data and invalid values set to NaN, but no rows will be deleted.
- `preprocessing = 'Fill_delete'`: Missing columns will be added, and all rows containing missing, duplicated, or invalid values will be deleted. This option is recommended before running genetic methods.
Defaults to `'Fill'`.

By default, and depending on the global preprocessing level (`'None'`, `'Fill'`, `'Fill_delete'`) chosen, the `preprocess_data` method of `genal.Geno` will run the following checks:
- Ensure the `CHR` (chromosome) and `POS` (genomic position) columns are integers.
- Ensure the `EA` (effect allele) and `NEA` (non-effect allele) columns are uppercase characters containing A, T, C, G letters. Multiallelic values are set to NaN.
- Validate the `P` (p-value) column for proper values.
- Check for no duplicated SNPs based on rsid.
- Determine if the `BETA` (effect) column contains beta estimates or odds ratios, and log-transform odds ratios if necessary.
- Create `SNP` column (containing rsids) using a reference panel if CHR and POS columns are present.
- Create `CHR` and/or `POS` column using a reference panel if `SNP` column is present.
- Create `NEA` (non-effect allele) column using a reference panel if `EA` (effect allele) column is present.
- Create the `SE` (standard-error) column if the `BETA` and `P` (p-value) columns are present.
- Create the `P` column if the `BETA` and `SE` columns are present.

If you do not wish to run certain steps, or wish to run only certain steps, you can use additional arguments. For more information, please refer to the `genal.Geno.preprocess_data` method in the API documentation.

In our case, the `SNP` column (for SNP identifier - rsid) was missing from our dataframe and has been added based on a 1000 genome reference panel:

    Using the EUR reference panel.
    The SNP column (rsID) has been created. 197511(2.787%) SNPs were not found in the reference data and their ID set to CHR:POS:EA.
    The BETA column looks like Beta estimates. Use effect_column='OR' if it is a column of Odds Ratios.

You can always check the data of a `genal.Geno` instance by accessing the 'data' attribute:

```python
SBP_Geno.data
```
|         | EA | NEA |   EAF |   BETA |     SE |        P | CHR |     POS |       SNP |
|---------|----|-----|-------|--------|--------|----------|-----|---------|-----------|
| 0       |  A |   G | 0.5660|  0.0523| 0.0303 | 0.083940 |  10 | 100000625 | rs7899632 |
| 1       |  A |   C | 0.7936|  0.0200| 0.0372 | 0.591100 |  10 | 100000645 | rs61875309 |


The `SNP` column with the rsids has been added based on the reference data.
You do not need to obtain the 1000 genome reference panel yourself, genal will download it the first time you use it. 

> **Note:**
> 
> By default, the reference panel used is the european (eur) one in build 37. You can specify another valid reference panel (afr, eas, sas, amr) with the reference_panel argument. For instance "AFR_38" for the african reference panel in build 38:

```python
SBP_Geno.preprocess_data(preprocessing = 'Fill_delete', reference_panel = "AFR_37")
```

You can also use a custom reference panel by specifying to the reference_panel argument a path to bed/bim/fam (plink v1.9 format) or pgen/pvar/psam files (plink v2.0 format), without the extension.

### Clumping <a name="paragraph3.3"></a>

Clumping, or C+T: Clumping + Thresholding, is the step at which we select the SNPs that will be used as our genetic instruments in future Polygenic Risk Scores and Mendelian Randomization analyses. The process involves identifying the SNPs that are strongly associated with our trait of interest (systolic blood pressure in this tutorial) and are independent from each other. This second step ensures that selected SNPs are not highly correlated, (i.e., they are not in high linkage disequilibrium). For this step, we again need to use a reference panel.

The SNP-data loaded in a `genal.Geno` instance can be clumped using the `genal.Geno.clump` method. It will return another `genal.Geno` instance containing only the clumped data:

```python
SBP_clumped = SBP_Geno.clump(p1 = 5e-8, r2 = 0.01, kb = 10000, reference_panel = "EUR_37")
```

You can specify the thresholds you want to use for the clumping with the following arguments:
- `p1`: P-value threshold during clumping. SNPs with a P-value higher than this value are excluded. Defaults to `5e-8`.
- `r2`: Linkage disequilibrium threshold for the independence check. Takes values between 0 and 1. Defaults to `0.01`.
- `kb`: Genomic window used for the independence check (the unit is thousands of base-pair positions). Defaults to `10000`.
- `reference_panel`: The reference population used to derive linkage disequilibrium values and select independent SNPs. Defaults to `EUR_37`.

### Polygenic Risk Scoring <a name="paragraph3.4"></a>

Computing a Polygenic Risk Score (PRS) can be done in one line with the `genal.Geno.prs` method:

```python
SBP_clumped.prs(name = "SBP_prs", path = "path/to/genetic/files")
```

The genetic files of the target population can be either contained in one triple of bed/bim/fam or pgen/pvar/psam files with information for all SNPs, or divided by chromosome (one bed/bim/fam or pgen/pvar/psam triple for chr 1, another for chr 2, etc...). In the latter case, provide the path by replacing the chromosome number by `$` and genal will extract the necessary SNPs from each chromosome and merge them before running the PRS. For instance, if the genetic files are named `Pop_chr1.bed`, `Pop_chr1.bim`, `Pop_chr1.fam`, `Pop_chr2.bed`, ..., you can use:

```python
SBP_clumped.prs(name = "SBP_prs", path = "Pop_chr$")
```

The `name` argument specifies the name of the .csv file that will be saved with the individual risk scores. 
The output of the `genal.Geno.prs` method will include how many SNPs were used to compute the risk score. It can happen that some of the SNPs are multiallelic in the genetic data (even if they are not multiallelic in our SNP data) and need to be excluded. It can also happen that some of the SNPs are missing from the genetic files of the target population (for instance if the data has not been imputed):

    CHR/POS columns present: SNPs searched based on genomic positions.
    Extracting SNPs for each chromosome...
    SNPs extracted for chr1.
    SNPs extracted for chr2.
    ...
    SNPs extracted for chr22.
    Merging SNPs extracted from each chromosome...
    Created bed/bim/fam fileset with extracted SNPs: tmp_GENAL/4f4ce6a7_allchr
    Extraction completed. 786(50.874%) SNPs were not extracted from the genetic data.
    Computing a weighted PRS using tmp_GENAL/4f4ce6a7_allchr data.
    The PRS computation was successful and used 759/1545 (49.126%) SNPs.
    PRS data saved to SBP_prs.csv

Here, we see that about half of the SNPs were not extracted from the data. In such cases, we may want to try and salvage some of these SNPs by looking for proxies (SNPs in high linkage disequilibrium, i.e. highly correlated SNPs). This can be done by specifying the `proxy = True`. argument:

```python
SBP_clumped.prs(name = "SBP_prs_proxy" ,path = "Pop_chr$", proxy = True, reference_panel = "eur", r2=0.8, kb=5000, window_snps=5000)
```

and the output is:

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
    ...
    SNPs extracted for chr22.
    Merging SNPs extracted from each chromosome...
    Created bed/bim/fam fileset with extracted SNPs: tmp_GENAL/4f4ce6a7_allchr
    Extraction completed. 208(13.524%) SNPs were not extracted from the genetic data.
    Computing a weighted PRS using tmp_GENAL/4f4ce6a7_allchr data.
    The PRS computation was successful and used 1330/1538 (86.476%) SNPs.
    PRS data saved to SBP_prs.csv

In our case, we have been able to find proxies for 578 of the 786 SNPs that were missing in the population genetic data (7 potential proxies have been removed because they were identical to SNPs already present in our data).

You can customize how the proxies are chosen with the following arguments:
- `reference_panel`: The reference population used to derive linkage disequilibrium values and find proxies. Defaults to `eur`.
- `kb`: Width of the genomic window to look for proxies (in thousands of base-pair positions). Defaults to `5000`.
- `r2`: Minimum linkage disequilibrium value with the original SNP for a proxy to be included. Defaults to `0.8`.
- `window_snps`: Width of the window to look for proxies (in number of SNPs). Defaults to `5000`.

> **Note:**
> 
> You can call the `genal.Geno.prs` method on any `Geno` instance (containing at least the EA, BETA, and either SNP or CHR/POS columns). The data does not need to be clumped, and there is no limit to the number of SNPs used to compute the scores.


### Mendelian Randomization <a name="paragraph3.5"></a>

To run MR, we need to load both our exposure and outcome SNP-level data in `genal.Geno` instances. In our case, the genetic instruments of the MR are the SNPs associated with blood pressure at genome-wide significant levels resulting from the clumping of the blood pressure GWAS. They are stored in our `SBP_clumped` `genal.Geno` instance which also include their association with the exposure trait (instrument-SBP estimates in the `BETA` column).

To get their association with the outcome trait (instrument-stroke estimates), we are going to use SNP-level data from a large GWAS of stroke performed by the GIGASTROKE consortium: [Link to study](https://www.nature.com/articles/s41586-022-05165-3). [Link to download](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104539/GCST90104539_buildGRCh37.tsv.gz):

```python
stroke_gwas = pd.read_csv("GCST90104539_buildGRCh37.tsv",sep="\t")
```

We load it in a `genal.Geno` instance:

```python
Stroke_Geno = genal.Geno(stroke_gwas, CHR = "chromosome", POS = "base_pair_location", EA = "effect_allele", NEA = "other_allele", BETA = "beta", SE = "standard_error", P = "p_value", EAF = "effect_allele_frequency", keep_columns = False)
```

We preprocess it as well to put it in the correct format and make sure there is no invalid values:

```python
Stroke_Geno.preprocess_data(preprocessing = 'Fill_delete')
```

Now, we need to extract our instruments (SNPs of the `SBP_clumped` data) from the outcome data to obtain their association with the outcome trait (stroke). It can be done by calling the `genal.Geno.query_outcome` method:

```python
SBP_clumped.query_outcome(Stroke_Geno, proxy = False)
```

Genal will print how many SNPs were successfully found and extracted from the outcome data:

    Outcome data successfully loaded from 'b352e412' geno instance.
    Identifying the exposure SNPs present in the outcome data...
    1541 SNPs out of 1545 are present in the outcome data.
    (Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute.
    

> **Note:**
>
> Here as well you have the option to use proxies for the instruments that are not present in the outcome data:
>
> ```python
> SBP_clumped.query_outcome(Stroke_geno, proxy = True, reference_panel = "EUR_37", kb = 5000, r2 = 0.8, window_snps = 5000)
> ```
> 
> And genal will print the number of missing instruments that have been proxied:
> 
>     Outcome data successfully loaded from 'b352e412' geno instance.
>     Identifying the exposure SNPs present in the outcome data...
>     1541 SNPs out of 1545 are present in the outcome data.
>     Searching proxies for 4 SNPs...
>     Using the EUR reference panel.
>     Found proxies for 4 SNPs.
>     (Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute.

After extracting the instruments from the outcome data, the `SBP_clumped` `genal.Geno` instance contains an `MR_data` attribute containing the instruments-exposure and instruments-outcome associations necessary to run MR. Running MR is now as simple as calling the `genal.Geno.MR` method of the SBP_clumped `genal.Geno` instance:

```python
SBP_clumped.MR(action = 2, exposure_name = "SBP", outcome_name = "Stroke_eur")
```

The `genal.Geno.MR` method prints a message specifying which SNPs have been excluded from the analysis (it depends on the action argument, as we will see):

    Action = 2: 42 SNPs excluded for being palindromic with intermediate allele frequencies: rs11817866, rs3802517, rs2788293, rs2274224, rs7310615, rs7953257, rs2024385, rs61912333, rs11632436, rs1012089, rs3851018, rs9899540, rs4617956, rs773432, rs11585169, rs7796, rs2487904, rs12321, rs73029563, rs4673238, rs3845811, rs2160236, rs10165271, rs9848170, rs2724535, rs6842486, rs4834792, rs990619, rs155364, rs480882, rs6875372, rs258951, rs1870735, rs1800795, rs12700814, rs1821002, rs3021500, rs28601761, rs7463212, rs907183, rs534523, rs520015 

It returns a dataframe containing the results for different MR methods:
| exposure | outcome    | method                                    | nSNP | b        | se       | pval          |
|----------|------------|-------------------------------------------|------|----------|----------|---------------|
| SBP      | Stroke_eur | Inverse-Variance Weighted                 | 1499 | 0.023049 | 0.001061 | 1.382645e-104 |
| SBP      | Stroke_eur | Inverse Variance Weighted (Fixed Effects) | 1499 | 0.023049 | 0.000754 | 4.390655e-205 |
| SBP      | Stroke_eur | Weighted Median                           | 1499 | 0.022365 | 0.001337 | 8.863203e-63  |
| SBP      | Stroke_eur | Simple mode                               | 1499 | 0.027125 | 0.007698 | 4.382993e-04  |
| SBP      | Stroke_eur | MR Egger                                  | 1499 | 0.027543 | 0.002849 | 1.723156e-21  |
| SBP      | Stroke_eur | Egger Intercept                           | 1499 | -0.001381| 0.000813 | 8.935529e-02  |

You can specify several arguments. We refer to the API for a full list, but the most important one is the `action` argument. It determines how palindromic SNPs are treated during the exposure-outcome harmonization step. Palindromic SNPs are SNPs where the nucleotide change reads the same forward and backward on complementary strands of DNA (for instance `EA = 'A'` and `NEA = 'T'`).

- `action = 1`: Palindromic SNPs are not treated (assumes all alleles are on the forward strand)
- `action = 2`: Uses effect allele frequencies to attempt to flip them (conservative, default)
- `action = 3`: Removes all palindromic SNPs (very conservative)

When choosing the option 2 or 3 (recommended), genal will print the list of palindromic SNPs that have been removed from the analysis.

By default, only some MR methods (inverse-variance weighted, weighted median, Simple mode, MR-Egger) are going to be run. But if you wish to run a different set of MR methods, you can pass a list of strings to the `methods` argument. The possible strings are:
- `IVW` for the classical Inverse-Variance Weighted method with random effects
- `IVW-RE` for the Inverse Variance Weighted method with Random Effects where the standard error is not corrected for under dispersion
- `IVW-FE` for the Inverse Variance Weighted with fixed effects
- `UWR` for the Unweighted Regression method
- `WM` for the Weighted Median method
- `WM-pen` for the penalised Weighted Median method
- `Simple-median` for the Simple Median method
- `Sign` for the Sign concordance test
- `Egger` for MR-Egger and the MR-Egger intercept
- `Egger-boot` for the bootstrapped version of MR-Egger and its intercept
- `Simple-mode` for the Simple mode method
- `Weighted-mode` for the Weighted mode method
- `all` to run all the above methods

For more fine-tuning, such as settings for the number of boostrapping iterations, please refer to the API: [https://genal.readthedocs.io/en/latest/modules.html#id4] (MR method).

If you want to visualize the obtained MR results, you can use the `genal.Geno.MR_plot` method that will plot each SNP in an `effect_on_exposure x effect_on_outcome` plane as well as lines corresponding to different MR methods:

```python
SBP_clumped.MR_plot(filename="MR_plot_SBP_AS")
```

![MR plot](docs/build/_images/MR_plot_SBP_AS.png)
You can select which MR methods to plot with the `methods` argument. Note that for an MR method to be plotted, they must be included in the latest `genal.Geno.MR` call of this `genal.Geno` instance.

To include the heterogeneity values (Cochran's Q) in the results, you can use the heterogeneity argument in the `genal.Geno.MR` call. Here, the heterogeneity for the inverse-variance weighted method:

```python
SBP_clumped.MR(action = 2, methods = ["Egger","IVW"], exposure_name = "SBP", outcome_name = "Stroke_eur", heterogeneity = True)
```

And that will give:
| exposure | outcome    | method                   | nSNP | b         | se        | pval          | Q           | Q_df | Q_pval       |
|----------|------------|--------------------------|------|-----------|-----------|---------------|-------------|------|--------------|
| SBP      | Stroke_eur | MR Egger                 | 1499 | 0.027543  | 0.002849  | 1.723156e-21  | 2959.965136 | 1497 | 1.253763e-98 |
| SBP      | Stroke_eur | Egger Intercept          | 1499 | -0.001381 | 0.000813  | 8.935529e-02  | 2959.965136 | 1497 | 1.253763e-98 |
| SBP      | Stroke_eur | Inverse-Variance Weighted| 1499 | 0.023049  | 0.001061  | 1.382645e-104 | 2965.678836 | 1498 | 4.280737e-99 |
    
To display the coefficients as odds ratios with confidence intervals for a binary outcome trait, you can use the `odds = True` argument:

```python
SBP_clumped.MR(action = 2, methods = ["Egger","IVW"], exposure_name = "SBP", outcome_name = "Stroke_eur", heterogeneity = True, odds = True)
```

As expected, many MR methods indicate that SBP is strongly associated with stroke, but there could be concerns for horizontal pleiotropy (instruments influencing the outcome through a different pathway than the one used as exposure) given the almost significant MR-Egger intercept p-value.
To investigate horizontal pleiotropy in more details, a very useful method is Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO). MR-PRESSO is a method designed to detect and correct for horizontal pleiotropy. It will identify which instruments are likely to be pleiotropic on their effect on the outcome, and it will rerun an inverse-variance weighted MR after excluding them. It can be run using the `genal.Geno.MRpresso` method:

```python
mod_table, GlobalTest, OutlierTest, BiasTest = SBP_clumped.MRpresso(action = 2, n_iterations = 30000)
```

As with the `genal.Geno.MR` method, the `action` argument determines how the pleiotropic SNPs will be treated. The output is a list containing:
- A table containing the original and outlier-corrected inverse variance-weighted results.
- The global test p-value indicating the presence of horizontal pleiotropy.
- A dataframe of p-values, one for each instrument, representing the likelihood that this instrument is pleiotropic (only relevant if the global test is significant).
- A dictionary containing the outputs of the distortion test. This test assesses whether the removal of the pleiotropic instruments has significantly altered the original MR estimate.
    - An array containing the indices of the pleiotropic SNPs.
    - The coefficient of the distortion test.
    - The p-value of the distortion test.

### SNP-association testing <a name="paragraph3.6"></a>

We may want to calibrate instrument-trait estimates in a specific population for which we have individual-level data (genetic files as well as phenotypic data). For instance, if the GWAS of SBP was done in a european population, we may want to adjust the estimates based on data coming from a population of a different ancestry. This can be done in 2 steps:
- Loading the phenotypic data in a dataframe and calling the `genal.Geno.set_phenotype` method
- Calling the `genal.Geno.association_test` method to run the association tests and update the estimates
Let's start by loading phenotypic data::

```python
df_pheno = pd.read_csv("path/to/trait/data")
```

> **Note:**
> 
>    One important point is to make sure that both the Family IDs (FID) and Individual IDs (IID) of the participants are identical in the phenotypic data and in the genetic data. 

Then, it is advised to make a copy of the `genal.Geno` instance containing our instruments as we are going to update their coefficients and to avoid any confusion:

```python
SBP_adjusted = SBP_clumped.copy()
```

We can then call the `genal.Geno.set_phenotype` method, specifying which column contains our trait of interest (for the association testing) and which column contains the individual IDs:

```python
SBP_adjusted.set_phenotype(df_pheno, PHENO = "htn", IID = "IID", FID = "FID")
```

At this point, genal will identify if the phenotype is binary or quantitative in order to choose the appropriate regression model. If the phenotype is binary, it will assume that the most frequent value is coding for control (and the other value for case), this can be changed with `alternate_control = True`:

    Detected a binary phenotype in the 'PHENO' column. Specify 'PHENO_type="quant"' if this is incorrect.
    Identified 0 as the control code in 'PHENO'. Set 'alternate_control=True' to inverse this interpretation.
    The phenotype data is stored in the .phenotype attribute.
    
We can then run the association tests, specifying the path to the genetic files in plink format, and any columns we may want to include as covariates in the regression tests (making sure that the covariates are all numerical):

```python
SBP_adjusted.association_test(covar=["age"], path = "path/to/genetic/files")
```

Genal will print information regarding the number of individuals used in the tests and the kind of tests performed. It is advised to make sure that these information are consistent with your data:

    CHR/POS columns present: SNPs searched based on genomic positions.
    Extracting SNPs for each chromosome...
    SNPs extracted for chr1.
    ...
    SNPs extracted for chr22.
    Merging SNPs extracted from each chromosome...
    Created bed/bim/fam fileset with extracted SNPs: tmp_GENAL/e415aab3_allchr
    39131 individuals are present in the genetic data and have a valid phenotype trait.
    Running single-SNP logistic regression tests on tmp_GENAL/e415aab3_allchr data with adjustment for: ['age'].
    The BETA, SE, P columns of the .data attribute have been updated.
    
The `BETA`, `SE`, and `P` columns of the `SBP_adjusted.data` attribute have been updated with the results of the association tests. 

### Lifting <a name="paragraph3.7"></a>

It is sometimes necessary to lift the SNP data to a different build. For instance, if the genetic data of our target population is in build 38 (hg38), but the GWAS summary statistics are in build 37 (hg19).
This can easily be done in genal using the `genal.Geno.lift` method:

```python
SBP_clumped.lift(start = "hg19", end = "hg38", replace = False)
```

This outputs a table with the lifted SBP instruments (stored in the `SBP_clumped` instance) from build 37 (hg19) to build 38 (hg38). We specified `replace = False` to not modify the `SBP_clumped.data` attribute, but we may want to modify it (before running a PRS in a population stored in build 38 for instance). 
Genal will download the appropriate chain files required for the lift, and it will be done in  python by default. However, if you plan to lift large datasets of SNPs (the whole summary statistics for instance), it may be useful to install the LiftOver executable that will run faster than the python version. It can be downloaded here: [https://genome-store.ucsc.edu/](https://genome-store.ucsc.edu/) You will need to create an account, scroll down to "LiftOver program", add it to your cart, and declare that you are a non-profit user.

You can specify the path of the LiftOver executable to the `liftover_path` argument:

```python
SBP_Geno.lift(start = "hg19", end = "hg38", replace = False, liftover_path = "path/to/liftover/exec")
```

### GWAS Catalog <a name="paragraph3.8"></a>

It is sometimes interesting to determine the traits associated with our SNPs. In Mendelian Randomization, for instance, we may want to exclude instruments that are associated with traits likely causing horizontal pleiotropy. For this purpose, we can use the `genal.Geno.query_gwas_catalog` method. This method will query the GWAS Catalog API to determine the list of traits associated with each of our SNPs and store the results in a list in the `ASSOC` column of the `.data` attribute:

```python
SBP_clumped.query_gwas_catalog(p_threshold=5e-8)
```
Which will output:

    Querying the GWAS Catalog and creating the ASSOC column. 
    Only associations with a p-value <= 5e-08 are reported. Use the p_threshold argument to change the threshold.
    To report the p-value for each association, use return_p=True.
    To report the study ID for each association, use return_study=True.
    The .data attribute will be modified. Use replace=False to leave it as is.
    100%|██████████| 1545/1545 [00:34<00:00, 44.86it/s]
    The ASSOC column has been successfully created.
    701 (45.37%) SNPs failed to query (not found in GWAS Catalog) and 7 (0.5%) SNPs timed out after 34.33 seconds. You can increase the timeout value with the timeout argument.
| EA  | NEA | EAF   | BETA   | SE     | CHR | POS        | SNP        | ASSOC                                                                 |
|-----|-----|-------|--------|--------|-----|------------|------------|------------------------------------------------------------------------|
| A   | G   | 0.1784| 0.2330 | 0.0402 | 10  | 102075479  | rs603424   | [eicosanoids measurement, decadienedioic acid (...]                     |
| A   | G   | 0.0706| -0.3873| 0.0626 | 10  | 102403682  | rs2996303  | FAILED_QUERY                                                           |
| T   | G   | 0.8872| 0.6846 | 0.0480 | 10  | 102553647  | rs1006545  | [diastolic blood pressure, systolic blood pressure...]                  |
| T   | G   | 0.6652| -0.2098| 0.0340 | 10  | 102558506  | rs12570050 | FAILED_QUERY                                                           |
| T   | C   | 0.3057| -0.2448| 0.0334 | 10  | 102603924  | rs4919502  | FAILED_QUERY                                                           |
| ... | ... | ...   | ...    | ...    | ... | ...        | ...        | ...                                                                    |                                          |
| T   | C   | 0.3514| 0.2203 | 0.0314 | 9   | 9350706    | rs1332813  | [diastolic blood pressure, systolic blood pressure...]                  |
| T   | C   | 0.6880| -0.1897| 0.0332 | 9   | 94201341   | rs10820855 | FAILED_QUERY                                                           |
| A   | T   | 0.3669| -0.1862| 0.0313 | 9   | 95201540   | rs7045409  | [protein measurement, pulse pressure measurement...]                   |

If you are also interested in the p-values of each SNP-trait association, or the ID of the study from which the association was reported, you can use the `return_p = True` and `return_study = True` arguments. Then, the `ASSOC` column will contain a list of tuples, where each tuple contains the trait name, the p-value, and the study ID:

```python
SBP_clumped.query_gwas_catalog(p_threshold=5e-8, return_p=True, return_study=True)
```

| EA  | NEA | EAF   | BETA   | SE     | CHR | POS        | SNP        | ASSOC                                                                 |
|-----|-----|-------|--------|--------|-----|------------|------------|------------------------------------------------------------------------|
| A   | G   | 0.1784| 0.2330 | 0.0402 | 10  | 102075479  | rs603424   | TIMEOUT                                                                |
| A   | G   | 0.0706| -0.3873| 0.0626 | 10  | 102403682  | rs2996303  | FAILED_QUERY                                                           |
| T   | G   | 0.8872| 0.6846 | 0.0480 | 10  | 102553647  | rs1006545  | [(heart rate response to exercise, 6e-12, GCST...                      |
| T   | G   | 0.6652| -0.2098| 0.0340 | 10  | 102558506  | rs12570050 | FAILED_QUERY                                                           |
| T   | C   | 0.3057| -0.2448| 0.0334 | 10  | 102603924  | rs4919502  | FAILED_QUERY                                                           |
| ... | ... | ...   | ...    | ...    | ... | ...        | ...        | ...                                                                    |                                                         |
| T   | C   | 0.3514| 0.2203 | 0.0314 | 9   | 9350706    | rs1332813  | [(diastolic blood pressure, 1e-12, GCST9031029...                      |
| T   | C   | 0.6880| -0.1897| 0.0332 | 9   | 94201341   | rs10820855 | FAILED_QUERY                                                           |
| A   | T   | 0.3669| -0.1862| 0.0313 | 9   | 95201540   | rs7045409  | [(systolic blood pressure, 9e-13, GCST006624),...                      |


> **Note:**
> 
> As you can see, many SNPs failed to be queried. This is normal as the GWAS Catalog is not exhaustive.