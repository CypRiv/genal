============
Installation
============

.. note::
    **Optional**: It is recommended to create a new environment to avoid dependencies conflicts. Here, we create a new conda environment called 'genal'.

    .. code-block:: bash

        conda create --name genal python=3.11
        conda activate genal

The genal package requires Python 3.11. Download and install it with pip: 

.. code-block:: bash

    pip install genal-python

And import it in a python environment with:

.. code-block:: python

    import genal

The main genal functionalities require a working installation of PLINK v1.9 (and not 2.0 as certain functionalities have not been updated yet) that can be downloaded here: https://www.cog-genomics.org/plink/ 
Once downloaded, the path to the plink 1.9 executable should be set with:

.. code-block:: python

    genal.set_plink(path="/path/to/plink/executable/file")

========
Tutorial
========

For the purpose of this tutorial, we are going to build a PRS of systolic blood pressure (SBP) and investigate the genetically-determined effect of SBP on the risk of stroke. We will use both summary statistics from Genome-Wide Association Studies (GWAS) and individual-level data from the UK Biobank as our test population. We are going to go through the following steps:

Table of contents
-----------------

a. `Data loading`_
b. `Data preprocessing`_
c. `Clumping`_
d. `Polygenic Risk Scoring`_
e. `Mendelian Randomization`_
f. `SNP-association testing`_
g. `Lifting`_
h. `GWAS Catalog`_


Data loading
============

We start this tutorial with publicly available summary statistics data from a large GWAS of systolic blood pressure (https://www.nature.com/articles/s41588-018-0205-x). After downloading and unzipping the summary statistics, we load the data into a pandas dataframe:

.. code-block:: python

    import pandas as pd
    sbp_gwas = pd.read_csv("Evangelou_30224653_SBP.txt", sep=" ")

The loaded dataframe:

.. code-block:: python

    >>> sbp_gwas.head(5)
             MarkerName Allele1 Allele2   Freq1  Effect  StdErr         P  TotalSampleSize  N_effective
    0  10:100000625:SNP       a       g  0.5660  0.0523  0.0303  0.083940           738170       736847
    1  10:100000645:SNP       a       c  0.7936  0.0200  0.0372  0.591100           738168       735018
    2  10:100003242:SNP       t       g  0.8831  0.1417  0.0469  0.002526           738168       733070
    3  10:100003304:SNP       a       g  0.9609  0.0245  0.0838  0.769800           737054       663809
    4  10:100003785:SNP       t       c  0.6406 -0.0680  0.0313  0.029870           738169       735681

We can now load this data into a :class:`~genal.Geno` object. The :class:`~genal.Geno` class is the central piece of the package. It is designed to store Single Nucleotide Polymorphisms (SNP) data and make it easy to preprocess and clean.

The :class:`~genal.Geno` takes as input a pandas dataframe where each row corresponds to a SNP, with columns describing the position and possibly the effect of the SNP for the given trait (SBP in our case). The following arguments can be passed to specify the column names:

* **CHR**: Column name for chromosome. Defaults to "CHR".
* **POS**: Column name for genomic position. Defaults to "POS".
* **SNP**: Column name for SNP identifier (rsid). Defaults to "SNP".
* **EA**: Column name for effect allele. Defaults to "EA".
* **NEA**: Column name for non-effect allele. Defaults to "NEA".
* **BETA**: Column name for effect estimate. Defaults to "BETA".
* **SE**: Column name for effect standard error. Defaults to "SE".
* **P**: Column name for effect p-value. Defaults to "P".
* **EAF**: Column name for effect allele frequency. Defaults to "EAF".

.. note::

   You do not need all columns to move forward, as not all columns are required by every function. Additionally, some columns can be imputed as we will see in the next paragraph.

In our case, and after inspecting the dataframe, we must first extract the chromosome and position information from the "MarkerName" column into two new columns "CHR" and "POS":

.. code-block:: python

    sbp_gwas[["CHR", "POS", "Filler"]] = sbp_gwas["MarkerName"].str.split(":", expand=True)

The resulting dataframe now has separate columns for the chromosome and genomic position:

.. code-block:: python

    >>> sbp_gwas.head(5)
             MarkerName Allele1 Allele2   Freq1  Effect  StdErr         P  TotalSampleSize  N_effective CHR        POS Filler
    0  10:100000625:SNP       a       g  0.5660  0.0523  0.0303  0.083940           738170       736847  10  100000625    SNP
    1  10:100000645:SNP       a       c  0.7936  0.0200  0.0372  0.591100           738168       735018  10  100000645    SNP
    2  10:100003242:SNP       t       g  0.8831  0.1417  0.0469  0.002526           738168       733070  10  100003242    SNP
    3  10:100003304:SNP       a       g  0.9609  0.0245  0.0838  0.769800           737054       663809  10  100003304    SNP
    4  10:100003785:SNP       t       c  0.6406 -0.0680  0.0313  0.029870           738169       735681  10  100003785    SNP

and it can now be loaded into a :class:`~genal.Geno` object:

.. code-block:: python

    import genal
    SBP_Geno = genal.Geno(sbp_gwas, CHR="CHR", POS="POS", EA="Allele1", NEA="Allele2", BETA="Effect", SE="StdErr", P="P", EAF="Freq1", keep_columns=False)

The last argument (``keep_columns = False``) indicates that we do not wish to keep the other (non-main) columns in the dataframe.

.. note::

   Make sure to read the readme file usually provided with the summary statistics to identify the correct columns. It is particularly important to correctly identify the allele that represents the effect allele.

Data preprocessing
===================

Now that we have loaded the data into a :class:`~genal.Geno` instance, we can begin cleaning and formatting it. Methods such as Polygenic Risk Scoring or Mendelian Randomization require the SNP data to be in a specific format. Additionally, raw summary statistics can sometimes contain missing or invalid values that need to be handled. Some columns may be missing from the data (such as the SNP rsid column or the non-effect allele column), and these columns can be created based on existing ones and a reference panel.

Genal can run all the basic cleaning and preprocessing steps in one command:

.. code-block:: python

    SBP_Geno.preprocess_data(preprocessing='Fill_delete')

The ``preprocessing`` argument specifies the global level of preprocessing applied to the data:

- ``preprocessing = 'None'``: The data won't be modified.
- ``preprocessing = 'Fill'``: Missing columns will be added based on reference data and invalid values set to NaN, but no rows will be deleted.
- ``preprocessing = 'Fill_delete'``: Missing columns will be added, and all rows containing missing, duplicated, or invalid values will be deleted. This option is recommended before running genetic methods.

Defaults to ``'Fill'``.

By default, and depending on the global preprocessing level (``'None'``, ``'Fill'``, ``'Fill_delete'``) chosen, the :meth:`~genal.Geno.preprocess_data` method of :class:`~genal.Geno` will run the following checks:

- Ensure the ``CHR`` (chromosome) and ``POS`` (genomic position) columns are integers.
- Ensure the ``EA`` (effect allele) and ``NEA`` (non-effect allele) columns are uppercase characters containing A, T, C, G letters. Multiallelic values are set to NaN.
- Validate the ``P`` (p-value) column for proper values.
- Check for no duplicated SNPs based on rsid.
- Determine if the ``BETA`` (effect) column contains beta estimates or odds ratios, and log-transform odds ratios if necessary.
- Create ``SNP`` column using a reference panel if ``CHR`` and ``POS`` columns are present.
- Create ``CHR`` and/or ``POS`` column using a reference panel if ``SNP`` column is present.
- Create ``NEA`` (non-effect allele) column using a reference panel if ``EA`` (effect allele) column is present.
- Create the ``SE`` (standard-error) column if the ``BETA`` and ``P`` (p-value) columns are present.
- Create the ``P`` column if the ``BETA`` and ``SE`` columns are present.

If you do not wish to run certain steps, or wish to run only certain steps, you can use additional arguments. For more information, please refer to the :meth:`~genal.Geno.preprocess_data` method in the API documentation.

In our case, the ``SNP`` column (for SNP identifier - rsid) was missing from our dataframe and has been added based on a 1000 genome reference panel::

    Using the EUR reference panel.
    The SNP column (rsID) has been created. 197511 (2.787%) SNPs were not found in the reference data and their ID set to CHR:POS:EA.
    The BETA column looks like Beta estimates. Use effect_column='OR' if it is a column of Odds Ratios.

You can always check the data of a ``genal.Geno`` instance by accessing the ``data`` attribute:

.. code-block:: python

    >>> SBP_Geno.data.head(5)
        EA NEA     EAF   BETA     SE        P  CHR       POS        SNP
    0    A   G  0.5660  0.0523  0.0303  0.083940   10  100000625  rs7899632
    1    A   C  0.7936  0.0200  0.0372  0.591100   10  100000645  rs61875309
    2    T   G  0.8831  0.1417  0.0469  0.002526   10  100003242  rs12258651
    3    A   G  0.9609  0.0245  0.0838  0.769800   10  100003304  rs72828461
    4    T   C  0.6406 -0.0680  0.0313  0.029870   10  100003785  rs1359508


And we see that the ``SNP`` column with the rsids has been added based on the reference data. You do not need to obtain the 1000 genome reference panel yourself, genal will download it the first time you use it. By default, the reference panel used is the European (EUR) one. You can specify another valid reference panel (AFR, EAS, SAS, AMR) with the ``reference_panel`` argument:

.. code-block:: python

    SBP_Geno.preprocess_data(preprocessing='Fill_delete', reference_panel="afr")

You can also use a custom reference panel by specifying the path to bed/bim/fam files (without the extension) in the ``reference_panel`` argument.


Clumping
--------

Clumping, or C+T: Clumping + Thresholding, is the step at which we select the SNPs that will be used as our genetic instruments in future Polygenic Risk Scores and Mendelian Randomization analyses. The process involves identifying the SNPs that are strongly associated with our trait of interest (systolic blood pressure in this tutorial) and are independent from each other. This second step ensures that selected SNPs are not highly correlated, (i.e., they are not in high linkage disequilibrium). For this step, we again need to use a reference panel.

The SNP-data loaded in a :class:`~genal.Geno` instance can be clumped using the :meth:`~genal.Geno.clump` method. It will return another :class:`~genal.Geno` instance containing only the clumped data:

.. code-block:: python

    SBP_clumped = SBP_Geno.clump(p1=5e-8, r2=0.1, kb=250, reference_panel="eur")

It will output the number of instruments obtained::

    Using the EUR reference panel.
    Warning: 760  top variant IDs missing
    1545 clumps formed from 73594 top variants.

You can specify the thresholds you want to use for the clumping with the following arguments:

* ``p1``: P-value threshold during clumping. SNPs with a P-value higher than this value are excluded. Defaults to ``5e-8``.
* ``r2``: Linkage disequilibrium threshold for the independence check. Takes values between 0 and 1. Defaults to ``0.1``.
* ``kb``: Genomic window used for the independence check (the unit is thousands of base-pair positions). Defaults to ``250``.
* ``reference_panel``: The reference population used to derive linkage disequilibrium values and select independent SNPs. Defaults to ``eur``.

Polygenic Risk Scoring
----------------------

Computing a Polygenic Risk Score (PRS) can be done in one line with the :meth:`~genal.Geno.prs` method:

.. code-block:: python

    SBP_clumped.prs(name="SBP_prs", path="path/to/genetic/files")

The genetic files of the target population can be either contained in one triple of bed/bim/fam files with information for all SNPs, or divided by chromosome (one bed/bim/fam triple for chr 1, another for chr 2, etc...). In the latter case, provide the path by replacing the chromosome number by ``$`` and genal will extract the necessary SNPs from each chromosome and merge them before running the PRS. For instance, if the genetic files are named ``Pop_chr1.bed``, ``Pop_chr1.bim``, ``Pop_chr1.fam``, ``Pop_chr2.bed``, ..., you can use:

.. code-block:: python

    SBP_clumped.prs(name="SBP_prs", path="Pop_chr$")

The ``name`` argument specifies the name of the .csv file that will be saved with the individual risk scores. 
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

Here, we see that about half of the SNPs were not extracted from the data. In such cases, we may want to try and salvage some of these SNPs by looking for proxies (SNPs in high linkage disequilibrium, i.e. highly correlated SNPs). This can be done by specifying the ``proxy = True`` argument:

.. code-block:: python

    SBP_clumped.prs(name="SBP_prs_proxy", path="Pop_chr$", proxy=True, reference_panel="eur", r2=0.8, kb=5000, window_snps=5000)

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

In our case, we have been able to find proxies for 578 of the 786 SNPs that were missing in the population genetic data (7 potential proxies have been removed because they were identical to SNPs already present in our data).

You can customize how the proxies are chosen with the following arguments:

* ``reference_panel``: The reference population used to derive linkage disequilibrium values and find proxies. Defaults to ``eur``.
* ``kb``: Width of the genomic window to look for proxies (in thousands of base-pair positions). Defaults to ``5000``.
* ``r2``: Minimum linkage disequilibrium value with the original SNP for a proxy to be included. Defaults to ``0.8``.
* ``window_snps``: Width of the window to look for proxies (in number of SNPs). Defaults to ``5000``.

.. note::
   You can call the :meth:`~genal.Geno.prs` method on any :class:`~genal.Geno` instance (containing at least the EA, BETA, and either SNP or CHR/POS columns). The data does not need to be clumped, and there is no limit to the number of SNPs used to compute the scores.

Mendelian Randomization
-----------------------

To run MR, we need to load both our exposure and outcome SNP-level data in :class:`~genal.Geno` instances. In our case, the genetic instruments of the MR are the SNPs associated with blood pressure at genome-wide significant levels resulting from the clumping of the blood pressure GWAS. They are stored in our ``SBP_clumped`` :class:`~genal.Geno` instance which also include their association with the exposure trait (instrument-SBP estimates in the ``BETA`` column).

To get their association with the outcome trait (instrument-stroke estimates), we are going to use SNP-level data from a large GWAS of stroke performed by the GIGASTROKE consortium (`Nature article <https://www.nature.com/articles/s41586-022-05165-3>`_):

.. code-block:: python

    stroke_gwas = pd.read_csv("GCST90104539_buildGRCh37.tsv", sep="\t")

We inspect it to determine the column names:

.. code-block:: python

    chromosome  base_pair_location  effect_allele_frequency   beta  standard_error  p_value  odds_ratio  ci_lower  ci_upper effect_allele other_allele
    0           5            29439275                    0.3569  0.0030         0.0070  0.6658   1.003005  0.989337  1.016861            T            C
    1           5            85928892                    0.0639 -0.0152         0.0137  0.2686   0.984915  0.958820  1.011720            T            C
    2          10           128341232                    0.4613  0.0025         0.0065  0.6998   1.002503  0.989812  1.015357            T            C
    3           3            62707519                    0.0536  0.0152         0.0152  0.3177   1.015316  0.985514  1.046019            T            C
    4           2            80464120                    0.9789  0.0057         0.0254  0.8223   1.005716  0.956874  1.057052            T            G

We load it in a :class:`~genal.Geno` instance:

.. code-block:: python

    Stroke_Geno = genal.Geno(stroke_gwas, CHR="chromosome", POS="base_pair_location", 
                             EA="effect_allele", NEA="other_allele", BETA="beta", 
                             SE="standard_error", P="p_value", 
                             EAF="effect_allele_frequency", keep_columns=False)

We preprocess it as well to put it in the correct format and make sure there are no invalid values:

.. code-block:: python

    Stroke_Geno.preprocess_data(preprocessing='Fill_delete')

Now, we need to extract our instruments (SNPs of the ``SBP_clumped`` data) from the outcome data to obtain their association with the outcome trait (stroke). It can be done by calling the :meth:`~genal.Geno.query_outcome` method:

.. code-block:: python

    SBP_clumped.query_outcome(Stroke_Geno, proxy=False)

Genal will print how many SNPs were successfully found and extracted from the outcome data::

    Outcome data successfully loaded from 'b352e412' geno instance.
    Identifying the exposure SNPs present in the outcome data...
    1541 SNPs out of 1545 are present in the outcome data.
    (Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute.
    
Here as well you have the option to use proxies for the instruments that are not present in the outcome data:

.. code-block:: python

    SBP_clumped.query_outcome(Stroke_geno, proxy=True, reference_panel="eur", 
                              kb=5000, r2=0.6, window_snps=5000)

And genal will print the number of missing instruments which have been proxied::

    Outcome data successfully loaded from 'b352e412' geno instance.
    Identifying the exposure SNPs present in the outcome data...
    1541 SNPs out of 1545 are present in the outcome data.
    Searching proxies for 4 SNPs...
    Using the EUR reference panel.
    Found proxies for 4 SNPs.
    (Exposure data, Outcome data, Outcome name) stored in the .MR_data attribute.

After extracting the instruments from the outcome data, the ``SBP_clumped`` :class:`~genal.Geno` instance contains an :attr:`~genal.Geno.MR` attribute containing the instruments-exposure and instruments-outcome associations necessary to run MR. Running MR is now as simple as calling the :meth:`~genal.Geno.MR` method of the SBP_clumped :class:`~genal.Geno` instance:

.. code-block:: python

    SBP_clumped.MR(action=2, exposure_name="SBP", outcome_name="Stroke_eur")

The :meth:`~genal.Geno.MR` method prints a message specifying which SNPs have been excluded from the analysis (it depends on the action argument, as we will see)::

    Action = 2: 42 SNPs excluded for being palindromic with intermediate allele frequencies: rs11817866, rs3802517, rs2788293, rs2274224, rs7310615, rs7953257, rs2024385, rs61912333, rs11632436, rs1012089, rs3851018, rs9899540, rs4617956, rs773432, rs11585169, rs7796, rs2487904, rs12321, rs73029563, rs4673238, rs3845811, rs2160236, rs10165271, rs9848170, rs2724535, rs6842486, rs4834792, rs990619, rs155364, rs480882, rs6875372, rs258951, rs1870735, rs1800795, rs12700814, rs1821002, rs3021500, rs28601761, rs7463212, rs907183, rs534523, rs520015 

It returns a dataframe containing the results for different MR methods:

+----------+------------+--------------------------------------------+------+----------+----------+---------------+
| exposure | outcome    | method                                     | nSNP | b        | se       | pval          |
+==========+============+============================================+======+==========+==========+===============+
| SBP      | Stroke_eur | Inverse-Variance Weighted                  | 1499 | 0.023049 | 0.001061 | 1.382645e-104 |
+----------+------------+--------------------------------------------+------+----------+----------+---------------+
| SBP      | Stroke_eur | Inverse Variance Weighted (Fixed Effects)  | 1499 | 0.023049 | 0.000754 | 4.390655e-205 |
+----------+------------+--------------------------------------------+------+----------+----------+---------------+
| SBP      | Stroke_eur | Weighted Median                            | 1499 | 0.022365 | 0.001337 | 8.863203e-63  |
+----------+------------+--------------------------------------------+------+----------+----------+---------------+
| SBP      | Stroke_eur | Simple mode                                | 1499 | 0.027125 | 0.007698 | 4.382993e-04  |
+----------+------------+--------------------------------------------+------+----------+----------+---------------+
| SBP      | Stroke_eur | MR Egger                                   | 1499 | 0.027543 | 0.002849 | 1.723156e-21  |
+----------+------------+--------------------------------------------+------+----------+----------+---------------+
| SBP      | Stroke_eur | Egger Intercept                            | 1499 | -0.001381| 0.000813 | 8.935529e-02  |
+----------+------------+--------------------------------------------+------+----------+----------+---------------+

You can specify several arguments. We refer to the API for a full list, but the most important one is the ``action`` argument. It determines how palindromic SNPs are treated during the exposure-outcome harmonization step. Palindromic SNPs are SNPs where the nucleotide change reads the same forward and backward on complementary strands of DNA (for instance ``EA = 'A'`` and ``NEA = 'T'``).

- ``action = 1``: Palindromic SNPs are not treated (assumes all alleles are on the forward strand)
- ``action = 2``: Uses effect allele frequencies to attempt to flip them (conservative, default)
- ``action = 3``: Removes all palindromic SNPs (very conservative)

If you choose the option 2 or 3 (recommended), genal will print the list of palindromic SNPs that have been removed from the analysis.

By default, only some MR methods (inverse-variance weighted, weighted median, Simple mode, MR-Egger) are going to be run. But if you wish to run a different set of MR methods, you can pass a list of strings to the ``methods`` argument. The possible strings are:

- ``IVW`` for the classical Inverse-Variance Weighted method with random effects
- ``IVW-RE`` for the Inverse Variance Weighted method with Random Effects where the standard error is not corrected for under dispersion
- ``IVW-FE`` for the Inverse Variance Weighted with fixed effects
- ``UWR`` for the Unweighted Regression method
- ``WM`` for the Weighted Median method
- ``WM-pen`` for the penalised Weighted Median method
- ``Simple-median`` for the Simple Median method
- ``Sign`` for the Sign concordance test
- ``Egger`` for MR-Egger and the MR-Egger intercept
- ``Egger-boot`` for the bootstrapped version of MR-Egger and its intercept
- ``Simple-mode`` for the Simple mode method
- ``Weighted-mode`` for the Weighted mode method
- ``all`` to run all the above methods

For more fine-tuning, such as settings for the number of boostrapping iterations, please refer to the API.

If you want to visualize the obtained MR results, you can use the :meth:`~genal.Geno.MR_plot` method that will plot each SNP in an ``effect_on_exposure x effect_on_outcome`` plane as well as lines corresponding to different MR methods:

.. code-block:: python

    SBP_clumped.MR_plot(filename="MR_plot_SBP_AS")

.. image:: Images/MR_plot_SBP_AS.png
   :alt: MR plot

You can select which MR methods you wish to plot with the ``methods`` argument. Note that for an MR method to be plotted, they must be included in the latest :meth:`~genal.Geno.MR` call of this :class:`~genal.Geno` instance.

If you wish to include the heterogeneity values (Cochran's Q) in the results, you can use the heterogeneity argument in the :meth:`~genal.Geno.MR` call. Here, the heterogeneity for the inverse-variance weighted method:

.. code-block:: python

    SBP_clumped.MR(action=2, methods=["Egger","IVW"], exposure_name="SBP", outcome_name="Stroke_eur", heterogeneity=True)

And that will give:

.. code-block:: python

      exposure     outcome                      method  nSNP        b       se          pval            Q  Q_df         Q_pval
    0      SBP  Stroke_eur                   MR Egger  1499  0.027543  0.002849  1.723156e-21  2959.965136  1497  1.253763e-98
    1      SBP  Stroke_eur            Egger Intercept  1499 -0.001381  0.000813  8.935529e-02  2959.965136  1497  1.253763e-98
    2      SBP  Stroke_eur  Inverse-Variance Weighted  1499  0.023049  0.001061  1.382645e-104 2965.678836  1498  4.280737e-99


    
As expected, many MR methods indicate that SBP is strongly associated with stroke, but there could be concerns for horizontal pleiotropy (instruments influencing the outcome through a different pathway than the one used as exposure) given the almost significant MR-Egger intercept p-value.

To investigate horizontal pleiotropy in more detail, a very useful method is Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO). 
MR-PRESSO is a method designed to detect and correct for horizontal pleiotropy. 
It will identify which instruments are likely to be pleiotropic on their effect on the outcome, and it will rerun an inverse-variance weighted MR after excluding them. 
It can be run using the :meth:`~genal.Geno.MRpresso` method:

.. code-block:: python

    mod_table, GlobalTest, OutlierTest, BiasTest = SBP_clumped.MRpresso(action=2, n_iterations=30000)

As with the :meth:`~genal.Geno.MR` method, the ``action`` argument determines how the pleiotropic SNPs will be treated. The output is a list containing:

- A table containing the original and outlier-corrected inverse variance-weighted results.
- The global test p-value indicating the presence of horizontal pleiotropy.
- A dataframe of p-values, one for each instrument, representing the likelihood that this instrument is pleiotropic (only relevant if the global test is significant).
- A dictionary containing the outputs of the distortion test. This test assesses whether the removal of the pleiotropic instruments has significantly altered the original MR estimate.
    - An array containing the indices of the pleiotropic SNPs.
    - The coefficient of the distortion test.
    - The p-value of the distortion test.

SNP-association testing
-----------------------

We may want to calibrate instrument-trait estimates in a specific population for which we have individual-level data (genetic files as well as phenotypic data). For instance, if the GWAS of SBP was done in a European population, we may want to adjust the estimates based on data coming from a population of a different ancestry. This can be done in 2 steps:

1. Loading the phenotypic data in a dataframe and calling the :meth:`~genal.Geno.set_phenotype` method
2. Calling the :meth:`~genal.Geno.association_test` method to run the association tests and update the estimates

Let's start by loading phenotypic data:

.. code-block:: python

    df_pheno = pd.read_csv("path/to/trait/data")

.. note::
   One important point is to make sure that the IDs of the participants are identical in the phenotypic data and in the genetic data.

Then, it is advised to make a copy of the :class:`~genal.Geno` instance containing our instruments as we are going to update their coefficients and to avoid any confusion:

.. code-block:: python

    SBP_adjusted = SBP_clumped.copy()

We can then call the :meth:`~genal.Geno.set_phenotype` method, specifying which column contains our trait of interest (for the association testing) and which column contains the individual IDs:

.. code-block:: python

    SBP_adjusted.set_phenotype(df_pheno, PHENO="htn", IID="IID")

At this point, genal will identify if the phenotype is binary or quantitative in order to choose the appropriate regression model. If the phenotype is binary, it will assume that the most frequent value is coding for control (and the other value for case), this can be changed with ``alternate_control=True``::

    Detected a binary phenotype in the 'PHENO' column. Specify 'PHENO_type="quant"' if this is incorrect.
    Identified 0 as the control code in 'PHENO'. Set 'alternate_control=True' to inverse this interpretation.
    The phenotype data is stored in the .phenotype attribute.
    
We can then run the association tests, specifying the path to the genetic files in plink format, and any columns we may want to include as covariates in the regression tests (making sure that the covariates are all numerical):

.. code-block:: python

    SBP_adjusted.association_test(covar=["age"], path="path/to/genetic/files")

Genal will print information regarding the number of individuals used in the tests and the kind of tests performed. It is advised to make sure that this information is consistent with your data::

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
    
The ``BETA``, ``SE``, and ``P`` columns of the ``SBP_adjusted.data`` attribute have been updated with the results of the association tests.

Lifting
-------

It is sometimes necessary to lift the SNP data to a different build. For instance, if the genetic data of our target population is in build 38 (hg38), but the GWAS summary statistics are in build 37 (hg19).
This can easily be done in genal using the :meth:`~genal.Geno.lift` method:

.. code-block:: python

    SBP_clumped.lift(start="hg19", end="hg38", replace=False)

This outputs a table with the lifted SBP instruments (stored in the ``SBP_clumped`` instance) from build 37 (hg19) to build 38 (hg38). We specified ``replace=False`` to not modify the ``SBP_clumped.data`` attribute, but we may want to modify it (before running a PRS in a population stored in build 38 for instance). 
Genal will download the appropriate chain files required for the lift, and it will be done in  python by default. However, if you plan to lift large datasets of SNPs (the whole summary statistics for instance), it may be useful to install the LiftOver executable that will run faster than the python version. It can be downloaded here: `<https://genome-store.ucsc.edu/>`_ You will need to create an account, scroll down to "LiftOver program", add it to your cart, and declare that you are a non-profit user.

You can specify the path of the LiftOver executable to the ``liftover_path`` argument:

.. code-block:: python

    SBP_Geno.lift(start="hg19", end="hg38", replace=False, liftover_path="path/to/liftover/exec")

GWAS Catalog
------------

It is sometimes interesting to determine the traits associated with our SNPs. In Mendelian Randomization, for instance, we may want to exclude instruments that are associated with traits likely causing horizontal pleiotropy. For this purpose, we can use the :meth:`~genal.Geno.query_gwas_catalog` method. This method will query the GWAS Catalog API to determine the list of traits associated with each of our SNPs and store the results in a list in the ``ASSOC`` column of the ``.data`` attribute:

.. code-block:: python

    SBP_clumped.query_gwas_catalog(p_threshold=5e-8)

Which will output::

        Querying the GWAS Catalog and creating the ASSOC column. 
        Only associations with a p-value <= 5e-08 are reported. Use the p_threshold argument to change the threshold.
        To report the p-value for each association, use return_p=True.
        To report the study ID for each association, use return_study=True.
        The .data attribute will be modified. Use replace=False to leave it as is.
        100%|██████████| 1545/1545 [00:34<00:00, 44.86it/s]
        The ASSOC column has been successfully created.
        701 (45.37%) SNPs failed to query (not found in GWAS Catalog) and 7 (0.5%) SNPs timed out after 34.33 seconds. You can increase the timeout value with the timeout argument.

And the :attr:`~genal.Geno.data` attribute now contains an `ASSOC` column::

        EA NEA    EAF    BETA     SE  CHR        POS         SNP                                               ASSOC
        0  A   G  0.1784  0.2330  0.0402   10  102075479    rs603424  [eicosanoids measurement, decadienedioic acid (...]
        1  A   G  0.0706 -0.3873  0.0626   10  102403682   rs2996303                                       FAILED_QUERY
        2  T   G  0.8872  0.6846  0.0480   10  102553647   rs1006545  [diastolic blood pressure, systolic blood pressure...]
        3  T   G  0.6652 -0.2098  0.0340   10  102558506  rs12570050                                       FAILED_QUERY
        4  T   C  0.3057 -0.2448  0.0334   10  102603924   rs4919502                                       FAILED_QUERY
        5  ... ...    ...    ...    ...  ...        ...         ...                                                ...
        6  T   C  0.3514  0.2203  0.0314    9   9350706    rs1332813  [diastolic blood pressure, systolic blood pressure...]
        7  T   C  0.6880 -0.1897  0.0332    9  94201341  rs10820855                                       FAILED_QUERY
        8  A   T  0.3669 -0.1862  0.0313    9  95201540   rs7045409  [protein measurement, pulse pressure measurement...]



If you are also interested in the p-values of each SNP-trait association, or the ID of the study from which the association was reported, you can use the ``return_p = True`` and ``return_study = True`` arguments. Then, the ``ASSOC`` column will contain a list of tuples, where each tuple contains the trait name, the p-value, and the study ID:

.. code-block:: python

    SBP_clumped.query_gwas_catalog(p_threshold=5e-8, return_p=True, return_study=True)

::

      EA NEA    EAF    BETA     SE  CHR        POS         SNP                                               ASSOC
    0  A   G  0.1784  0.2330  0.0402   10  102075479    rs603424                                            TIMEOUT
    1  A   G  0.0706 -0.3873  0.0626   10  102403682   rs2996303                                       FAILED_QUERY
    2  T   G  0.8872  0.6846  0.0480   10  102553647   rs1006545  [(heart rate response to exercise, 6e-12, GCST... 
    3  T   G  0.6652 -0.2098  0.0340   10  102558506  rs12570050                                       FAILED_QUERY
    4  T   C  0.3057 -0.2448  0.0334   10  102603924   rs4919502                                       FAILED_QUERY
    5  ... ...    ...    ...    ...  ...        ...         ...                                                ...
    6  T   C  0.3514  0.2203  0.0314    9   9350706    rs1332813  [(diastolic blood pressure, 1e-12, GCST9031029...
    7  T   C  0.6880 -0.1897  0.0332    9  94201341  rs10820855                                       FAILED_QUERY
    8  A   T  0.3669 -0.1862  0.0313    9  95201540   rs7045409  [(systolic blood pressure, 9e-13, GCST006624),...


.. note::
   As you can see, many SNPs failed to be queried. This is normal as the GWAS Catalog is not exhaustive.







































