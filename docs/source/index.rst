.. genal documentation master file, created by
   sphinx-quickstart on Thu Sep 14 14:04:16 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: Images/genal_logo.png
   :alt: genal_logo
   :width: 400px

genal: A Python Toolkit for Genetic Risk Scoring and Mendelian Randomization
============================================================================

:Author: Cyprien A. Rivier
:Date: |today|
:Version: "1.2"

Genal is a python module designed to make it easy to run genetic risk scores and mendelian randomization analyses. It integrates a collection of tools that facilitate the cleaning of single nucleotide polymorphism data (usually derived from Genome-Wide Association Studies) and enable the execution of key clinical population genetic workflows. The functionalities provided by genal include clumping, lifting, association testing, polygenic risk scoring, and Mendelian randomization analyses, all within a single Python module.

The module prioritizes user-friendliness and intuitive operation, aiming to reduce the complexity of data analysis for researchers. Despite its focus on simplicity, Genal does not sacrifice the depth of customization or the precision of analysis. Researchers can expect to maintain analytical rigour while benefiting from the streamlined experience.

Genal draws on concepts from well-established R packages such as TwoSampleMR, MR-Presso, MendelianRandomization, and gwasvcf, adapting their proven methodologies to the Python environment. This approach ensures that users have access to tried and tested techniques with the versatility of Python's data science tools. 

To install the latest release, type::

    pip install genal-python

Contents
--------

.. toctree::
   :maxdepth: 1
   
   Home <self>
   introduction.rst
   modules.rst
   api.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
--------

If you use genal in your work, please cite the following paper:

.. [Rivier.2024] *Genal: A Python Toolkit for Genetic Risk Scoring and Mendelian Randomization*
   Cyprien A. Rivier, Santiago Clocchiatti-Tuozzo, Shufan Huo, Victor Torres-Lopez, Daniela Renedo, Kevin N. Sheth, Guido J. Falcone, Julian N. Acosta. 
   Bioinformatics Advances. 2024 December; `10.1093/bioadv/vbae207 <https://doi.org/10.1093/bioadv/vbae207>`_.

References
----------

.. [Hemani.2018] *The MR-Base platform supports systematic causal inference across the human phenome.*
   Hemani G, Zheng J, Elsworth B, Wade KH, Baird D, Haberland V, Laurin C, Burgess S, Bowden J, Langdon R, Tan VY, Yarmolinsky J, Shihab HA, Timpson NJ, Evans DM, Relton C, Martin RM, Davey Smith G, Gaunt TR, Haycock PC, The MR-Base Collaboration
   eLife. 2018 May `10.7554/eLife.34408 <https://elifesciences.org/articles/34408>`_.
   PMID: `29846171 <https://pubmed.ncbi.nlm.nih.gov/29846171>`_.

.. [Verbanck.2018] *Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases.*
   Marie Verbanck, Chia-Yen Chen, Benjamin Neale, Ron Do.
   Nature Genetics 2018 May `10.1038/s41588-018-0099-7 <https://www.nature.com/articles/s41588-018-0099-7>`_.
   PMID: `29686387 <https://pubmed.ncbi.nlm.nih.gov/29686387/>`_.

.. [Lyon.2020] *The variant call format provides efficient and robust storage of GWAS summary statistics.*
   Matthew Lyon, Shea J Andrews, Ben Elsworth, Tom R Gaunt, Gibran Hemani, Edoardo Marcora.
   bioRxiv 2020 May 30 `2020.05.29.115824v1 <https://www.biorxiv.org/content/10.1101/2020.05.29.115824v1>`_.
   PMID: `33441155 <https://pubmed.ncbi.nlm.nih.gov/33441155/>`_.
   