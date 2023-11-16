<center><h1> genal: A Python Toolkit for Genetic Risk Scoring and Mendelian Randomization </h1></center>



**This project was develop by Cyprien A. Rivier**

# Table of contents
1. [Introduction](#introduction)
2. [Requirements for the GENAL module](#paragraph1)
3. [Installation and how to use GENAL](#paragraph2)
    1. [Installation](#subparagraph1)
4. [Tutorial and presentation of the main tools] 




## Introduction <a name="introduction"></a>
Genal is a python module designed to make it easy to run genetic risk scores and mendelian randomization analyses. It integrates a collection of tools that facilitate the cleaning of single nucleotide polymorphism data (usually derived from Genome-Wide Association Studies) and enable the execution of clinical population genetic workflows. The functionalities provided by genal include clumping, lifting, association testing, polygenic risk scoring, and Mendelian randomization analyses, all within a single Python module.

The module prioritizes user-friendliness and intuitive operation, aiming to reduce the complexity of data analysis for researchers. Despite its focus on simplicity, Genal does not sacrifice the depth of customization or the precision of analysis. Researchers can expect to maintain analytical rigour while benefiting from the streamlined experience.

Genal draws on concepts from well-established R packages such as TwoSampleMR, MR-Presso, MendelianRandomization, and gwasvcf, adapting their proven methodologies to the Python environment. This approach ensures that users have access to tried and tested techniques with the versatility of Python's data science tools. 

## Requirements for the GENAL module <a name="paragraph1"></a> 
***Python 3.7 or later***. https://www.python.org/ or  https://www.python.org/downloads/release/python-379/ <br> 


## Installation and How to use the GENAL module <a name="paragraph2"></a>

### Installation <a name="subparagraph1"></a>

Download and install the package with pip:    
```
pip install genal
```

The main genal functionalities require a working installation of PLINK v1.9 that can be downloaded here: https://www.cog-genomics.org/plink/ 
Once downloaded, the path to the plink executable can be set with:

```
genal.set_plink(path="/path/to/plink/executable/file")
```

### Tutorial <a name="subparagraph2"></a>
In a jupyter notebook or another python interface   