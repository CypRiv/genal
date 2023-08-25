<center><h1> GENAL: Python module to automate population genetic data processing and analysis </h1></center>



**This project was develop by Dr. Falcone Lab**
https://falconelab.org

# Table of contents
1. [Introduction](#introduction)
2. [Requirements for the GENAL module](#paragraph1)
3. [Installation and how to use GENAL](#paragraph2)
    1. [Installation](#subparagraph1)
4. [Tutorial and presentation of the main tools] 




## Introduction <a name="introduction"></a>
The focus of this project is to modularize the code used for the processing and analysis of GWAS-derived genetic data. Its main purpose is to save time by removing the need to manually adapt code. Secondary purposes include the prevention of errors during the data manipulation and processing steps and making genetic analysis accessible to researchers with little coding experience.    
The tools included in this module revolve around the data manipulation of GWAS-derived summary statistics, Polygenic Risk Scoring, and Mendelian Randomization.   
GWAS: Genome Wide Association Study

## Requirements for the GENAL module <a name="paragraph1"></a>
***Anaconda or Miniconda Environments*** https://docs.conda.io/en/latest/miniconda.html <br> 
***Python 3.7 or later***. https://www.python.org/ or  https://www.python.org/downloads/release/python-379/ <br> 
***R 4.1 or later***. https://www.r-project.org. 


## Installation and How to use the GENAL module <a name="paragraph2"></a>

### Installation <a name="subparagraph1"></a>
Optional: request enough memory. For instance, if you are working on a server using slurm:   
```
srun --pty -p interactive --mem=32G --constraint=oldest -c 4 bash
```
Navigate to the path where you wish to install GENAL. Download it there:    
```
git clone https://github.com/Vlaati/GENAL
```
Create the genal_env miniconda environment based on the .yml file:   
```
cd GENAL
module load miniconda
mamba env create -f genal_env.yml
mamba install -n genal_env bioconductor-biocgenerics bioconductor-s4vectors bioconductor-variantannotation r-processx bioconductor-genomicfeatures bioconductor-bsgenome
conda activate genal_env
```
Install the Mendelian Randomization R packages from MRCIEU:   
(Don't copy paste the "> ". When prompted whether to update the dependencies or not, always skip updates by selecting option 3: None).   
```
R
> IRkernel::installspec()
> remotes::install_github("MRCIEU/TwoSampleMR")
> remotes::install_github("MRCIEU/gwasvcf")
> remotes::install_github("MRCIEU/gwasglue")
> q()
```
Optional: install bash kernel if you plan to use it:   
```
conda install -c conda-forge bash_kernel 
python -m bash_kernel.install
```

### Using the module <a name="subparagraph2"></a>
In a jupyter notebook or another python interface   

