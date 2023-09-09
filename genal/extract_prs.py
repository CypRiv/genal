import pandas as pd
import numpy as np
import subprocess
from functools import partial
from concurrent.futures import ProcessPoolExecutor

from .tools import *


def prs_func(data, weighted, path, checks=[], ram=10000, name=""):
    """Compute a PRS
    weighted=True to perform a PRS weighted by the BETA column estimates. False for an unweighted PRS (equivalent to BETAs == 1)
    path = "". Can be used to provide a path to a bed/bim/fam set of genetic files to use for PRS calculation. If not provided, will use the genetic data extracted with the .extract_snps method.
    checks=list: list of column checks already performed on the data ("SNP", "P") to avoid repeating them.
    ram = ram: ram memory in MB to be used by plink
    """
    
    ## Check the path. If None, use data extracted with extract_snps
    if path is None:
        path = os.path.join('tmp_GENAL',f'{name}_allchr')
        if not check_bfiles(path):
            raise TypeError("No valid path found. You should either provide a path with .prs(path='...') or first call the .extract_snps(path='...') method.")
    else:
        path = os.path.splitext(path)[0] 
        if not check_bfiles(path):
            raise TypeError("The path provided does not lead to a valid set of bed/bim/fam genomic files.")
            
    ## Check mandatory columns
    for column in ["SNP","EA","BETA"]:
        if not(column in data.columns):
            raise ValueError(f"The column {column} is not found in the data!")
    data_prs = data.copy()
    if "SNP" not in checks:
        data_prs = check_snp_column(data_prs)
    if "EA" not in checks:
        data_prs = check_allele_columm(data_prs, "EA", keep_multi=False)
    ## Set the BETAs to 1 if weighted==False
    if weighted==False:
        data_prs["BETA"]=1
        print("Computing an unweighted PRS.")
    else:
        print("Computing a weighted PRS.")
        
    ## Write the data to_csv in the tmp folder and call plink on it
    data_prs=data_prs[["SNP","EA","BETA"]]
    data_prs_path = os.path.join("tmp_GENAL", "To_prs.txt")
    output_path = os.path.join("tmp_GENAL", f"prs_{name}")
    data_prs.to_csv(data_prs_path, sep="\t", index=False, header=True)
    output=subprocess.run(f"{get_plink19_path()} --memory {ram} --bfile {path} \
    --score {data_prs_path} 1 2 3 header --out tmp_GENAL/prs_{name}", shell=True, capture_output=True, text=True, check=True)

    ## Read the results file, change columns names
    prs_file = output_path + ".profile"
    if os.path.isfile(prs_file):
        print("The PRS computation was successfull!")
        df_score=pd.read_csv(prs_file, sep="\s+")
        df_score=df_score[["FID","IID","SCORE"]]
        return df_score
    else:
        print (output.stdout)
        raise ValueError(f"The PRS computation was not successfull. Check the {output_path + '.log'} file.")
    
    
#Test the non-split option

def extract_command_parallel(task_id, name, path, snp_list_path):
    bfile_path = path.replace('$',str(task_id))
    if not check_bfiles(bfile_path):
        return task_id
    output_path = os.path.join('tmp_GENAL',f'{name}_extract_chr{task_id}')
    command = f"{get_plink19_path()} --bfile {bfile_path} --extract {snp_list_path} --make-bed --out {output_path}"
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return
    
def extract_snps_func(snp_list, name, path=""): 
    
    config = read_config()
    if path is None:
        path = config["paths"]["geno_path"]
        print(f"Extracting from saved path: {path}")
    else:
        # Check path
        if path.count("$") > 1:
            raise TypeError("The path should contain at most 1 '$' sign. Use it to indicate the chromosome number if the data is split by chromosomes.")
        path = os.path.splitext(path)[0] #Delete file extension
        config["paths"]["geno_path"] = path # Add path to config.json file
        write_config(config)

    filetype = "split" if path.count("$") == 1 else "combined"
    
    #Define names and write SNP list
    snp_list_name=f"{name}_list.txt"
    bedlist_name=f"{name}_bedlist.txt"
    snp_list = snp_list.dropna()
    snp_list_path = os.path.join("tmp_GENAL",snp_list_name)
    snp_list.dropna().to_csv(snp_list_path, sep=" ", index=False, header=None)
    nrow=len(snp_list)
    
    #Extract job
    if filetype == "split":
        print("Extracting SNPs for each chromosome...")
        num_tasks = 22 
        partial_extract_command_parallel = partial(extract_command_parallel, name=name, path=path, snp_list_path=snp_list_path) #Wrapper function freezing the extract_command_parallel call
        with ProcessPoolExecutor() as executor:
            not_found = list(executor.map(partial_extract_command_parallel, range(1, num_tasks + 1)))
    else:
        print("Extracting SNPs...")
        output_path = os.path.join('tmp_GENAL',f'{name}_allchr')
        extract_command = f"{get_plink19_path()} --bfile {path} --extract {snp_list_path} --make-bed --out {output_path}"
        subprocess.run(extract_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #Merge job
    if filetype == "split":
        #Create list of chromosomes with at least 1 extracted SNP
        bedlist_path = os.path.join('tmp_GENAL', bedlist_name)
        create_bedlist(bedlist_path, os.path.join("tmp_GENAL", f"{name}_extract"), not_found) 
        print("Merging SNPs extracted from each chromosome...") 
        output_path = os.path.join('tmp_GENAL',f'{name}_allchr')
        merge_command = f"{get_plink19_path()} --merge-list {bedlist_path} --make-bed --out {output_path}"
        subprocess.run(merge_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Created merged bed/bim/fam fileset: {output_path}")
        
        def remove_multiallelic():
            snps_to_exclude=pd.read_csv(os.path.join("tmp_GENAL",f"{name}_allchr-merge.missnp"), header=None)
            for i in range(22):
                bim_path = os.path.join("tmp_GENAL",f"{name}_extract_chr{i+1}.bim")                       
                if os.path.isfile(bim_path):
                    bim=pd.read_csv(bim_path, sep="\t", header=None)
                    if len(set(bim[1]).intersection(set(snps_to_exclude[0]))) > 0:
                        bfile_path = os.path.join('tmp_GENAL', f'{name}_extract_chr{i+1}')
                        missnp_path = os.path.join('tmp_GENAL', f'{name}_allchr-merge.missnp')
                        output_path = os.path.join('tmp_GENAL', f'{name}_extract_chr{i+1}')
                        command = f"{get_plink19_path()} --bfile {bfile_path} --exclude {missnp_path} --make-bed --out {output_path}"
                        subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return
        
        ## Check for the most common error that can occur while merging: multiallelic variants. If the error is present, rerun the merge excluding them.            
        if "with 3+ alleles present" in open(os.path.join("tmp_GENAL",f"{name}_allchr.log")).read():
            print("Multiallelic variants detected: removing them before merging.")
            remove_multiallelic()
            print(f"Reattempting the merge after deletion of multiallelic variants.")
            subprocess.run(merge_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
        ## In rare cases, plink fails to identify all the multiallelic SNPs with the first pass and we have to exclude again.
            if "with 3+ alleles present" in open(os.path.join("tmp_GENAL",f"{name}_allchr.log")).read():
                print("Multiallelic variants still detected: removing them before merging.")
                remove_multiallelic()
                print(f"Reattempting the merge after deletion of multiallelic variants.")
                subprocess.run(merge_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    ## Report the number of SNPs not found in UKB data
    delta_nrow = nrow-int(subprocess.check_output(['wc', '-l', f"{os.path.join('tmp_GENAL',f'{name}_allchr.bim')}"]).split()[0])
    if delta_nrow > 0: print(f"{delta_nrow} SNPs were not extracted because they are not present in the data.")
    return


def create_bedlist(bedlist, output_name, not_found):
    with open(bedlist, "w+") as bedlist_file:
        for i in range(1,23):
            if i in not_found:
                print("bed/bim/fam files not found for chr{i}.")
            elif check_bfiles("{}_chr{}".format(output_name,i)):
                bedlist_file.write("{}_chr{}\n".format(output_name,i))
                print("SNPs extracted for chr{}.".format(i))
            else:
                print("No SNPs extracted for chr{}".format(i))












