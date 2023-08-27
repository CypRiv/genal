import os
import subprocess
import pandas as pd

def clump_data(data, plink19_path, reference_panel_path, kb=250, r2=0.1, p1=5e-8, p2=0.01, name="noname", ram = 10000):
        """ Clump the data in data and return it. The clumping is done with plink.
        data is a pandas dataframe with the standard GENO column names (at least SNP and P)
        name is the name used for the files created in the tmp_GENAL folder
        plink19_path: the path to plink 1.9
        ram: the amount of RAM in MB to be used by plink
        reference_panel_path = the path to the bed/bim/fam of the reference population
        kb=250 sets in thousands the window for the clumping
        r2=0.1 sets the linkage disequilibrium threshold 
        p1=5e-8 sets the p-value used during the clumping (the SNPs above this threshold are not considered)
        p2=0.01 sets the p-value used after the clumping to further filter the clumped SNPs (set p2<p1 to not use this
        """
        ## Check the existence of the required columns. Create the tmp file if necessary.
        for column in ["SNP","P"]:
            if not(column in data.columns):
                raise ValueError(f"The column {column} is not found in the data")
        if not os.path.exists("tmp_GENAL"):
            os.makedirs("tmp_GENAL")

        ## Write only the necessary column to file. Call plink with the correct arguments.
        data[["SNP","P"]].to_csv(f"tmp_GENAL/{name}_to_clump.txt", index=False, sep="\t")
        output=subprocess.run([f"{plink19_path} \
        --memory {ram} \
        --bfile {reference_panel_path} \
        --clump tmp_GENAL/{name}_to_clump.txt --clump-kb {kb} --clump-r2 {r2} --clump-p1 {p1} \
        --clump-p2 {p2} --out tmp_GENAL/{name}"], shell=True, capture_output=True, text=True, check=True)
        
        ## Print common outcomes. If no SNP were found: exit the function. 
        if "more top variant IDs missing" in output.stderr:
            n_missing=output.stderr.split('more top variant IDs missing')[0].split('\n')[-1]
            print(f"Warning: {n_missing} top variant IDs missing")
        if "No significant --clump results." in output.stderr:
            print("No SNPs remaining after clumping.")
            return
        clumped_message=output.stdout.split("--clump: ")[1].split("\n")[0]
        print(clumped_message)
        
        ## Get the list of clumped SNPs. Take the corresponding data subset from .data and attribute it to .data_clumped and return it too.
        with open(f"tmp_GENAL/{name}.list", 'wb') as f:
            subprocess.run([f"awk '{{print $3}}' tmp_GENAL/{name}.clumped"], shell=True, text=True, check=True, stdout=f)
        plink_clumped=pd.read_csv(f"tmp_GENAL/{name}.list",sep=" ")
        clumped_data=data[data["SNP"].isin(plink_clumped["SNP"])]
        clumped_data.reset_index(drop=True, inplace=True)
        return clumped_data