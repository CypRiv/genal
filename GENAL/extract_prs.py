

































def prs_regress(self, weighted=True,clumped=False,model="add",alternate_control=False,fastscore=False,error_1_code=2,
               step=5e-5,lower=5e-8,upper=0.5,maf=None,cpus=4,cpus_batch=8):
    """Run PRSice in UKB data with regression against the phenotype defined with the set_phenotype function and search for the ideal p-value threshold.
    The function launches a batchscript. The result can be retrieved from prs_regress_result function once the computation is done. prs_regress_result can also be used to monitor the progress of the batchscript.
    weighted=False will put all betas to 1 to create an unweighted PRS 
    By default, perform the search starting with the unclumped data, but it is possible to do it on the clumped data with clumped==True.
    model="add" is the genetic model used for regression: add for additive, dom for dominant, rec for recessive, het for heterozygous
    alternate_control=True if the control group is smaller than the case group (unlikely)
    fastscore=True to compute only at the main thresholds (and not in between)
    maf=None will threshold by minor allele frequency before performing the search
    error_1_code=2 allows to control the number of times PRSice should be rerun excluding duplicates. In many cases that will be 2 (default), but sometimes, only 1 time is necessary (then set error_1_code=1), or not at all (=0)
    step=5e-5 is the length of the step the software takes when iterating over the P theresholds, the lower the more precise but the more computationally intensive
    lower=5e-8 and upper=0.5 indicate the bounds of the lookup over the P thresholds. A common practice is to start by taking big steps on a big interval and then rerun the function with smaller steps on a smaller interval around the most significant threshold to obtain a precise result.
    cpus=4 is the number of threads used by PRSice locally
    cpus_batch=8 is the number of threads used by PRSice from the batch script

    """

    ##Check that a phenotype has been set with the set_phenotype function.
    if not(hasattr(self,"phenotype")):
        raise ValueError("You first need to set a phenotype with .set_phenotype(data,PHENO,PHENO_type,IID)!") 

    ## Set the BETAs to 1 if weighted==False
    data_to_prs=self.data.copy()
    if weighted==False:
        data_to_prs["BETA"]=1
        print("Computing an unweighted prs.")

    ## Check the mandatory columns
    for column in ["SNP","P","EA","NEA","BETA"]:
        if column not in data_to_prs.columns:
            raise ValueError("The column {column} is not found in the data!".format(column=column))

    ## For a binary trait: code it in 1/2 to comply with PRSice requirements and put the columns to Int64 format.
    phenotype_to_prs=self.phenotype.copy()
    if self.phenotype.type=="binary":
        phenotype_to_prs["PHENO"]+=1
        phenotype_to_prs["PHENO"]=phenotype_to_prs.PHENO.astype("Int64")

    ## Go the tmp folder, write the data and phenotype to_csv
    if not os.path.exists("tmp_GENAL"):
        os.makedirs("tmp_GENAL")
    data_to_prs=data_to_prs.loc[:,data_to_prs.columns.isin(["SNP","P","EA","NEA","BETA","EAF"])]
    data_to_prs.to_csv(f"tmp_GENAL/To_prs_{self.name}.txt",sep="\t",index=False,header=True)
    phenotype_to_prs.to_csv(f"tmp_GENAL/PHENO_{self.name}.txt",sep="\t",index=False,header=True)

    ## Call PRSice. Add arguments if binary trait, if fastscore, if maf threshold. 
    current_path=os.getcwd()+"/tmp_GENAL"
    command=f'Rscript {prsice_path}PRSice.R \
    --dir {prsice_path} \
    --prsice {prsice_path}PRSice_linux \
    --base {current_path}/To_prs_{self.name}.txt --A1 EA --A2 NEA --pvalue P --snp SNP --stat BETA \
    --target {ukb_geno_path}plinkfiltered_# \
    --type bed --out {current_path}/prs_regress_{self.name} --beta \
    --missing SET_ZERO --score avg \
    --seed 33 --no-clump --thread {cpus} --memory 90Gb \
    --bar-levels 5e-08,1e-05,5e-05,0.001,0.05,1 --lower {lower} --upper {upper} --interval {step} \
    --model {model} \
    --pheno {current_path}/PHENO_{self.name}.txt --pheno-col PHENO'

    if self.phenotype.type=="binary":
        command=command+" --binary-target T"
    else:
        command=command+" --binary-target F"
    if fastscore:
        command=command+" --fastscore"
    if maf!=None:
        if "EAF" not in data_to_prs.columns:
            print("A minor allele frequency column must be present in the data in order to use the maf threshold.")
        else:
            command=command+f" --base-maf EAF:{maf}"

    ## Handles a common PRSice error: duplicated SNP ID which requires to extract the SNP list from the first run
    if error_1_code>0:
        output=subprocess.run(command, shell=True,capture_output=True,text=True,check=False)
        if "duplicated SNP ID detected out of" in output.stderr:
            print("Rerunning PRSice after excluding the duplicated SNPs")
            command=command+f" --extract {current_path}/prs_regress_{self.name}.valid"
            if error_1_code>1:
                output2=subprocess.run(command, shell=True,capture_output=True,text=True,check=False)
                if "duplicated SNP ID detected out of" in output2.stderr:
                    print("Rerunning PRSice after excluding the duplicated SNPs (2)")

    ## Declare the name of the job, the name of the .sh file, cancel an existing job if it has the same name and launch the script.
    bash_name=f"prs_regress_{self.name}"
    bash_script_name=f"prs_regress_{self.name}.sh"

    if bash_name[:18] in subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout:
        print("Canceling a job with the same name. If it was not intentional, be careful not to call this function several times in a row. Call prs_regress_result to inspect the progression of the job.")
        subprocess.run(f"scancel --name={bash_name}",shell=True,check=True)

    #If we change the specs, don't forget to change the ram and threads in the PRSice command
    os.chdir("tmp_GENAL")
    create_bashscript(partition="pi_falcone", job_name=bash_name, bash_output_name=bash_name, ntasks=1, cpus=cpus_batch, mem_per_cpu=50000, bash_filename=bash_script_name, command=command)
    subprocess.run(["sbatch", f"{bash_script_name}"],check=True,text=True)
    os.chdir("..")

    return 

def prs_regress_result(self):
    """Function to check the status of the batchscript launched by prs_regress. 
    If it has ended properly, load the results in a PRS format.
    """
    ##Get our jobs in queue and check if the name corresponding to the prs regress job is listed.
    bash_name=f"prs_regress_{self.name}"
    summary_name=f"tmp_GENAL/prs_regress_{self.name}.summary"
    status=subprocess.run(["squeue","--me"],capture_output=True,text=True).stdout
    ##If the job is still in queue: print the time it has been running as well as the last line of the job status file.
    if bash_name[:18] in status:
        m=status.split(bash_name+" ")[1].split("R    ")[1].split(":")
        print(f"The job is still running. It has been running for: {m[0][-3:]+':'+m[1][:3]}")
        out=subprocess.run(f"tail -n 1 tmp_GENAL/{bash_name}",shell=True,capture_output=True,text=True,check=True)
        update_string=out.stdout.split('\n')[-1]
        print(f"Progression update: {update_string}")
    ##If not: if the summary file doesn't exist --> the job was not successful. If it exists --> job was successful, summary file is loaded and printed, result file is loaded and returned.
    else:
        if not os.path.isfile(summary_name):
            print("The prs_regress job has ended but was not successful. Refer to the .log file to see the error.")
        elif os.path.isfile(summary_name):
            print ("The job has ended and was successfull!")
            summary=pd.read_csv(summary_name,sep="\t")
            print (f"The p-value threshold selected is {summary.Threshold.values[0]}, it corresponds to {summary.Num_SNP.values[0]} SNPs, and it reached a regression P-value of {summary.P.values[0]:.3f} against the phenotype.")
            df_score=pd.read_csv(f"tmp_GENAL/prs_regress_{self.name}.best",sep=" ")
            df_score=df_score[["FID","IID","PRS"]].rename(columns={PRS:"SCORE"})
            return df_score