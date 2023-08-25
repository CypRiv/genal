#!/usr/bin/env Rscript

#Attention au référentiel pour les proxy! On a changé pour test.

library(argparse)

parser <- ArgumentParser()
parser$add_argument("--path", type="character", help="Path")
parser$add_argument("--name", type="character", help="Name")
parser$add_argument("--genal_mr_outcomes_path", type="character", help="Path to GENAL MR outcomes")
parser$add_argument("--plink19_path", type="character", help="Path to PLINK 1.9")
parser$add_argument("--ref_3kg_path", type="character", help="Path to 3KG reference")
parser$add_argument("--cpus", type="integer", help="Number of CPUs")
parser$add_argument("--pattern", type="character", help="File pattern")
parser$add_argument("--action", type="integer", help="Action parameter")

args <- parser$parse_args()

# Get the parameters in the argument list
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(doMC)
library(gwasvcf)

# Setup the different paths and names and create the results file
results_name <- paste0(args$path, "/tmp_GENAL/Results_MR_all_outcomes_", args$name, ".csv")
path_outcome <- args$genal_mr_outcomes_path
file_exposure <- paste0(args$path, "/tmp_GENAL/", args$name, "_clumped.txt")
Results=data.frame(matrix(ncol=9,nrow=0))
colnames(Results) <- c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval")
write.table(Results,results_name,sep=",",append=F,col.names=T)
files_outcome <- list.files(path = path_outcome,pattern=paste0(args$pattern,"$"))


## Convert the exposure to the correct format for MR
exposure=read.table(file_exposure,header=T,sep="\t")
clumped_dat <- format_data(exposure, 
                      type="exposure",
                      snp_col = "SNP", 
                      beta_col = "BETA", 
                      se_col = "SE", 
                      effect_allele_col = "EA", 
                      other_allele_col = "NEA",
                      eaf_col = "EAF",
                          pval_col="P")
gwasvcf::set_bcftools()
gwasvcf::set_plink(args$plink19_path)

# Declare the parallel loop
registerDoMC(args$cpus)

foreach (i=1:length(files_outcome))%dopar% {
    outcome=files_outcome[i]
    
    ## Try to query_gwas. This can fail especially if the exposure is from a recent GWAS because newly labeled SNPs are not present in the reference panel. If that is the case and the function fails, we try again by first removing the problematic SNPs from the exposure data.
    t=try(ht_outcome_dat <- gwasvcf::query_gwas(VariantAnnotation::readVcf(paste0(path_outcome,"/",outcome),"hg19"),rsid = clumped_dat$SNP, proxies = "yes",bfile=paste0(args$ref_3kg_path,"EUR")))
    if ("try-error" %in% class(t)){
        if (grepl("tag_r2",t[1])){
            print("Some of the SNPs in the exposure are not present in the reference panel, retrying after deleting them.")
            bfile_bim=read.table(paste0(args$ref_3kg_path,"EUR"))
            i=intersect(clumped_dat$SNP,bfile_bim$V2)
            clumped_dat_i=clumped_dat[clumped_dat$SNP %in% i,]
            ht_outcome_dat <- gwasvcf::query_gwas(VariantAnnotation::readVcf(paste0(path_outcome,"/",outcome),"hg19"),rsid = clumped_dat_i$SNP, proxies = "yes",bfile=paste0(args$ref_3kg_path,"EUR"))
            }
        else {
            print ("Error in query_gwas (different from the missing SNPs in reference panel).")
            }
    }
    ## Once the query is done, convert to the TwoSampleMR format, harmonize exposure and outcome, and run MR
    ht_outcome_dat <- gwasglue::gwasvcf_to_TwoSampleMR(ht_outcome_dat, "outcome")
    dat <- harmonise_data(clumped_dat, ht_outcome_dat,action=args$action)
    res <- mr(dat)
    
    ## Append the MR results to the results file
    if (nrow(res)>0){
        name_outcome=str_replace(outcome,".vcf.gz","")
        res$exposure=args$name
        res$outcome=name_outcome
        write.table(res,results_name,sep=",",append=T,col.names=F,row.names=F)
        }
    }
