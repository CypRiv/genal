MR <- function(data,path,action,ref_3kg_path,sensitivity,n){
    library(tidyverse)
    library(TwoSampleMR)
    library(ieugwasr)
    
    clumped_dat <- format_data(data, 
                          type="exposure",
                          snp_col = "SNP", 
                          beta_col = "BETA", 
                          se_col = "SE", 
                          effect_allele_col = "NEA", 
                          other_allele_col = "EA",
                          eaf_col = "EAF",
                           pval_col="P")
    gwasvcf::set_bcftools()
    gwasvcf::set_plink("/gpfs/ysm/project/falcone/gf272/software/plink2")
    ## Try to query_gwas. This can fail especially if the exposure is from a recent GWAS because newly labeled SNPs are not present in the reference panel. If that is the case and the function fails, we try again by first removing the problematic SNPs from the exposure data.
    t=try(ht_outcome_dat <- gwasvcf::query_gwas(VariantAnnotation::readVcf(path,"hg19"),rsid = clumped_dat$SNP, proxies = "yes",
                                      bfile=paste0(ref_3kg_path,"data_maf0.01_rs_ref")))
    if ("try-error" %in% class(t)){
        if (grepl("tag_r2",t[1])){
            print("Some of the SNPs in the exposure are not present in the reference panel, retrying after deleting them.")
            bfile_bim=read.table(paste0(ref_3kg_path,"data_maf0.01_rs_ref.bim"))
            i=intersect(clumped_dat$SNP,bfile_bim$V2)
            clumped_dat_i=clumped_dat[clumped_dat$SNP %in% i,]
            ht_outcome_dat <- gwasvcf::query_gwas(VariantAnnotation::readVcf(path,"hg19"),rsid = clumped_dat_i$SNP, proxies = "yes",bfile=paste0(ref_3kg_path,"data_maf0.01_rs_ref"))
            }
        else {
            print ("Error in query_gwas (different from the missing SNPs in reference panel).")
            }
    }
    ht_outcome_dat <- gwasglue::gwasvcf_to_TwoSampleMR(ht_outcome_dat, "outcome")
    dat <- harmonise_data(clumped_dat, ht_outcome_dat,action=action)
    res <- mr(dat)
    if (sensitivity==FALSE){return (res)}
    else {
        if (res[3,9] < 0.05){
            res_presso <- run_mr_presso(dat, NbDistribution = n)
            return (list(res,res_presso,dat))
        }
        else {
            res_presso=data.frame(matrix(ncol=0,nrow=0))
            return (list(res,res_presso,dat))
        }
    }
}



