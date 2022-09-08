tovcf <- function(data,path,name){
        library(tidyverse, quietly = T)
        library(gwasvcf, quietly = T)
    
    out <- create_vcf(chrom=data$CHR, 
                  pos=data$POS, 
                  nea=data$NEA,
                  ea=data$EA, 
                  snp=data$SNP,
                  effect=data$BETA, 
                  se=data$SE,
                  pval=data$P, ea_af=data$EAF,
                      name=name)
    writeVcf(out, file=path)
    
    return}



