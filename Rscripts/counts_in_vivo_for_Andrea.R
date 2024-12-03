#!/bioinfo/local/build/Centos/R/R-4.2.1/bin/Rscript

library(tidyverse)


counts_in_vivo <- readRDS("counts_in_vivo_jpp.rds")
head(counts_in_vivo)
counts_in_vivo_tpm <- read.csv("tablecounts_tpm_InVivo.csv")
colnames(counts_in_vivo_tpm)[1] <- "gene_id"
tpm_updated <- counts_in_vivo[,c("gene_id","gene_name")] %>% dplyr::inner_join(counts_in_vivo_tpm, by="gene_id") %>% as.data.frame() %>% dplyr::select(-"gene_id")
rownames(tpm_updated) <- make.names(tpm_updated$gene_name, unique=TRUE)
tpm_updated <- tpm_updated %>% dplyr::select(-"gene_name")
tpm_updated <- tpm_updated %>% dplyr::select(-1)
sample_plan_in_vivo <- readRDS("sample_plan_in_vivo.rds")

for (i in 1:ncol(tpm_updated)){
  colnames(tpm_updated)[i] <- sample_plan_in_vivo$new_name[which(sample_plan_in_vivo$`sample ID` %in% colnames(tpm_updated)[i])]
}


tpm_updated_MTX <- tpm_updated %>% dplyr::select(starts_with("MTX"))
tpm_updated_OXA <- tpm_updated %>% dplyr::select(starts_with("Oxa"))
tpm_updated_CIS <- tpm_updated %>% dplyr::select(starts_with("Cis")) 
tpm_updated_PBS <- tpm_updated %>% dplyr::select(starts_with("PBS"))
tpm_updated_PBS_MTX <- tpm_updated_PBS %>% dplyr::select(!contains("OCP"))
tpm_updated_PBS_OCP <- tpm_updated_PBS %>% dplyr::select(contains("OCP"))



andrea_function <- function(table_counts, treatment){
  for (i in 1:nrow(table_counts)){
    for (j in 1:ncol(table_counts)){
      table_counts[i,j] <- table_counts[i,j] - rowMeans(table_counts)[i]
    }
  }
  saveRDS(table_counts, paste("counts_in_vivo_andrea",treatment, sep="_"))
}


andrea_function(tpm_updated_MTX, "MTX")
andrea_function(tpm_updated_OXA, "OXA")
andrea_function(tpm_updated_CIS, "CIS")
andrea_function(tpm_updated_PBS_MTX, "PBS_MTX")
andrea_function(tpm_updated_PBS_OCP, "PBS_OCP")


