#!/bioinfo/local/build/Centos/R/R-4.2.1/bin/Rscript


library(clusterProfiler)
library(tidyverse)
run_gsea <- function(top_table_filename, contrast_interest, drug_directory, root_interest){
  

  top_table_contrast_interest <- read.csv(top_table_filename) %>% filter(contrast == contrast_interest) %>% mutate(rank = sign(logFC) * -log10(adj.P.Val)) %>% arrange(desc(rank))
  gene_list_interest <- top_table_contrast_interest$rank
  names(gene_list_interest) <- top_table_contrast_interest$genes
  gsea_contrast_interest <- gseGO(gene_list_interest, ont = "BP", OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", eps=0)
  gsea_contrast_interest_up_reg_pathways <- gsea_contrast_interest %>% as.data.frame() %>% filter(NES > 0 & p.adjust <= 0.05) %>% arrange(desc(NES))
  write.csv(gsea_contrast_interest_up_reg_pathways, paste0(root_interest,drug_directory,"/up_enriched_sets_gsea_",contrast_interest,".csv"), col.names = TRUE, row.names = TRUE)
  
  
  for (i in 1:dim(gsea_contrast_interest_up_reg_pathways)[1]) {
    
    gsea_contrast_interest_up_reg_pathways$core_enrichment[i] <- str_split(gsub("/"," ", gsea_contrast_interest_up_reg_pathways$core_enrichment[i]), " ")
    
  }
  
  
  
  leading_edge_ranked_genes <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(leading_edge_ranked_genes) <- c("leading_edge_genes_ranking_score", "leading_edge_genes_ranking_order")
  for (i in 1:dim(gsea_contrast_interest_up_reg_pathways)[1]) {
    leading_edge_ranked_genes <- rbind(leading_edge_ranked_genes, gene_list_interest[gsea_contrast_interest_up_reg_pathways$core_enrichment[i] %>% unlist()] %>% as.data.frame(row.names = names(gene_list_interest[gsea_contrast_interest_up_reg_pathways$core_enrichment[i] %>% unlist()])) %>% mutate(leading_edge_genes_ranking_order = c(1:length(gsea_contrast_interest_up_reg_pathways$core_enrichment[i] %>% unlist()))) %>%  dplyr::rename("leading_edge_genes_ranking_score" = ".") %>% mutate("up_pathways" = rep(gsea_contrast_interest_up_reg_pathways$Description[i], length(gsea_contrast_interest_up_reg_pathways$core_enrichment[i]))))
  }
  
  write.csv(leading_edge_ranked_genes, paste0(root_interest,drug_directory,"/leading_edge_genes_up_enriched_sets_gsea_",contrast_interest,".csv"), col.names = TRUE, row.names = TRUE)
  
}







# in vitro


# MTX vs Control
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_24h_vs_Control_24h", "MTX_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_18h_vs_Control_18h", "MTX_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_12h_vs_Control_12h", "MTX_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_6h_vs_Control_6h", "MTX_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_3h_vs_Control_3h", "MTX_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")


# OXA vs Control
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_24h_vs_Control_24h", "OXA_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_18h_vs_Control_18h", "OXA_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_12h_vs_Control_12h", "OXA_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_6h_vs_Control_6h", "OXA_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_3h_vs_Control_3h", "OXA_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")



# CIS vs Control
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_24h_vs_Control_24h", "CIS_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_18h_vs_Control_18h", "CIS_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_12h_vs_Control_12h", "CIS_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_6h_vs_Control_6h", "CIS_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")
run_gsea("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_3h_vs_Control_3h", "CIS_vs_Control", "Bulk_RNA_in_vitro/in_vitro_results/gsea_results/")




# in vivo

# MTX vs Control
run_gsea("top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3", "MTX_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")
run_gsea("top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10", "MTX_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")
run_gsea("top_table_in_vivo_for_Andrea.csv", "MTX_D22_vs_PBS_D22", "MTX_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")


# OXA vs Control
run_gsea("top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3", "OXA_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")
run_gsea("top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10", "OXA_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")
run_gsea("top_table_in_vivo_for_Andrea.csv", "OXA_D22_vs_PBS_D22", "OXA_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")


# CIS vs Control
run_gsea("top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3", "CIS_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")
run_gsea("top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10", "CIS_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")
run_gsea("top_table_in_vivo_for_Andrea.csv", "CIS_D22_vs_PBS_D22", "CIS_vs_Control", "Bulk_RNA_in_vivo/in_vivo_results/gsea_results/")








