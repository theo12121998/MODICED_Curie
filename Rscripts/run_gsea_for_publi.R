#!/bioinfo/local/build/Centos/R/R-4.2.1/bin/Rscript


# GSEA for publication


run_gsea_go <- function(top_table_filename, contrast_interest, setting, directory_interest){
  
  library(clusterProfiler)
  library(tidyverse)
  top_table_contrast_interest <- read.csv(top_table_filename) %>% filter(contrast == contrast_interest) %>% mutate(rank = sign(logFC) * -log10(adj.P.Val)) %>% arrange(desc(rank))
  gene_list_interest <- top_table_contrast_interest$rank
  names(gene_list_interest) <- top_table_contrast_interest$genes
  gsea_contrast_interest <- gseGO(gene_list_interest, ont = "ALL", OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", eps=0)
  for (i in 1:50){
    
    gseaplot(gsea_contrast_interest, geneSetID = i, by="runningScore", title = gsea_contrast_interest$Description[i]) + theme(plot.title = element_text(face = "bold", hjust=0.5)) + theme(axis.title.x = element_text(face = "")) + theme(axis.title.x = element_text(size = 15, face="italic")) + theme(axis.title.y = element_text(size=15, face="italic"))
    
    
    ggsave(paste0("figures_for_publication","/",setting,"/","GSEA/", directory_interest,"/gsea_plot_",i,".png"))
    
    
  }
}






# In vitro 

# MTX

run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_24h_vs_Control_24h","In_vitro", "MTX vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_18h_vs_Control_18h","In_vitro", "MTX vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_12h_vs_Control_12h","In_vitro", "MTX vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_6h_vs_Control_6h","In_vitro", "MTX vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_3h_vs_Control_3h","In_vitro", "MTX vs Control")



# OXA

run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_24h_vs_Control_24h","In_vitro", "OXA vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_18h_vs_Control_18h","In_vitro", "OXA vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_12h_vs_Control_12h","In_vitro", "OXA vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_6h_vs_Control_6h","In_vitro", "OXA vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_3h_vs_Control_3h","In_vitro", "OXA vs Control")



# MTX

run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_24h_vs_Control_24h","In_vitro", "CIS vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_18h_vs_Control_18h","In_vitro", "CIS vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_12h_vs_Control_12h","In_vitro", "CIS vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_6h_vs_Control_6h","In_vitro", "CIS vs Control")
run_gsea_go("Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_3h_vs_Control_3h","In_vitro", "CIS vs Control")






# In vivo



# MTX

run_gsea_go("top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3","In_vivo", "MTX vs Control")
run_gsea_go("top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10","In_vivo", "MTX vs Control")



# OXA

run_gsea_go("top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3","In_vivo", "OXA vs Control")
run_gsea_go("top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10","In_vivo", "OXA vs Control")



# CIS

run_gsea_go("top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3","In_vivo", "CIS vs Control")
run_gsea_go("top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10","In_vivo", "CIS vs Control")

