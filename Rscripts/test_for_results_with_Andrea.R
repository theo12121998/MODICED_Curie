


## Enrichment


## MTX


filename_1 <- "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_MTX(24h) vs PBS(24h)_all_diff_genes.csv"



enriched_terms_MTX_24h <- read.csv(filename_1)
head(enriched_terms_MTX_24h)

enriched_terms_MTX_24h_pull <- enriched_terms_MTX_24h %>% dplyr::pull(Description) 
write.csv(enriched_terms_MTX_24h_pull, "~/MODICED_link/enriched_terms_MTX_24h.csv")

# OXA

filename_2 <- "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_OXA(24h) vs PBS(24h)_all_diff_genes.csv" 
  
enriched_terms_OXA_24h <- read.csv(filename_2)
head(enriched_terms_OXA_24h)

enriched_terms_OXA_24h_pull <- enriched_terms_OXA_24h %>% dplyr::pull(Description) 
write.csv(enriched_terms_OXA_24h_pull, "~/MODICED_link/enriched_terms_OXA_24h.csv")



# CIS

filename_3 <- "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_CIS(24h) vs PBS(24h)_all_diff_genes.csv" 

enriched_terms_CIS_24h <- read.csv(filename_3)
head(enriched_terms_CIS_24h)

enriched_terms_CIS_24h_pull <- enriched_terms_CIS_24h %>% dplyr::pull(Description) 
write.csv(enriched_terms_CIS_24h_pull, "~/MODICED_link/enriched_terms_CIS_24h.csv")






## DEA




# In vitro

counts_in_vitro_df <- readRDS("~/MODICED_link/Bulk_RNA_in_vitro/preprocessing_in_vitro_RNA/counts_in_vitro_for_analysis.rds") %>% as.data.frame()
counts_in_vitro_df_new <- counts_in_vitro_df %>% rename("genes" = "gene_name")


  
top_table_in_vitro <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv")


MTX_24h_in_vitro_all_diff <- top_table_in_vitro %>% filter(contrast == "MTX_24h_vs_Control_24h") %>% filter(adj.P.Val <= 0.05 & abs(logFC) > 2) %>% mutate(group = if_else(logFC > 0, "up","down"))

OXA_24h_in_vitro_all_diff <- top_table_in_vitro %>% filter(contrast == "OXA_24h_vs_Control_24h") %>% filter(adj.P.Val <= 0.05 & abs(logFC) > 2) %>% mutate(group = if_else(logFC > 0, "up","down"))

CIS_24h_in_vitro_all_diff <- top_table_in_vitro %>% filter(contrast == "CIS_24h_vs_Control_24h") %>% filter(adj.P.Val <= 0.05 & abs(logFC) > 2) %>% mutate(group = if_else(logFC > 0, "up","down"))








library(ggvenn)
library(RColorBrewer)
list_all_diff <- list(CIS = CIS_24h_in_vitro_all_diff$genes, OXA = OXA_24h_in_vitro_all_diff$genes, MTX = MTX_24h_in_vitro_all_diff$genes)

venn <- ggvenn(list_all_diff, fill_color = c("green2","pink2","blue2"), show_percentage = FALSE, stroke_alpha = 0, stroke_size = 0) 
ggsave("~/MODICED_link/venn_diagram_in_vitro_24h.pdf", venn)
#+ geom_segment(x = 0, y = 0, xend = 0, yend = 1, lineend = "round", linejoin = "round", size = 1, arrow = arrow(length = unit(0.3, "inches")), colour = "#EC7014" )
                                                                                                                            



a <- inner_join(CIS_24h_in_vitro_all_diff, counts_in_vitro_df_new[,c("gene_type","genes")], by="genes")

b <- inner_join(MTX_24h_in_vitro_all_diff, counts_in_vitro_df_new[,c("gene_type","genes")], by="genes")

c <- inner_join(OXA_24h_in_vitro_all_diff, counts_in_vitro_df_new[,c("gene_type","genes")], by="genes")


table(b$gene_type)
table(c$gene_type)



#good_answer <- inner_join(dlist_filtered$genes, counts_in_vitro_df_new[,c("genes","gene_type")], relationship = "many-to-many")

table(good_answer$gene_type)

