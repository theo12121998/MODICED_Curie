
# Load libraries of interest
library(readODS)
library(tidyverse)
library(edgeR)
library(limma)
library(htmlwidgets)
library(plotly)



# Data wrangling (Import Gautier's list of genes)
not_all_na <- function(x) any(!is.na(x))
genes_ICD_core_model <- readODS::read_ods("ListGenes.ods")  %>% dplyr::select(where(not_all_na)) %>% dplyr::select(-"Node")


genes_for_dea <- c(genes_ICD_core_model$Genes, genes_ICD_core_model[,2])
for (i in 3:ncol(genes_ICD_core_model)){
  genes_for_dea <- c(genes_for_dea, genes_ICD_core_model[,i])
  
}
genes_for_dea_mouse <- genes_for_dea %>% str_to_title()
genes_for_dea_mouse <- genes_for_dea_mouse[!is.na(genes_for_dea_mouse)] %>% as.data.frame() %>% dplyr::rename("genes" = ".")



# Get counts data for Gautier's genes
sample_plan_in_vitro <- readRDS("~/MODICED_link/Bulk_RNA_in_vitro/preprocessing_in_vitro_RNA/sample_plan_in_vitro_for_analysis.rds")
attach(sample_plan_in_vitro)
counts_in_vitro_df <- readRDS("~/MODICED_link/Bulk_RNA_in_vitro/preprocessing_in_vitro_RNA/counts_in_vitro_for_analysis.rds") %>% as.data.frame()
counts_in_vitro_df_genes <- counts_in_vitro_df %>% dplyr::select(1:3)
counts_in_vitro_df_samples <- counts_in_vitro_df %>% dplyr::select(-c(1:3))
dlist <- DGEList(counts = counts_in_vitro_df_samples, samples = sample_plan_in_vitro$Sample_Name, genes = counts_in_vitro_df_genes$gene_name)
dsgn <- model.matrix(~ 0 + Treatment_hour, data=dlist$samples)
colnames(dsgn) <- sub("Treatment_hour","",colnames(dsgn))
keep <- filterByExpr(dlist$counts, dsgn)
dlist_filtered <- dlist[keep,]
dlist_filtered <- calcNormFactors(dlist_filtered, method="TMM")
voom <- voomWithQualityWeights(dlist_filtered, dsgn)
counts_for_gautier <- cbind(voom$E, voom$genes)
rownames(counts_for_gautier) <- make.names(counts_for_gautier$genes, unique=TRUE)
counts_for_gautier <- counts_for_gautier %>% dplyr::select(-"genes")
fit_dea <- lmFit(counts_for_gautier, dsgn)





# Get contrasts
contrasts_dea_3h <- makeContrasts(MP_3h = MTX_3h - Control_3h, OP_3h = OXA_3h - Control_3h, CP_3h = CIS_3h - Control_3h, CO_3h = CIS_3h - OXA_3h, CM_3h = CIS_3h - MTX_3h, levels = dsgn)
fit_contrast_3h <- contrasts.fit(fit_dea, contrasts_dea_3h)
fit_contrast_3h <- eBayes(fit_contrast_3h)

contrasts_dea_6h <- makeContrasts(MP_6h = MTX_6h - Control_6h, OP_6h = OXA_6h - Control_6h, CP_6h = CIS_6h - Control_6h, CO_6h = CIS_6h - OXA_6h, CM_6h = CIS_6h - MTX_6h, levels = dsgn)
fit_contrast_6h <- contrasts.fit(fit_dea, contrasts_dea_6h)
fit_contrast_6h <- eBayes(fit_contrast_6h)

contrasts_dea_12h <- makeContrasts(MP_12h = MTX_12h - Control_12h, OP_12h = OXA_12h - Control_12h, CP_12h = CIS_12h - Control_12h, CO_12h = CIS_12h - OXA_12h, CM_12h = CIS_12h - MTX_12h, levels = dsgn)
fit_contrast_12h <- contrasts.fit(fit_dea, contrasts_dea_12h)
fit_contrast_12h <- eBayes(fit_contrast_12h)

contrasts_dea_18h <- makeContrasts(MP_18h = MTX_18h - Control_18h, OP_18h = OXA_18h - Control_18h, CP_18h = CIS_18h - Control_18h, CO_18h = CIS_18h - OXA_18h, CM_18h = CIS_18h - MTX_18h, levels = dsgn)
fit_contrast_18h <- contrasts.fit(fit_dea, contrasts_dea_18h)
fit_contrast_18h <- eBayes(fit_contrast_18h)

contrasts_dea_24h <- makeContrasts(MP_24h = MTX_24h - Control_24h, OP_24h = OXA_24h - Control_24h, CP_24h = CIS_24h - Control_24h, CO_24h = CIS_24h - OXA_24h, CM_24h = CIS_24h - MTX_24h, levels = dsgn)
fit_contrast_24h <- contrasts.fit(fit_dea, contrasts_dea_24h)
fit_contrast_24h <- eBayes(fit_contrast_24h)







# Get results dea


get_results_dea <- function(fit_contrast, coef_interest, all_genes_Gautier_filename, genes_Gautier_up_reg_filename, genes_Gautier_down_reg_filename, all_diff_genes_Gautier_concat_filename, all_diff_genes_filename, label, volcano_plot_filename, volcano_plot_widget_filename){
  
  
  
  # Top table for Gautier's genes (+ Save it)
  all_genes_Gautier <- topTable(fit_contrast, sort.by= "P", adjust.method = "BH", number = Inf, coef = coef_interest) %>% rownames_to_column("gene_name") %>% mutate(Gautier_genes = if_else(gene_name %in% genes_for_dea_mouse$genes, "Yes","No")) %>% filter(Gautier_genes == "Yes") %>% dplyr::select(-"Gautier_genes")
 

  write.csv(all_genes_Gautier, all_genes_Gautier_filename, col.names = TRUE, row.names = FALSE)

  
  # Get up-reg Gautier's genes
  diff_genes_Gautier_up <- all_genes_Gautier %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% dplyr::arrange(adj.P.Val) %>% 
    dplyr::mutate(group = "up")
  
  write.csv(diff_genes_Gautier_up, genes_Gautier_up_reg_filename, col.names = TRUE, row.names = FALSE)
  
  
  # Get down-reg Gautier's genes
  diff_genes_Gautier_down <- all_genes_Gautier %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% dplyr::arrange(adj.P.Val) %>% 
    dplyr::mutate(group = "down")
  
  write.csv(diff_genes_Gautier_down, genes_Gautier_down_reg_filename, col.names = TRUE, row.names = FALSE)
  
  
  
  # Get up and down-reg Gautier's genes (concatenated)
  all_diff_genes_Gautier_concat <- rbind(diff_genes_Gautier_up, diff_genes_Gautier_down) 
  write.csv(all_diff_genes_Gautier_concat, all_diff_genes_Gautier_concat_filename, col.names = TRUE, row.names = FALSE)
  
  
  # Get all diff genes for later (upstream)
  all_diff_genes <- topTable(fit_contrast, sort.by= "P", adjust.method = "BH", number = Inf, coef = coef_interest) %>% rownames_to_column("gene_name") %>% mutate(Gautier_genes = if_else(gene_name %in% genes_for_dea_mouse$genes, "Yes","No")) %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0)
  write.csv(all_diff_genes, all_diff_genes_filename)
  
  
  
  # Volcano Plot - Colored by logFC/adj p-value (for Gautier's genes)
  volcano_plot <- all_genes_Gautier %>%  
    ggplot(aes(x=logFC, y=-log10(adj.P.Val), label=gene_name)) +
    geom_point(data = filter(all_genes_Gautier, logFC > 1 & adj.P.Val < 0.05), color="red", size=1.5) +
    geom_point(data = filter(all_genes_Gautier, logFC < 1 & adj.P.Val < 0.05), color="blue", size=1.5) +
    geom_point(data = filter(all_genes_Gautier, adj.P.Val < 0.05 & (logFC > 0 & logFC < 1)), color="orange", size=1.5) +
    geom_point(data = filter(all_genes_Gautier, adj.P.Val < 0.05 & (logFC < 0 & logFC >  -1)), color="lightblue", size=1.5) +
    geom_point(data = filter(all_genes_Gautier, adj.P.Val >= 0.05 & (logFC > -1 & logFC < 1)), color="darkgrey", size=1.5) +
    ggrepel::geom_label_repel(data = filter(all_genes_Gautier, logFC > 1 & adj.P.Val < 0.05), aes(label = gene_name), label.size = 0.01, label.padding = 0.04, max.overlaps = Inf) +
    ggrepel::geom_label_repel(data = filter(all_genes_Gautier, logFC < -1 & adj.P.Val < 0.05), aes(label = gene_name), label.size = 0.01, label.padding = 0.04, max.overlaps = Inf) +
    geom_vline(xintercept = c(-1,1), colour = "black", linetype = "dashed") +
    geom_hline(yintercept = (-log10(0.05)), colour = "black", linetype = "dashed") +
    ylab("-log10(adj.P.value)") +
    xlab("log2(FoldChange)") +
    theme_linedraw() +
    ggtitle(label) +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  ggsave(volcano_plot_filename, width = 9, height = 5)
  
  
  
  # Plotly (+ Save object)
  volcano_plot_widget <- plotly::ggplotly(volcano_plot)
  htmlwidgets::saveWidget(volcano_plot_widget, volcano_plot_widget_filename)


}





# Run and save results for DEA



# MTX vs PBS


# 3h

get_results_dea(fit_contrast_3h, "MP_3h", "Analysis_for_Gautier/DEA/MTX_vs_Control/3h/tables_csv/top_table_Gautier_genes_MTX_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/3h/tables_csv/up_reg_Gautier_genes_MTX_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/3h/tables_csv/down_reg_Gautier_genes_Gautier_MTX_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/3h/tables_csv/up_and_down_reg_Gautier_genes_concat_MTX_3h_vs_Control_3h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_3h_vs_Control_3h.csv","Volcano plot: MTX_3h_vs_Control_3h", "Analysis_for_Gautier/DEA/MTX_vs_Control/3h/Volcano_plots/Volcano_plot_MTX_3h_vs_Control_3h.pdf", "Analysis_for_Gautier/DEA/MTX_vs_Control/3h/Volcano_plots/Volcano_plot_MTX_3h_vs_Control_3h_widget.html")




# 6h

get_results_dea(fit_contrast_6h, "MP_6h", "Analysis_for_Gautier/DEA/MTX_vs_Control/6h/tables_csv/top_table_Gautier_genes_MTX_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/6h/tables_csv/up_reg_Gautier_genes_MTX_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/6h/tables_csv/down_reg_Gautier_genes_Gautier_MTX_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/6h/tables_csv/up_and_down_reg_Gautier_genes_concat_MTX_6h_vs_Control_6h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_6h_vs_Control_6h.csv","Volcano plot: MTX_6h_vs_Control_6h", "Analysis_for_Gautier/DEA/MTX_vs_Control/6h/Volcano_plots/Volcano_plot_MTX_6h_vs_Control_6h.pdf", "Analysis_for_Gautier/DEA/MTX_vs_Control/6h/Volcano_plots/Volcano_plot_MTX_6h_vs_Control_6h_widget.html")





# 12h

get_results_dea(fit_contrast_12h, "MP_12h", "Analysis_for_Gautier/DEA/MTX_vs_Control/12h/tables_csv/top_table_Gautier_genes_MTX_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/12h/tables_csv/up_reg_Gautier_genes_MTX_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/12h/tables_csv/down_reg_Gautier_genes_Gautier_MTX_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/12h/tables_csv/up_and_down_reg_Gautier_genes_concat_MTX_12h_vs_Control_12h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_12h_vs_Control_12h.csv","Volcano plot: MTX_12h_vs_Control_12h", "Analysis_for_Gautier/DEA/MTX_vs_Control/12h/Volcano_plots/Volcano_plot_MTX_12h_vs_Control_12h.pdf", "Analysis_for_Gautier/DEA/MTX_vs_Control/12h/Volcano_plots/Volcano_plot_MTX_12h_vs_Control_12h_widget.html")






# 18h

get_results_dea(fit_contrast_18h, "MP_18h", "Analysis_for_Gautier/DEA/MTX_vs_Control/18h/tables_csv/top_table_Gautier_genes_MTX_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/18h/tables_csv/up_reg_Gautier_genes_MTX_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/18h/tables_csv/down_reg_Gautier_genes_Gautier_MTX_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/18h/tables_csv/up_and_down_reg_Gautier_genes_concat_MTX_18h_vs_Control_18h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_18h_vs_Control_18h.csv","Volcano plot: MTX_18h_vs_Control_18h", "Analysis_for_Gautier/DEA/MTX_vs_Control/18h/Volcano_plots/Volcano_plot_MTX_18h_vs_Control_18h.pdf", "Analysis_for_Gautier/DEA/MTX_vs_Control/18h/Volcano_plots/Volcano_plot_MTX_18h_vs_Control_18h_widget.html")




# 24h

get_results_dea(fit_contrast_24h, "MP_24h", "Analysis_for_Gautier/DEA/MTX_vs_Control/24h/tables_csv/top_table_Gautier_genes_MTX_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/24h/tables_csv/up_reg_Gautier_genes_MTX_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/24h/tables_csv/down_reg_Gautier_genes_Gautier_MTX_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/MTX_vs_Control/24h/tables_csv/up_and_down_reg_Gautier_genes_concat_MTX_24h_vs_Control_24h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_24h_vs_Control_24h.csv","Volcano plot: MTX_24h_vs_Control_24h", "Analysis_for_Gautier/DEA/MTX_vs_Control/24h/Volcano_plots/Volcano_plot_MTX_24h_vs_Control_24h.pdf", "Analysis_for_Gautier/DEA/MTX_vs_Control/24h/Volcano_plots/Volcano_plot_MTX_24h_vs_Control_24h_widget.html")

  
  




# OXA vs PBS


# 3h

get_results_dea(fit_contrast_3h, "OP_3h", "Analysis_for_Gautier/DEA/OXA_vs_Control/3h/tables_csv/top_table_Gautier_genes_OXA_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/3h/tables_csv/up_reg_Gautier_genes_OXA_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/3h/tables_csv/down_reg_Gautier_genes_Gautier_OXA_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/3h/tables_csv/up_and_down_reg_Gautier_genes_concat_OXA_3h_vs_Control_3h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_3h_vs_Control_3h.csv","Volcano plot: OXA_3h_vs_Control_3h", "Analysis_for_Gautier/DEA/OXA_vs_Control/3h/Volcano_plots/Volcano_plot_OXA_3h_vs_Control_3h.pdf", "Analysis_for_Gautier/DEA/OXA_vs_Control/3h/Volcano_plots/Volcano_plot_OXA_3h_vs_Control_3h_widget.html")



# 6h

get_results_dea(fit_contrast_6h, "OP_6h", "Analysis_for_Gautier/DEA/OXA_vs_Control/6h/tables_csv/top_table_Gautier_genes_OXA_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/6h/tables_csv/up_reg_Gautier_genes_OXA_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/6h/tables_csv/down_reg_Gautier_genes_Gautier_OXA_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/6h/tables_csv/up_and_down_reg_Gautier_genes_concat_OXA_6h_vs_Control_6h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_6h_vs_Control_6h.csv","Volcano plot: OXA_6h_vs_Control_6h", "Analysis_for_Gautier/DEA/OXA_vs_Control/6h/Volcano_plots/Volcano_plot_OXA_6h_vs_Control_6h.pdf", "Analysis_for_Gautier/DEA/OXA_vs_Control/6h/Volcano_plots/Volcano_plot_OXA_6h_vs_Control_6h_widget.html")





# 12h

get_results_dea(fit_contrast_12h, "OP_12h", "Analysis_for_Gautier/DEA/OXA_vs_Control/12h/tables_csv/top_table_Gautier_genes_OXA_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/12h/tables_csv/up_reg_Gautier_genes_OXA_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/12h/tables_csv/down_reg_Gautier_genes_Gautier_OXA_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/12h/tables_csv/up_and_down_reg_Gautier_genes_concat_OXA_12h_vs_Control_12h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_12h_vs_Control_12h.csv","Volcano plot: OXA_12h_vs_Control_12h", "Analysis_for_Gautier/DEA/OXA_vs_Control/12h/Volcano_plots/Volcano_plot_OXA_12h_vs_Control_12h.pdf", "Analysis_for_Gautier/DEA/OXA_vs_Control/12h/Volcano_plots/Volcano_plot_OXA_12h_vs_Control_12h_widget.html")






# 18h

get_results_dea(fit_contrast_18h, "OP_18h", "Analysis_for_Gautier/DEA/OXA_vs_Control/18h/tables_csv/top_table_Gautier_genes_OXA_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/18h/tables_csv/up_reg_Gautier_genes_OXA_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/18h/tables_csv/down_reg_Gautier_genes_Gautier_OXA_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/18h/tables_csv/up_and_down_reg_Gautier_genes_concat_OXA_18h_vs_Control_18h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_18h_vs_Control_18h.csv","Volcano plot: OXA_18h_vs_Control_18h", "Analysis_for_Gautier/DEA/OXA_vs_Control/18h/Volcano_plots/Volcano_plot_OXA_18h_vs_Control_18h.pdf", "Analysis_for_Gautier/DEA/OXA_vs_Control/18h/Volcano_plots/Volcano_plot_OXA_18h_vs_Control_18h_widget.html")




# 24h

get_results_dea(fit_contrast_24h, "OP_24h", "Analysis_for_Gautier/DEA/OXA_vs_Control/24h/tables_csv/top_table_Gautier_genes_OXA_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/24h/tables_csv/up_reg_Gautier_genes_OXA_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/24h/tables_csv/down_reg_Gautier_genes_Gautier_OXA_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/OXA_vs_Control/24h/tables_csv/up_and_down_reg_Gautier_genes_concat_OXA_24h_vs_Control_24h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_24h_vs_Control_24h.csv","Volcano plot: OXA_24h_vs_Control_24h", "Analysis_for_Gautier/DEA/OXA_vs_Control/24h/Volcano_plots/Volcano_plot_OXA_24h_vs_Control_24h.pdf", "Analysis_for_Gautier/DEA/OXA_vs_Control/24h/Volcano_plots/Volcano_plot_OXA_24h_vs_Control_24h_widget.html")







# CIS vs PBS



# 3h

get_results_dea(fit_contrast_3h, "CP_3h", "Analysis_for_Gautier/DEA/CIS_vs_Control/3h/tables_csv/top_table_Gautier_genes_CIS_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/3h/tables_csv/up_reg_Gautier_genes_CIS_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/3h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_3h_vs_Control_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/3h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_3h_vs_Control_3h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_3h_vs_Control_3h.csv","Volcano plot: CIS_3h_vs_Control_3h", "Analysis_for_Gautier/DEA/CIS_vs_Control/3h/Volcano_plots/Volcano_plot_CIS_3h_vs_Control_3h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_Control/3h/Volcano_plots/Volcano_plot_CIS_3h_vs_Control_3h_widget.html")



# 6h

get_results_dea(fit_contrast_6h, "CP_6h", "Analysis_for_Gautier/DEA/CIS_vs_Control/6h/tables_csv/top_table_Gautier_genes_CIS_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/6h/tables_csv/up_reg_Gautier_genes_CIS_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/6h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_6h_vs_Control_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/6h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_6h_vs_Control_6h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_6h_vs_Control_6h.csv","Volcano plot: CIS_6h_vs_Control_6h", "Analysis_for_Gautier/DEA/CIS_vs_Control/6h/Volcano_plots/Volcano_plot_CIS_6h_vs_Control_6h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_Control/6h/Volcano_plots/Volcano_plot_CIS_6h_vs_Control_6h_widget.html")




# 12h

get_results_dea(fit_contrast_12h, "CP_12h", "Analysis_for_Gautier/DEA/CIS_vs_Control/12h/tables_csv/top_table_Gautier_genes_CIS_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/12h/tables_csv/up_reg_Gautier_genes_CIS_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/12h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_12h_vs_Control_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/12h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_12h_vs_Control_12h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_12h_vs_Control_12h.csv","Volcano plot: CIS_12h_vs_Control_12h", "Analysis_for_Gautier/DEA/CIS_vs_Control/12h/Volcano_plots/Volcano_plot_CIS_12h_vs_Control_12h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_Control/12h/Volcano_plots/Volcano_plot_CIS_12h_vs_Control_12h_widget.html")




# 18h

get_results_dea(fit_contrast_18h, "CP_18h", "Analysis_for_Gautier/DEA/CIS_vs_Control/18h/tables_csv/top_table_Gautier_genes_CIS_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/18h/tables_csv/up_reg_Gautier_genes_CIS_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/18h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_18h_vs_Control_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/18h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_18h_vs_Control_18h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_18h_vs_Control_18h.csv","Volcano plot: CIS_18h_vs_Control_18h", "Analysis_for_Gautier/DEA/CIS_vs_Control/18h/Volcano_plots/Volcano_plot_CIS_18h_vs_Control_18h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_Control/18h/Volcano_plots/Volcano_plot_CIS_18h_vs_Control_18h_widget.html")




# 24h

get_results_dea(fit_contrast_24h, "CP_24h", "Analysis_for_Gautier/DEA/CIS_vs_Control/24h/tables_csv/top_table_Gautier_genes_CIS_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/24h/tables_csv/up_reg_Gautier_genes_CIS_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/24h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_24h_vs_Control_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_Control/24h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_24h_vs_Control_24h.csv", "t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_24h_vs_Control_24h.csv","Volcano plot: CIS_24h_vs_Control_24h", "Analysis_for_Gautier/DEA/CIS_vs_Control/24h/Volcano_plots/Volcano_plot_CIS_24h_vs_Control_24h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_Control/24h/Volcano_plots/Volcano_plot_CIS_24h_vs_Control_24h_widget.html")







# CIS vs MTX



# 3h

get_results_dea(fit_contrast_3h, "CM_3h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/tables_csv/top_table_Gautier_genes_CIS_3h_vs_MTX_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/tables_csv/up_reg_Gautier_genes_CIS_3h_vs_MTX_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_3h_vs_MTX_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_3h_vs_MTX_3h.csv", "Volcano plot: CIS_3h_vs_MTX_3h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/Volcano_plots/Volcano_plot_CIS_3h_vs_MTX_3h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/Volcano_plots/Volcano_plot_CIS_3h_vs_MTX_3h_widget.html")



# 6h

get_results_dea(fit_contrast_6h, "CM_6h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/tables_csv/top_table_Gautier_genes_CIS_6h_vs_MTX_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/tables_csv/up_reg_Gautier_genes_CIS_6h_vs_MTX_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_6h_vs_MTX_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_6h_vs_MTX_6h.csv", "Volcano plot: CIS_6h_vs_MTX_6h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/Volcano_plots/Volcano_plot_CIS_6h_vs_MTX_6h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/Volcano_plots/Volcano_plot_CIS_6h_vs_MTX_6h_widget.html")




# 12h

get_results_dea(fit_contrast_12h, "CM_12h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/tables_csv/top_table_Gautier_genes_CIS_12h_vs_MTX_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/tables_csv/up_reg_Gautier_genes_CIS_12h_vs_MTX_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_12h_vs_MTX_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_12h_vs_MTX_12h.csv", "Volcano plot: CIS_12h_vs_MTX_12h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/Volcano_plots/Volcano_plot_CIS_12h_vs_MTX_12h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/Volcano_plots/Volcano_plot_CIS_12h_vs_MTX_12h_widget.html")




# 18h

get_results_dea(fit_contrast_18h, "CM_18h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/tables_csv/top_table_Gautier_genes_CIS_18h_vs_MTX_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/tables_csv/up_reg_Gautier_genes_CIS_18h_vs_MTX_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_18h_vs_MTX_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_18h_vs_MTX_18h.csv", "Volcano plot: CIS_18h_vs_MTX_18h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/Volcano_plots/Volcano_plot_CIS_18h_vs_MTX_18h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/Volcano_plots/Volcano_plot_CIS_18h_vs_MTX_18h_widget.html")



# 24h

get_results_dea(fit_contrast_24h, "CM_24h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/tables_csv/top_table_Gautier_genes_CIS_24h_vs_MTX_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/tables_csv/up_reg_Gautier_genes_CIS_24h_vs_MTX_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_24h_vs_MTX_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_24h_vs_MTX_24h.csv", "Volcano plot: CIS_24h_vs_MTX_24h", "Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/Volcano_plots/Volcano_plot_CIS_24h_vs_MTX_24h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/Volcano_plots/Volcano_plot_CIS_24h_vs_MTX_24h_widget.html")








# CIS vs OXA



# 3h

get_results_dea(fit_contrast_3h, "CO_3h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/tables_csv/top_table_Gautier_genes_CIS_3h_vs_OXA_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/tables_csv/up_reg_Gautier_genes_CIS_3h_vs_OXA_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_3h_vs_OXA_3h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_3h_vs_OXA_3h.csv", "Volcano plot: CIS_3h_vs_OXA_3h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/Volcano_plots/Volcano_plot_CIS_3h_vs_OXA_3h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/Volcano_plots/Volcano_plot_CIS_3h_vs_OXA_3h_widget.html")



# 6h

get_results_dea(fit_contrast_6h, "CO_6h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/tables_csv/top_table_Gautier_genes_CIS_6h_vs_OXA_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/tables_csv/up_reg_Gautier_genes_CIS_6h_vs_OXA_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_6h_vs_OXA_6h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_6h_vs_OXA_6h.csv", "Volcano plot: CIS_6h_vs_OXA_6h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/Volcano_plots/Volcano_plot_CIS_6h_vs_OXA_6h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/Volcano_plots/Volcano_plot_CIS_6h_vs_OXA_6h_widget.html")




# 12h

get_results_dea(fit_contrast_12h, "CO_12h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/tables_csv/top_table_Gautier_genes_CIS_12h_vs_OXA_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/tables_csv/up_reg_Gautier_genes_CIS_12h_vs_OXA_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_12h_vs_OXA_12h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_12h_vs_OXA_12h.csv", "Volcano plot: CIS_12h_vs_OXA_12h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/Volcano_plots/Volcano_plot_CIS_12h_vs_OXA_12h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/Volcano_plots/Volcano_plot_CIS_12h_vs_OXA_12h_widget.html")





# 18h

get_results_dea(fit_contrast_18h, "CO_18h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/tables_csv/top_table_Gautier_genes_CIS_18h_vs_OXA_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/tables_csv/up_reg_Gautier_genes_CIS_18h_vs_OXA_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_18h_vs_OXA_18h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_18h_vs_OXA_18h.csv", "Volcano plot: CIS_18h_vs_OXA_18h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/Volcano_plots/Volcano_plot_CIS_18h_vs_OXA_18h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/Volcano_plots/Volcano_plot_CIS_18h_vs_OXA_18h_widget.html")



# 24h

get_results_dea(fit_contrast_24h, "CO_24h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/tables_csv/top_table_Gautier_genes_CIS_24h_vs_OXA_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/tables_csv/up_reg_Gautier_genes_CIS_24h_vs_OXA_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_24h_vs_OXA_24h.csv", "Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/tables_csv/up_and_down_reg_Gautier_genes_concat_CIS_24h_vs_OXA_24h.csv", "Volcano plot: CIS_24h_vs_OXA_24h", "Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/Volcano_plots/Volcano_plot_CIS_24h_vs_OXA_24h.pdf", "Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/Volcano_plots/Volcano_plot_CIS_24h_vs_OXA_24h_widget.html")








gautier_genes_up_reg_CIS_MTX_24h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/tables_csv/up_reg_Gautier_genes_CIS_24h_vs_MTX_24h.csv")
gautier_genes_up_reg_CIS_OXA_24h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/tables_csv/up_reg_Gautier_genes_CIS_24h_vs_OXA_24h.csv")



gautier_genes_up_reg_CIS_MTX_24h_FC_2 <- gautier_genes_up_reg_CIS_MTX_24h[which(gautier_genes_up_reg_CIS_MTX_24h$logFC > 1 & gautier_genes_up_reg_CIS_MTX_24h$adj.P.Val <= 0.05),]


gautier_genes_up_reg_CIS_OXA_24h_FC_2 <- gautier_genes_up_reg_CIS_OXA_24h[which(gautier_genes_up_reg_CIS_OXA_24h$logFC > 1 & gautier_genes_up_reg_CIS_OXA_24h$adj.P.Val <= 0.05),]



intersect(gautier_genes_up_reg_CIS_MTX_24h_FC_2$gene_name, gautier_genes_up_reg_CIS_OXA_24h_FC_2$gene_name)



#bv <- list(gautier_genes_up_reg_CIS_MTX_24h_FC_2$gene_name, gautier_genes_up_reg_CIS_OXA_24h_FC_2$gene_name)
#Øvenn.diagram(bv)





gautier_genes_up_reg_CIS_MTX_18h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/tables_csv/up_reg_Gautier_genes_CIS_18h_vs_MTX_18h.csv")
gautier_genes_up_reg_CIS_OXA_18h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/tables_csv/up_reg_Gautier_genes_CIS_18h_vs_OXA_18h.csv")



gautier_genes_up_reg_CIS_MTX_18h_FC_2 <- gautier_genes_up_reg_CIS_MTX_18h[which(gautier_genes_up_reg_CIS_MTX_18h$logFC > 1 & gautier_genes_up_reg_CIS_MTX_18h$adj.P.Val <= 0.05),]


gautier_genes_up_reg_CIS_OXA_18h_FC_2 <- gautier_genes_up_reg_CIS_OXA_18h[which(gautier_genes_up_reg_CIS_OXA_18h$logFC > 1 & gautier_genes_up_reg_CIS_OXA_18h$adj.P.Val <= 0.05),]


intersect(gautier_genes_up_reg_CIS_MTX_18h_FC_2$gene_name, gautier_genes_up_reg_CIS_OXA_18h_FC_2$gene_name)



# 12h

gautier_genes_up_reg_CIS_MTX_12h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/tables_csv/up_reg_Gautier_genes_CIS_12h_vs_MTX_12h.csv")
gautier_genes_up_reg_CIS_OXA_12h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/tables_csv/up_reg_Gautier_genes_CIS_12h_vs_OXA_12h.csv")



gautier_genes_up_reg_CIS_MTX_12h_FC_2 <- gautier_genes_up_reg_CIS_MTX_12h[which(gautier_genes_up_reg_CIS_MTX_12h$logFC > 1 & gautier_genes_up_reg_CIS_MTX_12h$adj.P.Val <= 0.05),]


gautier_genes_up_reg_CIS_OXA_12h_FC_2 <- gautier_genes_up_reg_CIS_OXA_12h[which(gautier_genes_up_reg_CIS_OXA_12h$logFC > 1 & gautier_genes_up_reg_CIS_OXA_12h$adj.P.Val <= 0.05),]


intersect(gautier_genes_up_reg_CIS_MTX_12h_FC_2$gene_name, gautier_genes_up_reg_CIS_OXA_12h_FC_2$gene_name)



# 6h

gautier_genes_up_reg_CIS_MTX_6h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/tables_csv/up_reg_Gautier_genes_CIS_6h_vs_MTX_6h.csv")
gautier_genes_up_reg_CIS_OXA_6h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/tables_csv/up_reg_Gautier_genes_CIS_6h_vs_OXA_6h.csv")



gautier_genes_up_reg_CIS_MTX_6h_FC_2 <- gautier_genes_up_reg_CIS_MTX_6h[which(gautier_genes_up_reg_CIS_MTX_6h$logFC > 1 & gautier_genes_up_reg_CIS_MTX_6h$adj.P.Val <= 0.05),]


gautier_genes_up_reg_CIS_OXA_6h_FC_2 <- gautier_genes_up_reg_CIS_OXA_6h[which(gautier_genes_up_reg_CIS_OXA_6h$logFC > 1 & gautier_genes_up_reg_CIS_OXA_6h$adj.P.Val <= 0.05),]


intersect(gautier_genes_up_reg_CIS_MTX_6h_FC_2$gene_name, gautier_genes_up_reg_CIS_OXA_6h_FC_2$gene_name)




# 3h

gautier_genes_down_reg_CIS_MTX_3h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/3h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_3h_vs_MTX_3h.csv")
gautier_genes_down_reg_CIS_OXA_3h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/3h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_3h_vs_OXA_3h.csv")



gautier_genes_down_reg_CIS_MTX_3h_FC_2 <- gautier_genes_down_reg_CIS_MTX_3h[which(gautier_genes_down_reg_CIS_MTX_3h$logFC < -1 & gautier_genes_down_reg_CIS_MTX_3h$adj.P.Val <= 0.05),]


gautier_genes_down_reg_CIS_OXA_3h_FC_2 <- gautier_genes_down_reg_CIS_OXA_3h[which(gautier_genes_down_reg_CIS_OXA_3h$logFC < -1 & gautier_genes_down_reg_CIS_OXA_3h$adj.P.Val <= 0.05),]


intersect(gautier_genes_down_reg_CIS_MTX_3h_FC_2$gene_name, gautier_genes_down_reg_CIS_OXA_3h_FC_2$gene_name)









# 6h

gautier_genes_down_reg_CIS_MTX_6h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/6h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_6h_vs_MTX_6h.csv")
gautier_genes_down_reg_CIS_OXA_6h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/6h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_6h_vs_OXA_6h.csv")



gautier_genes_down_reg_CIS_MTX_6h_FC_2 <- gautier_genes_down_reg_CIS_MTX_6h[which(gautier_genes_down_reg_CIS_MTX_6h$logFC < -1 & gautier_genes_down_reg_CIS_MTX_6h$adj.P.Val <= 0.05),]


gautier_genes_down_reg_CIS_OXA_6h_FC_2 <- gautier_genes_down_reg_CIS_OXA_6h[which(gautier_genes_down_reg_CIS_OXA_6h$logFC < -1 & gautier_genes_down_reg_CIS_OXA_6h$adj.P.Val <= 0.05),]


intersect(gautier_genes_down_reg_CIS_MTX_6h_FC_2$gene_name, gautier_genes_down_reg_CIS_OXA_6h_FC_2$gene_name)





# 12h

gautier_genes_down_reg_CIS_MTX_12h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/12h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_12h_vs_MTX_12h.csv")
gautier_genes_down_reg_CIS_OXA_12h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/12h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_12h_vs_OXA_12h.csv")



gautier_genes_down_reg_CIS_MTX_12h_FC_2 <- gautier_genes_down_reg_CIS_MTX_12h[which(gautier_genes_down_reg_CIS_MTX_12h$logFC < -1 & gautier_genes_down_reg_CIS_MTX_12h$adj.P.Val <= 0.05),]


gautier_genes_down_reg_CIS_OXA_12h_FC_2 <- gautier_genes_down_reg_CIS_OXA_12h[which(gautier_genes_down_reg_CIS_OXA_12h$logFC < -1 & gautier_genes_down_reg_CIS_OXA_12h$adj.P.Val <= 0.05),]


down_12 <- intersect(gautier_genes_down_reg_CIS_MTX_12h_FC_2$gene_name, gautier_genes_down_reg_CIS_OXA_12h_FC_2$gene_name)





# 18h

gautier_genes_down_reg_CIS_MTX_18h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/18h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_18h_vs_MTX_18h.csv")
gautier_genes_down_reg_CIS_OXA_18h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/18h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_18h_vs_OXA_18h.csv")



gautier_genes_down_reg_CIS_MTX_18h_FC_2 <- gautier_genes_down_reg_CIS_MTX_18h[which(gautier_genes_down_reg_CIS_MTX_18h$logFC < -1 & gautier_genes_down_reg_CIS_MTX_18h$adj.P.Val <= 0.05),]


gautier_genes_down_reg_CIS_OXA_18h_FC_2 <- gautier_genes_down_reg_CIS_OXA_18h[which(gautier_genes_down_reg_CIS_OXA_18h$logFC < -1 & gautier_genes_down_reg_CIS_OXA_18h$adj.P.Val <= 0.05),]


down_18 <- intersect(gautier_genes_down_reg_CIS_MTX_18h_FC_2$gene_name, gautier_genes_down_reg_CIS_OXA_18h_FC_2$gene_name)





# 24h

gautier_genes_down_reg_CIS_MTX_24h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_MTX/24h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_24h_vs_MTX_24h.csv")
gautier_genes_down_reg_CIS_OXA_24h <- read.csv("Analysis_for_Gautier/DEA/CIS_vs_OXA/24h/tables_csv/down_reg_Gautier_genes_Gautier_CIS_24h_vs_OXA_24h.csv")



gautier_genes_down_reg_CIS_MTX_24h_FC_2 <- gautier_genes_down_reg_CIS_MTX_24h[which(gautier_genes_down_reg_CIS_MTX_24h$logFC < -1 & gautier_genes_down_reg_CIS_MTX_24h$adj.P.Val <= 0.05),]


gautier_genes_down_reg_CIS_OXA_24h_FC_2 <- gautier_genes_down_reg_CIS_OXA_24h[which(gautier_genes_down_reg_CIS_OXA_24h$logFC < -1 & gautier_genes_down_reg_CIS_OXA_24h$adj.P.Val <= 0.05),]


down_24 <- intersect(gautier_genes_down_reg_CIS_MTX_24h_FC_2$gene_name, gautier_genes_down_reg_CIS_OXA_24h_FC_2$gene_name)








