

# Load libraries of interest
library(ComplexHeatmap)
library(rRoma)
library(tidyverse)


# Import gmt file
tf_gmt <- ReadGMTFile("/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/m3.gtrd.v2023.2.Mm.symbols.gmt")


# Get tfs binding to at least one Gautier's list of genes
tfs_target_at_least_one_gene_gautier <- c()
n <- 0
for (i in 1:length(tf_gmt)){
  if (any(tf_gmt[[i]]$Genes %in% genes_for_dea_mouse$genes)){
    print(i)
    n <- n+1
    tfs_target_at_least_one_gene_gautier <- append(tfs_target_at_least_one_gene_gautier, tf_gmt[[i]]$Name)
  }
}



# Run rRoma
run_roma_all_timepoints <- rRoma.R(counts_for_gautier, tf_gmt)


# Retrieve module and sample matrix from rRoma output
get_roma_results_module_matrix <- run_roma_all_timepoints$ModuleMatrix %>% as.data.frame()
get_roma_results_sample_matrix <- run_roma_all_timepoints$SampleMatrix %>% as.data.frame()

get_roma_results_module_matrix <- read.csv("roma_module_matrix_for_Gautier.csv") %>% dplyr::rename(modules = "X")
get_roma_results_sample_matrix <- read.csv("roma_sample_matrix_for_Gautier.csv") %>% column_to_rownames("X")
colnames(get_roma_results_sample_matrix) <- sub("X","", colnames(get_roma_results_sample_matrix))




# Get the mean activity for each TF (i.e. average sample matrix across replicates for each drug_time condition)



# 3h

get_roma_results_sample_matrix_MTX_3h <- get_roma_results_sample_matrix %>% select(contains("3h") & contains("MTX"))

get_roma_results_sample_matrix_MTX_3h_mean <- rowMeans(get_roma_results_sample_matrix_MTX_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_3h" = ".")


get_roma_results_sample_matrix_OXA_3h <- get_roma_results_sample_matrix %>% select(contains("3h") & contains("OXA"))

get_roma_results_sample_matrix_OXA_3h_mean <- rowMeans(get_roma_results_sample_matrix_OXA_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_3h" = ".")


get_roma_results_sample_matrix_CIS_3h <- get_roma_results_sample_matrix %>% select(contains("3h") & contains("CIS"))

get_roma_results_sample_matrix_CIS_3h_mean <- rowMeans(get_roma_results_sample_matrix_CIS_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_3h" = ".")


get_roma_results_sample_matrix_Control_3h <- get_roma_results_sample_matrix %>% select(contains("3h") & contains("Control"))

get_roma_results_sample_matrix_Control_3h_mean <- rowMeans(get_roma_results_sample_matrix_Control_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_3h" = ".")



get_roma_results_sample_matrix_3h_mean <- cbind(get_roma_results_sample_matrix_MTX_3h_mean, get_roma_results_sample_matrix_OXA_3h_mean, get_roma_results_sample_matrix_CIS_3h_mean, get_roma_results_sample_matrix_Control_3h_mean)



# 6h

get_roma_results_sample_matrix_MTX_6h <- get_roma_results_sample_matrix %>% select(contains("6h") & contains("MTX"))

get_roma_results_sample_matrix_MTX_6h_mean <- rowMeans(get_roma_results_sample_matrix_MTX_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_6h" = ".")


get_roma_results_sample_matrix_OXA_6h <- get_roma_results_sample_matrix %>% select(contains("6h") & contains("OXA"))

get_roma_results_sample_matrix_OXA_6h_mean <- rowMeans(get_roma_results_sample_matrix_OXA_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_6h" = ".")


get_roma_results_sample_matrix_CIS_6h <- get_roma_results_sample_matrix %>% select(contains("6h") & contains("CIS"))

get_roma_results_sample_matrix_CIS_6h_mean <- rowMeans(get_roma_results_sample_matrix_CIS_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_6h" = ".")


get_roma_results_sample_matrix_Control_6h <- get_roma_results_sample_matrix %>% select(contains("6h") & contains("Control"))

get_roma_results_sample_matrix_Control_6h_mean <- rowMeans(get_roma_results_sample_matrix_Control_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_6h" = ".")



get_roma_results_sample_matrix_6h_mean <- cbind(get_roma_results_sample_matrix_MTX_6h_mean, get_roma_results_sample_matrix_OXA_6h_mean, get_roma_results_sample_matrix_CIS_6h_mean, get_roma_results_sample_matrix_Control_6h_mean)





# 12h

get_roma_results_sample_matrix_MTX_12h <- get_roma_results_sample_matrix %>% select(contains("12h") & contains("MTX"))

get_roma_results_sample_matrix_MTX_12h_mean <- rowMeans(get_roma_results_sample_matrix_MTX_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_12h" = ".")


get_roma_results_sample_matrix_OXA_12h <- get_roma_results_sample_matrix %>% select(contains("12h") & contains("OXA"))

get_roma_results_sample_matrix_OXA_12h_mean <- rowMeans(get_roma_results_sample_matrix_OXA_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_12h" = ".")


get_roma_results_sample_matrix_CIS_12h <- get_roma_results_sample_matrix %>% select(contains("12h") & contains("CIS"))

get_roma_results_sample_matrix_CIS_12h_mean <- rowMeans(get_roma_results_sample_matrix_CIS_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_12h" = ".")


get_roma_results_sample_matrix_Control_12h <- get_roma_results_sample_matrix %>% select(contains("12h") & contains("Control"))

get_roma_results_sample_matrix_Control_12h_mean <- rowMeans(get_roma_results_sample_matrix_Control_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_12h" = ".")



get_roma_results_sample_matrix_12h_mean <- cbind(get_roma_results_sample_matrix_MTX_12h_mean, get_roma_results_sample_matrix_OXA_12h_mean, get_roma_results_sample_matrix_CIS_12h_mean, get_roma_results_sample_matrix_Control_12h_mean)





# 18h

get_roma_results_sample_matrix_MTX_18h <- get_roma_results_sample_matrix %>% select(contains("18h") & contains("MTX"))

get_roma_results_sample_matrix_MTX_18h_mean <- rowMeans(get_roma_results_sample_matrix_MTX_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_18h" = ".")


get_roma_results_sample_matrix_OXA_18h <- get_roma_results_sample_matrix %>% select(contains("18h") & contains("OXA"))

get_roma_results_sample_matrix_OXA_18h_mean <- rowMeans(get_roma_results_sample_matrix_OXA_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_18h" = ".")


get_roma_results_sample_matrix_CIS_18h <- get_roma_results_sample_matrix %>% select(contains("18h") & contains("CIS"))

get_roma_results_sample_matrix_CIS_18h_mean <- rowMeans(get_roma_results_sample_matrix_CIS_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_18h" = ".")


get_roma_results_sample_matrix_Control_18h <- get_roma_results_sample_matrix %>% select(contains("18h") & contains("Control"))

get_roma_results_sample_matrix_Control_18h_mean <- rowMeans(get_roma_results_sample_matrix_Control_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_18h" = ".")



get_roma_results_sample_matrix_18h_mean <- cbind(get_roma_results_sample_matrix_MTX_18h_mean, get_roma_results_sample_matrix_OXA_18h_mean, get_roma_results_sample_matrix_CIS_18h_mean, get_roma_results_sample_matrix_Control_18h_mean)





# 24h

get_roma_results_sample_matrix_MTX_24h <- get_roma_results_sample_matrix %>% select(contains("24h") & contains("MTX"))

get_roma_results_sample_matrix_MTX_24h_mean <- rowMeans(get_roma_results_sample_matrix_MTX_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_24h" = ".")


get_roma_results_sample_matrix_OXA_24h <- get_roma_results_sample_matrix %>% select(contains("24h") & contains("OXA"))

get_roma_results_sample_matrix_OXA_24h_mean <- rowMeans(get_roma_results_sample_matrix_OXA_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_24h" = ".")


get_roma_results_sample_matrix_CIS_24h <- get_roma_results_sample_matrix %>% select(contains("24h") & contains("CIS"))

get_roma_results_sample_matrix_CIS_24h_mean <- rowMeans(get_roma_results_sample_matrix_CIS_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_24h" = ".")


get_roma_results_sample_matrix_Control_24h <- get_roma_results_sample_matrix %>% select(contains("24h") & contains("Control"))

get_roma_results_sample_matrix_Control_24h_mean <- rowMeans(get_roma_results_sample_matrix_Control_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_24h" = ".")



get_roma_results_sample_matrix_24h_mean <- cbind(get_roma_results_sample_matrix_MTX_24h_mean, get_roma_results_sample_matrix_OXA_24h_mean, get_roma_results_sample_matrix_CIS_24h_mean, get_roma_results_sample_matrix_Control_24h_mean)



get_roma_results_sample_matrix_all_timepoint_mean <- cbind(get_roma_results_sample_matrix_3h_mean,get_roma_results_sample_matrix_6h_mean, get_roma_results_sample_matrix_12h_mean, get_roma_results_sample_matrix_18h_mean, get_roma_results_sample_matrix_24h_mean)





overdispersed_modules_module_matrix <- get_roma_results_module_matrix %>% filter(ppv.L1 <= 0.05 & ppv.Median.Exp > 0.05) %>% arrange(ppv.L1) %>% filter(modules %in% tfs_target_at_least_one_gene_gautier)


shifted_modules_module_matrix <- get_roma_results_module_matrix %>% filter(ppv.Median.Exp <= 0.05) %>% arrange(ppv.Median.Exp) %>% filter(modules %in% tfs_target_at_least_one_gene_gautier) %>% dplyr::mutate(rank = seq(1,nrow(.)))


overdispersed_modules_sample_matrix <- get_roma_results_sample_matrix_all_timepoint_mean %>% rownames_to_column("modules") %>% filter(modules %in% overdispersed_modules_module_matrix$modules) %>% column_to_rownames("modules")
rownames(overdispersed_modules_sample_matrix) <- sub("_TARGET_GENES","",rownames(overdispersed_modules_sample_matrix)) %>% str_to_title()


shifted_modules_sample_matrix <- get_roma_results_sample_matrix_all_timepoint_mean %>% rownames_to_column("modules") %>% inner_join(.,shifted_modules_module_matrix[,c('modules','rank')], by="modules") %>% arrange(rank) %>% column_to_rownames("modules")
rownames(shifted_modules_sample_matrix) <- sub("_TARGET_GENES","",rownames(shifted_modules_sample_matrix)) %>% str_to_title()





# Visualization (Heatmap)


# Overdispersed TFs


annotation_columns_roma_mean <- data.frame(treatment = rep(rep(c("MTX","OXA","CIS", "Control"),c(1,1,1,1)), 5))
rownames(annotation_columns_roma_mean)[1:4] <- colnames(overdispersed_modules_sample_matrix)[1:4]
rownames(annotation_columns_roma_mean)[5:8] <- colnames(overdispersed_modules_sample_matrix)[5:8]
rownames(annotation_columns_roma_mean)[9:12] <- colnames(overdispersed_modules_sample_matrix)[9:12]
rownames(annotation_columns_roma_mean)[13:16] <- colnames(overdispersed_modules_sample_matrix)[13:16]
rownames(annotation_columns_roma_mean)[17:20] <- colnames(overdispersed_modules_sample_matrix)[17:20]


annot_colors_roma <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen", Control="darkgrey"))


heatmap_overdispersed_modules <- ComplexHeatmap::pheatmap(as.matrix(overdispersed_modules_sample_matrix), color= my_color,cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11, fontsize_col = 11, annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="28 overdispersed TFs", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")






# Shifted TFs


annotation_columns_roma_mean <- data.frame(treatment = rep(rep(c("MTX","OXA","CIS", "Control"),c(1,1,1,1)), 5))
rownames(annotation_columns_roma_mean)[1:4] <- colnames(shifted_modules_sample_matrix)[1:4]
rownames(annotation_columns_roma_mean)[5:8] <- colnames(shifted_modules_sample_matrix)[5:8]
rownames(annotation_columns_roma_mean)[9:12] <- colnames(shifted_modules_sample_matrix)[9:12]
rownames(annotation_columns_roma_mean)[13:16] <- colnames(shifted_modules_sample_matrix)[13:16]
rownames(annotation_columns_roma_mean)[17:20] <- colnames(shifted_modules_sample_matrix)[17:20]


annot_colors_roma <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen", Control="darkgrey"))


heatmap_top_50_shifted_modules <- ComplexHeatmap::pheatmap(as.matrix(shifted_modules_sample_matrix[1:50,-21]), color= my_color,cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11, fontsize_col = 11, annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="1-50 TFs", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



heatmap_last_50_shifted_modules <- ComplexHeatmap::pheatmap(as.matrix(shifted_modules_sample_matrix[51:99,-21]), color= my_color,cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11, fontsize_col = 11, annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="51 - 99 TFs", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



shifted_roma_patch <- patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(heatmap_top_50_shifted_modules)), grid.grabExpr(ComplexHeatmap::draw(heatmap_last_50_shifted_modules))), ) & patchwork::plot_annotation(tag_levels = "A", title="TF_activity") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face = "bold", hjust=0.5))






# Keep track of the targets


get_roma_results_module_matrix_for_targets <- get_roma_results_module_matrix %>% column_to_rownames("modules")
rownames(get_roma_results_module_matrix_for_targets) <- sub("_TARGET_GENES", "", rownames(get_roma_results_module_matrix_for_targets)) %>% str_to_title()




# Overdispersed TFs


for (i in 1:nrow(overdispersed_modules_sample_matrix)){
  
  tf_targets_all_info <- GetTopContrib(run_roma_all_timepoints, Selected = which(rownames(get_roma_results_module_matrix_for_targets) == rownames(overdispersed_modules_sample_matrix)[i]), nGenes = 3000, OrderType = "Abs", Mode = "Wei", Plot = TRUE)$Table %>% dplyr::select(c("Gene","Weight")) %>% as.data.frame() %>% dplyr::rename("target" = "Gene") %>% dplyr::arrange(desc(abs(Weight))) %>% dplyr::mutate(Gautier_genes = if_else(target %in% genes_for_dea_mouse$genes, "Yes", "No"))
  rownames(tf_targets_all_info) <- NULL
  
  if (rownames(overdispersed_modules_sample_matrix)[i] %in% collectri_df$source){
    collectri_df_tf_tracked <- collectri_df %>% filter(source == rownames(overdispersed_modules_sample_matrix)[i])
    tf_targets_all_info <- left_join(tf_targets_all_info, collectri_df_tf_tracked[,c("target","mor")], by = "target")
  }
  
  
  write.csv(tf_targets_all_info, paste0("Analysis_for_Gautier/Upstream_analysis/rRoma/TF_targets_info/Overdispersed_TFs/",rownames(overdispersed_modules_sample_matrix)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
}





# Shifted TFs 

for (i in 1:nrow(shifted_modules_sample_matrix)){
  
  tf_targets_all_info <- GetTopContrib(run_roma_all_timepoints, Selected = which(rownames(get_roma_results_module_matrix_for_targets) == rownames(shifted_modules_sample_matrix)[i]), nGenes = 3000, OrderType = "Abs", Mode = "Wei", Plot = TRUE)$Table %>% dplyr::select(c("Gene","Weight")) %>% as.data.frame() %>% dplyr::rename("target" = "Gene") %>% dplyr::arrange(desc(abs(Weight))) %>% dplyr::mutate(Gautier_genes = if_else(target %in% genes_for_dea_mouse$genes, "Yes", "No"))
  rownames(tf_targets_all_info) <- NULL
  
  if (rownames(shifted_modules_sample_matrix)[i] %in% collectri_df$source){
    collectri_df_tf_tracked <- collectri_df %>% filter(source == rownames(shifted_modules_sample_matrix)[i])
    tf_targets_all_info <- left_join(tf_targets, collectri_df_tf_tracked[,c("target","mor")], by = "target")
  }
  
  
  write.csv(tf_targets_all_info, paste0("Analysis_for_Gautier/Upstream_analysis/rRoma/TF_targets_info/Shifted_TFs/",rownames(shifted_modules_sample_matrix)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
}


