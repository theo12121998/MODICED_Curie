library(rRoma)
library(tidyverse)
library(edgeR)
library(ComplexHeatmap)
library(scales)
library(patchwork)
library(grid)



## ROMA


counts_in_vitro_df <- readRDS("~/MODICED_link/Bulk_RNA_in_vitro/preprocessing_in_vitro_RNA/counts_in_vitro_for_analysis.rds") %>% as.data.frame()
counts_in_vitro_df_genes <- counts_in_vitro_df %>% dplyr::select(1:3)
counts_in_vitro_df_samples <- counts_in_vitro_df %>% dplyr::select(-c(1:3))
sample_plan_in_vitro <- readRDS("~/MODICED_link/Bulk_RNA_in_vitro/preprocessing_in_vitro_RNA/sample_plan_in_vitro_for_analysis.rds")
attach(sample_plan_in_vitro)
dlist <- DGEList(counts = counts_in_vitro_df_samples, samples = sample_plan_in_vitro$Sample_Name, genes = counts_in_vitro_df_genes$gene_name)
dsgn <- model.matrix(~ 0 + Treatment_hour, data=dlist$samples)
colnames(dsgn) <- sub("Treatment_hour","",colnames(dsgn))
keep <- filterByExpr(dlist$counts, dsgn)
dlist_filtered <- dlist[keep,]



counts_for_roma <- cbind(dlist_filtered$counts, dlist_filtered$genes)
rownames(counts_for_roma) <- make.names(counts_for_roma$genes, unique=TRUE)
counts_for_roma <- counts_for_roma %>% select(-"genes")
counts_in_vitro_roma_all_timepoints <- cpm(counts_for_roma, log=TRUE) %>% as.data.frame()


tf_gmt <- ReadGMTFile("/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/m3.gtrd.v2023.2.Mm.symbols.gmt")
run_roma_tf_all_timepoints <- rRoma.R(counts_in_vitro_roma_all_timepoints, tf_gmt)
saveRDS(run_roma_tf_all_timepoints, "/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/roma_for_tf_analysis_all_timepoints.rds")


get_roma_results <- readRDS("/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/roma_for_tf_analysis_all_timepoints.rds")

get_roma_results_module_matrix <- get_roma_results$ModuleMatrix %>% as.data.frame()
get_roma_results_sample_matrix <- get_roma_results$SampleMatrix %>% as.data.frame()



p_cyto_interest <- c()
n_cyto_interest <- 0
for (i in 1:length(tf_gmt)){
  if ("Ccl5" %in% tf_gmt[[i]]$Genes | "Il6" %in% tf_gmt[[i]]$Genes | "Ccl4" %in% tf_gmt[[i]]$Genes | "Ccl20" %in% tf_gmt[[i]]$Genes | "Ccl9" %in% tf_gmt[[i]]$Genes){
    print(i)
    n_cyto_interest <- n_cyto_interest+1
    p_cyto_interest <- append(p_cyto_interest, tf_gmt[[i]]$Name)
  }
}
print(n_cyto_interest)
print(p_cyto_interest)


p_cyto_jo <- c()
n_cyto_jo <- 0
for (i in 1:length(tf_gmt)){
  if (any(tf_gmt[[i]]$Genes %in% jo_list_of_cyto$cytos)){
    print(i)
    n_cyto_jo <- n_cyto_jo+1
    p_cyto_jo <- append(p_cyto_jo, tf_gmt[[i]]$Name)
  }
}
print(n_cyto_jo)
print(p_cyto_jo)







# All roma








# Get the mean activity for each TF

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
saveRDS(get_roma_results_sample_matrix_all_timepoint_mean, "/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/roma_for_tf_analysis_all_timepoints_mean.rds" )















get_roma_results_sample_matrix_all_timepoint_mean <- readRDS("/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/roma_for_tf_analysis_all_timepoints_mean.rds")


get_roma_results <- readRDS("/data/users/tlassale/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/roma_for_tf_analysis_all_timepoints.rds")
get_roma_results_module_matrix <- get_roma_results$ModuleMatrix %>% as.data.frame()




# Get overdispersed modules
overdispersed_modules_arranged <- get_roma_results_module_matrix %>% filter(`ppv L1` <= 0.05) %>% rownames_to_column("modules") %>% arrange(`ppv L1`) %>% dplyr::mutate(rank = seq(1:nrow(.)))





# All TFs



get_roma_results_sample_matrix_all_timepoint_mean_overdispersed <- get_roma_results_sample_matrix_all_timepoint_mean %>% filter(rownames(.) %in% overdispersed_modules_arranged$modules)


get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50 <- get_roma_results_sample_matrix_all_timepoint_mean %>% filter(rownames(.) %in% overdispersed_modules_arranged$modules[1:50])




# Heatmap

palette_length <- 100
my_color <- viridis::viridis(palette_length)

rownames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50) <- sub("_TARGET_GENES","",rownames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50))


annotation_columns_roma_mean <- data.frame(treatment = rep(rep(c("MTX","OXA","CIS", "Control"),c(1,1,1,1)), 5))
rownames(annotation_columns_roma_mean)[1:4] <- colnames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50)[1:4]
rownames(annotation_columns_roma_mean)[5:8] <- colnames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50)[5:8]
rownames(annotation_columns_roma_mean)[9:12] <- colnames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50)[9:12]
rownames(annotation_columns_roma_mean)[13:16] <- colnames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50)[13:16]
rownames(annotation_columns_roma_mean)[17:20] <- colnames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50)[17:20]



annot_colors_roma <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen", Control="yellow"))




get_roma_results_sample_matrix_all_timepoint_mean_heatmap_top_50_overdispersed <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_top_50), cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 12 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="Activity of TFs that target at least one cytokine of interest", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")






# TFs with at least one cyto of interest

overdispersed_modules_arranged_cyto_interest <- overdispersed_modules_arranged  %>% filter(modules %in% p_cyto_interest) 
get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_interest <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed %>% filter(rownames(.) %in% overdispersed_modules_arranged_cyto_interest$modules)



# Heatmap


rownames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_interest) <- sub("_TARGET_GENES","",rownames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_interest))



get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_interest_heatmap <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_interest), color = my_color, cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="Overdispersed TFs (target at least one cytokine of interest)", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")







# TFs with at least one cyto of Jo

overdispersed_modules_arranged_jo <- overdispersed_modules_arranged  %>% filter(modules %in% p_cyto_jo)  
overdispersed_modules_arranged_jo_top_50 <- overdispersed_modules_arranged  %>% filter(modules %in% p_cyto_jo) %>% slice(1:50)
overdispersed_modules_arranged_jo_last_20 <- overdispersed_modules_arranged  %>% filter(modules %in% p_cyto_jo) %>% slice(51:70)



# top 50
get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed %>% filter(rownames(.) %in% overdispersed_modules_arranged_jo_top_50$modules)

# last 20
get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed %>% filter(rownames(.) %in% overdispersed_modules_arranged_jo_last_20$modules)



# Heatmaps



## top 50

rownames(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo) <- sub("_TARGET_GENES","",rownames(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo))


get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo_scaled_01 <- rescale(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo))


get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo_scaled <- scale(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo))

get_roma_results_sample_matrix_all_timepoint_mean_heatmap_top_50_overdispersed_cytos_jo <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo_scaled_01), cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="Top 50 overdispersed TFs that target at least one cytos of Jo", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



# standard scale
get_roma_results_sample_matrix_all_timepoint_mean_heatmap_top_50_overdispersed_cytos_jo_scaled <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo_scaled), cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="Top 50 overdispersed TFs that target at least one cytos of Jo", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



# not scaled
get_roma_results_sample_matrix_all_timepoint_mean_heatmap_top_50_overdispersed_cytos_jo <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo), main = "1 - 50 TFs", color = my_color, cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 11 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")




## last 20

rownames(get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo) <- sub("_TARGET_GENES","",rownames(get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo))


get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo_scaled_01 <- rescale(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo))


# not scaled
get_roma_results_sample_matrix_all_timepoint_mean_heatmap_last_20_overdispersed_cytos_jo <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo), main = "51 - 70 TFs",cluster_rows = TRUE, color = my_color,cluster_cols = TRUE, fontsize_row = 11 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")






tf_roma_patch <-  patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(get_roma_results_sample_matrix_all_timepoint_mean_heatmap_top_50_overdispersed_cytos_jo)), grid.grabExpr(ComplexHeatmap::draw(get_roma_results_sample_matrix_all_timepoint_mean_heatmap_last_20_overdispersed_cytos_jo)))) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13)) & patchwork::plot_annotation(tag_levels = "A", title = "Overdispersed TFs (target at least one cytokine of Jo's list)") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face="bold", hjust=0.5))







# Activity + expression of the overdispersed TFs of interest (target list jo) - CIS


# expression
overdispersed_top_50_tfs_jo <- get_roma_results_sample_matrix_all_timepoint_mean_top_50_overdispersed_cytos_jo %>% rownames_to_column("TFs") %>% pull("TFs")
overdispersed_top_50_tf_jo_mouse_symbols <- overdispersed_top_50_tfs_jo %>% str_to_title()
overdispersed_top_50_tf_jo_mouse_symbols_new <- gsub("_target_genes","",overdispersed_top_50_tf_jo_mouse_symbols)



overdispersed_last_20_tfs_jo <- get_roma_results_sample_matrix_all_timepoint_mean_last_20_overdispersed_cytos_jo %>% rownames_to_column("TFs") %>% pull("TFs")
overdispersed_last_20_tf_jo_mouse_symbols <- overdispersed_last_20_tfs_jo %>% str_to_title()
overdispersed_last_20_tf_jo_mouse_symbols_new <- gsub("_target_genes","",overdispersed_last_20_tf_jo_mouse_symbols)





tfs_first_14_interest <- counts_in_vitro_roma_all_timepoints %>% filter(rownames(.) %in% overdispersed_top_50_tf_jo_mouse_symbols_new) %>% filter(rownames(.) %in% c("Hoxc9","Tead2","Nfatc2","Zfp661","Pdx1","Prdm9","Sp5","Mta3","Zfp953","Dlx1","Myef2_Myef2l","Sall1","Zfp92","Smarca2"))
tfs_last_14_interest <- counts_in_vitro_roma_all_timepoints %>% filter(rownames(.) %in% overdispersed_top_50_tf_jo_mouse_symbols_new) %>% filter(rownames(.) %in% c("Nrob2","Pwwp2b","Gm14412","Zfp809","Msx1","Nfe2","Nacc1","Zfp449","Gcm2","Foxd3","Pcgf5","Tox2","Hdgfl2","Msgn1"))
tfs_last_13_interest <- counts_in_vitro_roma_all_timepoints %>% filter(rownames(.) %in% overdispersed_last_20_tf_jo_mouse_symbols_new) %>% filter(rownames(.) %in% c("Hells","Nfe2l1","Tbx4","Hdac4","Hlf","Zfp595","Nfatc1","Foxa3","Zfp429","Ruvbl1","Tbx20","Anp32e","Mllt3"))
tf_all_41_interest <- rbind(tfs_first_14_interest, tfs_last_14_interest, tfs_last_13_interest)



merged_41_tfs <- c()
for (i in 1:nrow(tf_all_41_interest)){
  counts_in_vitro_df_i_CIS <- tf_all_41_interest %>% dplyr::slice(i) %>% dplyr::select(contains("CIS"))
  counts_in_vitro_df_i_CIS_new <- counts_in_vitro_df_i_CIS[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(counts_in_vitro_df_i_CIS)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  counts_in_vitro_df_i_CIS_new$name <- rep(rownames(counts_in_vitro_df_i_CIS),5)
  counts_in_vitro_df_i_CIS_new$treatment <- rep("CIS",5)
  merged_41_tfs <- rbind(merged_41_tfs, counts_in_vitro_df_i_CIS_new)
}




get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_jo_scaled_01 <- rescale(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed %>% filter(rownames(.) %in% overdispersed_modules_arranged_jo$modules)))


# activation
activity_all_overdispersed_tfs <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_jo_scaled_01
rownames(activity_all_overdispersed_tfs) <- str_to_title(rownames(activity_all_overdispersed_tfs)) %>% gsub("_target_genes","",.)
activity_all_overdispersed_tfs_interest <- activity_all_overdispersed_tfs %>% as.data.frame() %>% filter(rownames(.) %in% merged_41_tfs$name) %>% select(contains("CIS"))
activity_all_overdispersed_tfs_interest_ranked <- activity_all_overdispersed_tfs_interest[unique(merged_41_tfs$name),] 
merged_41_tfs_for_activity <- merged_41_tfs %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv <- t(activity_all_overdispersed_tfs_interest_ranked) %>% as.data.frame()

mean_activity <- c()
for (i in 1:ncol(bv)){
  mean_activity <- rbind(mean_activity, bv[,i] %>% as.data.frame())
}

merged_41_tfs_for_activity$Mean_activity <- mean_activity$.




for (i in seq(1,140,by=5)){
  
  # expression
  merged_4_tfs_expression <- ggplot(merged_41_tfs[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[i,"name"], " in CIS ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  merged_4_tfs_activation <- ggplot(merged_41_tfs_for_activity[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[i,"name"]," in CIS ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(merged_4_tfs_expression, merged_4_tfs_activation) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  ggsave(paste0("tf_expression_activation/longitudinal_expression_activation_in_CIS_tf_",merged_41_tfs[i,"name"],".png"))
  
  
  
  
  
}







# Activity + expression of the overdispersed TFs of interest (target list jo) - OXA




# expression
merged_41_tfs_OXA <- c()
for (i in 1:nrow(tf_all_41_interest)){
  counts_in_vitro_df_i_OXA <- tf_all_41_interest %>% dplyr::slice(i) %>% dplyr::select(contains("OXA"))
  counts_in_vitro_df_i_OXA_new <- counts_in_vitro_df_i_OXA[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(counts_in_vitro_df_i_OXA)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  counts_in_vitro_df_i_OXA_new$name <- rep(rownames(counts_in_vitro_df_i_OXA),5)
  counts_in_vitro_df_i_OXA_new$treatment <- rep("OXA",5)
  merged_41_tfs_OXA <- rbind(merged_41_tfs_OXA, counts_in_vitro_df_i_OXA_new)
}





# activation
activity_all_overdispersed_tfs <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_jo_scaled_01
rownames(activity_all_overdispersed_tfs) <- str_to_title(rownames(activity_all_overdispersed_tfs))
activity_all_overdispersed_tfs_interest <- activity_all_overdispersed_tfs %>% as.data.frame() %>% filter(rownames(.) %in% merged_41_tfs_OXA$name) %>% select(contains("OXA"))
activity_all_overdispersed_tfs_interest_ranked <- activity_all_overdispersed_tfs_interest[unique(merged_41_tfs_OXA$name),] 
merged_41_tfs_OXA_for_activity <- merged_41_tfs_OXA %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv <- t(activity_all_overdispersed_tfs_interest_ranked) %>% as.data.frame()

mean_activity <- c()
for (i in 1:ncol(bv)){
  mean_activity <- rbind(mean_activity, bv[,i] %>% as.data.frame())
}

merged_41_tfs_OXA_for_activity$Mean_activity <- mean_activity$.





for (i in seq(1,140,by=5)){
  
  # expression
  merged_4_tfs_expression_OXA <- ggplot(merged_41_tfs_OXA[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[i,"name"], " in OXA ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  merged_4_tfs_activation_OXA <- ggplot(merged_41_tfs_OXA_for_activity[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[i,"name"]," in OXA ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(merged_4_tfs_expression_OXA, merged_4_tfs_activation_OXA) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  ggsave(paste0("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/OXA/longitudinal_expression_activation_in_OXA_tf_",merged_41_tfs_OXA[i,"name"],".png"))
  
  
  
  
  
}






# Activity + expression of the overdispersed TFs of interest (target list jo) - MTX




# expression
merged_41_tfs_MTX <- c()
for (i in 1:nrow(tf_all_41_interest)){
  counts_in_vitro_df_i_MTX <- tf_all_41_interest %>% dplyr::slice(i) %>% dplyr::select(contains("MTX"))
  counts_in_vitro_df_i_MTX_new <- counts_in_vitro_df_i_MTX[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(counts_in_vitro_df_i_MTX)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  counts_in_vitro_df_i_MTX_new$name <- rep(rownames(counts_in_vitro_df_i_MTX),5)
  counts_in_vitro_df_i_MTX_new$treatment <- rep("MTX",5)
  merged_41_tfs_MTX <- rbind(merged_41_tfs_MTX, counts_in_vitro_df_i_MTX_new)
}





# activation
activity_all_overdispersed_tfs <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_jo_scaled_01
rownames(activity_all_overdispersed_tfs) <- str_to_title(rownames(activity_all_overdispersed_tfs))
activity_all_overdispersed_tfs_interest <- activity_all_overdispersed_tfs %>% as.data.frame() %>% filter(rownames(.) %in% merged_41_tfs_MTX$name) %>% select(contains("MTX"))
activity_all_overdispersed_tfs_interest_ranked <- activity_all_overdispersed_tfs_interest[unique(merged_41_tfs_MTX$name),] 
merged_41_tfs_MTX_for_activity <- merged_41_tfs_MTX %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv <- t(activity_all_overdispersed_tfs_interest_ranked) %>% as.data.frame()

mean_activity <- c()
for (i in 1:ncol(bv)){
  mean_activity <- rbind(mean_activity, bv[,i] %>% as.data.frame())
}

merged_41_tfs_MTX_for_activity$Mean_activity <- mean_activity$.





for (i in seq(1,140,by=5)){
  
  # expression
  merged_4_tfs_expression_MTX <- ggplot(merged_41_tfs_MTX[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[i,"name"], " in MTX ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  merged_4_tfs_activation_MTX <- ggplot(merged_41_tfs_MTX_for_activity[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[i,"name"]," in MTX ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(merged_4_tfs_expression_MTX, merged_4_tfs_activation_MTX) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  ggsave(paste0("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/MTX/longitudinal_expression_activation_in_MTX_tf_",merged_41_tfs_MTX[i,"name"],".png"))
  
  
  
  

}


















# Overdispersed TFs with no cyto of Jo

overdispersed_modules_arranged_not_jo <- overdispersed_modules_arranged  %>% filter(!(modules %in% p_cyto_jo))




# all
get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed %>% filter(rownames(.) %in% overdispersed_modules_arranged_not_jo$modules)



# Heatmap 

rownames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo) <- sub("_TARGET_GENES","",rownames(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo))


#get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo_scaled_01 <- rescale(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo))


get_roma_results_sample_matrix_all_timepoint_mean_heatmap_overdispersed_cytos_not_jo <- ComplexHeatmap::pheatmap(as.matrix(get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo), color = my_color, cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 9 , annotation_col = annotation_columns_roma_mean, annotation_colors = annot_colors_roma, main="23 other overdispersed TFs (with no cytokines of Jo in targets)", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")










# Activity + expression of the overdispersed TFs of interest (no target list jo)


# expression
overdispersed_tfs_not_jo <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo %>% rownames_to_column("TFs") %>% select(TFs)
overdispersed_tf_not_jo_mouse_symbols <- overdispersed_tfs_not_jo$TFs %>% str_to_title() %>% as.data.frame() %>% dplyr::rename("TFs" = ".")
overdispersed_tf_not_jo_mouse_symbols[6,] <- "Gm20449"
overdispersed_tf_not_jo_mouse_symbols[7,] <- "Gm6710"
overdispersed_tf_not_jo_mouse_symbols[22,] <- "Zfp986"







tfs_interest <- counts_in_vitro_roma_all_timepoints %>% filter(rownames(.) %in% overdispersed_tf_not_jo_mouse_symbols$TFs) %>% filter(rownames(.) %in% c("Sox5","Dnajc2","Gli1","Lhx1","Gm14410","Zfp986","Gm6710","Klf6","Satb2","Gm20449","Smad6"))



# CIS
merged_tfs_interest <- c()
for (i in 1:nrow(tfs_interest)){
  counts_in_vitro_df_i_CIS <- tfs_interest %>% dplyr::slice(i) %>% dplyr::select(contains("CIS"))
  counts_in_vitro_df_i_CIS_new <- counts_in_vitro_df_i_CIS[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(counts_in_vitro_df_i_CIS)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  counts_in_vitro_df_i_CIS_new$name <- rep(rownames(counts_in_vitro_df_i_CIS),5)
  counts_in_vitro_df_i_CIS_new$treatment <- rep("CIS",5)
  merged_tfs_interest <- rbind(merged_tfs_interest, counts_in_vitro_df_i_CIS_new)
}




# activation
activity_all_overdispersed_tfs <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo_scaled_01
rownames(activity_all_overdispersed_tfs) <- str_to_title(rownames(activity_all_overdispersed_tfs))
rownames(activity_all_overdispersed_tfs)[6] <- "Gm20449"
rownames(activity_all_overdispersed_tfs)[7] <- "Gm6710"
rownames(activity_all_overdispersed_tfs)[22] <- "Zfp986"
activity_all_overdispersed_tfs_interest <- activity_all_overdispersed_tfs %>% as.data.frame() %>% filter(rownames(.) %in% merged_tfs_interest$name) %>% select(contains("CIS"))
activity_all_overdispersed_tfs_interest_ranked <- activity_all_overdispersed_tfs_interest[unique(merged_tfs_interest$name),] 
merged_tfs_interest_for_activity <- merged_tfs_interest %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv <- t(activity_all_overdispersed_tfs_interest_ranked) %>% as.data.frame()

mean_activity <- c()
for (i in 1:ncol(bv)){
  mean_activity <- rbind(mean_activity, bv[,i] %>% as.data.frame())
}

merged_tfs_interest_for_activity$Mean_activity <- mean_activity$.




for (i in seq(1,40,by=5)){
  
  # expression
  merged_4_tfs_expression <- ggplot(merged_tfs_interest[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",merged_tfs_interest[i,"name"], " in CIS ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  merged_4_tfs_activation <- ggplot(merged_tfs_interest_for_activity[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",merged_tfs_interest_for_activity[i,"name"]," in CIS ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(merged_4_tfs_expression, merged_4_tfs_activation) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  ggsave(paste0("tf_expression_activation/Overdispersed TFs (without cyto of jo)/longitudinal_expression_activation_in_CIS_tf_",merged_tfs_interest[i,"name"],".png"))
  
  
  
  
  
}





# MTX
merged_tfs_interest <- c()
for (i in 1:nrow(tfs_interest)){
  counts_in_vitro_df_i_MTX <- tfs_interest %>% dplyr::slice(i) %>% dplyr::select(contains("MTX"))
  counts_in_vitro_df_i_MTX_new <- counts_in_vitro_df_i_MTX[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(counts_in_vitro_df_i_MTX)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  counts_in_vitro_df_i_MTX_new$name <- rep(rownames(counts_in_vitro_df_i_MTX),5)
  counts_in_vitro_df_i_MTX_new$treatment <- rep("MTX",5)
  merged_tfs_interest <- rbind(merged_tfs_interest, counts_in_vitro_df_i_MTX_new)
}




# activation
activity_all_overdispersed_tfs <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo_scaled_01
rownames(activity_all_overdispersed_tfs) <- str_to_title(rownames(activity_all_overdispersed_tfs))
rownames(activity_all_overdispersed_tfs)[6] <- "Gm20449"
rownames(activity_all_overdispersed_tfs)[7] <- "Gm6710"
rownames(activity_all_overdispersed_tfs)[22] <- "Zfp986"
activity_all_overdispersed_tfs_interest <- activity_all_overdispersed_tfs %>% as.data.frame() %>% filter(rownames(.) %in% merged_tfs_interest$name) %>% select(contains("MTX"))
activity_all_overdispersed_tfs_interest_ranked <- activity_all_overdispersed_tfs_interest[unique(merged_tfs_interest$name),] 
merged_tfs_interest_for_activity <- merged_tfs_interest %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv <- t(activity_all_overdispersed_tfs_interest_ranked) %>% as.data.frame()

mean_activity <- c()
for (i in 1:ncol(bv)){
  mean_activity <- rbind(mean_activity, bv[,i] %>% as.data.frame())
}

merged_tfs_interest_for_activity$Mean_activity <- mean_activity$.




for (i in seq(1,40,by=5)){
  
  # expression
  merged_4_tfs_expression <- ggplot(merged_tfs_interest[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",merged_tfs_interest[i,"name"], " in MTX ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  merged_4_tfs_activation <- ggplot(merged_tfs_interest_for_activity[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",merged_tfs_interest_for_activity[i,"name"]," in MTX ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(merged_4_tfs_expression, merged_4_tfs_activation) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  ggsave(paste0("tf_expression_activation/Overdispersed TFs (without cyto of jo)/MTX/longitudinal_expression_activation_in_MTX_tf_",merged_tfs_interest[i,"name"],".png"))
  
  
  
  
  
}






# OXA
merged_tfs_interest <- c()
for (i in 1:nrow(tfs_interest)){
  counts_in_vitro_df_i_OXA <- tfs_interest %>% dplyr::slice(i) %>% dplyr::select(contains("OXA"))
  counts_in_vitro_df_i_OXA_new <- counts_in_vitro_df_i_OXA[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(counts_in_vitro_df_i_OXA)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  counts_in_vitro_df_i_OXA_new$name <- rep(rownames(counts_in_vitro_df_i_OXA),5)
  counts_in_vitro_df_i_OXA_new$treatment <- rep("OXA",5)
  merged_tfs_interest <- rbind(merged_tfs_interest, counts_in_vitro_df_i_OXA_new)
}




# activation
activity_all_overdispersed_tfs <- get_roma_results_sample_matrix_all_timepoint_mean_overdispersed_cytos_not_jo_scaled_01
rownames(activity_all_overdispersed_tfs) <- str_to_title(rownames(activity_all_overdispersed_tfs))
rownames(activity_all_overdispersed_tfs)[6] <- "Gm20449"
rownames(activity_all_overdispersed_tfs)[7] <- "Gm6710"
rownames(activity_all_overdispersed_tfs)[22] <- "Zfp986"
activity_all_overdispersed_tfs_interest <- activity_all_overdispersed_tfs %>% as.data.frame() %>% filter(rownames(.) %in% merged_tfs_interest$name) %>% select(contains("OXA"))
activity_all_overdispersed_tfs_interest_ranked <- activity_all_overdispersed_tfs_interest[unique(merged_tfs_interest$name),] 
merged_tfs_interest_for_activity <- merged_tfs_interest %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv <- t(activity_all_overdispersed_tfs_interest_ranked) %>% as.data.frame()

mean_activity <- c()
for (i in 1:ncol(bv)){
  mean_activity <- rbind(mean_activity, bv[,i] %>% as.data.frame())
}

merged_tfs_interest_for_activity$Mean_activity <- mean_activity$.




for (i in seq(1,40,by=5)){
  
  # expression
  merged_4_tfs_expression <- ggplot(merged_tfs_interest[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",merged_tfs_interest[i,"name"], " in OXA ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  merged_4_tfs_activation <- ggplot(merged_tfs_interest_for_activity[i:(i+4),], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",merged_tfs_interest_for_activity[i,"name"]," in OXA ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(merged_4_tfs_expression, merged_4_tfs_activation) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  ggsave(paste0("tf_expression_activation/Overdispersed TFs (without cyto of jo)/OXA/longitudinal_expression_activation_in_OXA_tf_",merged_tfs_interest[i,"name"],".png"))
  
  
  
  
  
}















# Prdm9 + Zfp661

merged_4_tfs_expression_1 <- ggplot(merged_41_tfs[31:35,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[31,"name"], " in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_Expression")



# activation
merged_4_tfs_activation_1 <- ggplot(merged_41_tfs_for_activity[31:35,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[31,"name"]," in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_activity_score")


n <- wrap_plots(merged_4_tfs_expression_1, merged_4_tfs_activation_1) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))







merged_4_tfs_expression_2 <- ggplot(merged_41_tfs[6:10,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[6,"name"], " in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_Expression")



# activation
merged_4_tfs_activation_2 <- ggplot(merged_41_tfs_for_activity[6:10,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[6,"name"]," in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_activity_score")


b <- wrap_plots(merged_4_tfs_expression_2, merged_4_tfs_activation_2) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))


wrap_plots(n,b, ncol=1) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))














# Dlx1 + Hoxc9 + Tead2 + Zfp449


merged_4_tfs_expression_3 <- ggplot(merged_41_tfs[1:5,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[1,"name"], " in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_Expression")



# activation
merged_4_tfs_activation_3 <- ggplot(merged_41_tfs_for_activity[1:5,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[1,"name"]," in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_activity_score")


u <- wrap_plots(merged_4_tfs_expression_3, merged_4_tfs_activation_3) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))






merged_4_tfs_expression_4 <- ggplot(merged_41_tfs[26:30,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[26,"name"], " in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_Expression")



# activation
merged_4_tfs_activation_4 <- ggplot(merged_41_tfs_for_activity[26:30,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[26,"name"]," in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_activity_score")


f <- wrap_plots(merged_4_tfs_expression_4, merged_4_tfs_activation_4) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))












merged_4_tfs_expression_5 <- ggplot(merged_41_tfs[16:20,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[16,"name"], " in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_Expression")



# activation
merged_4_tfs_activation_5 <- ggplot(merged_41_tfs_for_activity[16:20,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[16,"name"]," in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_activity_score")


g <- wrap_plots(merged_4_tfs_expression_5, merged_4_tfs_activation_5) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))







merged_4_tfs_expression_6 <- ggplot(merged_41_tfs[81:85,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal expression of ",merged_41_tfs[81,"name"], " in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_Expression")



# activation
merged_4_tfs_activation_6 <- ggplot(merged_41_tfs_for_activity[81:85,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
  geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
  ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
  geom_path(aes(color=name, group=name), size=1.1) +
  scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
  ggtitle(paste0("Longitudinal activity score of ",merged_41_tfs_for_activity[81,"name"]," in CIS ")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.title = element_text(face="bold")) +
  theme(axis.text = element_text(face="bold")) +
  labs(x = "Timepoint", y = "Mean_activity_score")


h <- wrap_plots(merged_4_tfs_expression_6, merged_4_tfs_activation_6) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))





wrap_plots(u,f,g,h, ncol=1) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))






















# Correlation (but now enough points, as only 5 timepoints)


mean_expression <- rbind(merged_41_tfs[31:35,]$Mean_Expression %>% as.data.frame(), merged_41_tfs[6:10,]$Mean_Expression %>% as.data.frame(), merged_41_tfs[1:5,]$Mean_Expression %>% as.data.frame(), merged_41_tfs[26:30,]$Mean_Expression %>% as.data.frame(), merged_41_tfs[26:30,]$Mean_Expression %>% as.data.frame(), merged_41_tfs[81:85,]$Mean_Expression %>% as.data.frame()) %>% dplyr::rename("mean_expression" = ".")

mean_activity <- rbind(merged_41_tfs_for_activity[31:35,]$Mean_activity %>% as.data.frame(), merged_41_tfs_for_activity[6:10,]$Mean_activity %>% as.data.frame(), merged_41_tfs_for_activity[1:5,]$Mean_activity %>% as.data.frame(), merged_41_tfs_for_activity[26:30,]$Mean_activity %>% as.data.frame(), merged_41_tfs_for_activity[26:30,]$Mean_activity %>% as.data.frame(), merged_41_tfs_for_activity[81:85,]$Mean_activity %>% as.data.frame()) %>% dplyr::rename("mean_activity" = ".")



mean_expression_activity <- cbind(mean_expression, mean_activity)


ggplot(mean_expression_activity[1:5,], aes(x=mean_expression, y=mean_activity)) + 
  geom_point() + geom_line()


cor(mean_expression[1:5,], mean_activity[1:5,])



ggplot(mean_expression_activity[6:10,], aes(x=mean_expression, y=mean_activity)) + 
  geom_point() + geom_line()


cor(mean_expression[6:10,], mean_activity[6:10,])



ggplot(mean_expression_activity[11:15,], aes(x=mean_expression, y=mean_activity)) + 
  geom_point() + geom_line()


cor(mean_expression[11:15,], mean_activity[11:15,])













# Keep track of targets of the overdispersed TFs


overdispersed_modules_arranged_jo_top_50_new <- sub("_TARGET_GENES","",overdispersed_modules_arranged_jo_top_50$modules)
overdispersed_modules_arranged_jo_top_50_new_final <- overdispersed_modules_arranged_jo_top_50_new %>% str_to_title()


overdispersed_modules_arranged_jo_last_20_new <- sub("_TARGET_GENES","",overdispersed_modules_arranged_jo_last_20$modules)
overdispersed_modules_arranged_jo_last_20_new_final <- overdispersed_modules_arranged_jo_last_20_new %>% str_to_title()


overdispersed_modules_arranged_not_jo_new <- sub("_TARGET_GENES","",overdispersed_modules_arranged_not_jo$modules)
overdispersed_modules_arranged_not_jo_new_final <- overdispersed_modules_arranged_not_jo_new %>% str_to_title()


overdispersed_modules_arranged_cyto_interest_new <- sub("_TARGET_GENES","",overdispersed_modules_arranged_cyto_interest$modules)
overdispersed_modules_arranged_cyto_interest_new_final  <- overdispersed_modules_arranged_cyto_interest_new %>% str_to_title()




# top 50 

for (i in 1:nrow(overdispersed_modules_arranged_jo_top_50)){
  
  tf_targets_all_info <- GetTopContrib(get_roma_results, Selected = which(rownames(get_roma_results_module_matrix) == overdispersed_modules_arranged_jo_top_50$modules[i]),
                nGenes = 3000, OrderType = "Abs", Mode = "Wei", Plot = TRUE)$Table %>% dplyr::select(c("Gene","Weight")) %>% as.data.frame() %>% dplyr::rename("target" = "Gene") %>% dplyr::arrange(desc(abs(Weight))) %>% dplyr::mutate(cyto_of_jo = if_else(target %in% jo_list_of_cyto$cytos, "Yes", "No"))
  rownames(tf_targets_all_info) <- NULL
  
  if (overdispersed_modules_arranged_jo_top_50_new_final[i] %in% collectri_df$source){
    collectri_df_tf_tracked <- collectri_df %>% filter(source == overdispersed_modules_arranged_jo_top_50_new_final[i])
    tf_targets_all_info <- left_join(tf_targets_all_info, collectri_df_tf_tracked[,c("target","mor")], by = "target")
  }
 

  #write.csv(tf_targets_all_info, paste0("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/",overdispersed_modules_arranged_jo_top_50$modules[i],".csv"), col.names = TRUE, row.names = FALSE)
  
  write.csv(tf_targets_all_info, paste0("tfs_roma_modiced_for_publi/targets_of_TFs/top_50/",overdispersed_modules_arranged_jo_top_50$modules[i],".csv"), col.names = TRUE, row.names = FALSE)
}





# last 20

for (i in 1:nrow(overdispersed_modules_arranged_jo_last_20)){
  
  tf_targets_all_info <- GetTopContrib(get_roma_results, Selected = which(rownames(get_roma_results_module_matrix) == overdispersed_modules_arranged_jo_last_20$modules[i]),
                                       nGenes = 3000, OrderType = "Abs", Mode = "Wei", Plot = TRUE)$Table %>% dplyr::select(c("Gene","Weight")) %>% as.data.frame() %>% dplyr::rename("target" = "Gene") %>% dplyr::arrange(desc(abs(Weight))) %>% dplyr::mutate(cyto_of_jo = if_else(target %in% jo_list_of_cyto$cytos, "Yes", "No"))
  rownames(tf_targets_all_info) <- NULL
  
  if (overdispersed_modules_arranged_jo_last_20_new_final[i] %in% collectri_df$source){
    collectri_df_tf_tracked <- collectri_df %>% filter(source == overdispersed_modules_arranged_jo_last_20_new_final[i])
    tf_targets_all_info <- left_join(tf_targets_all_info, collectri_df_tf_tracked[,c("target","mor")], by = "target")
  }
  
  
  #write.csv(tf_targets_all_info, paste0("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/",overdispersed_modules_arranged_jo_last_20$modules[i],".csv"), col.names = TRUE, row.names = FALSE)
  
  write.csv(tf_targets_all_info, paste0("tfs_roma_modiced_for_publi/targets_of_TFs/last_20/",overdispersed_modules_arranged_jo_last_20$modules[i],".csv"), col.names = TRUE, row.names = FALSE)
}



# cytos of interest

for (i in 1:nrow(overdispersed_modules_arranged_cyto_interest)){
  
  tf_targets_all_info <- GetTopContrib(get_roma_results, Selected = which(rownames(get_roma_results_module_matrix) == overdispersed_modules_arranged_cyto_interest$modules[i]),
                                       nGenes = 3000, OrderType = "Abs", Mode = "Wei", Plot = TRUE)$Table %>% dplyr::select(c("Gene","Weight")) %>% as.data.frame() %>% dplyr::rename("target" = "Gene") %>% dplyr::arrange(desc(abs(Weight))) %>% dplyr::mutate(cyto_of_interest = if_else(target %in% c("Il6","Ccl4","Ccl5","Ccl20","Ccl9"), "Yes", "No"))
  rownames(tf_targets_all_info) <- NULL
  
  if (overdispersed_modules_arranged_cyto_interest_new_final[i] %in% collectri_df$source){
    collectri_df_tf_tracked <- collectri_df %>% filter(source == overdispersed_modules_arranged_cyto_interest_new_final[i])
    tf_targets_all_info <- left_join(tf_targets_all_info, collectri_df_tf_tracked[,c("target","mor")], by = "target")
  }
  
  
  #write.csv(tf_targets_all_info, paste0("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/",overdispersed_modules_arranged_jo_cyto_interest$modules[i],".csv"), col.names = TRUE, row.names = FALSE)
  
  write.csv(tf_targets_all_info, paste0("tfs_roma_modiced_for_publi/targets_of_TFs (cytos of interest)/",overdispersed_modules_arranged_cyto_interest$modules[i],".csv"), col.names = TRUE, row.names = FALSE)
}






# Prdm9, Zfp661

prdm9_target_genes <- read.csv("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/PRDM9_TARGET_GENES.csv") %>% dplyr::mutate(TF_name = 'Prdm9') %>% dplyr::select(-"mor")

zfp661_target_genes <- read.csv("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/ZFP661_TARGET_GENES.csv") %>% dplyr::mutate(TF_name = 'Zfp661') %>% dplyr::select(-"mor")


prdm9_zfp661_target_genes <- inner_join(prdm9_target_genes, zfp661_target_genes, by="target")





tead2_target_genes <- read.csv("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/TEAD2_TARGET_GENES.csv") %>% dplyr::mutate(TF_name = 'Tead2') %>% dplyr::select(-"mor")

zfp449_target_genes <- read.csv("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/ZFP449_TARGET_GENES.csv") %>% dplyr::mutate(TF_name = 'Zfp449') %>% dplyr::select(-"mor")

dlx1_target_genes <- read.csv("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/DLX1_TARGET_GENES.csv") %>% dplyr::mutate(TF_name = 'Dlx1') %>% dplyr::select(-"mor")

hoxc9_target_genes <- read.csv("tf_expression_activation/Overdispersed TFs (at least one cyto of Jo)/TFs_target_genes/HOXC9_TARGET_GENES.csv") %>% dplyr::mutate(TF_name = 'Dlx1') %>% dplyr::select(-"mor")





tead2_zfp449_target_genes <- inner_join(tead2_target_genes, zfp449_target_genes, by="target")

inner_join(prdm9_zfp661_target_genes, tead2_zfp449_target_genes, by="target")

PlotGeneWeight(RomaData = get_roma_results, PlotGenes = 100, LogExpression = FALSE, Selected = which(rownames(get_roma_results_module_matrix) == "PRDM9_TARGET_GENES"))



