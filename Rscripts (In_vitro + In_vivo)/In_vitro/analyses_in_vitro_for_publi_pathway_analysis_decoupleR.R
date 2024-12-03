library(tidyverse)
library(edgeR)
library(limma)
library(decoupleR)
library(pheatmap)



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
counts_for_progeny <- cbind(voom$E, voom$genes)
rownames(counts_for_progeny) <- make.names(counts_for_progeny$genes, unique=TRUE)
counts_for_progeny <- counts_for_progeny %>% dplyr::select(-"genes")


net <- get_progeny(organism = 'mouse', top = 500)


progeny_in_vitro <- run_mlm(mat=counts_for_progeny, net=net, .source='source', .target='target',
                            .mor='weight', minsize = 5)


progeny_in_vitro_mat_scaled_by_pathway <- progeny_in_vitro %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>% scale()


progeny_in_vitro_mat_scaled_by_sample <- progeny_in_vitro %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>% t() %>% scale()
 


progeny_in_vitro_mat_scaled_by_pathway_transposed <- progeny_in_vitro_mat_scaled_by_pathway %>% t() %>% as.data.frame()



palette_length = 100
my_color <- viridis::viridis(palette_length)



# Plot
#progeny_in_vitro_24h_matrix_scale_per_pathway <- progeny_in_vitro_mat_scaled_by_pathway %>% as.data.frame() %>%  dplyr::slice(which(grepl("24h", rownames(.)))) %>% as.matrix()







# 3h


progeny_activity_MTX_3h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("3h") & contains("MTX"))

progeny_activity_MTX_3h_mean <- rowMeans(progeny_activity_MTX_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_3h" = ".")


progeny_activity_OXA_3h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("3h") & contains("OXA"))

progeny_activity_OXA_3h_mean <- rowMeans(progeny_activity_OXA_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_3h" = ".")



progeny_activity_CIS_3h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("3h") & contains("CIS"))

progeny_activity_CIS_3h_mean <- rowMeans(progeny_activity_CIS_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_3h" = ".")


progeny_activity_Control_3h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("3h") & contains("Control"))

progeny_activity_Control_3h_mean <- rowMeans(progeny_activity_Control_3h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_3h" = ".")



progeny_activity_3h_mean <- cbind(progeny_activity_MTX_3h_mean, progeny_activity_OXA_3h_mean, progeny_activity_CIS_3h_mean, progeny_activity_Control_3h_mean)



# 6h


progeny_activity_MTX_6h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("6h") & contains("MTX"))

progeny_activity_MTX_6h_mean <- rowMeans(progeny_activity_MTX_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_6h" = ".")


progeny_activity_OXA_6h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("6h") & contains("OXA"))

progeny_activity_OXA_6h_mean <- rowMeans(progeny_activity_OXA_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_6h" = ".")



progeny_activity_CIS_6h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("6h") & contains("CIS"))

progeny_activity_CIS_6h_mean <- rowMeans(progeny_activity_CIS_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_6h" = ".")


progeny_activity_Control_6h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("6h") & contains("Control"))

progeny_activity_Control_6h_mean <- rowMeans(progeny_activity_Control_6h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_6h" = ".")



progeny_activity_6h_mean <- cbind(progeny_activity_MTX_6h_mean, progeny_activity_OXA_6h_mean, progeny_activity_CIS_6h_mean, progeny_activity_Control_6h_mean)






# 12h


progeny_activity_MTX_12h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("12h") & contains("MTX"))

progeny_activity_MTX_12h_mean <- rowMeans(progeny_activity_MTX_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_12h" = ".")


progeny_activity_OXA_12h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("12h") & contains("OXA"))

progeny_activity_OXA_12h_mean <- rowMeans(progeny_activity_OXA_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_12h" = ".")



progeny_activity_CIS_12h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("12h") & contains("CIS"))

progeny_activity_CIS_12h_mean <- rowMeans(progeny_activity_CIS_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_12h" = ".")


progeny_activity_Control_12h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("12h") & contains("Control"))

progeny_activity_Control_12h_mean <- rowMeans(progeny_activity_Control_12h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_12h" = ".")



progeny_activity_12h_mean <- cbind(progeny_activity_MTX_12h_mean, progeny_activity_OXA_12h_mean, progeny_activity_CIS_12h_mean, progeny_activity_Control_12h_mean)





# 18h


progeny_activity_MTX_18h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("18h") & contains("MTX"))

progeny_activity_MTX_18h_mean <- rowMeans(progeny_activity_MTX_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_18h" = ".")


progeny_activity_OXA_18h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("18h") & contains("OXA"))

progeny_activity_OXA_18h_mean <- rowMeans(progeny_activity_OXA_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_18h" = ".")



progeny_activity_CIS_18h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("18h") & contains("CIS"))

progeny_activity_CIS_18h_mean <- rowMeans(progeny_activity_CIS_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_18h" = ".")


progeny_activity_Control_18h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("18h") & contains("Control"))

progeny_activity_Control_18h_mean <- rowMeans(progeny_activity_Control_18h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_18h" = ".")



progeny_activity_18h_mean <- cbind(progeny_activity_MTX_18h_mean, progeny_activity_OXA_18h_mean, progeny_activity_CIS_18h_mean, progeny_activity_Control_18h_mean)




# 24h


progeny_activity_MTX_24h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("24h") & contains("MTX"))

progeny_activity_MTX_24h_mean <- rowMeans(progeny_activity_MTX_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_24h" = ".")


progeny_activity_OXA_24h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("24h") & contains("OXA"))

progeny_activity_OXA_24h_mean <- rowMeans(progeny_activity_OXA_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("OXA_24h" = ".")



progeny_activity_CIS_24h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("24h") & contains("CIS"))

progeny_activity_CIS_24h_mean <- rowMeans(progeny_activity_CIS_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("CIS_24h" = ".")


progeny_activity_Control_24h <- progeny_in_vitro_mat_scaled_by_pathway_transposed %>% select(contains("24h") & contains("Control"))

progeny_activity_Control_24h_mean <- rowMeans(progeny_activity_Control_24h[,1:3]) %>% as.data.frame() %>% dplyr::rename("Control_24h" = ".")



progeny_activity_24h_mean <- cbind(progeny_activity_MTX_24h_mean, progeny_activity_OXA_24h_mean, progeny_activity_CIS_24h_mean, progeny_activity_Control_24h_mean)





progeny_activity_all_timepoints <- cbind(progeny_activity_3h_mean, progeny_activity_6h_mean, progeny_activity_12h_mean, progeny_activity_18h_mean, progeny_activity_24h_mean) %>% as.matrix()



annotation_columns_progeny_mean <- data.frame(treatment = rep(rep(c("MTX","OXA","CIS", "Control"),c(1,1,1,1)), 5))
rownames(annotation_columns_progeny_mean)[1:4] <- colnames(progeny_activity_all_timepoints)[1:4]
rownames(annotation_columns_progeny_mean)[5:8] <- colnames(progeny_activity_all_timepoints)[5:8]
rownames(annotation_columns_progeny_mean)[9:12] <- colnames(progeny_activity_all_timepoints)[9:12]
rownames(annotation_columns_progeny_mean)[13:16] <- colnames(progeny_activity_all_timepoints)[13:16]
rownames(annotation_columns_progeny_mean)[17:20] <- colnames(progeny_activity_all_timepoints)[17:20]



annot_colors_progeny <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen", Control="darkgrey"))




pathway_activity_progeny_heatmap_expression <- ComplexHeatmap::pheatmap(progeny_activity_all_timepoints, color = my_color, annotation_col = annotation_columns_progeny_mean, annotation_colors = annot_colors_progeny, fontsize_row = 11.5, fontsize_col = 11.5, cluster_rows = TRUE, cluster_cols = TRUE, main="Pathway activity (PROGENy) in all the conditions - gene expression values", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white") %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")




# With t-values




# MTX vs Control

t_values_diff_genes_MTX_Control_3h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/3h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_3h) <- make.names(t_values_diff_genes_MTX_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_3h <- t_values_diff_genes_MTX_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_6h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/6h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_6h) <- make.names(t_values_diff_genes_MTX_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_6h <- t_values_diff_genes_MTX_Control_6h %>% dplyr::select(-"gene_name")



t_values_diff_genes_MTX_Control_12h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/12h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_12h) <- make.names(t_values_diff_genes_MTX_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_12h <- t_values_diff_genes_MTX_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_18h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/18h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_18h) <- make.names(t_values_diff_genes_MTX_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_18h <- t_values_diff_genes_MTX_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_24h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/24h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_24h) <- make.names(t_values_diff_genes_MTX_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_24h <- t_values_diff_genes_MTX_Control_24h %>% dplyr::select(-"gene_name")




# OXA vs Control

t_values_diff_genes_OXA_Control_3h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/3h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_3h) <- make.names(t_values_diff_genes_OXA_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_3h <- t_values_diff_genes_OXA_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_6h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/6h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_6h) <- make.names(t_values_diff_genes_OXA_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_6h <- t_values_diff_genes_OXA_Control_6h %>% dplyr::select(-"gene_name")



t_values_diff_genes_OXA_Control_12h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/12h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_12h) <- make.names(t_values_diff_genes_OXA_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_12h <- t_values_diff_genes_OXA_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_18h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/18h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_18h) <- make.names(t_values_diff_genes_OXA_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_18h <- t_values_diff_genes_OXA_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_24h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/24h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_24h) <- make.names(t_values_diff_genes_OXA_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_24h <- t_values_diff_genes_OXA_Control_24h %>% dplyr::select(-"gene_name")





# CIS vs Control

t_values_diff_genes_CIS_Control_3h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/3h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_3h) <- make.names(t_values_diff_genes_CIS_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_3h <- t_values_diff_genes_CIS_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_6h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/6h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_6h) <- make.names(t_values_diff_genes_CIS_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_6h <- t_values_diff_genes_CIS_Control_6h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_12h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/12h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_12h) <- make.names(t_values_diff_genes_CIS_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_12h <- t_values_diff_genes_CIS_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_18h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/18h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_18h) <- make.names(t_values_diff_genes_CIS_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_18h <- t_values_diff_genes_CIS_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_24h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/24h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_24h) <- make.names(t_values_diff_genes_CIS_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_24h <- t_values_diff_genes_CIS_Control_24h %>% dplyr::select(-"gene_name")







## run decoupleR - MTX



### 3h

run_decoupler_MTX_3h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_3h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_3h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_3h, seed = 2)


# Get results
run_decoupler_MTX_3h_consensus_pivot <- run_decoupler_MTX_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_MTX_3h <- t(run_decoupler_MTX_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_3h" = "t")




### 6h

run_decoupler_MTX_6h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_6h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_6h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_6h, seed = 2)


# Get results
run_decoupler_MTX_6h_consensus_pivot <- run_decoupler_MTX_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_6h <- t(run_decoupler_MTX_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_6h" = "t")




### 12h

run_decoupler_MTX_12h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_12h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_12h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_12h, seed = 2)


# Get results
run_decoupler_MTX_12h_consensus_pivot <- run_decoupler_MTX_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_12h <- t(run_decoupler_MTX_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_12h" = "t")



### 18h

run_decoupler_MTX_18h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_18h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_18h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_18h, seed = 2)


# Get results
run_decoupler_MTX_18h_consensus_pivot <- run_decoupler_MTX_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_18h <- t(run_decoupler_MTX_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_18h" = "t")



### 24h

run_decoupler_MTX_24h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_24h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_24h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_24h, seed = 2)


# Get results
run_decoupler_MTX_24h_consensus_pivot <- run_decoupler_MTX_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_MTX_24h <- t(run_decoupler_MTX_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_24h" = "t")




## run decoupleR - OXA


### 3h

run_decoupler_OXA_3h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_3h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_3h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_3h, seed = 2)


# Get results
run_decoupler_OXA_3h_consensus_pivot <- run_decoupler_OXA_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_OXA_3h <- t(run_decoupler_OXA_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_3h" = "t")




### 6h

run_decoupler_OXA_6h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_6h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_6h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_6h, seed = 2)


# Get results
run_decoupler_OXA_6h_consensus_pivot <- run_decoupler_OXA_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_OXA_6h <- t(run_decoupler_OXA_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_6h" = "t")




### 12h

run_decoupler_OXA_12h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_12h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_12h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_12h, seed = 2)


# Get results

run_decoupler_OXA_12h_consensus_pivot <- run_decoupler_OXA_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_OXA_12h <- t(run_decoupler_OXA_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_12h" = "t")




### 18h

run_decoupler_OXA_18h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_18h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_18h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_18h, seed = 2)


# Get results
run_decoupler_OXA_18h_consensus_pivot <- run_decoupler_OXA_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

tf_target_cytos_interest_OXA_18h <- t(run_decoupler_OXA_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_18h" = "t")



### 24h

run_decoupler_OXA_24h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_24h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_24h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_24h, seed = 2)


# Get results
run_decoupler_OXA_24h_consensus_pivot <- run_decoupler_OXA_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_OXA_24h <- t(run_decoupler_OXA_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_24h" = "t")





## run decoupleR - CIS


### 3h

run_decoupler_CIS_3h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_3h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_3h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_3h, seed = 2)


# Get results
run_decoupler_CIS_3h_consensus_pivot <- run_decoupler_CIS_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_CIS_3h <- t(run_decoupler_CIS_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_3h" = "t")




### 6h


run_decoupler_CIS_6h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_6h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_6h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_6h, seed = 2)


# Get results
run_decoupler_CIS_6h_consensus_pivot <- run_decoupler_CIS_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_CIS_6h <- t(run_decoupler_CIS_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_6h" = "t")




### 12h

run_decoupler_CIS_12h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_12h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_12h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_12h, seed = 2)


# Get results
run_decoupler_CIS_12h_consensus_pivot <- run_decoupler_CIS_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_CIS_12h <- t(run_decoupler_CIS_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_12h" = "t")




### 18h

run_decoupler_CIS_18h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_18h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_18h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_18h, seed = 2)


# Get results
run_decoupler_CIS_18h_consensus_pivot <- run_decoupler_CIS_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_CIS_18h <- t(run_decoupler_CIS_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_18h" = "t")



### 24h

run_decoupler_CIS_24h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_24h, net = net, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_24h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_24h, seed = 2)


# Get results
run_decoupler_CIS_24h_consensus_pivot <- run_decoupler_CIS_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() 


tf_target_cytos_interest_CIS_24h <- t(run_decoupler_CIS_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_24h" = "t")





# Join 

# MTX
tf_target_cytos_interest_MTX_3h_6h <- inner_join(tf_target_cytos_interest_MTX_3h, tf_target_cytos_interest_MTX_6h, by="TFs")
tf_target_cytos_interest_MTX_3h_6h_12h <- inner_join(tf_target_cytos_interest_MTX_3h_6h, tf_target_cytos_interest_MTX_12h, by="TFs")
tf_target_cytos_interest_MTX_3h_6h_12h_18h <- inner_join(tf_target_cytos_interest_MTX_3h_6h_12h, tf_target_cytos_interest_MTX_18h, by="TFs")
tf_target_cytos_interest_MTX_all_timepoints <- inner_join(tf_target_cytos_interest_MTX_3h_6h_12h_18h, tf_target_cytos_interest_MTX_24h, by="TFs")


# OXA
tf_target_cytos_interest_OXA_3h_6h <- inner_join(tf_target_cytos_interest_OXA_3h, tf_target_cytos_interest_OXA_6h, by="TFs")
tf_target_cytos_interest_OXA_3h_6h_12h <- inner_join(tf_target_cytos_interest_OXA_3h_6h, tf_target_cytos_interest_OXA_12h, by="TFs")
tf_target_cytos_interest_OXA_3h_6h_12h_18h <- inner_join(tf_target_cytos_interest_OXA_3h_6h_12h, tf_target_cytos_interest_OXA_18h, by="TFs")
tf_target_cytos_interest_OXA_all_timepoints <- inner_join(tf_target_cytos_interest_OXA_3h_6h_12h_18h, tf_target_cytos_interest_OXA_24h, by="TFs")


# CIS
tf_target_cytos_interest_CIS_3h_6h <- inner_join(tf_target_cytos_interest_CIS_3h, tf_target_cytos_interest_CIS_6h, by="TFs")
tf_target_cytos_interest_CIS_3h_6h_12h <- inner_join(tf_target_cytos_interest_CIS_3h_6h, tf_target_cytos_interest_CIS_12h, by="TFs")
tf_target_cytos_interest_CIS_3h_6h_12h_18h <- inner_join(tf_target_cytos_interest_CIS_3h_6h_12h, tf_target_cytos_interest_CIS_18h, by="TFs")
tf_target_cytos_interest_CIS_all_timepoints <- inner_join(tf_target_cytos_interest_CIS_3h_6h_12h_18h, tf_target_cytos_interest_CIS_24h, by="TFs")


# All conditions
tf_target_cytos_interest_MTX_OXA_all_timepoints <- inner_join(tf_target_cytos_interest_MTX_all_timepoints, 
                                                              tf_target_cytos_interest_OXA_all_timepoints, by="TFs")


tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints <- inner_join(tf_target_cytos_interest_MTX_OXA_all_timepoints, tf_target_cytos_interest_CIS_all_timepoints, by = "TFs") %>% column_to_rownames("TFs")

pathway_progeny_results_t_values <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints %>% t() %>% scale() %>% t()







palette_length <- 100 

my_color <- viridis::viridis(palette_length)


annotation_columns_treatment <- data.frame(treatment = rep(c("MTX","OXA","CIS"),c(5,5,5)))
rownames(annotation_columns_treatment)[1:5] <- colnames(pathway_progeny_results)[1:5]
rownames(annotation_columns_treatment)[6:10] <- colnames(pathway_progeny_results)[6:10]
rownames(annotation_columns_treatment)[11:15] <- colnames(pathway_progeny_results)[11:15]

annot_colors_all <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"))



pathway_activity_progeny_heatmap_t_values <- ComplexHeatmap::pheatmap(pathway_progeny_results_t_values, color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11.5, fontsize_col = 11.5, cluster_rows = TRUE, cluster_cols = TRUE, main="Pathway activity (PROGENy) in all the conditions - t-values", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white") %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



pathway_activity_progeny_heatmap_expression <- ComplexHeatmap::pheatmap(progeny_activity_all_timepoints,color = my_color, annotation_col = annotation_columns_progeny_mean, annotation_colors = annot_colors_progeny, fontsize_row = 11.5, fontsize_col = 11.5, cluster_rows = TRUE, cluster_cols = TRUE, main="Pathway activity (PROGENy) in all the conditions - gene expression values", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white") %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")





patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(pathway_activity_progeny_heatmap_expression)), grid.grabExpr(ComplexHeatmap::draw(pathway_activity_progeny_heatmap_t_values)))) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))





# Targets


# MTX vs Control

t_values_diff_genes_MTX_Control_3h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/3h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_3h) <- make.names(t_values_diff_genes_MTX_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_3h <- t_values_diff_genes_MTX_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_6h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/6h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_6h) <- make.names(t_values_diff_genes_MTX_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_6h <- t_values_diff_genes_MTX_Control_6h %>% dplyr::select(-"gene_name")



t_values_diff_genes_MTX_Control_12h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/12h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_12h) <- make.names(t_values_diff_genes_MTX_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_12h <- t_values_diff_genes_MTX_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_18h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/18h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_18h) <- make.names(t_values_diff_genes_MTX_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_18h <- t_values_diff_genes_MTX_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_24h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_MTX_vs_Control/24h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_24h) <- make.names(t_values_diff_genes_MTX_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_24h <- t_values_diff_genes_MTX_Control_24h %>% dplyr::select(-"gene_name")




# OXA vs Control

t_values_diff_genes_OXA_Control_3h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/3h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_3h) <- make.names(t_values_diff_genes_OXA_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_3h <- t_values_diff_genes_OXA_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_6h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/6h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_6h) <- make.names(t_values_diff_genes_OXA_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_6h <- t_values_diff_genes_OXA_Control_6h %>% dplyr::select(-"gene_name")



t_values_diff_genes_OXA_Control_12h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/12h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_12h) <- make.names(t_values_diff_genes_OXA_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_12h <- t_values_diff_genes_OXA_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_18h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/18h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_18h) <- make.names(t_values_diff_genes_OXA_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_18h <- t_values_diff_genes_OXA_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_24h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_OXA_vs_Control/24h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_24h) <- make.names(t_values_diff_genes_OXA_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_24h <- t_values_diff_genes_OXA_Control_24h %>% dplyr::select(-"gene_name")





# CIS vs Control

t_values_diff_genes_CIS_Control_3h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/3h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_3h) <- make.names(t_values_diff_genes_CIS_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_3h <- t_values_diff_genes_CIS_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_6h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/6h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_6h) <- make.names(t_values_diff_genes_CIS_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_6h <- t_values_diff_genes_CIS_Control_6h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_12h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/12h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_12h) <- make.names(t_values_diff_genes_CIS_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_12h <- t_values_diff_genes_CIS_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_18h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/18h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_18h) <- make.names(t_values_diff_genes_CIS_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_18h <- t_values_diff_genes_CIS_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_24h <- read.csv("~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/Analysis_CIS_vs_Control/24h/diff_genes.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_24h) <- make.names(t_values_diff_genes_CIS_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_24h <- t_values_diff_genes_CIS_Control_24h %>% dplyr::select(-"gene_name")





# MTX vs Control
t_values_diff_genes_MTX_Control_3h_df <-  t_values_diff_genes_MTX_Control_3h %>% rownames_to_column("target")
t_values_diff_genes_MTX_Control_6h_df <-  t_values_diff_genes_MTX_Control_6h %>% rownames_to_column("target")
t_values_diff_genes_MTX_Control_12h_df <-  t_values_diff_genes_MTX_Control_12h %>% rownames_to_column("target")
t_values_diff_genes_MTX_Control_18h_df <-  t_values_diff_genes_MTX_Control_18h %>% rownames_to_column("target")
t_values_diff_genes_MTX_Control_24h_df <-  t_values_diff_genes_MTX_Control_24h %>% rownames_to_column("target")

# OXA vs Control
t_values_diff_genes_OXA_Control_3h_df <-  t_values_diff_genes_OXA_Control_3h %>% rownames_to_column("target")
t_values_diff_genes_OXA_Control_6h_df <-  t_values_diff_genes_OXA_Control_6h %>% rownames_to_column("target")
t_values_diff_genes_OXA_Control_12h_df <-  t_values_diff_genes_OXA_Control_12h %>% rownames_to_column("target")
t_values_diff_genes_OXA_Control_18h_df <-  t_values_diff_genes_OXA_Control_18h %>% rownames_to_column("target")
t_values_diff_genes_OXA_Control_24h_df <-  t_values_diff_genes_OXA_Control_24h %>% rownames_to_column("target")


# CIS vs Control
t_values_diff_genes_CIS_Control_3h_df <-  t_values_diff_genes_CIS_Control_3h %>% rownames_to_column("target")
t_values_diff_genes_CIS_Control_6h_df <-  t_values_diff_genes_CIS_Control_6h %>% rownames_to_column("target")
t_values_diff_genes_CIS_Control_12h_df <-  t_values_diff_genes_CIS_Control_12h %>% rownames_to_column("target")
t_values_diff_genes_CIS_Control_18h_df <-  t_values_diff_genes_CIS_Control_18h %>% rownames_to_column("target")
t_values_diff_genes_CIS_Control_24h_df <-  t_values_diff_genes_CIS_Control_24h %>% rownames_to_column("target")



net_for_targets <- net %>% mutate(target_type = if_else(weight > 0, "activator","inhibitor"))

for (i in 1:nrow(pathway_progeny_results_t_values)){
  
  
  
  # MTX
  
  
  progeny_targets_MTX_3h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  progeny_targets_MTX_3h_activators <- progeny_targets_MTX_3h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_3h_inhibitors <- progeny_targets_MTX_3h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_3h_new <- rbind(progeny_targets_MTX_3h_activators, progeny_targets_MTX_3h_inhibitors)
  write.csv(progeny_targets_MTX_3h_new, paste0("figures_for_publication/In_vitro/PROGENy/MTX vs Control/3h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_MTX_6h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  progeny_targets_MTX_6h_activators <- progeny_targets_MTX_6h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_6h_inhibitors <- progeny_targets_MTX_6h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_6h_new <- rbind(progeny_targets_MTX_6h_activators, progeny_targets_MTX_6h_inhibitors)
  write.csv(progeny_targets_MTX_6h_new, paste0("figures_for_publication/In_vitro/PROGENy/MTX vs Control/6h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_MTX_12h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  progeny_targets_MTX_12h_activators <- progeny_targets_MTX_12h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_12h_inhibitors <- progeny_targets_MTX_12h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_12h_new <- rbind(progeny_targets_MTX_12h_activators, progeny_targets_MTX_12h_inhibitors)
  write.csv(progeny_targets_MTX_12h_new, paste0("figures_for_publication/In_vitro/PROGENy/MTX vs Control/12h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  progeny_targets_MTX_18h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  progeny_targets_MTX_18h_activators <- progeny_targets_MTX_18h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_18h_inhibitors <- progeny_targets_MTX_18h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_18h_new <- rbind(progeny_targets_MTX_18h_activators, progeny_targets_MTX_18h_inhibitors)
  write.csv(progeny_targets_MTX_18h_new, paste0("figures_for_publication/In_vitro/PROGENy/MTX vs Control/18h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_MTX_24h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  progeny_targets_MTX_24h_activators <- progeny_targets_MTX_24h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_24h_inhibitors <- progeny_targets_MTX_24h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_MTX_24h_new <- rbind(progeny_targets_MTX_24h_activators, progeny_targets_MTX_24h_inhibitors)
  write.csv(progeny_targets_MTX_24h_new, paste0("figures_for_publication/In_vitro/PROGENy/MTX vs Control/24h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  # OXA
  
  
  progeny_targets_OXA_3h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  progeny_targets_OXA_3h_activators <- progeny_targets_OXA_3h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_3h_inhibitors <- progeny_targets_OXA_3h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_3h_new <- rbind(progeny_targets_OXA_3h_activators, progeny_targets_OXA_3h_inhibitors)
  write.csv(progeny_targets_OXA_3h_new, paste0("figures_for_publication/In_vitro/PROGENy/OXA vs Control/3h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_OXA_6h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  progeny_targets_OXA_6h_activators <- progeny_targets_OXA_6h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_6h_inhibitors <- progeny_targets_OXA_6h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_6h_new <- rbind(progeny_targets_OXA_6h_activators, progeny_targets_OXA_6h_inhibitors)
  write.csv(progeny_targets_OXA_6h_new, paste0("figures_for_publication/In_vitro/PROGENy/OXA vs Control/6h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_OXA_12h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  progeny_targets_OXA_12h_activators <- progeny_targets_OXA_12h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_12h_inhibitors <- progeny_targets_OXA_12h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_12h_new <- rbind(progeny_targets_OXA_12h_activators, progeny_targets_OXA_12h_inhibitors)
  write.csv(progeny_targets_OXA_12h_new, paste0("figures_for_publication/In_vitro/PROGENy/OXA vs Control/12h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  progeny_targets_OXA_18h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  progeny_targets_OXA_18h_activators <- progeny_targets_OXA_18h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_18h_inhibitors <- progeny_targets_OXA_18h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_18h_new <- rbind(progeny_targets_OXA_18h_activators, progeny_targets_OXA_18h_inhibitors)
  write.csv(progeny_targets_OXA_18h_new, paste0("figures_for_publication/In_vitro/PROGENy/OXA vs Control/18h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_OXA_24h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  progeny_targets_OXA_24h_activators <- progeny_targets_OXA_24h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_24h_inhibitors <- progeny_targets_OXA_24h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_OXA_24h_new <- rbind(progeny_targets_OXA_24h_activators, progeny_targets_OXA_24h_inhibitors)
  write.csv(progeny_targets_OXA_24h_new, paste0("figures_for_publication/In_vitro/PROGENy/OXA vs Control/24h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  # CIS
  
  
  progeny_targets_CIS_3h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  progeny_targets_CIS_3h_activators <- progeny_targets_CIS_3h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_3h_inhibitors <- progeny_targets_CIS_3h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_3h_new <- rbind(progeny_targets_CIS_3h_activators, progeny_targets_CIS_3h_inhibitors)
  write.csv(progeny_targets_CIS_3h_new, paste0("figures_for_publication/In_vitro/PROGENy/CIS vs Control/3h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_CIS_6h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  progeny_targets_CIS_6h_activators <- progeny_targets_CIS_6h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_6h_inhibitors <- progeny_targets_CIS_6h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_6h_new <- rbind(progeny_targets_CIS_6h_activators, progeny_targets_CIS_6h_inhibitors)
  write.csv(progeny_targets_CIS_6h_new, paste0("figures_for_publication/In_vitro/PROGENy/CIS vs Control/6h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_CIS_12h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  progeny_targets_CIS_12h_activators <- progeny_targets_CIS_12h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_12h_inhibitors <- progeny_targets_CIS_12h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_12h_new <- rbind(progeny_targets_CIS_12h_activators, progeny_targets_CIS_12h_inhibitors)
  write.csv(progeny_targets_CIS_12h_new, paste0("figures_for_publication/In_vitro/PROGENy/CIS vs Control/12h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  progeny_targets_CIS_18h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  progeny_targets_CIS_18h_activators <- progeny_targets_CIS_18h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_18h_inhibitors <- progeny_targets_CIS_18h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_18h_new <- rbind(progeny_targets_CIS_18h_activators, progeny_targets_CIS_18h_inhibitors)
  write.csv(progeny_targets_CIS_18h_new, paste0("figures_for_publication/In_vitro/PROGENy/CIS vs Control/18h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  progeny_targets_CIS_24h <- net_for_targets %>% filter(source == rownames(pathway_progeny_results_t_values)[i]) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  progeny_targets_CIS_24h_activators <- progeny_targets_CIS_24h %>% filter(target_type == "activator") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_24h_inhibitors <- progeny_targets_CIS_24h %>% filter(target_type == "inhibitor") %>% dplyr::arrange(desc(abs(t)))
  progeny_targets_CIS_24h_new <- rbind(progeny_targets_CIS_24h_activators, progeny_targets_CIS_24h_inhibitors)
  write.csv(progeny_targets_CIS_24h_new, paste0("figures_for_publication/In_vitro/PROGENy/CIS vs Control/24h/",rownames(pathway_progeny_results_t_values)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}




