library(decoupleR)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(scales)


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





# Get collectri with all the tfs
collectri_df <- get_collectri(organism="mouse", split_complexes = FALSE) %>% as.data.frame()





# TFs that include at least one of the cytokines of interest (Il6, Ccl4, Ccl5, Ccl20, Ccl9)



collectri_df_Ccl5 <- collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_Ccl5 = any(target == "Ccl5")) %>% as.data.frame() %>% filter(target_Ccl5 == "TRUE")
collectri_df_for_Ccl5 <- collectri_df %>% dplyr::filter(source %in% collectri_df_Ccl5$source)



collectri_df_Il6 <- collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_Il6 = any(target == "Il6")) %>% as.data.frame() %>% filter(target_Il6 == "TRUE")
collectri_df_for_Il6 <- collectri_df %>% dplyr::filter(source %in% collectri_df_Il6$source)



collectri_df_Ccl4 <- collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_Ccl4 = any(target == "Ccl4")) %>% as.data.frame() %>% filter(target_Ccl4 == "TRUE")
collectri_df_for_Ccl4 <- collectri_df %>% dplyr::filter(source %in% collectri_df_Ccl4$source)



collectri_df_Ccl20 <- collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_Ccl20 = any(target == "Ccl20")) %>% as.data.frame() %>% filter(target_Ccl20 == "TRUE")
collectri_df_for_Ccl20 <- collectri_df %>% dplyr::filter(source %in% collectri_df_Ccl20$source)



collectri_df_Ccl9 <- collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_Ccl9 = any(target == "Ccl9")) %>% as.data.frame() %>% filter(target_Ccl9 == "TRUE")
collectri_df_for_Ccl9 <- collectri_df %>% dplyr::filter(source %in% collectri_df_Ccl9$source)


# tfs that target at least one cyto of interest
collectri_df_at_least_one_cytokine <-  collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_cyto_interest = any(target == "Ccl4" | target == "Ccl5" | target == "Ccl20" | target == "Ccl9" | target == "Il6")) %>% as.data.frame() %>% filter(target_cyto_interest == "TRUE")
collectri_df_target_at_least_one_df <- collectri_df %>% filter(source %in% collectri_df_at_least_one_cytokine$source)


# tf that target all cyto of interest
collectri_target_Il6_Ccl4 <- inner_join(collectri_df_Il6, collectri_df_Ccl4, by="source")
collectri_target_Il6_Ccl4_Ccl5 <- inner_join(collectri_target_Il6_Ccl4, collectri_df_Ccl5, by="source")
collectri_target_Il6_Ccl4_Ccl5_Ccl20 <- inner_join(collectri_target_Il6_Ccl4_Ccl5, collectri_df_Ccl20, by="source")












# TFs that target cytos of interest



## run decoupleR - MTX


### 3h

run_decoupler_MTX_3h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_3h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_3h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_3h, seed = 2)


# Get results
run_decoupler_MTX_3h_consensus_pivot <- run_decoupler_MTX_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_3h <- t(run_decoupler_MTX_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_3h" = "t")




### 6h

run_decoupler_MTX_6h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_6h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_6h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_6h, seed = 2)


# Get results
run_decoupler_MTX_6h_consensus_pivot <- run_decoupler_MTX_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_6h <- t(run_decoupler_MTX_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_6h" = "t")




### 12h

run_decoupler_MTX_12h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_12h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_12h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_12h, seed = 2)


# Get results
run_decoupler_MTX_12h_consensus_pivot <- run_decoupler_MTX_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_12h <- t(run_decoupler_MTX_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_12h" = "t")



### 18h

run_decoupler_MTX_18h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_18h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_18h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_18h, seed = 2)


# Get results
run_decoupler_MTX_18h_consensus_pivot <- run_decoupler_MTX_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_MTX_18h <- t(run_decoupler_MTX_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_18h" = "t")



### 24h

run_decoupler_MTX_24h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_24h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
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

run_decoupler_OXA_3h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_3h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_3h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_3h, seed = 2)


# Get results
run_decoupler_OXA_3h_consensus_pivot <- run_decoupler_OXA_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_OXA_3h <- t(run_decoupler_OXA_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_3h" = "t")




### 6h

run_decoupler_OXA_6h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_6h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_6h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_6h, seed = 2)


# Get results
run_decoupler_OXA_6h_consensus_pivot <- run_decoupler_OXA_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_OXA_6h <- t(run_decoupler_OXA_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_6h" = "t")




### 12h

run_decoupler_OXA_12h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_12h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_12h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_12h, seed = 2)


# Get results

run_decoupler_OXA_12h_consensus_pivot <- run_decoupler_OXA_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_OXA_12h <- t(run_decoupler_OXA_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_12h" = "t")




### 18h

run_decoupler_OXA_18h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_18h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_18h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_18h, seed = 2)


# Get results
run_decoupler_OXA_18h_consensus_pivot <- run_decoupler_OXA_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_OXA_18h <- t(run_decoupler_OXA_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_18h" = "t")



### 24h

run_decoupler_OXA_24h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_24h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
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

run_decoupler_CIS_3h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_3h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_3h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_3h, seed = 2)


# Get results
run_decoupler_CIS_3h_consensus_pivot <- run_decoupler_CIS_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_CIS_3h <- t(run_decoupler_CIS_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_3h" = "t")




### 6h


run_decoupler_CIS_6h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_6h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_6h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_6h, seed = 2)


# Get results
run_decoupler_CIS_6h_consensus_pivot <- run_decoupler_CIS_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_CIS_6h <- t(run_decoupler_CIS_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_6h" = "t")




### 12h

run_decoupler_CIS_12h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_12h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_12h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_12h, seed = 2)


# Get results
run_decoupler_CIS_12h_consensus_pivot <- run_decoupler_CIS_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_CIS_12h <- t(run_decoupler_CIS_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_12h" = "t")




### 18h

run_decoupler_CIS_18h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_18h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_18h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_18h, seed = 2)


# Get results
run_decoupler_CIS_18h_consensus_pivot <- run_decoupler_CIS_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_interest_CIS_18h <- t(run_decoupler_CIS_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_18h" = "t")



### 24h

run_decoupler_CIS_24h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_24h, net = collectri_df_target_at_least_one_df, .source = "source", .target = "target", consensus_score = FALSE)
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


write.csv(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints, "tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints.csv")



# scaled by TFs

tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tf <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints %>% t() %>% scale() %>% t() %>% as.data.frame()

tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tf %>% rowwise() %>%
  mutate(standard_deviation = c_across(everything()) %>% sd()) %>%
  ungroup() %>% as.data.frame()

tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final <- cbind(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tf, tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd$standard_deviation) %>% dplyr::rename("sd" = "tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd$standard_deviation") %>% dplyr::arrange(desc(sd)) %>% as.matrix()





# Visualization




palette_length <- 100 

my_color <- viridis::viridis(palette_length)


annotation_columns_treatment <- data.frame(treatment = rep(c("MTX","OXA","CIS"),c(5,5,5)))
rownames(annotation_columns_treatment)[1:5] <- colnames(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints)[1:5]
rownames(annotation_columns_treatment)[6:10] <- colnames(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints)[6:10]
rownames(annotation_columns_treatment)[11:15] <- colnames(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints)[11:15]

annot_colors_all <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"))




tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_sd <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints %>% rowwise() %>%
  mutate(standard_deviation = c_across(everything()) %>% sd()) %>%
  ungroup() %>% as.data.frame()


tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final <- cbind(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints, tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_sd$standard_deviation) %>% dplyr::rename("sd" = "tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_sd$standard_deviation") %>% dplyr::arrange(desc(sd)) %>% as.matrix()


#tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled <- scales::rescale(as.matrix(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final[,-16]))

tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap <- scale(as.matrix(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final[,-16]))

tf_activity_target_al_least_one_cyto_interest_heatmap <- ComplexHeatmap::pheatmap(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap[1:50,], color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11.5, fontsize_col = 11.5, cluster_rows = TRUE, cluster_cols = TRUE, main="TF activity (targetting at least one cyto of interest) in all the conditions - top 50", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white") %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")





tf_activity_target_al_least_one_cytos_interest_heatmap_scaled_tf <- ComplexHeatmap::pheatmap(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16], color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, main="TF activity (targetting at least one cyto of interest) in all the conditions - top 50", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")





tf_activity_target_all_cytos_interest_heatmap_scaled <- ComplexHeatmap::pheatmap(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[collectri_target_Il6_Ccl4_Ccl5_Ccl20$source,-16], color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11.5, fontsize_col = 11.5, cluster_rows = TRUE, cluster_cols = TRUE, main="TF activity (targetting all the cytos of interest) in all the conditions", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 15) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")











# Snai1


# expression
tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap %>% as.data.frame() %>% filter(rownames(.) == "Snai1") %>% dplyr::select(contains("CIS"))


Snai1_expression_CIS <- counts_in_vitro_roma_all_timepoints %>% filter(rownames(.) == "Snai1") %>% select(contains("CIS"))


Snai1_expression_CIS_new <- Snai1_expression_CIS[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(Snai1_expression_CIS)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
Snai1_expression_CIS_new$name <- rep(rownames(Snai1_expression_CIS),5)
Snai1_expression_CIS_new$treatment <- rep("CIS",5)



# activity

activity_Snai1_CIS <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap %>% as.data.frame() %>% filter(rownames(.) == "Snai1") %>% dplyr::select(contains("CIS")) %>% as.matrix()
Snai1_activity_CIS_new <- Snai1_expression_CIS_new %>% dplyr::rename("Mean_activity" = "Mean_Expression")
bv_2 <- t(activity_Snai1_CIS) %>% as.data.frame()

Snai1_activity_CIS_new$Mean_activity <- bv_2$Snai1


# Plot

  # expression
  Snai1_expression_plot_CIS <- ggplot(Snai1_expression_CIS_new[1:5,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",Snai1_expression_CIS_new[1,"name"], " in CIS ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  Snai1_activity_plot_CIS <- ggplot(Snai1_activity_CIS_new[1:5,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",Snai1_activity_CIS_new[1,"name"]," in CIS ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(Snai1_expression_plot_CIS, Snai1_activity_plot_CIS) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))
  
  #ggsave(paste0("tf_expression_activation/Overdispersed TFs (without cyto of jo)/longitudinal_expression_activation_in_CIS_tf_",merged_tfs_interest[i,"name"],".png"))
  
  
  
  
  


  Hoxb13_expression_OXA <- counts_in_vitro_roma_all_timepoints %>% filter(rownames(.) == "Hoxb13") %>% select(contains("OXA"))
  
  
  Hoxb13_expression_OXA_new <- Hoxb13_expression_OXA[,c(1:3,7:9,13:15,4:6,10:12)] %>% as.data.frame() %>% t() %>% as.data.frame() %>%  dplyr::rename("Mean_Expression" = rownames(Hoxb13_expression_OXA)) %>% dplyr::mutate(timepoint = c(rep(3,3),rep(6,3),rep(12,3),rep(18,3),rep(24,3))) %>% group_by(timepoint) %>% summarise(Mean_Expression = mean(Mean_Expression)) %>% as.data.frame()
  Hoxb13_expression_OXA_new$name <- rep(rownames(Hoxb13_expression_OXA),5)
  Hoxb13_expression_OXA_new$treatment <- rep("OXA",5)
  
  
  
  # activity
  
  activity_Hoxb13_OXA <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap %>% as.data.frame() %>% filter(rownames(.) == "Hoxb13") %>% dplyr::select(contains("OXA")) %>% as.matrix()
  Hoxb13_activity_OXA_new <- Hoxb13_expression_OXA_new %>% dplyr::rename("Mean_activity" = "Mean_Expression")
  bv_2 <- t(activity_Hoxb13_OXA) %>% as.data.frame()
  
  Hoxb13_activity_OXA_new$Mean_activity <- bv_2$Hoxb13
  
  
  # Plot
  
  # expression
  Hoxb13_expression_plot_OXA <- ggplot(Hoxb13_expression_OXA_new[1:5,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_Expression)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal expression of ",Hoxb13_expression_OXA_new[1,"name"], " in OXA ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_Expression")
  
  
  
  # activation
  Hoxb13_activity_plot_OXA <- ggplot(Hoxb13_activity_OXA_new[1:5,], mapping = aes(x = factor(timepoint, levels=c("3","6","12","18","24")), y = Mean_activity)) +
    geom_point(size = 2, shape=21, aes(color = name, fill=name)) +
    ggrepel::geom_label_repel(aes(label=name, fill=name), label.size = 0.15, label.padding = 0.10) +
    geom_path(aes(color=name, group=name), size=1.1) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    scale_color_brewer(palette = "Set2", guide = guide_legend(override.aes = aes(label = "")), name = "TF_name") +
    ggtitle(paste0("Longitudinal activity score of ",Hoxb13_activity_OXA_new[1,"name"]," in OXA ")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold")) +
    theme(axis.title = element_text(face="bold")) +
    theme(axis.text = element_text(face="bold")) +
    labs(x = "Timepoint", y = "Mean_activity_score")
  
  
  wrap_plots(Hoxb13_expression_plot_OXA, Hoxb13_activity_plot_OXA) & patchwork::plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size=13))


  
  
  
  
write.csv(tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap, "figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/tf_targets_at_least_one_cytos.csv", row.names = TRUE, col.names = TRUE)


write.csv(tf_all_target_cytos_interest_top_50, "figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/tf_targets_at_least_one_cytos_top_50.csv", row.names = TRUE, col.names = TRUE)
  
  
  
  
  
  
  # Retrieve TFs with differentially expressed targets (+ with mode of reg) 
  
  
  tf_all_target_cytos_interest_top_50 <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap[1:50,]
  
  
  
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
  
  
  
  for (i in 1:nrow(tf_all_target_cytos_interest_top_50)){
    
    
    
    # MTX
    
    
    collectri_df_tf_interest_MTX_3h <-  collectri_df_target_at_least_one_df %>% filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
    collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors)
    write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/MTX vs Control/3h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    collectri_df_tf_interest_MTX_6h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
    collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors)
    write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/MTX vs Control/6h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    collectri_df_tf_interest_MTX_12h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
    collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors)
    write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/MTX vs Control/12h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    collectri_df_tf_interest_MTX_18h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
    collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors)
    write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/MTX vs Control/18h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    collectri_df_tf_interest_MTX_24h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
    collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors)
    write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/MTX vs Control/24h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    
    
    
    
    # OXA
    
    
    collectri_df_tf_interest_OXA_3h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
    collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors)
    write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/OXA vs Control/3h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    collectri_df_tf_interest_OXA_6h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
    collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors)
    write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/OXA vs Control/6h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    collectri_df_tf_interest_OXA_12h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
    collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors)
    write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/OXA vs Control/12h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    collectri_df_tf_interest_OXA_18h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
    collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors)
    write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/OXA vs Control/18h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    collectri_df_tf_interest_OXA_24h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
    collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors)
    write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/OXA vs Control/24h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    
    
    # CIS
    
    
    collectri_df_tf_interest_CIS_3h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
    collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors)
    write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/CIS vs Control/3h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    collectri_df_tf_interest_CIS_6h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
    collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_activators)
    write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/CIS vs Control/6h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    collectri_df_tf_interest_CIS_12h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
    collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors)
    write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/CIS vs Control/12h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    collectri_df_tf_interest_CIS_18h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
    collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors)
    write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/CIS vs Control/18h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    collectri_df_tf_interest_CIS_24h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
    collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
    collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors)
    write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cyto_interest/CIS vs Control/24h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
    
    
    
    
  }
  
  
  
  
  
  
  
# with scaled TFs (target info)
  
tf_all_target_cytos_interest_top_50 <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16]
  

  
  
for (i in 1:nrow(tf_target_cytos_interest_top_50)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_df_target_at_least_one_df %>% filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors)
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/MTX vs Control/3h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors)
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/MTX vs Control/6h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_MTX_12h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors)
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/MTX vs Control/12h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_MTX_18h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors)
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/MTX vs Control/18h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors)
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/MTX vs Control/24h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors)
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/OXA vs Control/3h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors)
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/OXA vs Control/6h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_OXA_12h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors)
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/OXA vs Control/12h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_OXA_18h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors)
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/OXA vs Control/18h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors)
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/OXA vs Control/24h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors)
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/CIS vs Control/3h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_activators)
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/CIS vs Control/6h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_CIS_12h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors)
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/CIS vs Control/12h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_CIS_18h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors)
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/CIS vs Control/18h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <- collectri_df_target_at_least_one_df %>%  filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors)
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_interest_new/CIS vs Control/24h/",rownames(tf_all_target_cytos_interest_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
}



  
# Get activator/repressor info

#tf_all_target_cytos_interest_top_50 <- tf_target_cytos_interest_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16]

top_50_tfs_targets_interest <- collectri_df %>% filter(source %in% rownames(tf_all_target_cytos_interest_top_50)) %>% filter(target == "Il6" | target == "Ccl4" | target == "Ccl5" | target == "Ccl20" | target == "Ccl9") %>% group_by(source, target) %>% as.data.frame()

for (i in 1:nrow(top_50_tfs_targets_interest)){
  if (top_50_tfs_targets_interest$mor[i] == 1){
    top_50_tfs_targets_interest$state[i] <- paste0("activator_",top_50_tfs_targets_interest$target[i])
  }
  else{
    top_50_tfs_targets_interest$state[i] <- paste0("repressor_",top_50_tfs_targets_interest$target[i])
  }
  
}
  
top_50_for_heatmap <- data.frame()
for (i in 1:nrow(tf_all_target_cytos_interest_top_50)){
  tf_targets_interest <- top_50_tfs_targets_interest %>% filter(source == rownames(tf_all_target_cytos_interest_top_50)[i]) %>% pivot_wider(., names_from = "target", values_from = "state", id_cols = "source") %>% unite(state_heatmap, all_of(2:ncol(.)), sep = " + ", na.rm = TRUE, remove = FALSE) %>% dplyr::select(1:2) %>% as.data.frame()
  top_50_for_heatmap <- rbind(tf_targets_interest, top_50_for_heatmap)
  
}


tf_all_target_cytos_interest_top_50_for_heatmap <- tf_all_target_cytos_interest_top_50 %>% as.data.frame() %>% rownames_to_column("source") %>% inner_join(., top_50_for_heatmap, by="source") %>% column_to_rownames("source") %>% as.data.frame()





tf_all_target_cytos_interest_top_50_for_heatmap_1 <- tf_all_target_cytos_interest_top_50 %>% as.data.frame() %>% rownames_to_column("source") %>% inner_join(., top_50_for_heatmap, by="source") %>% column_to_rownames("source") %>% as.data.frame() %>% filter(state_heatmap %in% c("activator_Il6","repressor_Il6", "activator_Ccl20","repressor_Ccl20","activator_Ccl4","repressor_Ccl4", "activator_Ccl5", "repressor_Ccl5"))


annotation_columns_treatment <- data.frame(treatment = rep(c("MTX","OXA","CIS"),c(5,5,5)))
rownames(annotation_columns_treatment)[1:5] <- colnames(tf_all_target_cytos_interest_top_50_for_heatmap_1)[1:5]
rownames(annotation_columns_treatment)[6:10] <- colnames(tf_all_target_cytos_interest_top_50_for_heatmap_1)[6:10]
rownames(annotation_columns_treatment)[11:15] <- colnames(tf_all_target_cytos_interest_top_50_for_heatmap_1)[11:15]


annotation_row_tfs_interest <- data.frame(tf_status = tf_all_target_cytos_interest_top_50_for_heatmap_1$state_heatmap)
rownames(annotation_row_tfs_interest) <- rownames(tf_all_target_cytos_interest_top_50_for_heatmap_1)


annot_colors_all <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"), tf_status = c(activator_Ccl20 = "green", activator_Il6 = "blue2", repressor_Il6 = "red", not_cyto_interest = "purple", repressor_Ccl4 = "grey", activator_Ccl5 = "orange"))



tf_activity_target_at_least_one_cytos_interest_heatmap_1 <- ComplexHeatmap::pheatmap(as.matrix(tf_all_target_cytos_interest_top_50_for_heatmap_1[,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, annotation_row = annotation_row_tfs_interest, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")




tf_all_target_cytos_interest_top_50_for_heatmap_2 <- tf_all_target_cytos_interest_top_50 %>% as.data.frame() %>% rownames_to_column("source") %>% inner_join(., top_50_for_heatmap, by="source") %>% column_to_rownames("source") %>% as.data.frame() %>% filter(!(state_heatmap %in% c("activator_Il6","repressor_Il6", "activator_Ccl20","repressor_Ccl20","activator_Ccl4","repressor_Ccl4", "activator_Ccl5", "repressor_Ccl5")))



annot_colors_all_multiple_cytos <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"))


annotation_row_tfs_interest_2 <- data.frame(tf_status = tf_all_target_cytos_interest_top_50_for_heatmap_2$state_heatmap)
rownames(annotation_row_tfs_interest_2) <- rownames(tf_all_target_cytos_interest_top_50_for_heatmap_2)


tf_activity_target_at_least_one_cytos_interest_heatmap_2 <- ComplexHeatmap::pheatmap(as.matrix(tf_all_target_cytos_interest_top_50_for_heatmap_2[,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all_multiple_cytos, annotation_row = annotation_row_tfs_interest_2, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(tf_activity_target_at_least_one_cytos_interest_heatmap_1)), grid.grabExpr(ComplexHeatmap::draw(tf_activity_target_at_least_one_cytos_interest_heatmap_2)))) & patchwork::plot_annotation(tag_levels = "A", title = "TF activity (targetting at least one cyto of interest) - Top 50 variable ") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face="bold", hjust=0.5))


