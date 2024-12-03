library(decoupleR)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(scales)
library(nichenetr)
library(patchwork)

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



# Load Jo list
jo_list_of_cyto <- read.csv("~/JONATHAN_LIST 1.csv", header = FALSE) %>% dplyr::rename("cytos" = "V1") %>% as.data.frame()
human_symbols <- jo_list_of_cyto$cytos
mouse_symbols <- human_symbols %>% convert_human_to_mouse_symbols()
jo_list_of_cyto <- mouse_symbols %>% as.data.frame() %>% drop_na() %>% dplyr::rename("cytos" = ".")
jo_list_of_cyto[159, "cytos"] <- "Ccl9"



# tfs that target at least one cyto of jo
collectri_df_at_least_one_cytokine_of_Jo <-  collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_cyto_jo = any(target %in% jo_list_of_cyto$cytos)) %>% as.data.frame() %>% filter(target_cyto_jo == "TRUE")
collectri_df_target_at_least_one_cyto_jo <- collectri_df %>% filter(source %in% collectri_df_at_least_one_cytokine_of_Jo$source)












# TFs that target cytos of Jo



## run decoupleR - MTX


### 3h

run_decoupler_MTX_3h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_3h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_3h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_3h, seed = 2)


# Get results
run_decoupler_MTX_3h_consensus_pivot <- run_decoupler_MTX_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_MTX_3h <- t(run_decoupler_MTX_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_3h" = "t")




### 6h

run_decoupler_MTX_6h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_6h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_6h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_6h, seed = 2)


# Get results
run_decoupler_MTX_6h_consensus_pivot <- run_decoupler_MTX_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_MTX_6h <- t(run_decoupler_MTX_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_6h" = "t")




### 12h

run_decoupler_MTX_12h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_12h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_12h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_12h, seed = 2)


# Get results
run_decoupler_MTX_12h_consensus_pivot <- run_decoupler_MTX_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_MTX_12h <- t(run_decoupler_MTX_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_12h" = "t")



### 18h

run_decoupler_MTX_18h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_18h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_18h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_18h, seed = 2)


# Get results
run_decoupler_MTX_18h_consensus_pivot <- run_decoupler_MTX_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_MTX_18h <- t(run_decoupler_MTX_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_18h" = "t")



### 24h

run_decoupler_MTX_24h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_24h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_MTX_24h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_24h, seed = 2)


# Get results
run_decoupler_MTX_24h_consensus_pivot <- run_decoupler_MTX_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_MTX_24h <- t(run_decoupler_MTX_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_24h" = "t")




## run decoupleR - OXA


### 3h

run_decoupler_OXA_3h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_3h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_3h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_3h, seed = 2)


# Get results
run_decoupler_OXA_3h_consensus_pivot <- run_decoupler_OXA_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_OXA_3h <- t(run_decoupler_OXA_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_3h" = "t")




### 6h

run_decoupler_OXA_6h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_6h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_6h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_6h, seed = 2)


# Get results
run_decoupler_OXA_6h_consensus_pivot <- run_decoupler_OXA_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_OXA_6h <- t(run_decoupler_OXA_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_6h" = "t")




### 12h

run_decoupler_OXA_12h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_12h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_12h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_12h, seed = 2)


# Get results

run_decoupler_OXA_12h_consensus_pivot <- run_decoupler_OXA_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_OXA_12h <- t(run_decoupler_OXA_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_12h" = "t")




### 18h

run_decoupler_OXA_18h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_18h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_18h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_18h, seed = 2)


# Get results
run_decoupler_OXA_18h_consensus_pivot <- run_decoupler_OXA_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_OXA_18h <- t(run_decoupler_OXA_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_18h" = "t")



### 24h

run_decoupler_OXA_24h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_24h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_OXA_24h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_24h, seed = 2)


# Get results
run_decoupler_OXA_24h_consensus_pivot <- run_decoupler_OXA_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_OXA_24h <- t(run_decoupler_OXA_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_24h" = "t")





## run decoupleR - CIS


### 3h

run_decoupler_CIS_3h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_3h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_3h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_3h, seed = 2)


# Get results
run_decoupler_CIS_3h_consensus_pivot <- run_decoupler_CIS_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_CIS_3h <- t(run_decoupler_CIS_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_3h" = "t")




### 6h


run_decoupler_CIS_6h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_6h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_6h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_6h, seed = 2)


# Get results
run_decoupler_CIS_6h_consensus_pivot <- run_decoupler_CIS_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_CIS_6h <- t(run_decoupler_CIS_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_6h" = "t")




### 12h

run_decoupler_CIS_12h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_12h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_12h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_12h, seed = 2)


# Get results
run_decoupler_CIS_12h_consensus_pivot <- run_decoupler_CIS_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_CIS_12h <- t(run_decoupler_CIS_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_12h" = "t")




### 18h

run_decoupler_CIS_18h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_18h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_18h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_18h, seed = 2)


# Get results
run_decoupler_CIS_18h_consensus_pivot <- run_decoupler_CIS_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_CIS_18h <- t(run_decoupler_CIS_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_18h" = "t")



### 24h

run_decoupler_CIS_24h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_24h, net = collectri_df_target_at_least_one_cyto_jo, .source = "source", .target = "target", consensus_score = FALSE)
run_decoupler_CIS_24h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_24h, seed = 2)


# Get results
run_decoupler_CIS_24h_consensus_pivot <- run_decoupler_CIS_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_target_cytos_jo_CIS_24h <- t(run_decoupler_CIS_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_24h" = "t")








# Join 

# MTX
tf_target_cytos_jo_MTX_3h_6h <- inner_join(tf_target_cytos_jo_MTX_3h, tf_target_cytos_jo_MTX_6h, by="TFs")
tf_target_cytos_jo_MTX_3h_6h_12h <- inner_join(tf_target_cytos_jo_MTX_3h_6h, tf_target_cytos_jo_MTX_12h, by="TFs")
tf_target_cytos_jo_MTX_3h_6h_12h_18h <- inner_join(tf_target_cytos_jo_MTX_3h_6h_12h, tf_target_cytos_jo_MTX_18h, by="TFs")
tf_target_cytos_jo_MTX_all_timepoints <- inner_join(tf_target_cytos_jo_MTX_3h_6h_12h_18h, tf_target_cytos_jo_MTX_24h, by="TFs")


# OXA
tf_target_cytos_jo_OXA_3h_6h <- inner_join(tf_target_cytos_jo_OXA_3h, tf_target_cytos_jo_OXA_6h, by="TFs")
tf_target_cytos_jo_OXA_3h_6h_12h <- inner_join(tf_target_cytos_jo_OXA_3h_6h, tf_target_cytos_jo_OXA_12h, by="TFs")
tf_target_cytos_jo_OXA_3h_6h_12h_18h <- inner_join(tf_target_cytos_jo_OXA_3h_6h_12h, tf_target_cytos_jo_OXA_18h, by="TFs")
tf_target_cytos_jo_OXA_all_timepoints <- inner_join(tf_target_cytos_jo_OXA_3h_6h_12h_18h, tf_target_cytos_jo_OXA_24h, by="TFs")


# CIS
tf_target_cytos_jo_CIS_3h_6h <- inner_join(tf_target_cytos_jo_CIS_3h, tf_target_cytos_jo_CIS_6h, by="TFs")
tf_target_cytos_jo_CIS_3h_6h_12h <- inner_join(tf_target_cytos_jo_CIS_3h_6h, tf_target_cytos_jo_CIS_12h, by="TFs")
tf_target_cytos_jo_CIS_3h_6h_12h_18h <- inner_join(tf_target_cytos_jo_CIS_3h_6h_12h, tf_target_cytos_jo_CIS_18h, by="TFs")
tf_target_cytos_jo_CIS_all_timepoints <- inner_join(tf_target_cytos_jo_CIS_3h_6h_12h_18h, tf_target_cytos_jo_CIS_24h, by="TFs")




# All conditions
tf_target_cytos_jo_MTX_OXA_all_timepoints <- inner_join(tf_target_cytos_jo_MTX_all_timepoints, 
                                                              tf_target_cytos_jo_OXA_all_timepoints, by="TFs")


tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints <- inner_join(tf_target_cytos_jo_MTX_OXA_all_timepoints, tf_target_cytos_jo_CIS_all_timepoints, by = "TFs") %>% column_to_rownames("TFs")

#write.csv(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints, "tf_targets_cytos_jo_all_timepoints.csv")


# scaled by TFs
tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tf <- tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints %>% t() %>% scale() %>% t() %>% as.data.frame()

tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd <- tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tf %>% rowwise() %>%
  mutate(standard_deviation = c_across(everything()) %>% sd()) %>%
  ungroup() %>% as.data.frame()

tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final <- cbind(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tf, tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd$standard_deviation) %>% dplyr::rename("sd" = "tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd$standard_deviation") %>% dplyr::arrange(desc(sd)) %>% as.matrix()



# Visualization




palette_length <- 100 

my_color <- viridis::viridis(palette_length)


annotation_columns_treatment <- data.frame(treatment = rep(c("MTX","OXA","CIS"),c(5,5,5)))
rownames(annotation_columns_treatment)[1:5] <- colnames(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints)[1:5]
rownames(annotation_columns_treatment)[6:10] <- colnames(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints)[6:10]
rownames(annotation_columns_treatment)[11:15] <- colnames(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints)[11:15]

annot_colors_all <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"))




tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_sd <- tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints %>% rowwise() %>%
  mutate(standard_deviation = c_across(everything()) %>% sd()) %>%
  ungroup() %>% as.data.frame()


tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final <- cbind(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints, tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_sd$standard_deviation) %>% dplyr::rename("sd" = "tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_sd$standard_deviation") %>% dplyr::arrange(desc(sd)) %>% as.matrix()



#tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final_scaled <- scales::rescale(as.matrix(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final[,-16]))



tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap <- scale(as.matrix(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final[,-16]))
tf_activity_target_al_least_one_cytos_jo_heatmap <- ComplexHeatmap::pheatmap(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap[1:50,], color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, main="TF activity (targetting at least one cyto of Jo) in all the conditions - top 50", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")


tf_activity_target_at_least_one_cytos_jo_heatmap_scaled_tf <- ComplexHeatmap::pheatmap(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16], color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, main="TF activity (targetting at least one cyto of Jo) in all the conditions - top 50", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")




patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(tf_activity_all_targets_heatmap)), grid.grabExpr(ComplexHeatmap::draw(tf_activity_target_al_least_one_cyto_jo_heatmap)),grid.grabExpr(ComplexHeatmap::draw(tf_activity_target_all_cytos_interest_heatmap)), grid.grabExpr(ComplexHeatmap::draw(tf_activity_target_al_least_one_cytos_jo_heatmap))))





write.csv(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap, "figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/tf_targets_at_least_one_cytos_jo.csv", row.names = TRUE, col.names = TRUE)


write.csv(tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap[1:50,], "figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/top_50_tf_targets_at_least_one_cytos_jo.csv", row.names = TRUE, col.names = TRUE)






# Retrieve TFs with differentially expressed targets (+ with mode of reg) 


tf_all_target_cytos_jo_top_50 <- tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_final_scaled_heatmap[1:50,]



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



for (i in 1:nrow(tf_all_target_cytos_jo_top_50)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_df_target_at_least_one_cyto_jo %>% filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors)
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/MTX vs Control/3h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors)
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/MTX vs Control/6h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_MTX_12h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors)
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/MTX vs Control/12h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_MTX_18h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors)
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/MTX vs Control/18h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors)
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/MTX vs Control/24h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors)
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/OXA vs Control/3h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors)
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/OXA vs Control/6h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_OXA_12h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors)
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/OXA vs Control/12h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_OXA_18h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors)
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/OXA vs Control/18h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors)
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/OXA vs Control/24h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors)
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/CIS vs Control/3h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_activators)
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/CIS vs Control/6h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_CIS_12h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors)
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/CIS vs Control/12h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_CIS_18h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors)
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/CIS vs Control/18h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors)
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_target_cytos_jo/CIS vs Control/24h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}








# with scaled TFs


tf_target_cytos_jo_top_50 <- tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16]



for (i in 1:nrow(tf_target_cytos_jo_top_50)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_df_target_at_least_one_cyto_jo %>% filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors)
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/MTX vs Control/3h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors)
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/MTX vs Control/6h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_MTX_12h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors)
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/MTX vs Control/12h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_MTX_18h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors)
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/MTX vs Control/18h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors)
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/MTX vs Control/24h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors)
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/OXA vs Control/3h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors)
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/OXA vs Control/6h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_OXA_12h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors)
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/OXA vs Control/12h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_OXA_18h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors)
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/OXA vs Control/18h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors)
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/OXA vs Control/24h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors)
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/CIS vs Control/3h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_activators)
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/CIS vs Control/6h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_CIS_12h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors)
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/CIS vs Control/12h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  collectri_df_tf_interest_CIS_18h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors)
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/CIS vs Control/18h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <- collectri_df_target_at_least_one_cyto_jo %>%  filter(source == rownames(tf_all_target_cytos_jo_top_50)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors)
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("figures_for_publication/In_vitro/TF_analysis/list_of_tf_targets_cytos_jo_new/CIS vs Control/24h/",rownames(tf_all_target_cytos_jo_top_50)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}


# Get activator/repressor info

#tf_target_cytos_jo_top_50 <- tf_target_cytos_jo_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16]


top_50_tfs_targets_jo <- collectri_df %>% filter(source %in% rownames(tf_target_cytos_jo_top_50)) %>% filter(target == "Il6" | target == "Ccl4" | target == "Ccl5" | target == "Ccl20" | target == "Ccl9") %>% group_by(source, target) %>% as.data.frame()


tf_target_cytos_jo_top_50_variable_with_activators_repressors_interest <- tf_target_cytos_jo_top_50 %>% as.data.frame() %>% filter(rownames(.) %in% unique(top_50_tfs_targets_jo$source)) 




for (i in 1:nrow(top_50_tfs_targets_jo)){
  if (top_50_tfs_targets_jo$mor[i] == 1){
    top_50_tfs_targets_jo$state[i] <- paste0("activator_",top_50_tfs_targets_jo$target[i])
  }
  else{
    top_50_tfs_targets_jo$state[i] <- paste0("repressor_",top_50_tfs_targets_jo$target[i])
  }
  
}

top_50_for_heatmap <- data.frame()
for (i in 1:nrow(tf_target_cytos_jo_top_50_variable_with_activators_repressors_interest)){
  tf_targets_jo <- top_50_tfs_targets_jo %>% filter(source == rownames(tf_target_cytos_jo_top_50_variable_with_activators_repressors_interest)[i]) %>% pivot_wider(., names_from = "target", values_from = "state", id_cols = "source") %>% unite(state_heatmap, all_of(2:ncol(.)), sep = " + ", na.rm = TRUE, remove = FALSE) %>% dplyr::select(1:2) %>% as.data.frame()
  top_50_for_heatmap <- rbind(tf_targets_jo, top_50_for_heatmap)
  
}


tf_target_cytos_jo_18_out_of_50 <- tf_target_cytos_jo_top_50_variable_with_activators_repressors_interest %>% as.data.frame() %>% rownames_to_column("source") %>% inner_join(., top_50_for_heatmap, by="source") %>% column_to_rownames("source") %>% as.data.frame()


tf_target_cytos_jo_32_out_of_50_not_interest <- tf_target_cytos_jo_top_50 %>% as.data.frame() %>% filter(!(rownames(.) %in% unique(top_50_tfs_targets_jo$source))) %>% mutate(state_heatmap = rep("not_cyto_interest", nrow(.)))

tf_targets_cytos_jo_with_and_without_cytos_interest <- rbind(tf_target_cytos_jo_18_out_of_50, tf_target_cytos_jo_32_out_of_50_not_interest)



tf_targets_cytos_jo_with_and_without_cytos_interest_filtered <- tf_targets_cytos_jo_with_and_without_cytos_interest %>% filter(state_heatmap %in% c("activator_Il6","repressor_Il6", "activator_Ccl20","repressor_Ccl20","activator_Ccl4","repressor_Ccl4", "activator_Ccl5", "repressor_Ccl5", "not_cyto_interest"))
                                                                                                                    

annotation_columns_treatment <- data.frame(treatment = rep(c("MTX","OXA","CIS"),c(5,5,5)))
rownames(annotation_columns_treatment)[1:5] <- colnames(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered)[1:5]
rownames(annotation_columns_treatment)[6:10] <- colnames(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered)[6:10]
rownames(annotation_columns_treatment)[11:15] <- colnames(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered)[11:15]





annotation_row_tfs_jo <- data.frame(tf_status = tf_targets_cytos_jo_with_and_without_cytos_interest_filtered$state_heatmap)
rownames(annotation_row_tfs_jo) <- rownames(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered)


annot_colors_all <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"), tf_status = c(activator_Ccl20 = "green", activator_Il6 = "blue2", repressor_Il6 = "red", not_cyto_interest = "purple"))


tf_activity_cytos_jo_1 <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered[,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, annotation_row = annotation_row_tfs_jo, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")




tf_targets_cytos_jo_with_and_without_cytos_interest_filtered_with_more_than_one_cyto_interest <- tf_targets_cytos_jo_with_and_without_cytos_interest %>% filter(!(state_heatmap %in% c("activator_Il6","repressor_Il6", "activator_Ccl20","repressor_Ccl20","activator_Ccl4","repressor_Ccl4", "activator_Ccl5", "repressor_Ccl5", "not_cyto_interest")))


annotation_row_tfs_jo_with_more_cyto_target <- data.frame(tf_status = tf_targets_cytos_jo_with_and_without_cytos_interest_filtered_with_more_than_one_cyto_interest$state_heatmap)
rownames(annotation_row_tfs_jo_with_more_cyto_target) <- rownames(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered_with_more_than_one_cyto_interest)


annot_colors_all_for_cyto_interest <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"))

tf_activity_cytos_jo_2 <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_cytos_jo_with_and_without_cytos_interest_filtered_with_more_than_one_cyto_interest[,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all_for_cyto_interest, annotation_row = annotation_row_tfs_jo_with_more_cyto_target, fontsize_row = 11, fontsize_col = 12, cluster_rows = TRUE, cluster_cols = TRUE, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")


patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(tf_activity_cytos_jo_1)), grid.grabExpr(ComplexHeatmap::draw(tf_activity_cytos_jo_2)))) & patchwork::plot_annotation(tag_levels = "A", title = "TF activity (targetting at least one cyto of Jo) - Top 50 variable ") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face="bold", hjust=0.5))



#"TF activity (targetting at least one cyto of interest) in all the conditions - top 50"
