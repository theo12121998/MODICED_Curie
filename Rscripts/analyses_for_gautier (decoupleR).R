

# Load libraries of interest
library(ComplexHeatmap)
library(decoupleR)
library(tidyverse)




# Import t-values of DEGs


# MTX vs Control

t_values_diff_genes_MTX_Control_3h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_3h_vs_Control_3h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_3h) <- make.names(t_values_diff_genes_MTX_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_3h <- t_values_diff_genes_MTX_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_6h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_6h_vs_Control_6h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_6h) <- make.names(t_values_diff_genes_MTX_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_6h <- t_values_diff_genes_MTX_Control_6h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_12h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_12h_vs_Control_12h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_12h) <- make.names(t_values_diff_genes_MTX_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_12h <- t_values_diff_genes_MTX_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_18h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_18h_vs_Control_18h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_18h) <- make.names(t_values_diff_genes_MTX_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_18h <- t_values_diff_genes_MTX_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_MTX_Control_24h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_MTX_24h_vs_Control_24h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_MTX_Control_24h) <- make.names(t_values_diff_genes_MTX_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_MTX_Control_24h <- t_values_diff_genes_MTX_Control_24h %>% dplyr::select(-"gene_name")




# OXA vs Control

t_values_diff_genes_OXA_Control_3h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_3h_vs_Control_3h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_3h) <- make.names(t_values_diff_genes_OXA_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_3h <- t_values_diff_genes_OXA_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_6h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_6h_vs_Control_6h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_6h) <- make.names(t_values_diff_genes_OXA_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_6h <- t_values_diff_genes_OXA_Control_6h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_12h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_12h_vs_Control_12h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_12h) <- make.names(t_values_diff_genes_OXA_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_12h <- t_values_diff_genes_OXA_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_18h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_18h_vs_Control_18h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_18h) <- make.names(t_values_diff_genes_OXA_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_18h <- t_values_diff_genes_OXA_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_OXA_Control_24h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_OXA_24h_vs_Control_24h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_OXA_Control_24h) <- make.names(t_values_diff_genes_OXA_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_OXA_Control_24h <- t_values_diff_genes_OXA_Control_24h %>% dplyr::select(-"gene_name")





# CIS vs Control

t_values_diff_genes_CIS_Control_3h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_3h_vs_Control_3h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_3h) <- make.names(t_values_diff_genes_CIS_Control_3h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_3h <- t_values_diff_genes_CIS_Control_3h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_6h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_6h_vs_Control_6h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_6h) <- make.names(t_values_diff_genes_CIS_Control_6h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_6h <- t_values_diff_genes_CIS_Control_6h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_12h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_12h_vs_Control_12h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_12h) <- make.names(t_values_diff_genes_CIS_Control_12h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_12h <- t_values_diff_genes_CIS_Control_12h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_18h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_18h_vs_Control_18h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_18h) <- make.names(t_values_diff_genes_CIS_Control_18h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_18h <- t_values_diff_genes_CIS_Control_18h %>% dplyr::select(-"gene_name")


t_values_diff_genes_CIS_Control_24h <- read.csv("t_values_all_diff_genes_for_decoupleR_gautier/all_diff_genes_CIS_24h_vs_Control_24h.csv") %>% select(c("t","gene_name"))
rownames(t_values_diff_genes_CIS_Control_24h) <- make.names(t_values_diff_genes_CIS_Control_24h$gene_name, unique=TRUE)
t_values_diff_genes_CIS_Control_24h <- t_values_diff_genes_CIS_Control_24h %>% dplyr::select(-"gene_name")









# Get collectri with all the tfs
collectri_df <- get_collectri(organism="mouse", split_complexes = FALSE) %>% as.data.frame()
gautier_target_genes <- collectri_df %>% dplyr::group_by(source) %>% dplyr::summarise(target_gautier_genes = any(target %in% genes_for_dea_mouse$genes)) %>% as.data.frame() %>% filter(target_gautier_genes == "TRUE")
collectri_filtered_gautier_target_genes <- collectri_df %>% filter(source %in% gautier_target_genes$source)






## run decoupleR - MTX


### 3h

run_decoupler_MTX_3h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_3h, net = collectri_filtered_gautier_target_genes , .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_MTX_3h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_3h, seed = 2)


# Get results
run_decoupler_MTX_3h_consensus_pivot <- run_decoupler_MTX_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_MTX_3h <- t(run_decoupler_MTX_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_3h" = "t")




### 6h

run_decoupler_MTX_6h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_6h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_MTX_6h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_6h, seed = 2)


# Get results
run_decoupler_MTX_6h_consensus_pivot <- run_decoupler_MTX_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_MTX_6h <- t(run_decoupler_MTX_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_6h" = "t")




### 12h

run_decoupler_MTX_12h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_12h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_MTX_12h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_12h, seed = 2)


# Get results
run_decoupler_MTX_12h_consensus_pivot <- run_decoupler_MTX_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_MTX_12h <- t(run_decoupler_MTX_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_12h" = "t")



### 18h

run_decoupler_MTX_18h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_18h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_MTX_18h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_18h, seed = 2)


# Get results
run_decoupler_MTX_18h_consensus_pivot <- run_decoupler_MTX_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_MTX_18h <- t(run_decoupler_MTX_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_18h" = "t")



### 24h

run_decoupler_MTX_24h <- decoupleR::decouple(mat = t_values_diff_genes_MTX_Control_24h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_MTX_24h_consensus <- decoupleR::run_consensus(run_decoupler_MTX_24h, seed = 2)


# Get results
run_decoupler_MTX_24h_consensus_pivot <- run_decoupler_MTX_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_MTX_24h <- t(run_decoupler_MTX_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("MTX_24h" = "t")




## run decoupleR - OXA


### 3h

run_decoupler_OXA_3h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_3h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_OXA_3h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_3h, seed = 2)


# Get results
run_decoupler_OXA_3h_consensus_pivot <- run_decoupler_OXA_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_OXA_3h <- t(run_decoupler_OXA_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_3h" = "t")




### 6h

run_decoupler_OXA_6h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_6h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_OXA_6h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_6h, seed = 2)


# Get results
run_decoupler_OXA_6h_consensus_pivot <- run_decoupler_OXA_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_OXA_6h <- t(run_decoupler_OXA_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_6h" = "t")




### 12h

run_decoupler_OXA_12h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_12h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_OXA_12h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_12h, seed = 2)


# Get results

run_decoupler_OXA_12h_consensus_pivot <- run_decoupler_OXA_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_OXA_12h <- t(run_decoupler_OXA_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_12h" = "t")




### 18h

run_decoupler_OXA_18h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_18h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_OXA_18h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_18h, seed = 2)


# Get results
run_decoupler_OXA_18h_consensus_pivot <- run_decoupler_OXA_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_OXA_18h <- t(run_decoupler_OXA_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_18h" = "t")



### 24h

run_decoupler_OXA_24h <- decoupleR::decouple(mat = t_values_diff_genes_OXA_Control_24h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_OXA_24h_consensus <- decoupleR::run_consensus(run_decoupler_OXA_24h, seed = 2)


# Get results
run_decoupler_OXA_24h_consensus_pivot <- run_decoupler_OXA_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_OXA_24h <- t(run_decoupler_OXA_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("OXA_24h" = "t")





## run decoupleR - CIS


### 3h

run_decoupler_CIS_3h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_3h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_CIS_3h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_3h, seed = 2)


# Get results
run_decoupler_CIS_3h_consensus_pivot <- run_decoupler_CIS_3h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_CIS_3h <- t(run_decoupler_CIS_3h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_3h" = "t")




### 6h


run_decoupler_CIS_6h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_6h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_CIS_6h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_6h, seed = 2)


# Get results
run_decoupler_CIS_6h_consensus_pivot <- run_decoupler_CIS_6h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_CIS_6h <- t(run_decoupler_CIS_6h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_6h" = "t")




### 12h

run_decoupler_CIS_12h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_12h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_CIS_12h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_12h, seed = 2)


# Get results
run_decoupler_CIS_12h_consensus_pivot <- run_decoupler_CIS_12h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_CIS_12h <- t(run_decoupler_CIS_12h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_12h" = "t")




### 18h

run_decoupler_CIS_18h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_18h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_CIS_18h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_18h, seed = 2)


# Get results
run_decoupler_CIS_18h_consensus_pivot <- run_decoupler_CIS_18h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_CIS_18h <- t(run_decoupler_CIS_18h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_18h" = "t")



### 24h

run_decoupler_CIS_24h <- decoupleR::decouple(mat = t_values_diff_genes_CIS_Control_24h, net = collectri_filtered_gautier_target_genes, .source = "source", .target = "target", consensus_score = FALSE, minsize = 7)
run_decoupler_CIS_24h_consensus <- decoupleR::run_consensus(run_decoupler_CIS_24h, seed = 2)


# Get results
run_decoupler_CIS_24h_consensus_pivot <- run_decoupler_CIS_24h_consensus %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


tf_all_targets_CIS_24h <- t(run_decoupler_CIS_24h_consensus_pivot) %>% as.data.frame() %>% rownames_to_column("TFs") %>% dplyr::rename("CIS_24h" = "t")








# Join 

# MTX
tf_all_targets_MTX_3h_6h <- inner_join(tf_all_targets_MTX_3h, tf_all_targets_MTX_6h, by="TFs")
tf_all_targets_MTX_3h_6h_12h <- inner_join(tf_all_targets_MTX_3h_6h, tf_all_targets_MTX_12h, by="TFs")
tf_all_targets_MTX_3h_6h_12h_18h <- inner_join(tf_all_targets_MTX_3h_6h_12h, tf_all_targets_MTX_18h, by="TFs")
tf_all_targets_MTX_all_timepoints <- inner_join(tf_all_targets_MTX_3h_6h_12h_18h, tf_all_targets_MTX_24h, by="TFs")


# OXA
tf_all_targets_OXA_3h_6h <- inner_join(tf_all_targets_OXA_3h, tf_all_targets_OXA_6h, by="TFs")
tf_all_targets_OXA_3h_6h_12h <- inner_join(tf_all_targets_OXA_3h_6h, tf_all_targets_OXA_12h, by="TFs")
tf_all_targets_OXA_3h_6h_12h_18h <- inner_join(tf_all_targets_OXA_3h_6h_12h, tf_all_targets_OXA_18h, by="TFs")
tf_all_targets_OXA_all_timepoints <- inner_join(tf_all_targets_OXA_3h_6h_12h_18h, tf_all_targets_OXA_24h, by="TFs")


# CIS
tf_all_targets_CIS_3h_6h <- inner_join(tf_all_targets_CIS_3h, tf_all_targets_CIS_6h, by="TFs")
tf_all_targets_CIS_3h_6h_12h <- inner_join(tf_all_targets_CIS_3h_6h, tf_all_targets_CIS_12h, by="TFs")
tf_all_targets_CIS_3h_6h_12h_18h <- inner_join(tf_all_targets_CIS_3h_6h_12h, tf_all_targets_CIS_18h, by="TFs")
tf_all_targets_CIS_all_timepoints <- inner_join(tf_all_targets_CIS_3h_6h_12h_18h, tf_all_targets_CIS_24h, by="TFs")


# All conditions
tf_all_targets_MTX_OXA_all_timepoints <- inner_join(tf_all_targets_MTX_all_timepoints, 
                                                    tf_all_targets_OXA_all_timepoints, by="TFs")


# non-scaled
tf_targets_MTX_OXA_CIS_all_timepoints <- inner_join(tf_all_targets_MTX_OXA_all_timepoints, tf_all_targets_CIS_all_timepoints, by = "TFs") %>% column_to_rownames("TFs")


# scaled
tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tf <- tf_targets_MTX_OXA_CIS_all_timepoints %>% t() %>% scale() %>% t() %>% as.data.frame()


# Get standard deviation of TFs
tf_targets_genes_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd <- tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tf %>% rowwise() %>%
  mutate(standard_deviation = c_across(everything()) %>% sd()) %>%
  ungroup() %>% as.data.frame()

tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final <- cbind(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tf, tf_targets_genes_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd$standard_deviation) %>% dplyr::rename("sd" = "tf_targets_genes_MTX_OXA_CIS_all_timepoints_scaled_by_tf_sd$standard_deviation") %>% dplyr::arrange(desc(sd)) %>% as.matrix()




# Visualization of the TF activity (Heatmap)



palette_length <- 100 

my_color <- viridis::viridis(palette_length)


annotation_columns_treatment <- data.frame(treatment = rep(c("MTX","OXA","CIS"),c(5,5,5)))
rownames(annotation_columns_treatment)[1:5] <- colnames(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tf)[1:5]
rownames(annotation_columns_treatment)[6:10] <- colnames(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tf)[6:10]
rownames(annotation_columns_treatment)[11:15] <- colnames(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tf)[11:15]

annot_colors_all <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen"))


top_50_tfs <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 11, cluster_rows = TRUE, cluster_cols = TRUE, main="1 - 50 TFs",fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_p = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")


top_100_tfs <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[51:100,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 11, cluster_rows = TRUE, cluster_cols = TRUE, main="51 - 100 TFs", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_p = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")




top_150_tfs <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[101:150,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 11, cluster_rows = TRUE, cluster_cols = TRUE, main="101 - 150 TFs",fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_p = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



top_200_tfs <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[151:200,-16]), color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 11, cluster_rows = TRUE, cluster_cols = TRUE, main="151 - 200 TFs", fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_p = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")



last_70_tfs <- ComplexHeatmap::pheatmap(as.matrix(tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final)[201:270,-16], color = my_color, annotation_col = annotation_columns_treatment, annotation_colors = annot_colors_all, fontsize_row = 11, fontsize_col = 11, cluster_rows = TRUE, cluster_cols = TRUE, main= "201 - 270 TFs",fontface = 3,show_colnames = TRUE, heatmap_legend_param = list(title = "activity_score", title_p = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), border_color = "white", cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")





top_50_100_tfs_patch <- patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(top_50_tfs)), grid.grabExpr(ComplexHeatmap::draw(top_100_tfs))), nrow = 1) & patchwork::plot_annotation(tag_levels = "A", title = "TF activity") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face="bold", hjust=0.5))


top_150_200_patch <- patchwork::wrap_plots(list(grid.grabExpr(ComplexHeatmap::draw(top_150_tfs)), grid.grabExpr(ComplexHeatmap::draw(top_200_tfs))), nrow = 1) & patchwork::plot_annotation(tag_levels = "A", title = "TF activity") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face="bold", hjust=0.5))


top_200_tfs_patch <- patchwork::wrap_plots(n,i, ncol=1) & patchwork::plot_annotation(tag_levels = "A", title = "TF activity") & theme(plot.tag = element_text(face = "bold", size=13)) & theme(plot.title = element_text(face="bold", hjust=0.5))











# Get differentially expressed targets for all the TFs


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





# 1 - 50 TFs

top_50_tf_targets_at_least_one_gautier_genes <- tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[1:50,-16]


for (i in 1:nrow(top_50_tf_targets_at_least_one_gautier_genes)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/MTX_vs_Control/3h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/MTX_vs_Control/6h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/MTX_vs_Control/12h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_MTX_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/MTX_vs_Control/18h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/MTX_vs_Control/24h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/OXA_vs_Control/3h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/OXA_vs_Control/6h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/OXA_vs_Control/12h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_OXA_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/OXA_vs_Control/18h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/OXA_vs_Control/24h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/CIS_vs_Control/3h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/CIS_vs_Control/6h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/CIS_vs_Control/12h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_CIS_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/CIS_vs_Control/18h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_50_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/Top 50 variable/CIS_vs_Control/24h/",rownames(top_50_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}



# 51 - 100 TFs

top_100_tf_targets_at_least_one_gautier_genes <- tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[51:100,-16]


for (i in 1:nrow(top_100_tf_targets_at_least_one_gautier_genes)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/MTX_vs_Control/3h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/MTX_vs_Control/6h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/MTX_vs_Control/12h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_MTX_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/MTX_vs_Control/18h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/MTX_vs_Control/24h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/OXA_vs_Control/3h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/OXA_vs_Control/6h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/OXA_vs_Control/12h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_OXA_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/OXA_vs_Control/18h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/OXA_vs_Control/24h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/CIS_vs_Control/3h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/CIS_vs_Control/6h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/CIS_vs_Control/12h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_CIS_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/CIS_vs_Control/18h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_100_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/51 - 100 TFs/CIS_vs_Control/24h/",rownames(top_100_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}








# 101 - 150 TFs

top_150_tf_targets_at_least_one_gautier_genes <- tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[101:150,-16]


for (i in 1:nrow(top_150_tf_targets_at_least_one_gautier_genes)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/MTX_vs_Control/3h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/MTX_vs_Control/6h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/MTX_vs_Control/12h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_MTX_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/MTX_vs_Control/18h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/MTX_vs_Control/24h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/OXA_vs_Control/3h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/OXA_vs_Control/6h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/OXA_vs_Control/12h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_OXA_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/OXA_vs_Control/18h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/OXA_vs_Control/24h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/CIS_vs_Control/3h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/CIS_vs_Control/6h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/CIS_vs_Control/12h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_CIS_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/CIS_vs_Control/18h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_150_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/101 - 150 TFs/CIS_vs_Control/24h/",rownames(top_150_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}




# 151 - 200 TFs




top_200_tf_targets_at_least_one_gautier_genes <- tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[151:200,-16]


for (i in 1:nrow(top_200_tf_targets_at_least_one_gautier_genes)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/MTX_vs_Control/3h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/MTX_vs_Control/6h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/MTX_vs_Control/12h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_MTX_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/MTX_vs_Control/18h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/MTX_vs_Control/24h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/OXA_vs_Control/3h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/OXA_vs_Control/6h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/OXA_vs_Control/12h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_OXA_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/OXA_vs_Control/18h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/OXA_vs_Control/24h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/CIS_vs_Control/3h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/CIS_vs_Control/6h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/CIS_vs_Control/12h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_CIS_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/CIS_vs_Control/18h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_200_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/151 - 200 TFs/CIS_vs_Control/24h/",rownames(top_200_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}



# 201 - 270 TFs



top_270_tf_targets_at_least_one_gautier_genes <- tf_targets_MTX_OXA_CIS_all_timepoints_scaled_by_tfs_final[201:270,-16]


for (i in 1:nrow(top_270_tf_targets_at_least_one_gautier_genes)){
  
  
  
  # MTX
  
  
  collectri_df_tf_interest_MTX_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_3h_df, by="target")
  collectri_df_tf_interest_MTX_3h_activators <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_repressors <- collectri_df_tf_interest_MTX_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_3h_new <- rbind(collectri_df_tf_interest_MTX_3h_activators, collectri_df_tf_interest_MTX_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/MTX_vs_Control/3h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_6h_df, by="target")
  collectri_df_tf_interest_MTX_6h_activators <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_repressors <- collectri_df_tf_interest_MTX_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_6h_new <- rbind(collectri_df_tf_interest_MTX_6h_activators, collectri_df_tf_interest_MTX_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/MTX_vs_Control/6h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_12h_df, by="target")
  collectri_df_tf_interest_MTX_12h_activators <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_repressors <- collectri_df_tf_interest_MTX_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_12h_new <- rbind(collectri_df_tf_interest_MTX_12h_activators, collectri_df_tf_interest_MTX_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/MTX_vs_Control/12h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_MTX_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_18h_df, by="target")
  collectri_df_tf_interest_MTX_18h_activators <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_repressors <- collectri_df_tf_interest_MTX_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_18h_new <- rbind(collectri_df_tf_interest_MTX_18h_activators, collectri_df_tf_interest_MTX_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/MTX_vs_Control/18h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_MTX_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_MTX_Control_24h_df, by="target")
  collectri_df_tf_interest_MTX_24h_activators <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_repressors <- collectri_df_tf_interest_MTX_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_MTX_24h_new <- rbind(collectri_df_tf_interest_MTX_24h_activators, collectri_df_tf_interest_MTX_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_MTX_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/MTX_vs_Control/24h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  
  # OXA
  
  
  collectri_df_tf_interest_OXA_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_3h_df, by="target")
  collectri_df_tf_interest_OXA_3h_activators <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_repressors <- collectri_df_tf_interest_OXA_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_3h_new <- rbind(collectri_df_tf_interest_OXA_3h_activators, collectri_df_tf_interest_OXA_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/OXA_vs_Control/3h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_6h_df, by="target")
  collectri_df_tf_interest_OXA_6h_activators <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_repressors <- collectri_df_tf_interest_OXA_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_6h_new <- rbind(collectri_df_tf_interest_OXA_6h_activators, collectri_df_tf_interest_OXA_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/OXA_vs_Control/6h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_12h_df, by="target")
  collectri_df_tf_interest_OXA_12h_activators <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_repressors <- collectri_df_tf_interest_OXA_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_12h_new <- rbind(collectri_df_tf_interest_OXA_12h_activators, collectri_df_tf_interest_OXA_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/OXA_vs_Control/12h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_OXA_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_18h_df, by="target")
  collectri_df_tf_interest_OXA_18h_activators <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_repressors <- collectri_df_tf_interest_OXA_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_18h_new <- rbind(collectri_df_tf_interest_OXA_18h_activators, collectri_df_tf_interest_OXA_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/OXA_vs_Control/18h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_OXA_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_OXA_Control_24h_df, by="target")
  collectri_df_tf_interest_OXA_24h_activators <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_repressors <- collectri_df_tf_interest_OXA_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_OXA_24h_new <- rbind(collectri_df_tf_interest_OXA_24h_activators, collectri_df_tf_interest_OXA_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_OXA_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/OXA_vs_Control/24h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  
  # CIS
  
  
  collectri_df_tf_interest_CIS_3h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_3h_df, by="target")
  collectri_df_tf_interest_CIS_3h_activators <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_repressors <- collectri_df_tf_interest_CIS_3h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_3h_new <- rbind(collectri_df_tf_interest_CIS_3h_activators, collectri_df_tf_interest_CIS_3h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_3h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/CIS_vs_Control/3h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_6h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_6h_df, by="target")
  collectri_df_tf_interest_CIS_6h_activators <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_repressors <- collectri_df_tf_interest_CIS_6h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_6h_new <- rbind(collectri_df_tf_interest_CIS_6h_activators, collectri_df_tf_interest_CIS_6h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_6h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/CIS_vs_Control/6h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_12h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_12h_df, by="target")
  collectri_df_tf_interest_CIS_12h_activators <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_repressors <- collectri_df_tf_interest_CIS_12h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_12h_new <- rbind(collectri_df_tf_interest_CIS_12h_activators, collectri_df_tf_interest_CIS_12h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_12h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/CIS_vs_Control/12h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
  collectri_df_tf_interest_CIS_18h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_18h_df, by="target")
  collectri_df_tf_interest_CIS_18h_activators <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_repressors <- collectri_df_tf_interest_CIS_18h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_18h_new <- rbind(collectri_df_tf_interest_CIS_18h_activators, collectri_df_tf_interest_CIS_18h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_18h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/CIS_vs_Control/18h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  collectri_df_tf_interest_CIS_24h <-  collectri_filtered_gautier_target_genes %>% filter(source == rownames(top_270_tf_targets_at_least_one_gautier_genes)[i]) %>% dplyr::mutate(mor = recode(mor, '1' = 'activator', '-1' = 'repressor')) %>% inner_join(.,t_values_diff_genes_CIS_Control_24h_df, by="target")
  collectri_df_tf_interest_CIS_24h_activators <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "activator") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_repressors <- collectri_df_tf_interest_CIS_24h %>% filter(mor == "repressor") %>% dplyr::arrange(desc(abs(t)))
  collectri_df_tf_interest_CIS_24h_new <- rbind(collectri_df_tf_interest_CIS_24h_activators, collectri_df_tf_interest_CIS_24h_repressors) %>% mutate(Gautier_gene = if_else(target %in% genes_for_dea_mouse$genes,"Yes","No"))
  write.csv(collectri_df_tf_interest_CIS_24h_new, paste0("Analysis_for_Gautier/Upstream_analysis/DecoupleR/Differentially_expressed_targets_of_the_TFs/201 - 270 TFs/CIS_vs_Control/24h/",rownames(top_270_tf_targets_at_least_one_gautier_genes)[i],"_targets_all_info.csv"), col.names = TRUE, row.names = FALSE)
  
  
  
  
}









