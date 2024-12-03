
# Intersect overdispersed MTX
get_plots_intersect_overdispersed_MTX_reactome <- function(){
  
  library(rRoma)
  library(tidyverse)
  all_overdispersed_reactome_MTX_3h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/MTX/all_overdispersed_MTX_3h.csv")
  all_overdispersed_reactome_MTX_6h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/MTX/all_overdispersed_MTX_6h.csv")
  all_overdispersed_reactome_MTX_12h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/MTX/all_overdispersed_MTX_12h.csv")
  all_overdispersed_reactome_MTX_18h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/MTX/all_overdispersed_MTX_18h.csv")
  all_overdispersed_reactome_MTX_24h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/MTX/all_overdispersed_MTX_24h.csv")
  
  intersect_overdispersed_reactome_MTX_3h_6h_12h <- intersect(intersect(all_overdispersed_reactome_MTX_3h$overdispersed_modules, all_overdispersed_reactome_MTX_6h$overdispersed_modules), all_overdispersed_reactome_MTX_12h$overdispersed_modules)
  intersect_overdispersed_reactome_MTX_all <- intersect(intersect_overdispersed_reactome_MTX_3h_6h_12h, intersect(all_overdispersed_reactome_MTX_18h$overdispersed_modules, all_overdispersed_reactome_MTX_24h$overdispersed_modules))
  
  roma_output_MTX_3h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/MTX/roma_output_MTX_3h.rds") 
  roma_output_MTX_6h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/MTX/roma_output_MTX_6h.rds") 
  roma_output_MTX_12h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/MTX/roma_output_MTX_12h.rds") 
  roma_output_MTX_18h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/MTX/roma_output_MTX_18h.rds") 
  roma_output_MTX_24h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/MTX/roma_output_MTX_24h.rds")
  
  overdispersed.modules_of_interest_MTX_3h <- which(roma_output_MTX_3h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_MTX_3h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_MTX_6h <- which(roma_output_MTX_6h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_MTX_6h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_MTX_12h <- which(roma_output_MTX_12h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_MTX_12h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_MTX_18h <- which(roma_output_MTX_18h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_MTX_18h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_MTX_24h <- which(roma_output_MTX_24h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_MTX_24h$ModuleMatrix[,"ppv L1"] <= 0.05)
  
  
  roma_output_MTX_all_time <- roma_output_MTX_3h
  
  roma_output_MTX_3h_df <- roma_output_MTX_3h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_6h_df <- roma_output_MTX_6h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_12h_df <- roma_output_MTX_12h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_18h_df <- roma_output_MTX_18h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_24h_df <- roma_output_MTX_24h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  
  roma_output_MTX_3h_6h_df <- inner_join(roma_output_MTX_3h_df, roma_output_MTX_6h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_MTX_3h_6h_12h_df <- inner_join(roma_output_MTX_3h_6h_df, roma_output_MTX_12h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_MTX_3h_6h_12h_18h_df <- inner_join(roma_output_MTX_3h_6h_12h_df, roma_output_MTX_18h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_all_df <- inner_join(roma_output_MTX_3h_6h_12h_18h_df, roma_output_MTX_24h_df[,c("modules", "ppv Median Exp")], by="modules") %>% column_to_rownames("modules") %>% as.matrix() 
  
  
  roma_output_MTX_3h_df_s <- roma_output_MTX_3h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_6h_df_s <- roma_output_MTX_6h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_12h_df_s <- roma_output_MTX_12h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_18h_df_s <- roma_output_MTX_18h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_MTX_24h_df_s <- roma_output_MTX_24h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  
  
  roma_output_MTX_3h_6h_df_s <- inner_join(roma_output_MTX_3h_df_s, roma_output_MTX_6h_df_s, by="modules")
  roma_output_MTX_3h_6h_12h_df_s <- inner_join(roma_output_MTX_3h_6h_df_s, roma_output_MTX_12h_df_s, by="modules")
  roma_output_MTX_3h_6h_12h_18h_df_s <- inner_join(roma_output_MTX_3h_6h_12h_df_s, roma_output_MTX_18h_df_s, by="modules")
  roma_output_all_df_s <- inner_join(roma_output_MTX_3h_6h_12h_18h_df_s, roma_output_MTX_24h_df_s, by="modules")
  
  roma_output_all_df_s <- roma_output_all_df_s %>% filter(modules %in% rownames(roma_output_all_df)) %>% column_to_rownames("modules") %>% as.matrix()
  
  
  
  roma_output_MTX_all_time$ModuleMatrix <- roma_output_all_df
  roma_output_MTX_all_time$SampleMatrix <- roma_output_all_df_s
  
  overdispersed <- which(roma_output_MTX_all_time$ModuleMatrix[, "ppv Median Exp"] > 0.05 & roma_output_MTX_all_time$ModuleMatrix[,"ppv L1"] <= 0.05)
  
  
  
  signaling_MTX <- intersect_overdispersed_reactome_MTX_all[grepl("SIGNALING|APOPTOSIS|DNA", intersect_overdispersed_reactome_MTX_all)]
  
  
  for (i in 1:length(signaling_MTX)){
    
    dir.create(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]))
    overdispersed_3h_MTX <- all_overdispersed_reactome_MTX_3h %>% filter(overdispersed_modules == signaling_MTX[i]) %>% dplyr::pull("L1")
    overdispersed_6h_MTX <- all_overdispersed_reactome_MTX_6h %>% filter(overdispersed_modules == signaling_MTX[i]) %>% dplyr::pull("L1")
    overdispersed_12h_MTX <- all_overdispersed_reactome_MTX_12h %>% filter(overdispersed_modules == signaling_MTX[i]) %>% dplyr::pull("L1")
    overdispersed_18h_MTX <- all_overdispersed_reactome_MTX_18h %>% filter(overdispersed_modules == signaling_MTX[i]) %>% dplyr::pull("L1")
    overdispersed_24h_MTX <- all_overdispersed_reactome_MTX_24h %>% filter(overdispersed_modules == signaling_MTX[i]) %>% dplyr::pull("L1")
    
    time <- c("3","6","12","18","24")
    time_ordered <- factor(time, order=TRUE, levels = c("3", "6", "12", "18", "24"))
    L1_explained <- c(overdispersed_3h_MTX, overdispersed_6h_MTX, overdispersed_12h_MTX, overdispersed_18h_MTX, overdispersed_24h_MTX)
    df <- data.frame(time_ordered, L1_explained)
    
    
    ggplot(df, mapping = aes(x=time_ordered, y=L1_explained, group=1)) + geom_point() + geom_line() + labs(title = signaling_MTX[i]) + theme(plot.title = element_text(face = "bold", hjust=0.5))
    ggsave(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","L1_explained.png"))
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","activity_scores_MTX_3h.png"), width=1000, height=502)
    Plot.Genesets.Samples(roma_output_MTX_3h, Selected = overdispersed.modules_of_interest_MTX_3h[which(names(overdispersed.modules_of_interest_MTX_3h) == paste("REACTOME", signaling_MTX[i], sep="_"))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_MTX[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","activity_scores_MTX_6h.png"))
    Plot.Genesets.Samples(roma_output_MTX_6h, Selected = overdispersed.modules_of_interest_MTX_6h[which(names(overdispersed.modules_of_interest_MTX_6h) == paste0("REACTOME_",signaling_MTX[i]))], Transpose = TRUE, cluster_cols = TRUE,HMTite = str_remove(signaling_MTX[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","activity_scores_MTX_12h.png"))
    Plot.Genesets.Samples(roma_output_MTX_12h, Selected = overdispersed.modules_of_interest_MTX_12h[which(names(overdispersed.modules_of_interest_MTX_12h) == paste0("REACTOME_",signaling_MTX[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_MTX[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","activity_scores_MTX_18h.png"))
    Plot.Genesets.Samples(roma_output_MTX_18h, Selected = overdispersed.modules_of_interest_MTX_18h[which(names(overdispersed.modules_of_interest_MTX_18h) == paste0("REACTOME_",signaling_MTX[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_MTX[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","activity_scores_MTX_24h.png"))
    Plot.Genesets.Samples(roma_output_MTX_24h, Selected = overdispersed.modules_of_interest_MTX_24h[which(names(overdispersed.modules_of_interest_MTX_24h) == paste0("REACTOME_",signaling_MTX[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_MTX[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/MTX/",signaling_MTX[i]),"/","activity_scores_MTX_all_time.png"), width = 700, height=800)
    Plot.Genesets.Samples(roma_output_MTX_all_time, Selected = overdispersed[which(names(overdispersed) == paste0("REACTOME_",signaling_MTX[i]))], Transpose = TRUE, HMTite = str_remove(signaling_MTX[i], "REACTOME_"), cluster_cols = TRUE)
    dev.off()
    
    
    
    
  }
  
}



# Intersect overdispersed OXA
get_plots_intersect_overdispersed_OXA_reactome <- function(){
  
  library(rRoma)
  library(tidyverse)
  all_overdispersed_reactome_OXA_3h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/OXA/all_overdispersed_OXA_3h.csv")
  all_overdispersed_reactome_OXA_6h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/OXA/all_overdispersed_OXA_6h.csv")
  all_overdispersed_reactome_OXA_12h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/OXA/all_overdispersed_OXA_12h.csv")
  all_overdispersed_reactome_OXA_18h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/OXA/all_overdispersed_OXA_18h.csv")
  all_overdispersed_reactome_OXA_24h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/OXA/all_overdispersed_OXA_24h.csv")
  
  intersect_overdispersed_reactome_OXA_3h_6h_12h <- intersect(intersect(all_overdispersed_reactome_OXA_3h$overdispersed_modules, all_overdispersed_reactome_OXA_6h$overdispersed_modules), all_overdispersed_reactome_OXA_12h$overdispersed_modules)
  intersect_overdispersed_reactome_OXA_all <- intersect(intersect_overdispersed_reactome_OXA_3h_6h_12h, intersect(all_overdispersed_reactome_OXA_18h$overdispersed_modules, all_overdispersed_reactome_OXA_24h$overdispersed_modules))
  
  roma_output_OXA_3h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/OXA/roma_output_OXA_3h.rds") 
  roma_output_OXA_6h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/OXA/roma_output_OXA_6h.rds") 
  roma_output_OXA_12h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/OXA/roma_output_OXA_12h.rds") 
  roma_output_OXA_18h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/OXA/roma_output_OXA_18h.rds") 
  roma_output_OXA_24h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/OXA/roma_output_OXA_24h.rds")
  
  overdispersed.modules_of_interest_OXA_3h <- which(roma_output_OXA_3h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_OXA_3h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_OXA_6h <- which(roma_output_OXA_6h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_OXA_6h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_OXA_12h <- which(roma_output_OXA_12h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_OXA_12h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_OXA_18h <- which(roma_output_OXA_18h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_OXA_18h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_OXA_24h <- which(roma_output_OXA_24h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_OXA_24h$ModuleMatrix[,"ppv L1"] <= 0.05)
  
  
  roma_output_OXA_all_time <- roma_output_OXA_3h
  
  roma_output_OXA_3h_df <- roma_output_OXA_3h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_6h_df <- roma_output_OXA_6h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_12h_df <- roma_output_OXA_12h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_18h_df <- roma_output_OXA_18h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_24h_df <- roma_output_OXA_24h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  
  roma_output_OXA_3h_6h_df <- inner_join(roma_output_OXA_3h_df, roma_output_OXA_6h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_OXA_3h_6h_12h_df <- inner_join(roma_output_OXA_3h_6h_df, roma_output_OXA_12h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_OXA_3h_6h_12h_18h_df <- inner_join(roma_output_OXA_3h_6h_12h_df, roma_output_OXA_18h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_all_df <- inner_join(roma_output_OXA_3h_6h_12h_18h_df, roma_output_OXA_24h_df[,c("modules", "ppv Median Exp")], by="modules") %>% column_to_rownames("modules") %>% as.matrix() 
  
  
  roma_output_OXA_3h_df_s <- roma_output_OXA_3h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_6h_df_s <- roma_output_OXA_6h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_12h_df_s <- roma_output_OXA_12h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_18h_df_s <- roma_output_OXA_18h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_OXA_24h_df_s <- roma_output_OXA_24h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  
  
  roma_output_OXA_3h_6h_df_s <- inner_join(roma_output_OXA_3h_df_s, roma_output_OXA_6h_df_s, by="modules")
  roma_output_OXA_3h_6h_12h_df_s <- inner_join(roma_output_OXA_3h_6h_df_s, roma_output_OXA_12h_df_s, by="modules")
  roma_output_OXA_3h_6h_12h_18h_df_s <- inner_join(roma_output_OXA_3h_6h_12h_df_s, roma_output_OXA_18h_df_s, by="modules")
  roma_output_all_df_s <- inner_join(roma_output_OXA_3h_6h_12h_18h_df_s, roma_output_OXA_24h_df_s, by="modules")
  
  roma_output_all_df_s <- roma_output_all_df_s %>% filter(modules %in% rownames(roma_output_all_df)) %>% column_to_rownames("modules") %>% as.matrix()
  
  
  
  roma_output_OXA_all_time$ModuleMatrix <- roma_output_all_df
  roma_output_OXA_all_time$SampleMatrix <- roma_output_all_df_s
  
  overdispersed <- which(roma_output_OXA_all_time$ModuleMatrix[, "ppv Median Exp"] > 0.05 & roma_output_OXA_all_time$ModuleMatrix[,"ppv L1"] <= 0.05)
  
  
  
  signaling_OXA <- intersect_overdispersed_reactome_OXA_all[grepl("SIGNALING|APOPTOSIS|DNA", intersect_overdispersed_reactome_OXA_all)]
  
  
  for (i in 1:length(signaling_OXA)){
    
    dir.create(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]))
    overdispersed_3h_OXA <- all_overdispersed_reactome_OXA_3h %>% filter(overdispersed_modules == signaling_OXA[i]) %>% dplyr::pull("L1")
    overdispersed_6h_OXA <- all_overdispersed_reactome_OXA_6h %>% filter(overdispersed_modules == signaling_OXA[i]) %>% dplyr::pull("L1")
    overdispersed_12h_OXA <- all_overdispersed_reactome_OXA_12h %>% filter(overdispersed_modules == signaling_OXA[i]) %>% dplyr::pull("L1")
    overdispersed_18h_OXA <- all_overdispersed_reactome_OXA_18h %>% filter(overdispersed_modules == signaling_OXA[i]) %>% dplyr::pull("L1")
    overdispersed_24h_OXA <- all_overdispersed_reactome_OXA_24h %>% filter(overdispersed_modules == signaling_OXA[i]) %>% dplyr::pull("L1")
    
    time <- c("3","6","12","18","24")
    time_ordered <- factor(time, order=TRUE, levels = c("3", "6", "12", "18", "24"))
    L1_explained <- c(overdispersed_3h_OXA, overdispersed_6h_OXA, overdispersed_12h_OXA, overdispersed_18h_OXA, overdispersed_24h_OXA)
    df <- data.frame(time_ordered, L1_explained)
    
    
    ggplot(df, mapping = aes(x=time_ordered, y=L1_explained, group=1)) + geom_point() + geom_line() + labs(title = signaling_OXA[i]) + theme(plot.title = element_text(face = "bold", hjust=0.5))
    ggsave(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","L1_explained.png"))
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","activity_scores_OXA_3h.png"), width=1000, height=502)
    Plot.Genesets.Samples(roma_output_OXA_3h, Selected = overdispersed.modules_of_interest_OXA_3h[which(names(overdispersed.modules_of_interest_OXA_3h) == paste("REACTOME", signaling_OXA[i], sep="_"))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_OXA[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","activity_scores_OXA_6h.png"))
    Plot.Genesets.Samples(roma_output_OXA_6h, Selected = overdispersed.modules_of_interest_OXA_6h[which(names(overdispersed.modules_of_interest_OXA_6h) == paste0("REACTOME_",signaling_OXA[i]))], Transpose = TRUE, cluster_cols = TRUE,HMTite = str_remove(signaling_OXA[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","activity_scores_OXA_12h.png"))
    Plot.Genesets.Samples(roma_output_OXA_12h, Selected = overdispersed.modules_of_interest_OXA_12h[which(names(overdispersed.modules_of_interest_OXA_12h) == paste0("REACTOME_",signaling_OXA[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_OXA[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","activity_scores_OXA_18h.png"))
    Plot.Genesets.Samples(roma_output_OXA_18h, Selected = overdispersed.modules_of_interest_OXA_18h[which(names(overdispersed.modules_of_interest_OXA_18h) == paste0("REACTOME_",signaling_OXA[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_OXA[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","activity_scores_OXA_24h.png"))
    Plot.Genesets.Samples(roma_output_OXA_24h, Selected = overdispersed.modules_of_interest_OXA_24h[which(names(overdispersed.modules_of_interest_OXA_24h) == paste0("REACTOME_",signaling_OXA[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_OXA[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/OXA/",signaling_OXA[i]),"/","activity_scores_OXA_all_time.png"), width = 700, height=800)
    Plot.Genesets.Samples(roma_output_OXA_all_time, Selected = overdispersed[which(names(overdispersed) == paste0("REACTOME_",signaling_OXA[i]))], Transpose = TRUE, HMTite = str_remove(signaling_OXA[i], "REACTOME_"), cluster_cols = TRUE)
    dev.off()
    
    
    
    
  }
  
}



# Intersect overdispersed CIS
get_plots_intersect_overdispersed_CIS_reactome <- function(){
  
  library(rRoma)
  library(tidyverse)
  all_overdispersed_reactome_CIS_3h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/CIS/all_overdispersed_CIS_3h.csv")
  all_overdispersed_reactome_CIS_6h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/CIS/all_overdispersed_CIS_6h.csv")
  all_overdispersed_reactome_CIS_12h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/CIS/all_overdispersed_CIS_12h.csv")
  all_overdispersed_reactome_CIS_18h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/CIS/all_overdispersed_CIS_18h.csv")
  all_overdispersed_reactome_CIS_24h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/CIS/all_overdispersed_CIS_24h.csv")
  
  intersect_overdispersed_reactome_CIS_3h_6h_12h <- intersect(intersect(all_overdispersed_reactome_CIS_3h$overdispersed_modules, all_overdispersed_reactome_CIS_6h$overdispersed_modules), all_overdispersed_reactome_CIS_12h$overdispersed_modules)
  intersect_overdispersed_reactome_CIS_all <- intersect(intersect_overdispersed_reactome_CIS_3h_6h_12h, intersect(all_overdispersed_reactome_CIS_18h$overdispersed_modules, all_overdispersed_reactome_CIS_24h$overdispersed_modules))
  
  roma_output_CIS_3h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/CIS/roma_output_CIS_3h.rds") 
  roma_output_CIS_6h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/CIS/roma_output_CIS_6h.rds") 
  roma_output_CIS_12h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/CIS/roma_output_CIS_12h.rds") 
  roma_output_CIS_18h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/CIS/roma_output_CIS_18h.rds") 
  roma_output_CIS_24h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/CIS/roma_output_CIS_24h.rds")
  
  overdispersed.modules_of_interest_CIS_3h <- which(roma_output_CIS_3h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_CIS_3h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_CIS_6h <- which(roma_output_CIS_6h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_CIS_6h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_CIS_12h <- which(roma_output_CIS_12h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_CIS_12h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_CIS_18h <- which(roma_output_CIS_18h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_CIS_18h$ModuleMatrix[,"ppv L1"] <= 0.05)
  overdispersed.modules_of_interest_CIS_24h <- which(roma_output_CIS_24h$ModuleMatrix[,"ppv Median Exp"] > 0.05 & roma_output_CIS_24h$ModuleMatrix[,"ppv L1"] <= 0.05)
  
  
  roma_output_CIS_all_time <- roma_output_CIS_3h
  
  roma_output_CIS_3h_df <- roma_output_CIS_3h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_6h_df <- roma_output_CIS_6h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_12h_df <- roma_output_CIS_12h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_18h_df <- roma_output_CIS_18h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_24h_df <- roma_output_CIS_24h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  
  roma_output_CIS_3h_6h_df <- inner_join(roma_output_CIS_3h_df, roma_output_CIS_6h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_CIS_3h_6h_12h_df <- inner_join(roma_output_CIS_3h_6h_df, roma_output_CIS_12h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_CIS_3h_6h_12h_18h_df <- inner_join(roma_output_CIS_3h_6h_12h_df, roma_output_CIS_18h_df[,c("modules", "ppv Median Exp")], by="modules")
  roma_output_all_df <- inner_join(roma_output_CIS_3h_6h_12h_18h_df, roma_output_CIS_24h_df[,c("modules", "ppv Median Exp")], by="modules") %>% column_to_rownames("modules") %>% as.matrix() 
  
  
  roma_output_CIS_3h_df_s <- roma_output_CIS_3h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_6h_df_s <- roma_output_CIS_6h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_12h_df_s <- roma_output_CIS_12h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_18h_df_s <- roma_output_CIS_18h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  roma_output_CIS_24h_df_s <- roma_output_CIS_24h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
  
  
  roma_output_CIS_3h_6h_df_s <- inner_join(roma_output_CIS_3h_df_s, roma_output_CIS_6h_df_s, by="modules")
  roma_output_CIS_3h_6h_12h_df_s <- inner_join(roma_output_CIS_3h_6h_df_s, roma_output_CIS_12h_df_s, by="modules")
  roma_output_CIS_3h_6h_12h_18h_df_s <- inner_join(roma_output_CIS_3h_6h_12h_df_s, roma_output_CIS_18h_df_s, by="modules")
  roma_output_all_df_s <- inner_join(roma_output_CIS_3h_6h_12h_18h_df_s, roma_output_CIS_24h_df_s, by="modules")
  
  roma_output_all_df_s <- roma_output_all_df_s %>% filter(modules %in% rownames(roma_output_all_df)) %>% column_to_rownames("modules") %>% as.matrix()
  
  
  
  roma_output_CIS_all_time$ModuleMatrix <- roma_output_all_df
  roma_output_CIS_all_time$SampleMatrix <- roma_output_all_df_s
  
  overdispersed <- which(roma_output_CIS_all_time$ModuleMatrix[, "ppv Median Exp"] > 0.05 & roma_output_CIS_all_time$ModuleMatrix[,"ppv L1"] <= 0.05)
  
  
  
  signaling_CIS <- intersect_overdispersed_reactome_CIS_all[grepl("SIGNALING|APOPTOSIS|DNA", intersect_overdispersed_reactome_CIS_all)]
  
  
  for (i in 1:length(signaling_CIS)){
    
    dir.create(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]))
    overdispersed_3h_CIS <- all_overdispersed_reactome_CIS_3h %>% filter(overdispersed_modules == signaling_CIS[i]) %>% dplyr::pull("L1")
    overdispersed_6h_CIS <- all_overdispersed_reactome_CIS_6h %>% filter(overdispersed_modules == signaling_CIS[i]) %>% dplyr::pull("L1")
    overdispersed_12h_CIS <- all_overdispersed_reactome_CIS_12h %>% filter(overdispersed_modules == signaling_CIS[i]) %>% dplyr::pull("L1")
    overdispersed_18h_CIS <- all_overdispersed_reactome_CIS_18h %>% filter(overdispersed_modules == signaling_CIS[i]) %>% dplyr::pull("L1")
    overdispersed_24h_CIS <- all_overdispersed_reactome_CIS_24h %>% filter(overdispersed_modules == signaling_CIS[i]) %>% dplyr::pull("L1")
    
    time <- c("3","6","12","18","24")
    time_ordered <- factor(time, order=TRUE, levels = c("3", "6", "12", "18", "24"))
    L1_explained <- c(overdispersed_3h_CIS, overdispersed_6h_CIS, overdispersed_12h_CIS, overdispersed_18h_CIS, overdispersed_24h_CIS)
    df <- data.frame(time_ordered, L1_explained)
    
    
    ggplot(df, mapping = aes(x=time_ordered, y=L1_explained, group=1)) + geom_point() + geom_line() + labs(title = signaling_CIS[i]) + theme(plot.title = element_text(face = "bold", hjust=0.5))
    ggsave(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","L1_explained.png"))
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","activity_scores_CIS_3h.png"), width=1000, height=502)
    Plot.Genesets.Samples(roma_output_CIS_3h, Selected = overdispersed.modules_of_interest_CIS_3h[which(names(overdispersed.modules_of_interest_CIS_3h) == paste("REACTOME", signaling_CIS[i], sep="_"))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_CIS[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","activity_scores_CIS_6h.png"))
    Plot.Genesets.Samples(roma_output_CIS_6h, Selected = overdispersed.modules_of_interest_CIS_6h[which(names(overdispersed.modules_of_interest_CIS_6h) == paste0("REACTOME_",signaling_CIS[i]))], Transpose = TRUE, cluster_cols = TRUE,HMTite = str_remove(signaling_CIS[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","activity_scores_CIS_12h.png"))
    Plot.Genesets.Samples(roma_output_CIS_12h, Selected = overdispersed.modules_of_interest_CIS_12h[which(names(overdispersed.modules_of_interest_CIS_12h) == paste0("REACTOME_",signaling_CIS[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_CIS[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","activity_scores_CIS_18h.png"))
    Plot.Genesets.Samples(roma_output_CIS_18h, Selected = overdispersed.modules_of_interest_CIS_18h[which(names(overdispersed.modules_of_interest_CIS_18h) == paste0("REACTOME_",signaling_CIS[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_CIS[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","activity_scores_CIS_24h.png"))
    Plot.Genesets.Samples(roma_output_CIS_24h, Selected = overdispersed.modules_of_interest_CIS_24h[which(names(overdispersed.modules_of_interest_CIS_24h) == paste0("REACTOME_",signaling_CIS[i]))], Transpose = TRUE, cluster_cols = TRUE, HMTite = str_remove(signaling_CIS[i], "REACTOME_"))
    dev.off()
    
    png(paste0(paste0("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/Intersection_overdispersed_pathways/CIS/",signaling_CIS[i]),"/","activity_scores_CIS_all_time.png"), width = 700, height=800)
    Plot.Genesets.Samples(roma_output_CIS_all_time, Selected = overdispersed[which(names(overdispersed) == paste0("REACTOME_",signaling_CIS[i]))], Transpose = TRUE, HMTite = str_remove(signaling_CIS[i], "REACTOME_"), cluster_cols = TRUE)
    dev.off()
    
    
    
    
  }
  
}










#Go envoyer


roma_output_MTX_3h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/MTX/roma_output_MTX_3h.rds") 
roma_output_MTX_6h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/MTX/roma_output_MTX_6h.rds") 
roma_output_MTX_12h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/MTX/roma_output_MTX_12h.rds") 
roma_output_MTX_18h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/MTX/roma_output_MTX_18h.rds") 
roma_output_MTX_24h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/MTX/roma_output_MTX_24h.rds")





roma_output_MTX_3h_df <- roma_output_MTX_3h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_3h_df$modules <- sub("REACTOME_", "", roma_output_MTX_3h_df$modules)


roma_output_MTX_6h_df <- roma_output_MTX_6h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 
roma_output_MTX_6h_df$modules <- sub("REACTOME_", "", roma_output_MTX_6h_df$modules)


roma_output_MTX_12h_df <- roma_output_MTX_12h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_12h_df$modules <- sub("REACTOME_", "", roma_output_MTX_12h_df$modules)


roma_output_MTX_18h_df <- roma_output_MTX_18h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_18h_df$modules <- sub("REACTOME_", "", roma_output_MTX_18h_df$modules)


roma_output_MTX_24h_df <- roma_output_MTX_24h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_24h_df$modules <- sub("REACTOME_", "", roma_output_MTX_24h_df$modules)



roma_output_MTX_3h_6h_df_new_MTX <- inner_join(roma_output_MTX_3h_df, roma_output_MTX_6h_df, by="modules")
roma_output_MTX_3h_6h_12h_df_new_MTX <- inner_join(roma_output_MTX_3h_6h_df_new_MTX, roma_output_MTX_12h_df, by="modules")
roma_output_MTX_3h_6h_12h_18h_df_new_MTX <- inner_join(roma_output_MTX_3h_6h_12h_df_new_MTX, roma_output_MTX_18h_df, by="modules")
roma_output_all_df_new_MTX <- inner_join(roma_output_MTX_3h_6h_12h_18h_df_new_MTX, roma_output_MTX_24h_df, by="modules")



all_overdispersed_reactome_MTX_3h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/MTX/all_overdispersed_MTX_3h.csv")
all_overdispersed_reactome_MTX_6h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/MTX/all_overdispersed_MTX_6h.csv")
all_overdispersed_reactome_MTX_12h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/MTX/all_overdispersed_MTX_12h.csv")
all_overdispersed_reactome_MTX_18h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/MTX/all_overdispersed_MTX_18h.csv")
all_overdispersed_reactome_MTX_24h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/MTX/all_overdispersed_MTX_24h.csv")


b <-  which(roma_output_all_df_new_MTX$modules %in% unique(c(all_overdispersed_reactome_MTX_3h$overdispersed_modules,all_overdispersed_reactome_MTX_6h$overdispersed_modules, all_overdispersed_reactome_MTX_12h$overdispersed_modules,all_overdispersed_reactome_MTX_18h$overdispersed_modules,all_overdispersed_reactome_MTX_24h$overdispersed_modules)))


n <- roma_output_all_df_new_MTX$modules[b]
roma_output_all_df_new_MTX <- roma_output_all_df_new_MTX %>% dplyr::filter(modules %in% n) %>% column_to_rownames("modules") %>% as.matrix()


roma_output_MTX_3h_df_s_new_MTX <- roma_output_MTX_3h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_6h_df_s_new_MTX <- roma_output_MTX_6h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_12h_df_s_new_MTX <- roma_output_MTX_12h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_18h_df_s_new_MTX <- roma_output_MTX_18h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_MTX_24h_df_s_new_MTX <- roma_output_MTX_24h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")


roma_output_MTX_3h_6h_df_s_new_MTX <- inner_join(roma_output_MTX_3h_df_s_new_MTX, roma_output_MTX_6h_df_s_new_MTX, by="modules")
roma_output_MTX_3h_6h_12h_df_s_new_MTX <- inner_join(roma_output_MTX_3h_6h_df_s_new_MTX, roma_output_MTX_12h_df_s_new_MTX, by="modules")
roma_output_MTX_3h_6h_12h_18h_df_s_new_MTX <- inner_join(roma_output_MTX_3h_6h_12h_df_s_new_MTX, roma_output_MTX_18h_df_s_new_MTX, by="modules")
roma_output_all_df_s_new_MTX <- inner_join(roma_output_MTX_3h_6h_12h_18h_df_s_new_MTX, roma_output_MTX_24h_df_s_new_MTX, by="modules")

roma_output_all_df_s_new_MTX$modules <- sub("REACTOME_", "", roma_output_all_df_s_new_MTX$modules)
roma_output_all_df_s_new_MTX <- roma_output_all_df_s_new_MTX %>% filter(modules %in% n) %>% column_to_rownames("modules") %>% as.matrix()


roma_output_for_MTX_OXA_CIS_1  <- roma_output_MTX_3h
roma_output_for_MTX_OXA_CIS_1$ModuleMatrix <- roma_output_all_df_new_MTX
roma_output_for_MTX_OXA_CIS_1$SampleMatrix <- roma_output_all_df_s_new_MTX


Plot.Genesets.Samples(roma_output_for_MTX_OXA_CIS, Selected = which(roma_output_for_MTX_OXA_CIS$ModuleMatrix[, "ppv L1.x"] <= 0.05)[1:50], cluster_cols = TRUE)








roma_output_OXA_3h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/OXA/roma_output_OXA_3h.rds") 
roma_output_OXA_6h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/OXA/roma_output_OXA_6h.rds") 
roma_output_OXA_12h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/OXA/roma_output_OXA_12h.rds") 
roma_output_OXA_18h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/OXA/roma_output_OXA_18h.rds") 
roma_output_OXA_24h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/OXA/roma_output_OXA_24h.rds")





roma_output_OXA_3h_df <- roma_output_OXA_3h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_3h_df$modules <- sub("REACTOME_", "", roma_output_OXA_3h_df$modules)


roma_output_OXA_6h_df <- roma_output_OXA_6h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 
roma_output_OXA_6h_df$modules <- sub("REACTOME_", "", roma_output_OXA_6h_df$modules)


roma_output_OXA_12h_df <- roma_output_OXA_12h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_12h_df$modules <- sub("REACTOME_", "", roma_output_OXA_12h_df$modules)


roma_output_OXA_18h_df <- roma_output_OXA_18h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_18h_df$modules <- sub("REACTOME_", "", roma_output_OXA_18h_df$modules)


roma_output_OXA_24h_df <- roma_output_OXA_24h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_24h_df$modules <- sub("REACTOME_", "", roma_output_OXA_24h_df$modules)



roma_output_OXA_3h_6h_df_new_OXA <- inner_join(roma_output_OXA_3h_df, roma_output_OXA_6h_df, by="modules")
roma_output_OXA_3h_6h_12h_df_new_OXA <- inner_join(roma_output_OXA_3h_6h_df_new_OXA, roma_output_OXA_12h_df, by="modules")
roma_output_OXA_3h_6h_12h_18h_df_new_OXA <- inner_join(roma_output_OXA_3h_6h_12h_df_new_OXA, roma_output_OXA_18h_df, by="modules")
roma_output_all_df_new_OXA <- inner_join(roma_output_OXA_3h_6h_12h_18h_df_new_OXA, roma_output_OXA_24h_df, by="modules")



all_overdispersed_reactome_OXA_3h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/OXA/all_overdispersed_OXA_3h.csv")
all_overdispersed_reactome_OXA_6h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/OXA/all_overdispersed_OXA_6h.csv")
all_overdispersed_reactome_OXA_12h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/OXA/all_overdispersed_OXA_12h.csv")
all_overdispersed_reactome_OXA_18h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/OXA/all_overdispersed_OXA_18h.csv")
all_overdispersed_reactome_OXA_24h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/OXA/all_overdispersed_OXA_24h.csv")


b <-  which(roma_output_all_df_new_OXA$modules %in% unique(c(all_overdispersed_reactome_OXA_3h$overdispersed_modules,all_overdispersed_reactome_OXA_6h$overdispersed_modules, all_overdispersed_reactome_OXA_12h$overdispersed_modules,all_overdispersed_reactome_OXA_18h$overdispersed_modules,all_overdispersed_reactome_OXA_24h$overdispersed_modules)))


n <- roma_output_all_df_new_OXA$modules[b]
roma_output_all_df_new_OXA <- roma_output_all_df_new_OXA %>% dplyr::filter(modules %in% n) %>% column_to_rownames("modules") %>% as.matrix()


roma_output_OXA_3h_df_s_new_OXA <- roma_output_OXA_3h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_6h_df_s_new_OXA <- roma_output_OXA_6h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_12h_df_s_new_OXA <- roma_output_OXA_12h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_18h_df_s_new_OXA <- roma_output_OXA_18h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_OXA_24h_df_s_new_OXA <- roma_output_OXA_24h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")


roma_output_OXA_3h_6h_df_s_new_OXA <- inner_join(roma_output_OXA_3h_df_s_new_OXA, roma_output_OXA_6h_df_s_new_OXA, by="modules")
roma_output_OXA_3h_6h_12h_df_s_new_OXA <- inner_join(roma_output_OXA_3h_6h_df_s_new_OXA, roma_output_OXA_12h_df_s_new_OXA, by="modules")
roma_output_OXA_3h_6h_12h_18h_df_s_new_OXA <- inner_join(roma_output_OXA_3h_6h_12h_df_s_new_OXA, roma_output_OXA_18h_df_s_new_OXA, by="modules")
roma_output_all_df_s_new_OXA <- inner_join(roma_output_OXA_3h_6h_12h_18h_df_s_new_OXA, roma_output_OXA_24h_df_s_new_OXA, by="modules")

roma_output_all_df_s_new_OXA$modules <- sub("REACTOME_", "", roma_output_all_df_s_new_OXA$modules)
roma_output_all_df_s_new_OXA <- roma_output_all_df_s_new_OXA %>% filter(modules %in% n) %>% column_to_rownames("modules") %>% as.matrix()


roma_output_for_MTX_OXA_CIS_2  <- roma_output_OXA_3h
roma_output_for_MTX_OXA_CIS_2$ModuleMatrix <- roma_output_all_df_new_OXA
roma_output_for_MTX_OXA_CIS_2$SampleMatrix <- roma_output_all_df_s_new_OXA


Plot.Genesets.Samples(roma_output_for_MTX_OXA_CIS_2, Selected = which(roma_output_for_MTX_OXA_CIS_2$ModuleMatrix[, "ppv L1.x"] <= 0.05)[1:50], cluster_cols = TRUE)






roma_output_CIS_3h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/CIS/roma_output_CIS_3h.rds") 
roma_output_CIS_6h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/CIS/roma_output_CIS_6h.rds") 
roma_output_CIS_12h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/CIS/roma_output_CIS_12h.rds") 
roma_output_CIS_18h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/CIS/roma_output_CIS_18h.rds") 
roma_output_CIS_24h <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/CIS/roma_output_CIS_24h.rds")





roma_output_CIS_3h_df <- roma_output_CIS_3h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_3h_df$modules <- sub("REACTOME_", "", roma_output_CIS_3h_df$modules)


roma_output_CIS_6h_df <- roma_output_CIS_6h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 
roma_output_CIS_6h_df$modules <- sub("REACTOME_", "", roma_output_CIS_6h_df$modules)


roma_output_CIS_12h_df <- roma_output_CIS_12h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_12h_df$modules <- sub("REACTOME_", "", roma_output_CIS_12h_df$modules)


roma_output_CIS_18h_df <- roma_output_CIS_18h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_18h_df$modules <- sub("REACTOME_", "", roma_output_CIS_18h_df$modules)


roma_output_CIS_24h_df <- roma_output_CIS_24h$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_24h_df$modules <- sub("REACTOME_", "", roma_output_CIS_24h_df$modules)



roma_output_CIS_3h_6h_df_new <- inner_join(roma_output_CIS_3h_df, roma_output_CIS_6h_df, by="modules")
roma_output_CIS_3h_6h_12h_df_new <- inner_join(roma_output_CIS_3h_6h_df_new, roma_output_CIS_12h_df, by="modules")
roma_output_CIS_3h_6h_12h_18h_df_new <- inner_join(roma_output_CIS_3h_6h_12h_df_new, roma_output_CIS_18h_df, by="modules")
roma_output_all_df_new_CIS <- inner_join(roma_output_CIS_3h_6h_12h_18h_df_new, roma_output_CIS_24h_df, by="modules")



all_overdispersed_reactome_CIS_3h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/3h/CIS/all_overdispersed_CIS_3h.csv")
all_overdispersed_reactome_CIS_6h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/6h/CIS/all_overdispersed_CIS_6h.csv")
all_overdispersed_reactome_CIS_12h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/12h/CIS/all_overdispersed_CIS_12h.csv")
all_overdispersed_reactome_CIS_18h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/18h/CIS/all_overdispersed_CIS_18h.csv")
all_overdispersed_reactome_CIS_24h <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/24h/CIS/all_overdispersed_CIS_24h.csv")


b <-  which(roma_output_all_df_new_CIS$modules %in% unique(c(all_overdispersed_reactome_CIS_3h$overdispersed_modules,all_overdispersed_reactome_CIS_6h$overdispersed_modules, all_overdispersed_reactome_CIS_12h$overdispersed_modules,all_overdispersed_reactome_CIS_18h$overdispersed_modules,all_overdispersed_reactome_CIS_24h$overdispersed_modules)))


n <- roma_output_all_df_new_CIS$modules[b]
roma_output_all_df_new_CIS <- roma_output_all_df_new_CIS %>% dplyr::filter(modules %in% n) %>% column_to_rownames("modules") %>% as.matrix()


roma_output_CIS_3h_df_s_new <- roma_output_CIS_3h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_6h_df_s_new <- roma_output_CIS_6h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_12h_df_s_new <- roma_output_CIS_12h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_18h_df_s_new <- roma_output_CIS_18h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")
roma_output_CIS_24h_df_s_new <- roma_output_CIS_24h$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules")


roma_output_CIS_3h_6h_df_s_new <- inner_join(roma_output_CIS_3h_df_s_new, roma_output_CIS_6h_df_s_new, by="modules")
roma_output_CIS_3h_6h_12h_df_s_new <- inner_join(roma_output_CIS_3h_6h_df_s_new, roma_output_CIS_12h_df_s_new, by="modules")
roma_output_CIS_3h_6h_12h_18h_df_s_new <- inner_join(roma_output_CIS_3h_6h_12h_df_s_new, roma_output_CIS_18h_df_s_new, by="modules")
roma_output_all_df_s_new_CIS <- inner_join(roma_output_CIS_3h_6h_12h_18h_df_s_new, roma_output_CIS_24h_df_s_new, by="modules")

roma_output_all_df_s_new_CIS$modules <- sub("REACTOME_", "", roma_output_all_df_s_new_CIS$modules)
roma_output_all_df_s_new_CIS <- roma_output_all_df_s_new_CIS %>% filter(modules %in% n) %>% column_to_rownames("modules") %>% as.matrix()


roma_output_for_MTX_OXA_CIS_3  <- roma_output_CIS_3h
roma_output_for_MTX_OXA_CIS_3$ModuleMatrix <- roma_output_all_df_new_CIS
roma_output_for_MTX_OXA_CIS_3$SampleMatrix <- roma_output_all_df_s_new_CIS


Plot.Genesets.Samples(roma_output_for_MTX_OXA_CIS_3, Selected = which(roma_output_for_MTX_OXA_CIS_3$ModuleMatrix[, "ppv L1.x"] <= 0.05)[1:50], cluster_cols = TRUE)


intersect_MTX_OXA_CIS_overdispersed_at_least_one <- intersect(rownames(roma_output_for_MTX_OXA_CIS_1$ModuleMatrix), intersect(rownames(roma_output_for_MTX_OXA_CIS_2$ModuleMatrix), rownames(roma_output_for_MTX_OXA_CIS_3$ModuleMatrix)))

roma_output_final_MTX <- roma_output_MTX_3h
roma_output_final_MTX$ModuleMatrix <- roma_output_all_df_new_MTX
roma_output_final_MTX$SampleMatrix <- roma_output_all_df_s_new_MTX


roma_output_final_OXA <- roma_output_OXA_3h
roma_output_final_OXA$ModuleMatrix <- roma_output_all_df_new_OXA
roma_output_final_OXA$SampleMatrix <- roma_output_all_df_s_new_OXA

roma_output_final_CIS <- roma_output_CIS_3h
roma_output_final_CIS$ModuleMatrix <- roma_output_all_df_new_CIS
roma_output_final_CIS$SampleMatrix <- roma_output_all_df_s_new_CIS



roma_output_final_MTX$ModuleMatrix <- roma_output_final_MTX$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") %>% filter(modules %in% intersect_MTX_OXA_CIS_overdispersed_at_least_one) %>% column_to_rownames("modules") %>% as.matrix()
saveRDS(roma_output_final_MTX$ModuleMatrix, file="/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_MTX/roma_output_final_MTX_ModuleMatrix.rds")

roma_output_final_MTX$SampleMatrix <- roma_output_final_MTX$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules") %>% filter(modules %in% intersect_MTX_OXA_CIS_overdispersed_at_least_one) %>% column_to_rownames("modules") %>% as.matrix()
saveRDS(roma_output_final_MTX$SampleMatrix, file="/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_MTX/roma_output_final_MTX_SampleMatrix.rds")


roma_output_final_OXA$ModuleMatrix <- roma_output_final_OXA$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") %>% filter(modules %in% intersect_MTX_OXA_CIS_overdispersed_at_least_one) %>% column_to_rownames("modules") %>% as.matrix()
saveRDS(roma_output_final_OXA$ModuleMatrix, file="/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_OXA/roma_output_final_OXA_ModuleMatrix.rds")

roma_output_final_OXA$SampleMatrix <- roma_output_final_OXA$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules") %>% filter(modules %in% intersect_MTX_OXA_CIS_overdispersed_at_least_one) %>% column_to_rownames("modules") %>% as.matrix()
saveRDS(roma_output_final_OXA$SampleMatrix, file="/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_OXA/roma_output_final_OXA_SampleMatrix.rds")

roma_output_final_CIS$ModuleMatrix <- roma_output_final_CIS$ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") %>% filter(modules %in% intersect_MTX_OXA_CIS_overdispersed_at_least_one) %>% column_to_rownames("modules") %>% as.matrix()
saveRDS(roma_output_final_CIS$ModuleMatrix, file="/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_CIS/roma_output_final_CIS_ModuleMatrix.rds")

roma_output_final_CIS$SampleMatrix <- roma_output_final_CIS$SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules") %>% filter(modules %in% intersect_MTX_OXA_CIS_overdispersed_at_least_one) %>% column_to_rownames("modules") %>% as.matrix()
saveRDS(roma_output_final_CIS$SampleMatrix, file="/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_CIS/roma_output_final_CIS_SampleMatrix.rds")




#Plot.Genesets.Samples(roma_output_final_CIS, Selected = which(roma_output_final_CIS$ModuleMatrix[, "ppv L1.x"] >= 0)[1:50])











# final roma analysis

roma_output_final_MTX_ModuleMatrix <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_MTX/roma_output_final_MTX_ModuleMatrix.rds")
roma_output_final_MTX_SampleMatrix <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_MTX/roma_output_final_MTX_SampleMatrix.rds")

roma_output_final_OXA_ModuleMatrix <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_OXA/roma_output_final_OXA_ModuleMatrix.rds")
roma_output_final_OXA_SampleMatrix <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_OXA/roma_output_final_OXA_SampleMatrix.rds")

roma_output_final_CIS_ModuleMatrix <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_CIS/roma_output_final_CIS_ModuleMatrix.rds")
roma_output_final_CIS_SampleMatrix <- readRDS("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_CIS/roma_output_final_CIS_SampleMatrix.rds")


roma_output_final_MTX_ModuleMatrix_df <- roma_output_final_MTX_ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 
roma_output_final_MTX_SampleMatrix_df <- roma_output_final_MTX_SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 

for (i in 1:length(rownames(roma_output_final_MTX_ModuleMatrix))){
  if (roma_output_final_MTX_ModuleMatrix_df[i, "ppv L1.x"] > 0.05){
    
    roma_output_final_MTX_SampleMatrix_df[i, 2:7] <- 0
    
  }
  
  
  if (roma_output_final_MTX_ModuleMatrix_df[i, "ppv L1.y"] > 0.05){
    
    roma_output_final_MTX_SampleMatrix_df[i, 8:13] <- 0
    
  }
  
  
  if (roma_output_final_MTX_ModuleMatrix_df[i, "ppv L1.x.x"] > 0.05){
    
    roma_output_final_MTX_SampleMatrix_df[i, 14:19] <- 0
    
  }
  
  
  if (roma_output_final_MTX_ModuleMatrix_df[i, "ppv L1.y.y"] > 0.05){
    
    roma_output_final_MTX_SampleMatrix_df[i, 20:25] <- 0
    
  }
  
  
  if (roma_output_final_MTX_ModuleMatrix_df[i, "ppv L1"] > 0.05){
    
    roma_output_final_MTX_SampleMatrix_df[i, 26:31] <- 0
    
  }
  
  
}





roma_output_final_OXA_ModuleMatrix_df <- roma_output_final_OXA_ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 
roma_output_final_OXA_SampleMatrix_df <- roma_output_final_OXA_SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 

for (i in 1:length(rownames(roma_output_final_OXA_ModuleMatrix))){
  if (roma_output_final_OXA_ModuleMatrix_df[i, "ppv L1.x"] > 0.05){
    
    roma_output_final_OXA_SampleMatrix_df[i, 2:7] <- 0
    
  }
  
  
  if (roma_output_final_OXA_ModuleMatrix_df[i, "ppv L1.y"] > 0.05){
    
    roma_output_final_OXA_SampleMatrix_df[i, 8:13] <- 0
    
  }
  
  
  if (roma_output_final_OXA_ModuleMatrix_df[i, "ppv L1.x.x"] > 0.05){
    
    roma_output_final_OXA_SampleMatrix_df[i, 14:19] <- 0
    
  }
  
  
  if (roma_output_final_OXA_ModuleMatrix_df[i, "ppv L1.y.y"] > 0.05){
    
    roma_output_final_OXA_SampleMatrix_df[i, 20:25] <- 0
    
  }
  
  
  if (roma_output_final_OXA_ModuleMatrix_df[i, "ppv L1"] > 0.05){
    
    roma_output_final_OXA_SampleMatrix_df[i, 26:31] <- 0
    
  }
  
  
}







roma_output_final_CIS_ModuleMatrix_df <- roma_output_final_CIS_ModuleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 
roma_output_final_CIS_SampleMatrix_df <- roma_output_final_CIS_SampleMatrix %>% as.data.frame() %>% rownames_to_column("modules") 



for (i in 1:length(rownames(roma_output_final_CIS_ModuleMatrix))){
  if (roma_output_final_CIS_ModuleMatrix_df[i, "ppv L1.x"] > 0.05){
    
    roma_output_final_CIS_SampleMatrix_df[i, 2:7] <- 0
    
  }
  
  
  if (roma_output_final_CIS_ModuleMatrix_df[i, "ppv L1.y"] > 0.05){
    
    roma_output_final_CIS_SampleMatrix_df[i, 8:13] <- 0
    
  }
  
  
  if (roma_output_final_CIS_ModuleMatrix_df[i, "ppv L1.x.x"] > 0.05){
    
    roma_output_final_CIS_SampleMatrix_df[i, 14:19] <- 0
    
  }
  
  
  if (roma_output_final_CIS_ModuleMatrix_df[i, "ppv L1.y.y"] > 0.05){
    
    roma_output_final_CIS_SampleMatrix_df[i, 20:25] <- 0
    
  }
  
  
  if (roma_output_final_CIS_ModuleMatrix_df[i, "ppv L1"] > 0.05){
    
    roma_output_final_CIS_SampleMatrix_df[i, 26:31] <- 0
    
  }
  
  
}



roma_output_final_MTX_ModuleMatrix_mat <- roma_output_final_MTX_ModuleMatrix_df %>% column_to_rownames("modules") %>% as.matrix()
roma_output_final_MTX_SampleMatrix_mat <- roma_output_final_MTX_SampleMatrix_df %>% column_to_rownames("modules") %>% as.matrix()


roma_output_final_OXA_ModuleMatrix_mat <- roma_output_final_OXA_ModuleMatrix_df %>% column_to_rownames("modules") %>% as.matrix()
roma_output_final_OXA_SampleMatrix_mat <- roma_output_final_OXA_SampleMatrix_df %>% column_to_rownames("modules") %>% as.matrix()


roma_output_final_CIS_ModuleMatrix_mat <- roma_output_final_CIS_ModuleMatrix_df %>% column_to_rownames("modules") %>% as.matrix()
roma_output_final_CIS_SampleMatrix_mat <- roma_output_final_CIS_SampleMatrix_df %>% column_to_rownames("modules") %>% as.matrix()


roma_output_final_MTX <- roma_output_MTX_3h
roma_output_final_MTX$ModuleMatrix <- roma_output_final_MTX_ModuleMatrix_mat
roma_output_final_MTX$SampleMatrix <- roma_output_final_MTX_SampleMatrix_mat

roma_output_final_OXA <- roma_output_OXA_3h
roma_output_final_OXA$ModuleMatrix <- roma_output_final_OXA_ModuleMatrix_mat
roma_output_final_OXA$SampleMatrix <- roma_output_final_OXA_SampleMatrix_mat


roma_output_final_CIS <- roma_output_CIS_3h
roma_output_final_CIS$ModuleMatrix <- roma_output_final_CIS_ModuleMatrix_mat
roma_output_final_CIS$SampleMatrix <- roma_output_final_CIS_SampleMatrix_mat


roma_output_final_MTX_treated <- roma_output_MTX_3h
roma_output_final_MTX_treated$ModuleMatrix <- roma_output_final_MTX_ModuleMatrix_mat
roma_output_final_MTX_treated$SampleMatrix <- roma_output_final_MTX_SampleMatrix_mat %>% as.data.frame() %>% dplyr::select(!contains("Control")) %>% as.matrix()
library(gsubfn)
colnames(roma_output_final_MTX_treated$SampleMatrix)<- strapplyc(X = colnames(roma_output_final_MTX_treated$SampleMatrix),
pattern = ".+\\_MTX.+\\_(.+h)", simplify = TRUE)


roma_output_final_MTX_control <- roma_output_MTX_3h
roma_output_final_MTX_control$ModuleMatrix <- roma_output_final_MTX_ModuleMatrix_mat
roma_output_final_MTX_control$SampleMatrix <- roma_output_final_MTX_SampleMatrix_mat %>% as.data.frame() %>% dplyr::select(contains("Control")) %>% as.matrix()
colnames(roma_output_final_MTX_control$SampleMatrix)<- strapplyc(X = colnames(roma_output_final_MTX_control$SampleMatrix),
                                                                 pattern = ".+\\_Control.+\\_(.+h)", simplify = TRUE)




roma_output_final_OXA_treated <- roma_output_OXA_3h
roma_output_final_OXA_treated$ModuleMatrix <- roma_output_final_OXA_ModuleMatrix_mat
roma_output_final_OXA_treated$SampleMatrix <- roma_output_final_OXA_SampleMatrix_mat %>% as.data.frame() %>% dplyr::select(!contains("Control")) %>% as.matrix()
library(gsubfn)
colnames(roma_output_final_OXA_treated$SampleMatrix)<- strapplyc(X = colnames(roma_output_final_OXA_treated$SampleMatrix),
                                                                 pattern = ".+\\_OXA.+\\_(.+h)", simplify = TRUE)

roma_output_final_OXA_control <- roma_output_OXA_3h
roma_output_final_OXA_control$ModuleMatrix <- roma_output_final_OXA_ModuleMatrix_mat
roma_output_final_OXA_control$SampleMatrix <- roma_output_final_OXA_SampleMatrix_mat %>% as.data.frame() %>% dplyr::select(contains("Control")) %>% as.matrix()
colnames(roma_output_final_OXA_control$SampleMatrix)<- strapplyc(X = colnames(roma_output_final_OXA_control$SampleMatrix),
                                                                 pattern = ".+\\_Control.+\\_(.+h)", simplify = TRUE)




roma_output_final_CIS_treated <- roma_output_CIS_3h
roma_output_final_CIS_treated$ModuleMatrix <- roma_output_final_CIS_ModuleMatrix_mat
roma_output_final_CIS_treated$SampleMatrix <- roma_output_final_CIS_SampleMatrix_mat %>% as.data.frame() %>% dplyr::select(!contains("Control")) %>% as.matrix()
library(gsubfn)
colnames(roma_output_final_CIS_treated$SampleMatrix)<- strapplyc(X = colnames(roma_output_final_CIS_treated$SampleMatrix),
                                                                 pattern = ".+\\_CIS.+\\_(.+h)", simplify = TRUE)

roma_output_final_CIS_control <- roma_output_CIS_3h
roma_output_final_CIS_control$ModuleMatrix <- roma_output_final_CIS_ModuleMatrix_mat
roma_output_final_CIS_control$SampleMatrix <- roma_output_final_CIS_SampleMatrix_mat %>% as.data.frame() %>% dplyr::select(contains("Control")) %>% as.matrix()
colnames(roma_output_final_CIS_control$SampleMatrix)<- strapplyc(X = colnames(roma_output_final_CIS_control$SampleMatrix),
                                                                 pattern = ".+\\_Control.+\\_(.+h)", simplify = TRUE)





dir_MTX <- "/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_MTX/"
for (i in seq(1,330, by=30)){
  from_to <- paste0(i,"_","to","_",i+29,"_","pathways_overdispersed_at_least_one_time_MTX.png")
  path_title <- paste0(dir_MTX,from_to)
  png(file = path_title, width=1500, height=700)
  v <- rownames(roma_output_final_OXA_treated$SampleMatrix)
  roma_output_final_MTX_treated$SampleMatrix[v,]
  pheatmap::pheatmap(roma_output_final_MTX_treated$SampleMatrix[rownames(roma_output_final_MTX_treated$SampleMatrix)[i:(i+29)],], cluster_cols = FALSE, cluster_rows=FALSE, main = paste(i,"to",i+29,"Overdispersed_at_least_one_time_point_MTX", sep="_"), color=colorRamps::blue2red(20))
  
  
  #Plot.Genesets.Samples(roma_output_final_MTX_treated, Selected = which(roma_output_final_MTX_treated$ModuleMatrix[, "ppv L1.x"] >= 0)[1:30], HMTite = paste(i+1,"to",i+30,"Overdispersed_at_least_one_time_point_MTX", sep="_"))
  dev.off()
  
}


dir_OXA <- "/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_OXA/"
for (i in seq(1,330, by=30)){
  from_to <- paste0(i,"_","to","_",i+29,"_","pathways_overdispersed_at_least_one_time_OXA.png")
  path_title <- paste0(dir_OXA,from_to)
  png(file = path_title, width=1500, height=700)
  pheatmap::pheatmap(roma_output_final_OXA_treated$SampleMatrix[rownames(roma_output_final_OXA_treated$SampleMatrix)[i:(i+29)],], cluster_cols = FALSE, cluster_rows=FALSE, main = paste(i,"to",i+29,"Overdispersed_at_least_one_time_point_OXA", sep="_"), color=colorRamps::blue2red(20))
  
  
  #Plot.Genesets.Samples(roma_output_final_OXA_treated, Selected = which(roma_output_final_OXA_treated$ModuleMatrix[, "ppv L1.x"] >= 0)[1:30], HMTite = paste(i+1,"to",i+30,"Overdispersed_at_least_one_time_point_OXA", sep="_"))
  dev.off()
  
}


dir_CIS <- "/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/pathways_overdispersed_at_least_one_time_CIS/"
for (i in seq(1,330, by=30)){
  from_to <- paste0(i,"_","to","_",i+29,"_","pathways_overdispersed_at_least_one_time_CIS.png")
  path_title <- paste0(dir_CIS,from_to)
  png(file = path_title, width=1500, height=700)
  pheatmap::pheatmap(roma_output_final_CIS_treated$SampleMatrix[rownames(roma_output_final_CIS_treated$SampleMatrix)[i:(i+29)],], cluster_cols = FALSE, cluster_rows=FALSE, main = paste(i,"to",i+29,"Overdispersed_at_least_one_time_point_CIS", sep="_"), color=colorRamps::blue2red(20))
  
  
  #Plot.Genesets.Samples(roma_output_final_CIS_treated, Selected = which(roma_output_final_CIS_treated$ModuleMatrix[, "ppv L1.x"] >= 0)[1:30], HMTite = paste(i+1,"to",i+30,"Overdispersed_at_least_one_time_point_CIS", sep="_"))
  dev.off()
  
}



png(paste0(dir_MTX,"331_to_332_pathways_overdispersed_at_least_one_time_MTX.png"), width=1500, height=700)
pheatmap::pheatmap(roma_output_final_MTX_treated$SampleMatrix[rownames(roma_output_final_MTX_treated$SampleMatrix)[331:332],], cluster_cols = FALSE, cluster_rows=FALSE, main = "331_to_332_Overdispersed_at_least_one_time_point_MTX", sep="_", color=colorRamps::blue2red(20))
dev.off()


png(paste0(dir_OXA,"331_to_332_pathways_overdispersed_at_least_one_time_OXA.png"), width=1500, height=700)
pheatmap::pheatmap(roma_output_final_OXA_treated$SampleMatrix[rownames(roma_output_final_OXA_treated$SampleMatrix)[331:332],], cluster_cols = FALSE, cluster_rows=FALSE, main = "331_to_332_Overdispersed_at_least_one_time_point_OXA", sep="_")
dev.off()


png(paste0(dir_CIS,"331_to_332_pathways_overdispersed_at_least_one_time_CIS.png"), width=1500, height=700)
pheatmap::pheatmap(roma_output_final_CIS_treated$SampleMatrix[rownames(roma_output_final_CIS_treated$SampleMatrix)[331:332],], cluster_cols = FALSE, cluster_rows=FALSE, main = "331_to_332_Overdispersed_at_least_one_time_point_CIS", sep="_", color=colorRamps::blue2red(20))
dev.off()




roma_output_final_MTX_OXA_CIS_treated <- cbind(roma_output_final_MTX_treated$SampleMatrix, roma_output_final_OXA_treated$SampleMatrix, roma_output_final_CIS_treated$SampleMatrix )


chimio <- c("MTX", "OXA", "CIS")
time <- c("3h", "6h", "12h", "18h", "24h")
annot_col <- data.frame(Time = rep(time, each=3, times=3), Treatment = rep(chimio, times = c(15,15,15)))
rownames(annot_col) <- colnames(roma_output_final_MTX_OXA_CIS_treated)


pheatmap::pheatmap(roma_output_final_MTX_OXA_CIS_treated[rownames(roma_output_final_MTX_OXA_CIS_treated)[1:10],], cluster_cols = FALSE, cluster_rows=FALSE, main = "331_to_332_Overdispersed_at_least_one_time_point_MTX", sep="_", color=colorRamps::blue2red(20), annotation_col = annot_col, show_colnames = FALSE, gp=gpar(fontsize=10))
              

#length(intersect_MTX_OXA_CIS_overdispersed_at_least_one)
#write.csv(intersect_MTX_OXA_CIS_overdispersed_at_least_one, "/data/users/tlassale/MODICED_link/in_vitro_results#/roma_results/Reactome_gmt_results/List_pathways_overdispersed_at_least_one_time_in_MTX_OXA_CIS_tumors/list_pathways_overdispersed_at_least_one_time_in_MTX_OXA_CIS_tumors.csv", quote = FALSE)





#time <- c("3h", "6h", "12h", "18h", "24h")
#annot_col <- data.frame(Time = rep(time, each=3, times=3), Treatment = rep("MTX", times = 15))
#rownames(annot_col) <- colnames(roma_output_final_MTX_OXA_CIS_treated)


#pheatmap::pheatmap(roma_output_final_MTX_treated$SampleMatrix[rownames(roma_output_final_MTX_treated$SampleMatrix)[1:30],], cluster_cols = FALSE, cluster_rows=FALSE, main = "331_to_332_Overdispersed_at_least_one_time_point_MTX", sep="_", color=colorRamps::blue2red(20), annotation_col = annot_col, show_colnames = FALSE, gp=gpar(fontsize=10))


list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS <- read.csv("/data/users/tlassale/MODICED_link/in_vitro_results/roma_results/Reactome_gmt_results/List_pathways_overdispersed_at_least_one_time_in_MTX_OXA_CIS_tumors/list_pathways_overdispersed_at_least_one_time_in_MTX_OXA_CIS_tumors.csv") %>% dplyr::rename("pathways" = "x") %>% dplyr::select(-1)
list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS <- list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS[grep("SIGNALING", list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS$pathways),]


pheatmap::pheatmap(roma_output_final_MTX_OXA_CIS_treated[rownames(roma_output_final_MTX_OXA_CIS_treated) %in% list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS,], cluster_cols = FALSE, cluster_rows=FALSE, main = "331_to_332_Overdispersed_at_least_one_time_point_MTX", sep="_", color=colorRamps::blue2red(20), annotation_col = annot_col, show_colnames = FALSE, gp=gpar(fontsize=10))



time <- c("3h", "6h", "12h", "18h", "24h")
annot_col <- data.frame(Time = rep(time, each=3, times=3), Treatment = rep("MTX", times = 15))
rownames(annot_col) <- colnames(roma_output_final_MTX_OXA_CIS_treated)



png(paste0(dir_MTX,"Signaling_pathways_overdispersed_at_least_one_time_point_MTX.png"), width=1500, height=700)
pheatmap::pheatmap(roma_output_final_MTX_OXA_CIS_treated[rownames(roma_output_final_MTX_OXA_CIS_treated) %in% list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS,1:15], cluster_cols = FALSE, cluster_rows=FALSE, main = "Signaling_pathways_overdispersed_at_least_one_time_point_MTX", sep="_", color=colorRamps::blue2red(20), annotation_col = annot_col, show_colnames = FALSE, gp=gpar(fontsize=10))
dev.off()

png(paste0(dir_OXA,"Signaling_pathways_overdispersed_at_least_one_time_point_OXA.png"), width=1500, height=700)
pheatmap::pheatmap(roma_output_final_MTX_OXA_CIS_treated[rownames(roma_output_final_MTX_OXA_CIS_treated) %in% list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS,16:30], cluster_cols = FALSE, cluster_rows=FALSE, main = "Signaling_pathways_overdispersed_at_least_one_time_point_OXA", sep="_", color=colorRamps::blue2red(20), annotation_col = annot_col, show_colnames = FALSE, gp=gpar(fontsize=10))
dev.off()

png(paste0(dir_CIS,"Signaling_pathways_overdispersed_at_least_one_time_point_CIS.png"), width=1500, height=700)
pheatmap::pheatmap(roma_output_final_MTX_OXA_CIS_treated[rownames(roma_output_final_MTX_OXA_CIS_treated) %in% list_pathways_overdispersed_at_least_one_time_MTX_OXA_CIS,31:45], cluster_cols = FALSE, cluster_rows=FALSE, main = "Signaling_pathways_overdispersed_at_least_one_time_point_CIS", sep="_", color=colorRamps::blue2red(20), annotation_col = annot_col, show_colnames = FALSE, gp=gpar(fontsize=10))
dev.off()
