library(tidyverse)
library(ComplexHeatmap)
library(scales)
library(edgeR)
library(rtracklayer)
library(patchwork)
library(FSA)
library(rstatix)
library(ggpubr)


top_table_in_vivo <- read.csv("Bulk_RNA_in_vitro/in_vitro_results/top_table_in_vivo_for_Andrea.csv")

counts_in_vivo <- read.csv("Bulk_RNA_in_vivo/preprocessing_in_vivo_RNA/results_pipeline/counts/tablecounts_raw.csv")




gtf <- import("/data/annotations/pipelines/Mouse/mm10/gtf/gencode.vM22.annotation.gtf")
gtf_annot <- gtf %>% as_tibble() %>% filter(type == "gene") %>%  mutate(gene_id = str_remove(gene_id, "\\..*"))


in.dir <- "Bulk_RNA_in_vivo/preprocessing_in_vivo_RNA/results_pipeline/"
annot_gene <- read_delim(file.path(in.dir,"counts/tableannot.csv"), delim = ",") %>%
  dplyr::select(gene_id, gene_name) %>%
  full_join(gtf_annot[, c("gene_id","gene_type")], by="gene_id")


sample_plan <- read.csv("Bulk_RNA_in_vivo/preprocessing_in_vivo_RNA/results_pipeline/samples/sample_plan.csv", header=FALSE)
colnames(sample_plan) <- c("sample ID", "sample Name", "Path_to_R1", "Path_to_R2")
rownames(sample_plan) <- NULL
sample_plan[grepl("Ct", sample_plan$`sample ID`) | grepl("CT", sample_plan$`sample ID`), "treatment"] <- "PBS"
sample_plan[1:3, "treatment"] <- "MTX"
sample_plan[4:6, "treatment"] <- "PBS"
sample_plan[!grepl("PBS", sample_plan$`sample ID`) & grepl("MTX", sample_plan$`sample ID`), "treatment"] <- "MTX"
sample_plan[grepl("Cis", sample_plan$`sample ID`),"treatment"] <- "Cis"
sample_plan[grepl("Oxa", sample_plan$`sample ID`),"treatment"] <- "Oxa"
sample_plan[grepl("D3", sample_plan$`sample ID`), "Time_point"] <- 3
sample_plan[grepl("D10", sample_plan$`sample ID`), "Time_point"] <- 10
sample_plan[grepl("D22", sample_plan$`sample ID`) | grepl("D24", sample_plan$`sample ID`) , "Time_point"] <- 22

sample_annot <- read.csv("Bulk_RNA_in_vivo/preprocessing_in_vivo_RNA/results_pipeline/samples/sample_annot.csv", check.names = FALSE)
sample_annot$`Sample ID` <- make.names(sample_annot$`Sample ID`)
sample_annot <- sample_annot %>% dplyr::rename("sample ID" = "Sample ID")
sample_plan$`sample ID`<- make.names(sample_plan$`sample ID`)
all_sample <- full_join(sample_plan, sample_annot, by = "sample ID")
all_sample$ID <- stringr::str_remove(stringr::str_remove_all(all_sample$`sample ID`, "_D3|_D10|_D22|.D22"), "MTX_PBS_|OCP_")

sample_plan <- tidyr::unite(all_sample, col = "new_name", Treatment, Day, Sheet, ID, remove = FALSE, sep = "_")




counts <- read_delim(file.path(in.dir,"counts/tablecounts_raw.csv"), delim = ",") %>%
  dplyr::rename("gene_id" = "...1") %>% 
  inner_join(., annot_gene, by = "gene_id") %>% 
  relocate(gene_name, .after = gene_id) %>%
  relocate(gene_type, .after = gene_name)

colnames(counts) <- make.names(colnames(counts))
counts <- counts[,c("gene_id","gene_name", "gene_type", sample_plan$`sample ID`)]
colnames(counts)[4:ncol(counts)] <- sample_plan$new_name




# **3. Differential Expression Analysis** 

# 1st dataset (counts_table + sample_plan) - MTX vs PBS
counts_MP <- counts %>% 
  dplyr::select(gene_id, gene_name, gene_type, contains("MTX") & !contains("CT1_D24")) %>% 
  as.data.frame() %>% 
  column_to_rownames("gene_id")

sample_plan_MP <- sample_plan %>% filter(new_name %in% colnames(counts_MP))
sample_plan_MP[,"Treatment"] <- paste(sample_plan_MP$treatment, sample_plan_MP$Day, sep ="_")


# 2nd dataset (counts_table + sample_plan) - Oxa/Cis vs PBS                          
counts_OCP <- counts %>% 
  dplyr::select(gene_id, gene_name, gene_type, contains("OCP") & !contains("CT1_D24")) %>% 
  as.data.frame() %>% 
  column_to_rownames("gene_id")

sample_plan_OCP <- sample_plan %>% filter(new_name %in% colnames(counts_OCP))
sample_plan_OCP[,"Treatment"] <- paste(sample_plan_OCP$treatment, sample_plan_OCP$Day, sep ="_")





y_raw_MP <- DGEList(counts = dplyr::select(counts_MP, -c("gene_name","gene_type")) , samples = sample_plan_MP, genes = counts_MP$gene_name)
design_MP <- model.matrix(~ 0 + Treatment, data=y_raw_MP$samples)
colnames(design_MP) <- sub("Treatment","", colnames(design_MP))
keep <- filterByExpr(y_raw_MP,design_MP)
y_MP <- y_raw_MP[keep,]
y_MP <- calcNormFactors(y_MP, method="TMM")
v_MP <- voomWithQualityWeights(y_MP, design_MP)
fit_MP <- lmFit(v_MP, design_MP)


y_raw_OCP <- DGEList(counts = dplyr::select(counts_OCP, -c("gene_name","gene_type")) , samples = sample_plan_OCP, genes = counts_OCP$gene_name)
design_OCP <- model.matrix(~ 0 + Treatment, data=y_raw_OCP$samples)
colnames(design_OCP) <- sub("Treatment","", colnames(design_OCP))
keep <- filterByExpr(y_raw_OCP,design_OCP)
y_OCP <- y_raw_OCP[keep,]
y_OCP <- calcNormFactors(y_OCP, method="TMM")
v_OCP <- voomWithQualityWeights(y_OCP, design_OCP)
fit_OCP <- lmFit(v_OCP, design_OCP)






counts_MTX_Control <- cbind(v_MP$E, v_MP$genes) %>% rownames_to_column("gene_id") %>% dplyr::relocate(genes, .after = gene_id)
rownames(counts_MTX_Control) <- make.names(counts_MTX_Control$genes, unique=TRUE)


counts_OCP_Control <- cbind(v_OCP$E, v_OCP$genes) %>% rownames_to_column("gene_id") %>% dplyr::relocate(genes, .after = gene_id) 
rownames(counts_OCP_Control) <- make.names(counts_OCP_Control$genes, unique=TRUE)




counts_MTX_Control_for_tregs <- counts_MTX_Control %>% filter(rownames(.) %in% c("Cd4", "Il2","Cd274","Pdcd1","Ctla4","Foxp3","Cd25","Il10", "Tgfb1", "Cd80","Cd86","Il2")) %>% as.data.frame() %>% dplyr::select(-c("gene_id","genes"))

counts_OCP_Control_for_tregs <- counts_OCP_Control %>% filter(rownames(.) %in% c("Cd4", "Il2","Cd274","Pdcd1","Ctla4","Foxp3","Cd25","Il10", "Tgfb1", "Cd80", "Cd86","Il2")) %>% as.data.frame()




counts_MTX_control_for_tregs_MTX_D3 <- rowMeans(counts_MTX_Control_for_tregs[,3:5]) %>% as.data.frame() %>% dplyr::rename("MTX_D3" = ".")

counts_MTX_control_for_tregs_PBS_D3 <- rowMeans(counts_MTX_Control_for_tregs[,6:8]) %>% as.data.frame() %>% dplyr::rename("PBS_MTX_D3" = ".")


counts_MTX_control_for_tregs_MTX_D10 <- rowMeans(counts_MTX_Control_for_tregs[,13:16]) %>% as.data.frame() %>% dplyr::rename("MTX_D10" = ".")

counts_MTX_control_for_tregs_PBS_D10 <- rowMeans(counts_MTX_Control_for_tregs[,9:12]) %>% as.data.frame() %>% dplyr::rename("PBS_MTX_D10" = ".")



counts_MTX_control_for_tregs_PBS_D22 <- rowMeans(counts_MTX_Control_for_tregs[,17:21]) %>% as.data.frame() %>% dplyr::rename("PBS_MTX_D22" = ".")

counts_MTX_control_for_tregs_MTX_D22 <- rowMeans(counts_MTX_Control_for_tregs[,22:26]) %>% as.data.frame() %>% dplyr::rename("MTX_D22" = ".")


counts_MTX_Control_final <- cbind(counts_MTX_control_for_tregs_MTX_D3, counts_MTX_control_for_tregs_PBS_D3, counts_MTX_control_for_tregs_MTX_D10, 
                            counts_MTX_control_for_tregs_PBS_D10, counts_MTX_control_for_tregs_MTX_D22, 
                            counts_MTX_control_for_tregs_PBS_D22)








counts_OCP_control_for_tregs_OXA_D3 <- rowMeans(counts_OCP_Control_for_tregs[,11:14]) %>% as.data.frame() %>% dplyr::rename("OXA_D3" = ".")

counts_OCP_control_for_tregs_PBS_D3 <- rowMeans(counts_OCP_Control_for_tregs[,7:10]) %>% as.data.frame() %>% dplyr::rename("PBS_OCP_D3" = ".")

counts_OCP_control_for_tregs_CIS_D3 <- rowMeans(counts_OCP_Control_for_tregs[,3:6]) %>% as.data.frame() %>% dplyr::rename("CIS_D3" = ".")





counts_OCP_control_for_tregs_OXA_D10 <- rowMeans(counts_OCP_Control_for_tregs[,23:26]) %>% as.data.frame() %>% dplyr::rename("OXA_D10" = ".")

counts_OCP_control_for_tregs_PBS_D10 <- rowMeans(counts_OCP_Control_for_tregs[,19:22]) %>% as.data.frame() %>% dplyr::rename("PBS_OCP_D10" = ".")

counts_OCP_control_for_tregs_CIS_D10 <- rowMeans(counts_OCP_Control_for_tregs[,15:18]) %>% as.data.frame() %>% dplyr::rename("CIS_D10" = ".")





counts_OCP_control_for_tregs_OXA_D22 <- rowMeans(counts_OCP_Control_for_tregs[,35:38]) %>% as.data.frame() %>% dplyr::rename("OXA_D22" = ".")

counts_OCP_control_for_tregs_PBS_D22 <- rowMeans(counts_OCP_Control_for_tregs[,31:34]) %>% as.data.frame() %>% dplyr::rename("PBS_OCP_D22" = ".")

counts_OCP_control_for_tregs_CIS_D22 <- rowMeans(counts_OCP_Control_for_tregs[,27:30]) %>% as.data.frame() %>% dplyr::rename("CIS_D22" = ".")






counts_OCP_Control_final <- cbind(counts_OCP_control_for_tregs_OXA_D3, counts_OCP_control_for_tregs_PBS_D3, counts_OCP_control_for_tregs_CIS_D3, 
                                  counts_OCP_control_for_tregs_OXA_D10, counts_OCP_control_for_tregs_PBS_D10, counts_OCP_control_for_tregs_CIS_D10
                                  ,counts_OCP_control_for_tregs_OXA_D22, counts_OCP_control_for_tregs_PBS_D22, counts_OCP_control_for_tregs_CIS_D22)






counts_MTX_control_final <- counts_MTX_Control_final %>% rownames_to_column("gene_name")
counts_OCP_Control_final <- counts_OCP_Control_final %>% rownames_to_column("gene_name")

counts_MTX_OCP_control_final_merged <- inner_join(counts_MTX_control_final, counts_OCP_Control_final, by="gene_name") %>% column_to_rownames("gene_name") %>% dplyr::select(!contains("PBS"))
counts_MTX_OCP_control_final_merged_scaled <- scale(as.matrix(counts_MTX_OCP_control_final_merged))



palette_length <- 100 

my_color <- viridis::viridis(palette_length)



annotation_columns_in_vivo <- data.frame(treatment = c("MTX","MTX","MTX","OXA","CIS","OXA","CIS","OXA","CIS"))
rownames(annotation_columns_in_vivo) <- colnames(counts_MTX_OCP_control_final_merged_scaled)


annotation_rows_in_vivo <- data.frame(signal = c(rep("checkpoint", 2), rep("imm_suppress",1), rep("Treg_marker",1), rep("imm_suppress",1), rep("checkpoint", 2),rep("checkpoint", 1), rep("Treg_marker", 1)))
rownames(annotation_rows_in_vivo) <- rownames(counts_MTX_OCP_control_final_merged_scaled)


annot_colors_vivo <- list(treatment=c(MTX="royalblue", OXA="firebrick", CIS = "forestgreen", PBS_MTX = "yellow", PBS_OCP = "grey"), signal = c(checkpoint = "red", imm_suppress= "blue", Treg_marker = "purple"))



ComplexHeatmap::pheatmap(as.matrix(counts_MTX_OCP_control_final_merged_scaled),color = my_color, cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 12, main="Treg marker+ Immunosuppressive cytokines + immune checkpoint expression", annotation_row =  annotation_rows_in_vivo, annotation_col = annotation_columns_in_vivo, annotation_colors = annot_colors_vivo, fontface = 3, show_colnames = TRUE, heatmap_legend_param = list(title = "scaled_mean_expression", title_gp = gpar(fontsize = 10, fontface="bold"), title_position= "topcenter",legend_height = grid::unit(35,"mm"), grid_width = grid::unit(5,"mm"), labels_gp = gpar(col = "black", font = 2), legend_direction = "horizontal", legend_width = unit(5,"cm")), cellwidth = 20) %>% ComplexHeatmap::draw(heatmap_legend_side = "bottom", annotation_legend_side = "right")










# Histogram for Andrea for the different markers


counts_MP_andrea <- cbind(y_MP$counts, y_MP$genes) %>% rownames_to_column("gene_id") %>% dplyr::relocate(genes, .after = gene_id) %>% dplyr::select(-"gene_id")
rownames(counts_MP_andrea) <- make.names(counts_MP_andrea$genes, unique=TRUE)
counts_MP_andrea_final <- counts_MP_andrea %>% select(-"genes")
counts_MP_andrea_cpm <- cpm(counts_MP_andrea_final, log=FALSE) %>% as.data.frame()


counts_OCP_andrea <- cbind(y_OCP$counts, y_OCP$genes) %>% rownames_to_column("gene_id") %>% dplyr::relocate(genes, .after = gene_id) %>% dplyr::select(-"gene_id")
rownames(counts_OCP_andrea) <- make.names(counts_OCP_andrea$genes, unique=TRUE)
counts_OCP_andrea_final <- counts_OCP_andrea %>% select(-"genes")
counts_OCP_andrea_cpm <- cpm(counts_OCP_andrea_final, log=FALSE) %>% as.data.frame()





counts_MTX_Control_for_tregs <- counts_MP_andrea_cpm %>% filter(rownames(.) %in% c("Cd4", "Il2","Cd274","Pdcd1","Ctla4","Foxp3","Cd25","Il10", "Tgfb1", "Cd80","Cd86","Il2"))

counts_OCP_Control_for_tregs <- counts_OCP_andrea_cpm %>% filter(rownames(.) %in% c("Cd4", "Il2","Cd274","Pdcd1","Ctla4","Foxp3","Cd25","Il10", "Tgfb1", "Cd80", "Cd86","Il2"))




counts_MTX_control_for_tregs_MTX_D3 <- rowMeans(counts_MTX_Control_for_tregs[,1:3]) %>% as.data.frame() %>% dplyr::rename("MTX_D3" = ".")

counts_MTX_control_for_tregs_PBS_D3 <- rowMeans(counts_MTX_Control_for_tregs[,4:6]) %>% as.data.frame() %>% dplyr::rename("PBS_MTX_D3" = ".")


counts_MTX_control_for_tregs_MTX_D10 <- rowMeans(counts_MTX_Control_for_tregs[,11:14]) %>% as.data.frame() %>% dplyr::rename("MTX_D10" = ".")

counts_MTX_control_for_tregs_PBS_D10 <- rowMeans(counts_MTX_Control_for_tregs[,7:10]) %>% as.data.frame() %>% dplyr::rename("PBS_MTX_D10" = ".")



counts_MTX_control_for_tregs_PBS_D22 <- rowMeans(counts_MTX_Control_for_tregs[,15:19]) %>% as.data.frame() %>% dplyr::rename("PBS_MTX_D22" = ".")

counts_MTX_control_for_tregs_MTX_D22 <- rowMeans(counts_MTX_Control_for_tregs[,20:24]) %>% as.data.frame() %>% dplyr::rename("MTX_D22" = ".")


counts_MTX_Control_final <- cbind(counts_MTX_control_for_tregs_MTX_D3, counts_MTX_control_for_tregs_PBS_D3, counts_MTX_control_for_tregs_MTX_D10, counts_MTX_control_for_tregs_PBS_D10, counts_MTX_control_for_tregs_MTX_D22, counts_MTX_control_for_tregs_PBS_D22)








counts_OCP_control_for_tregs_CIS_D3 <- rowMeans(counts_OCP_Control_for_tregs[,1:4]) %>% as.data.frame() %>% dplyr::rename("CIS_D3" = ".")

counts_OCP_control_for_tregs_PBS_D3 <- rowMeans(counts_OCP_Control_for_tregs[,5:8]) %>% as.data.frame() %>% dplyr::rename("PBS_OCP_D3" = ".")

counts_OCP_control_for_tregs_OXA_D3 <- rowMeans(counts_OCP_Control_for_tregs[,9:12]) %>% as.data.frame() %>% dplyr::rename("OXA_D3" = ".")





counts_OCP_control_for_tregs_CIS_D10 <- rowMeans(counts_OCP_Control_for_tregs[,13:16]) %>% as.data.frame() %>% dplyr::rename("CIS_D10" = ".")

counts_OCP_control_for_tregs_PBS_D10 <- rowMeans(counts_OCP_Control_for_tregs[,17:20]) %>% as.data.frame() %>% dplyr::rename("PBS_OCP_D10" = ".")

counts_OCP_control_for_tregs_OXA_D10 <- rowMeans(counts_OCP_Control_for_tregs[,21:24]) %>% as.data.frame() %>% dplyr::rename("OXA_D10" = ".")



counts_OCP_control_for_tregs_CIS_D22 <- rowMeans(counts_OCP_Control_for_tregs[,25:28]) %>% as.data.frame() %>% dplyr::rename("CIS_D22" = ".")

counts_OCP_control_for_tregs_PBS_D22 <- rowMeans(counts_OCP_Control_for_tregs[,29:32]) %>% as.data.frame() %>% dplyr::rename("PBS_OCP_D22" = ".")

counts_OCP_control_for_tregs_OXA_D22 <- rowMeans(counts_OCP_Control_for_tregs[,33:36]) %>% as.data.frame() %>% dplyr::rename("OXA_D22" = ".")



counts_OCP_Control_final <- cbind(counts_OCP_control_for_tregs_OXA_D3, counts_OCP_control_for_tregs_PBS_D3, counts_OCP_control_for_tregs_CIS_D3, counts_OCP_control_for_tregs_OXA_D10, counts_OCP_control_for_tregs_PBS_D10, counts_OCP_control_for_tregs_CIS_D10, counts_OCP_control_for_tregs_OXA_D22, counts_OCP_control_for_tregs_PBS_D22, counts_OCP_control_for_tregs_CIS_D22)



counts_MTX_control_final <- counts_MTX_Control_final %>% rownames_to_column("gene_name")
counts_OCP_Control_final <- counts_OCP_Control_final %>% rownames_to_column("gene_name")


counts_MTX_OCP_control_final_merged <- inner_join(counts_MTX_control_final, counts_OCP_Control_final, by="gene_name") %>% column_to_rownames("gene_name")
write.csv(counts_MTX_OCP_control_final_merged, "counts_in_vivo_andrea_tregs_immune_markers.csv")

counts_MTX_OCP_control_final_merged <- read.csv("counts_in_vivo_andrea_tregs_immune_markers.csv") %>% column_to_rownames("X")









# Ctla4, Pdcd1, Il10


signal_1 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(1) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[1])



signal_2 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(2) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[2])


signal_3 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(3) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[3])






signal_1_plot <- ggplot(signal_1, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[1]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))
                                                                                                                                                                                                                                                                                                                                      

signal_2_plot <- ggplot(signal_2, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[2]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))   



signal_3_plot <- ggplot(signal_3, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[3]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))   



signals_1_3 <- wrap_plots(signal_1_plot, signal_2_plot, signal_3_plot, ncol=1)








# Cd4, Tgfb1, Cd86


signal_4 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(4) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[4])



signal_5 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(5) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[5])


signal_6 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(6) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[6])






signal_4_plot <- ggplot(signal_4, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[4]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))


signal_5_plot <- ggplot(signal_5, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[5]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))   


signal_6_plot <- ggplot(signal_6, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[6]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))   



signals_4_6 <- wrap_plots(signal_4_plot, signal_5_plot, signal_6_plot, ncol=1)






# Cd80, Cd274, Foxp3


signal_7 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(7) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[7])



signal_8 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(8) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[8])


signal_9 <- counts_MTX_OCP_control_final_merged %>% dplyr::slice(9) %>% t() %>% as.data.frame() %>% dplyr::mutate(c(rep("D3",2), rep("D10",2), rep("D22",2), rep("D3",3), rep("D10",3), rep("D22",3))) %>% dplyr::rename("Day" = "c(...)") %>% dplyr::mutate(treatment = c("MTX","PBS_MTX","MTX","PBS_MTX","MTX","PBS_MTX","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS","OXA","PBS_OCP","CIS")) %>% dplyr::rename("mean_expression (cpm)" = rownames(counts_MTX_OCP_control_final_merged)[9])






signal_7_plot <- ggplot(signal_7, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[7]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))


signal_8_plot <- ggplot(signal_8, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[8]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))   


signal_9_plot <- ggplot(signal_9, mapping=aes(x = treatment, y = `mean_expression (cpm)`, fill = treatment)) + geom_col() + facet_wrap(~ factor(Day, levels = c("D3","D10","D22")), nrow = 1) + theme_bw() + scale_fill_manual(values = c('red','blue','green', 'darkgrey','yellow')) + labs(title = rownames(counts_MTX_OCP_control_final_merged)[9]) + theme(plot.title = element_text(face="bold", hjust = 0.5)) + theme(axis.title.x = element_text(face = "bold")) + theme(axis.title.y = element_text(face = "bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.text = element_text(size=10)) + theme(strip.text = element_text(size = 13, face="bold")) + theme(axis.text.x = element_text(face="bold", size=8))   



signals_7_9 <- wrap_plots(signal_7_plot, signal_8_plot, signal_9_plot, ncol=1)










# Boxplot + stats




# Ctla4

counts_MTX_Control_for_tregs_for_boxplot_Ctla4 <- counts_MTX_Control_for_tregs %>% dplyr::slice(1) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[1])

counts_OCP_Control_for_tregs_for_boxplot_Ctla4 <- counts_OCP_Control_for_tregs %>% dplyr::slice(1) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[1])    

counts_for_boxplot_Ctla4 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Ctla4, counts_OCP_Control_for_tregs_for_boxplot_Ctla4)



# D3

counts_for_boxplot_Ctla4_D3 <- counts_for_boxplot_Ctla4 %>% filter(Day == "D3")
kruskal_test_D3_Ctla4 <- counts_for_boxplot_Ctla4_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Ctla4 <- counts_for_boxplot_Ctla4_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")



counts_for_boxplot_Ctla4_D3$Treatment <- factor(counts_for_boxplot_Ctla4_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Ctla4_D3_boxplot_Ctla4 <- ggboxplot(counts_for_boxplot_Ctla4_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Ctla4, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Ctla4, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Ctla4_D3_boxplot_Ctla4 <- ggpar(Ctla4_D3_boxplot_Ctla4, legend = "none") 


# D10

counts_for_boxplot_Ctla4_D10 <- counts_for_boxplot_Ctla4 %>% filter(Day == "D10")
kruskal_test_D10_Ctla4 <- counts_for_boxplot_Ctla4_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Ctla4 <- counts_for_boxplot_Ctla4_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Ctla4_D10$Treatment <- factor(counts_for_boxplot_Ctla4_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Ctla4_D10_boxplot_Ctla4 <- ggboxplot(counts_for_boxplot_Ctla4_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Ctla4, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Ctla4, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Ctla4_D10_boxplot_Ctla4 <- ggpar(Ctla4_D10_boxplot_Ctla4, legend = "none") 



# D22

counts_for_boxplot_Ctla4_D22 <- counts_for_boxplot_Ctla4 %>% filter(Day == "D22")
kruskal_test_D22_Ctla4 <- counts_for_boxplot_Ctla4_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Ctla4 <- counts_for_boxplot_Ctla4_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Ctla4_D22$Treatment <- factor(counts_for_boxplot_Ctla4_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Ctla4_D22_boxplot_Ctla4 <- ggboxplot(counts_for_boxplot_Ctla4_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Ctla4, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Ctla4, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Ctla4_all_timepoints <- patchwork::wrap_plots(Ctla4_D3_boxplot_Ctla4, Ctla4_D10_boxplot_Ctla4, Ctla4_D22_boxplot_Ctla4, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[1]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))

bn <- ggpubr::ggarrange(Ctla4_D3_boxplot_Ctla4, Ctla4_D10_boxplot_Ctla4, Ctla4_D22_boxplot_Ctla4, nrow=1, labels = "Ctla4")






# Pdcd1

counts_MTX_Control_for_tregs_for_boxplot_Pdcd1 <- counts_MTX_Control_for_tregs %>% dplyr::slice(2) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[2])

counts_OCP_Control_for_tregs_for_boxplot_Pdcd1 <- counts_OCP_Control_for_tregs %>% dplyr::slice(2) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[2])    

counts_for_boxplot_Pdcd1 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Pdcd1, counts_OCP_Control_for_tregs_for_boxplot_Pdcd1)



# D3

counts_for_boxplot_Pdcd1_D3 <- counts_for_boxplot_Pdcd1 %>% filter(Day == "D3")
kruskal_test_D3_Pdcd1 <- counts_for_boxplot_Pdcd1_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Pdcd1 <- counts_for_boxplot_Pdcd1_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Pdcd1_D3$Treatment <- factor(counts_for_boxplot_Pdcd1_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Pdcd1_D3_boxplot_Pdcd1 <- ggboxplot(counts_for_boxplot_Pdcd1_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Pdcd1, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Pdcd1, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Pdcd1_D3_boxplot_Pdcd1 <- ggpar(Pdcd1_D3_boxplot_Pdcd1, legend = "none") 


# D10

counts_for_boxplot_Pdcd1_D10 <- counts_for_boxplot_Pdcd1 %>% filter(Day == "D10")
kruskal_test_D10_Pdcd1 <- counts_for_boxplot_Pdcd1_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Pdcd1 <- counts_for_boxplot_Pdcd1_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Pdcd1_D10$Treatment <- factor(counts_for_boxplot_Pdcd1_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Pdcd1_D10_boxplot_Pdcd1 <- ggboxplot(counts_for_boxplot_Pdcd1_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Pdcd1, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Pdcd1, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Pdcd1_D10_boxplot_Pdcd1 <- ggpar(Pdcd1_D10_boxplot_Pdcd1, legend = "none") 



# D22

counts_for_boxplot_Pdcd1_D22 <- counts_for_boxplot_Pdcd1 %>% filter(Day == "D22")
kruskal_test_D22_Pdcd1 <- counts_for_boxplot_Pdcd1_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Pdcd1 <- counts_for_boxplot_Pdcd1_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Pdcd1_D22$Treatment <- factor(counts_for_boxplot_Pdcd1_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Pdcd1_D22_boxplot_Pdcd1 <- ggboxplot(counts_for_boxplot_Pdcd1_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Pdcd1, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Pdcd1, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Pdcd1_all_timepoints <- patchwork::wrap_plots(Pdcd1_D3_boxplot_Pdcd1, Pdcd1_D10_boxplot_Pdcd1, Pdcd1_D22_boxplot_Pdcd1, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[2]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))




bn_2 <- ggpubr::ggarrange(Pdcd1_D3_boxplot_Pdcd1, Pdcd1_D10_boxplot_Pdcd1, Pdcd1_D22_boxplot_Pdcd1, nrow=1, labels = "Pdcd1")






# Il10

counts_MTX_Control_for_tregs_for_boxplot_Il10 <- counts_MTX_Control_for_tregs %>% dplyr::slice(3) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[3])

counts_OCP_Control_for_tregs_for_boxplot_Il10 <- counts_OCP_Control_for_tregs %>% dplyr::slice(3) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[3])    

counts_for_boxplot_Il10 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Il10, counts_OCP_Control_for_tregs_for_boxplot_Il10)





# D3

counts_for_boxplot_Il10_D3 <- counts_for_boxplot_Il10 %>% filter(Day == "D3")
kruskal_test_D3_Il10 <- counts_for_boxplot_Il10_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Il10 <- counts_for_boxplot_Il10_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Il10_D3$Treatment <- factor(counts_for_boxplot_Il10_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Il10_D3_boxplot_Il10 <- ggboxplot(counts_for_boxplot_Il10_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Il10, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Il10, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Il10_D3_boxplot_Il10 <- ggpar(Il10_D3_boxplot_Il10, legend = "none") 



# D10

counts_for_boxplot_Il10_D10 <- counts_for_boxplot_Il10 %>% filter(Day == "D10")
kruskal_test_D10_Il10 <- counts_for_boxplot_Il10_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Il10 <- counts_for_boxplot_Il10_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Il10_D10$Treatment <- factor(counts_for_boxplot_Il10_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Il10_D10_boxplot_Il10 <- ggboxplot(counts_for_boxplot_Il10_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Il10, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Il10, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Il10_D10_boxplot_Il10 <- ggpar(Il10_D10_boxplot_Il10, legend = "none") 



# D22

counts_for_boxplot_Il10_D22 <- counts_for_boxplot_Il10 %>% filter(Day == "D22")
kruskal_test_D22_Il10 <- counts_for_boxplot_Il10_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Il10 <- counts_for_boxplot_Il10_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Il10_D22$Treatment <- factor(counts_for_boxplot_Il10_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Il10_D22_boxplot_Il10 <- ggboxplot(counts_for_boxplot_Il10_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Il10, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Il10, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))

 




Il10_all_timepoints <- patchwork::wrap_plots(Il10_D3_boxplot_Il10, Il10_D10_boxplot_Il10, Il10_D22_boxplot_Il10, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[3]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))





bn_3 <- ggpubr::ggarrange(Il10_D3_boxplot_Il10, Il10_D10_boxplot_Il10, Il10_D22_boxplot_Il10, nrow=1, labels = "Il10")


p <- ggarrange(bn,bn_2,bn_3, ncol = 1, nrow=3)



patchwork::wrap_plots(Ctla4_all_timepoints, Pdcd1_all_timepoints, Il10_all_timepoints, ncol = 1)














# Cd4, Tgfb1, Cd86




# Cd4

counts_MTX_Control_for_tregs_for_boxplot_Cd4 <- counts_MTX_Control_for_tregs %>% dplyr::slice(4) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[4])

counts_OCP_Control_for_tregs_for_boxplot_Cd4 <- counts_OCP_Control_for_tregs %>% dplyr::slice(4) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[4])    

counts_for_boxplot_Cd4 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Cd4, counts_OCP_Control_for_tregs_for_boxplot_Cd4)



# D3

counts_for_boxplot_Cd4_D3 <- counts_for_boxplot_Cd4 %>% filter(Day == "D3")
kruskal_test_D3_Cd4 <- counts_for_boxplot_Cd4_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Cd4 <- counts_for_boxplot_Cd4_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")



counts_for_boxplot_Cd4_D3$Treatment <- factor(counts_for_boxplot_Cd4_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd4_D3_boxplot_Cd4 <- ggboxplot(counts_for_boxplot_Cd4_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Cd4, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Cd4, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd4_D3_boxplot_Cd4 <- ggpar(Cd4_D3_boxplot_Cd4, legend = "none") 


# D10

counts_for_boxplot_Cd4_D10 <- counts_for_boxplot_Cd4 %>% filter(Day == "D10")
kruskal_test_D10_Cd4 <- counts_for_boxplot_Cd4_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Cd4 <- counts_for_boxplot_Cd4_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Cd4_D10$Treatment <- factor(counts_for_boxplot_Cd4_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd4_D10_boxplot_Cd4 <- ggboxplot(counts_for_boxplot_Cd4_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Cd4, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Cd4, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd4_D10_boxplot_Cd4 <- ggpar(Cd4_D10_boxplot_Cd4, legend = "none") 



# D22

counts_for_boxplot_Cd4_D22 <- counts_for_boxplot_Cd4 %>% filter(Day == "D22")
kruskal_test_D22_Cd4 <- counts_for_boxplot_Cd4_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Cd4 <- counts_for_boxplot_Cd4_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Cd4_D22$Treatment <- factor(counts_for_boxplot_Cd4_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd4_D22_boxplot_Cd4 <- ggboxplot(counts_for_boxplot_Cd4_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Cd4, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Cd4, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Cd4_all_timepoints <- patchwork::wrap_plots(Cd4_D3_boxplot_Cd4, Cd4_D10_boxplot_Cd4, Cd4_D22_boxplot_Cd4, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[4]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))

cd4_arranged <- ggpubr::ggarrange(Cd4_D3_boxplot_Cd4, Cd4_D10_boxplot_Cd4, Cd4_D22_boxplot_Cd4, nrow=1, labels = "Cd4")






# Tgfb1

counts_MTX_Control_for_tregs_for_boxplot_Tgfb1 <- counts_MTX_Control_for_tregs %>% dplyr::slice(5) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[5])

counts_OCP_Control_for_tregs_for_boxplot_Tgfb1 <- counts_OCP_Control_for_tregs %>% dplyr::slice(5) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[5])    

counts_for_boxplot_Tgfb1 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Tgfb1, counts_OCP_Control_for_tregs_for_boxplot_Tgfb1)



# D3

counts_for_boxplot_Tgfb1_D3 <- counts_for_boxplot_Tgfb1 %>% filter(Day == "D3")
kruskal_test_D3_Tgfb1 <- counts_for_boxplot_Tgfb1_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Tgfb1 <- counts_for_boxplot_Tgfb1_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Tgfb1_D3$Treatment <- factor(counts_for_boxplot_Tgfb1_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Tgfb1_D3_boxplot_Tgfb1 <- ggboxplot(counts_for_boxplot_Tgfb1_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Tgfb1, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Tgfb1, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Tgfb1_D3_boxplot_Tgfb1 <- ggpar(Tgfb1_D3_boxplot_Tgfb1, legend = "none") 



# D10

counts_for_boxplot_Tgfb1_D10 <- counts_for_boxplot_Tgfb1 %>% filter(Day == "D10")
kruskal_test_D10_Tgfb1 <- counts_for_boxplot_Tgfb1_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Tgfb1 <- counts_for_boxplot_Tgfb1_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Tgfb1_D10$Treatment <- factor(counts_for_boxplot_Tgfb1_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Tgfb1_D10_boxplot_Tgfb1 <- ggboxplot(counts_for_boxplot_Tgfb1_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Tgfb1, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Tgfb1, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Tgfb1_D10_boxplot_Tgfb1 <- ggpar(Tgfb1_D10_boxplot_Tgfb1, legend = "none") 



# D22

counts_for_boxplot_Tgfb1_D22 <- counts_for_boxplot_Tgfb1 %>% filter(Day == "D22")
kruskal_test_D22_Tgfb1 <- counts_for_boxplot_Tgfb1_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Tgfb1 <- counts_for_boxplot_Tgfb1_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Tgfb1_D22$Treatment <- factor(counts_for_boxplot_Tgfb1_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Tgfb1_D22_boxplot_Tgfb1 <- ggboxplot(counts_for_boxplot_Tgfb1_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Tgfb1, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Tgfb1, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Tgfb1_all_timepoints <- patchwork::wrap_plots(Tgfb1_D3_boxplot_Tgfb1, Tgfb1_D10_boxplot_Tgfb1, Tgfb1_D22_boxplot_Tgfb1, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[5]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))




tgfb1_arranged <- ggpubr::ggarrange(Tgfb1_D3_boxplot_Tgfb1, Tgfb1_D10_boxplot_Tgfb1, Tgfb1_D22_boxplot_Tgfb1, nrow=1, labels = "Tgfb1")






# Cd86


counts_OCP_Control_for_tregs_for_boxplot_Cd86 <- counts_OCP_Control_for_tregs %>% dplyr::slice(6) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[6])    

counts_for_boxplot_Cd86 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Cd86, counts_OCP_Control_for_tregs_for_boxplot_Cd86)





# D3

counts_for_boxplot_Cd86_D3 <- counts_for_boxplot_Cd86 %>% filter(Day == "D3")
kruskal_test_D3_Cd86 <- counts_for_boxplot_Cd86_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Cd86 <- counts_for_boxplot_Cd86_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Cd86_D3$Treatment <- factor(counts_for_boxplot_Cd86_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd86_D3_boxplot_Cd86 <- ggboxplot(counts_for_boxplot_Cd86_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Cd86, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Cd86, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd86_D3_boxplot_Cd86 <- ggpar(Cd86_D3_boxplot_Cd86, legend = "none") 



# D10

counts_for_boxplot_Cd86_D10 <- counts_for_boxplot_Cd86 %>% filter(Day == "D10")
kruskal_test_D10_Cd86 <- counts_for_boxplot_Cd86_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Cd86 <- counts_for_boxplot_Cd86_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Cd86_D10$Treatment <- factor(counts_for_boxplot_Cd86_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd86_D10_boxplot_Cd86 <- ggboxplot(counts_for_boxplot_Cd86_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Cd86, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Cd86, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd86_D10_boxplot_Cd86 <- ggpar(Cd86_D10_boxplot_Cd86, legend = "none") 



# D22

counts_for_boxplot_Cd86_D22 <- counts_for_boxplot_Cd86 %>% filter(Day == "D22")
kruskal_test_D22_Cd86 <- counts_for_boxplot_Cd86_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Cd86 <- counts_for_boxplot_Cd86_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Cd86_D22$Treatment <- factor(counts_for_boxplot_Cd86_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd86_D22_boxplot_Cd86 <- ggboxplot(counts_for_boxplot_Cd86_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Cd86, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Cd86, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Cd86_all_timepoints <- patchwork::wrap_plots(Cd86_D3_boxplot_Cd86, Cd86_D10_boxplot_Cd86, Cd86_D22_boxplot_Cd86, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[6]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))





cd86_arranged <- ggpubr::ggarrange(Cd86_D3_boxplot_Cd86, Cd86_D10_boxplot_Cd86, Cd86_D22_boxplot_Cd86, nrow=1, labels = "Cd86")


cd4_tgfb1_cd86_arranged <- ggarrange(cd4_arranged, tgfb1_arranged, cd86_arranged, ncol = 1, nrow=3)







# Cd80, Cd274, Foxp3






# Cd80

counts_MTX_Control_for_tregs_for_boxplot_Cd80 <- counts_MTX_Control_for_tregs %>% dplyr::slice(7) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[7])

counts_OCP_Control_for_tregs_for_boxplot_Cd80 <- counts_OCP_Control_for_tregs %>% dplyr::slice(7) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[7])    

counts_for_boxplot_Cd80 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Cd80, counts_OCP_Control_for_tregs_for_boxplot_Cd80)



# D3

counts_for_boxplot_Cd80_D3 <- counts_for_boxplot_Cd80 %>% filter(Day == "D3")
kruskal_test_D3_Cd80 <- counts_for_boxplot_Cd80_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Cd80 <- counts_for_boxplot_Cd80_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")



counts_for_boxplot_Cd80_D3$Treatment <- factor(counts_for_boxplot_Cd80_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd80_D3_boxplot_Cd80 <- ggboxplot(counts_for_boxplot_Cd80_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Cd80, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Cd80, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd80_D3_boxplot_Cd80 <- ggpar(Cd80_D3_boxplot_Cd80, legend = "none") 


# D10

counts_for_boxplot_Cd80_D10 <- counts_for_boxplot_Cd80 %>% filter(Day == "D10")
kruskal_test_D10_Cd80 <- counts_for_boxplot_Cd80_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Cd80 <- counts_for_boxplot_Cd80_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Cd80_D10$Treatment <- factor(counts_for_boxplot_Cd80_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd80_D10_boxplot_Cd80 <- ggboxplot(counts_for_boxplot_Cd80_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Cd80, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Cd80, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd80_D10_boxplot_Cd80 <- ggpar(Cd80_D10_boxplot_Cd80, legend = "none") 



# D22

counts_for_boxplot_Cd80_D22 <- counts_for_boxplot_Cd80 %>% filter(Day == "D22")
kruskal_test_D22_Cd80 <- counts_for_boxplot_Cd80_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Cd80 <- counts_for_boxplot_Cd80_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Cd80_D22$Treatment <- factor(counts_for_boxplot_Cd80_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd80_D22_boxplot_Cd80 <- ggboxplot(counts_for_boxplot_Cd80_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Cd80, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Cd80, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Cd80_all_timepoints <- patchwork::wrap_plots(Cd80_D3_boxplot_Cd80, Cd80_D10_boxplot_Cd80, Cd80_D22_boxplot_Cd80, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))

Cd80_arranged <- ggpubr::ggarrange(Cd80_D3_boxplot_Cd80, Cd80_D10_boxplot_Cd80, Cd80_D22_boxplot_Cd80, nrow=1, labels = "Cd80")






# Cd274

counts_MTX_Control_for_tregs_for_boxplot_Cd274 <- counts_MTX_Control_for_tregs %>% dplyr::slice(8) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[8])

counts_OCP_Control_for_tregs_for_boxplot_Cd274 <- counts_OCP_Control_for_tregs %>% dplyr::slice(8) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[8])    

counts_for_boxplot_Cd274 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Cd274, counts_OCP_Control_for_tregs_for_boxplot_Cd274)



# D3

counts_for_boxplot_Cd274_D3 <- counts_for_boxplot_Cd274 %>% filter(Day == "D3")
kruskal_test_D3_Cd274 <- counts_for_boxplot_Cd274_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Cd274 <- counts_for_boxplot_Cd274_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")

counts_for_boxplot_Cd274_D3$Treatment <- factor(counts_for_boxplot_Cd274_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd274_D3_boxplot_Cd274 <- ggboxplot(counts_for_boxplot_Cd274_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Cd274, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Cd274, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd274_D3_boxplot_Cd274 <- ggpar(Cd274_D3_boxplot_Cd274, legend = "none") 



# D10

counts_for_boxplot_Cd274_D10 <- counts_for_boxplot_Cd274 %>% filter(Day == "D10")
kruskal_test_D10_Cd274 <- counts_for_boxplot_Cd274_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Cd274 <- counts_for_boxplot_Cd274_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Cd274_D10$Treatment <- factor(counts_for_boxplot_Cd274_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd274_D10_boxplot_Cd274 <- ggboxplot(counts_for_boxplot_Cd274_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Cd274, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Cd274, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Cd274_D10_boxplot_Cd274 <- ggpar(Cd274_D10_boxplot_Cd274, legend = "none") 



# D22

counts_for_boxplot_Cd274_D22 <- counts_for_boxplot_Cd274 %>% filter(Day == "D22")
kruskal_test_D22_Cd274 <- counts_for_boxplot_Cd274_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Cd274 <- counts_for_boxplot_Cd274_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Cd274_D22$Treatment <- factor(counts_for_boxplot_Cd274_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Cd274_D22_boxplot_Cd274 <- ggboxplot(counts_for_boxplot_Cd274_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Cd274, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Cd274, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Cd274_all_timepoints <- patchwork::wrap_plots(Cd274_D3_boxplot_Cd274, Cd274_D10_boxplot_Cd274, Cd274_D22_boxplot_Cd274, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[8]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))




Cd274_arranged <- ggpubr::ggarrange(Cd274_D3_boxplot_Cd274, Cd274_D10_boxplot_Cd274, Cd274_D22_boxplot_Cd274, nrow=1, labels = "Cd274")








# Foxp3

counts_MTX_Control_for_tregs_for_boxplot_Foxp3 <- counts_MTX_Control_for_tregs %>% dplyr::slice(9) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",6), rep("D10",8), rep("D22",10))) %>% dplyr::mutate(Treatment = c(rep("MTX",3), rep("PBS",7), rep("MTX",4), rep("PBS",5), rep("MTX",5))) %>% dplyr::rename("Expression" = rownames(counts_MTX_Control_for_tregs)[9])

counts_OCP_Control_for_tregs_for_boxplot_Foxp3 <- counts_OCP_Control_for_tregs %>% dplyr::slice(9) %>% t() %>% as.data.frame() %>% dplyr::mutate(Day = c(rep("D3",12), rep("D10",12), rep("D22",12))) %>% dplyr::mutate(Treatment = c(rep("CIS",4), rep("PBS",4), rep("OXA",4), rep("CIS",4), rep("PBS",4),rep("OXA",4), rep("CIS",4), rep("PBS",4), rep("OXA",4))) %>% dplyr::rename("Expression" = rownames(counts_OCP_Control_for_tregs)[9])    

counts_for_boxplot_Foxp3 <- rbind(counts_MTX_Control_for_tregs_for_boxplot_Foxp3, counts_OCP_Control_for_tregs_for_boxplot_Foxp3)





# D3

counts_for_boxplot_Foxp3_D3 <- counts_for_boxplot_Foxp3 %>% filter(Day == "D3")
kruskal_test_D3_Foxp3 <- counts_for_boxplot_Foxp3_D3 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D3_Foxp3 <- counts_for_boxplot_Foxp3_D3 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Foxp3_D3$Treatment <- factor(counts_for_boxplot_Foxp3_D3$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Foxp3_D3_boxplot_Foxp3 <- ggboxplot(counts_for_boxplot_Foxp3_D3, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D3_Foxp3, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D3_Foxp3, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D3") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Foxp3_D3_boxplot_Foxp3 <- ggpar(Foxp3_D3_boxplot_Foxp3, legend = "none") 



# D10

counts_for_boxplot_Foxp3_D10 <- counts_for_boxplot_Foxp3 %>% filter(Day == "D10")
kruskal_test_D10_Foxp3 <- counts_for_boxplot_Foxp3_D10 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D10_Foxp3 <- counts_for_boxplot_Foxp3_D10 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Foxp3_D10$Treatment <- factor(counts_for_boxplot_Foxp3_D10$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Foxp3_D10_boxplot_Foxp3 <- ggboxplot(counts_for_boxplot_Foxp3_D10, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D10_Foxp3, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D10_Foxp3, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D10") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))


Foxp3_D10_boxplot_Foxp3 <- ggpar(Foxp3_D10_boxplot_Foxp3, legend = "none") 



# D22

counts_for_boxplot_Foxp3_D22 <- counts_for_boxplot_Foxp3 %>% filter(Day == "D22")
kruskal_test_D22_Foxp3 <- counts_for_boxplot_Foxp3_D22 %>% kruskal_test(Expression ~ Treatment)
dunn_test_D22_Foxp3 <- counts_for_boxplot_Foxp3_D22 %>% dunn_test(Expression ~ Treatment, p.adjust.method = "bonferroni") %>% add_xy_position(x = "Treatment")


counts_for_boxplot_Foxp3_D22$Treatment <- factor(counts_for_boxplot_Foxp3_D22$Treatment, levels = c("PBS","CIS","OXA","MTX"))
Foxp3_D22_boxplot_Foxp3 <- ggboxplot(counts_for_boxplot_Foxp3_D22, x = "Treatment", y = "Expression", fill="Treatment", palette = c(MTX = 'royalblue', PBS = 'darkgrey', CIS = 'forestgreen', OXA = '#DC143C')) +
  stat_pvalue_manual(dunn_test_D22_Foxp3, hide.ns = TRUE) +
  labs(subtitle = get_test_label(kruskal_test_D22_Foxp3, detailed = TRUE)) + theme_pubr(legend = "right") + ggtitle("D22") + ylab("Expression (cpm)") + theme(plot.title = element_text(hjust= 0.5)) +  theme(axis.title = element_text(size=10)) + theme(plot.title = element_text(face="bold"))






Foxp3_all_timepoints <- patchwork::wrap_plots(Foxp3_D3_boxplot_Foxp3, Foxp3_D10_boxplot_Foxp3, Foxp3_D22_boxplot_Foxp3, nrow=1) & patchwork::plot_annotation(title = rownames(counts_MTX_Control_for_tregs)[9]) & theme(plot.title = element_text(face = "bold", size=13, hjust = 0.5))





Foxp3_arranged <- ggpubr::ggarrange(Foxp3_D3_boxplot_Foxp3, Foxp3_D10_boxplot_Foxp3, Foxp3_D22_boxplot_Foxp3, nrow=1, labels = "Foxp3")


Cd80_Cd274_Foxp3_arranged <- ggarrange(Cd80_arranged, Cd274_arranged, Foxp3_arranged, ncol = 1, nrow=3)
















