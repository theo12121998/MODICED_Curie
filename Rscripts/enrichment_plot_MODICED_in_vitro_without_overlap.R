




library(patchwork)
library(tidyverse)
library(clusterProfiler)



enrich_plot_go_for_publi <- function(ordered_go_terms_overlap_filename, treatment, timepoint, dotplot_filename, top_table_filename, contrast_interest, setting, dotplot_gsea_filename, nmax){
  
  
  
  ordered_go_terms_overlap <- read.csv(ordered_go_terms_overlap_filename)
  
  ordered_go_terms_by_prop_overlap <- ordered_go_terms_overlap %>% arrange(desc(prop_overlap))
  
  if (!is.null(ordered_go_terms_overlap)){
    
    dotplot <- ggplot2::ggplot(data = ordered_go_terms_by_prop_overlap[1:min(nmax, nrow(ordered_go_terms_by_prop_overlap)),], 
                               mapping = aes(x = prop_overlap,
                                             y = reorder(Biological_Processes,prop_overlap),
                                             color = p.adjust,
                                             size = prop_overlap)) + 
      geom_point() + 
      scale_color_gradient(low = "red", high = "blue") + 
      theme_bw() +
      xlab("Proportion of Overlap (in %)") +
      ylab("GO (BP/MF/CC) ontology") +
      ggtitle(paste0(treatment," ","(",timepoint,")")) +
      theme(axis.text= element_text(face = "bold")) +
      theme(axis.title.x = element_text(size = 12, face="italic"))  +
      theme(axis.title.y = element_text(size = 12, face="italic"))  +
      theme(plot.title = element_text(size=12, face="bold", hjust = 0.5)) 
    
    
    
    
    if (nrow(ordered_go_terms_overlap) >= 3){
      
      top_table_contrast_interest <- read.csv(top_table_filename) %>% filter(contrast == contrast_interest) %>% mutate(rank = sign(logFC) * -log10(adj.P.Val)) %>% arrange(desc(rank))
      gene_list_interest <- top_table_contrast_interest$rank
      names(gene_list_interest) <- top_table_contrast_interest$genes
      gsea_contrast_interest <- gseGO(gene_list_interest, ont = "ALL", OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", pvalueCutoff = 1, minGSSize = 1)
      
      
      
      gse_plots <- lapply(ordered_go_terms_by_prop_overlap$Description[1:3], function(x){ 
        id <- which(gsea_contrast_interest$Description == x)
        gseaplot(gsea_contrast_interest, geneSetID = id, by="runningScore", title = x) + theme(plot.title = element_text(face = "bold", hjust=0.5, size=10)) + theme(axis.title.x = element_text(face = "")) + theme(axis.title.x = element_text(size = 11, face="italic")) + theme(axis.title.y = element_text(size=11, face="italic"))
      })
      
      
      
      wrap <- wrap_plots(list(dotplot , wrap_plots(gse_plots, ncol = 1)), widths = c(1.5,2)) + plot_annotation(title = paste("EnrichGO: ",treatment,"(",timepoint,")"," ","(",setting,")", sep=""), theme = theme(plot.title = element_text(face= "bold", hjust=0.5)))
      ggsave(dotplot_gsea_filename, wrap, units = "in", width = 13, height = 6.46, create.dir=TRUE)
      
      
    } else {
      ggsave(dotplot_filename, dotplot, width=8, height=7, dpi=400, create.dir=TRUE)
    }
    
  }
}



# In vitro


## MTX

# 24h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_MTX(24h) vs PBS(24h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "24h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/24h/dotplot_in_vitro_MTX_24h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_24h_vs_Control_24h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/24h/dotplot_gsea_in_vitro_MTX_24h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}


# 18h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_MTX(18h) vs PBS(18h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "18h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/18h/dotplot_in_vitro_MTX_18h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_18h_vs_Control_18h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/18h/dotplot_gsea_in_vitro_MTX_18h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 12h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_MTX(12h) vs PBS(12h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "12h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/12h/dotplot_in_vitro_MTX_12h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_12h_vs_Control_12h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/12h/dotplot_gsea_in_vitro_MTX_12h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 6h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_MTX(6h) vs PBS(6h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "6h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/6h/dotplot_in_vitro_MTX_6h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_6h_vs_Control_6h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/6h/dotplot_gsea_in_vitro_MTX_6h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 3h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_MTX(3h) vs PBS(3h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "3h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/3h/dotplot_in_vitro_MTX_3h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "MTX_3h_vs_Control_3h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/MTX/3h/dotplot_gsea_in_vitro_MTX_3h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}










## OXA

# 24h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_OXA(24h) vs PBS(24h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "24h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/24h/dotplot_in_vitro_OXA_24h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_24h_vs_Control_24h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/24h/dotplot_gsea_in_vitro_OXA_24h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}


# 18h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_OXA(18h) vs PBS(18h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "18h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/18h/dotplot_in_vitro_OXA_18h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_18h_vs_Control_18h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/18h/dotplot_gsea_in_vitro_OXA_18h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 12h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_OXA(12h) vs PBS(12h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "12h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/12h/dotplot_in_vitro_OXA_12h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_12h_vs_Control_12h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/12h/dotplot_gsea_in_vitro_OXA_12h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 6h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_OXA(6h) vs PBS(6h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "6h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/6h/dotplot_in_vitro_OXA_6h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_6h_vs_Control_6h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/6h/dotplot_gsea_in_vitro_OXA_6h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 3h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_OXA(3h) vs PBS(3h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "3h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/3h/dotplot_in_vitro_OXA_3h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "OXA_3h_vs_Control_3h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/OXA/3h/dotplot_gsea_in_vitro_OXA_3h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}






## CIS

# 24h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_CIS(24h) vs PBS(24h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "24h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/24h/dotplot_in_vitro_CIS_24h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_24h_vs_Control_24h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/24h/dotplot_gsea_in_vitro_CIS_24h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}


# 18h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_CIS(18h) vs PBS(18h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "18h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/18h/dotplot_in_vitro_CIS_18h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_18h_vs_Control_18h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/18h/dotplot_gsea_in_vitro_CIS_18h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 12h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_CIS(12h) vs PBS(12h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "12h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/12h/dotplot_in_vitro_CIS_12h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_12h_vs_Control_12h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/12h/dotplot_gsea_in_vitro_CIS_12h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 6h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_CIS(6h) vs PBS(6h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "6h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/6h/dotplot_in_vitro_CIS_6h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_6h_vs_Control_6h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/6h/dotplot_gsea_in_vitro_CIS_6h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 3h

filename <- "Bulk_RNA_in_vitro/in_vitro_results/final_enriched_go_terms_CIS(3h) vs PBS(3h)_all_diff_genes.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "3h",  "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/3h/dotplot_in_vitro_CIS_3h_all_diff_genes.pdf", "~/MODICED_link/Bulk_RNA_in_vitro/in_vitro_results/top_table_drug_control_in_vitro.csv", "CIS_3h_vs_Control_3h", "In vitro", "~/MODICED_link/figures_for_publication/In_vitro/all_diff_genes_ORA_GSEA/CIS/3h/dotplot_gsea_in_vitro_CIS_3h_all_diff_genes.pdf", nmax=20)
} else{
  
  print("coucou")
}


