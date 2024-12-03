




library(patchwork)
library(tidyverse)
library(clusterProfiler)



enrich_plot_go_for_publi <- function(ordered_go_terms_overlap_filename, treatment, timepoint_1, timepoint_2, dotplot_overlap_filename, top_table_filename, contrast_interest, setting, dotplot_gsea_filename, nmax){
  
  
  
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
      ggtitle(paste0(treatment," ","(",timepoint_1,")")) +
      theme(axis.text= element_text(face = "bold")) +
      theme(axis.title.x = element_text(size = 12, face="italic"))  +
      theme(axis.title.y = element_text(size = 12, face="italic"))  +
      theme(plot.title = element_text(size=12, face="bold", hjust = 0.5)) 
    
    
    
    
    if (nrow(ordered_go_terms_overlap) >= 3){
      
      top_table_contrast_interest <- read.csv(top_table_filename) %>% filter(contrast == contrast_interest) %>% mutate(rank = sign(logFC) * -log10(adj.P.Val)) %>% arrange(desc(rank))
      gene_list_interest <- top_table_contrast_interest$rank
      names(gene_list_interest) <- top_table_contrast_interest$genes
      gsea_contrast_interest <- gseGO(gene_list_interest, ont = "ALL", OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", eps=0, pvalueCutoff = 1, minGSSize = 1)
      
      
      
      gse_plots <- lapply(ordered_go_terms_by_prop_overlap$Description[1:3], function(x){ 
        id <- which(gsea_contrast_interest$Description == x)
        gseaplot(gsea_contrast_interest, geneSetID = id, by="runningScore", title = x) + theme(plot.title = element_text(face = "bold", hjust=0.5, size=10)) + theme(axis.title.x = element_text(face = "")) + theme(axis.title.x = element_text(size = 11, face="italic")) + theme(axis.title.y = element_text(size=11, face="italic"))
      })
      
      
      
      wrap <- wrap_plots(list(dotplot , wrap_plots(gse_plots, ncol = 1)), widths = c(1.5,2)) + plot_annotation(title = paste("EnrichGO: ",treatment,"(",timepoint_1,")", "/", treatment,"(",timepoint_2,")"," ", "-"," ","overlap"," ","(",setting,")", sep=""), theme = theme(plot.title = element_text(face= "bold", hjust=0.5)))
      ggsave(dotplot_gsea_filename, wrap, units = "in", width = 13, height = 6.46, create.dir = TRUE)
      
      
    } else {
      ggsave(dotplot_overlap_filename, dotplot, width=8, height=7, dpi=400, create.dir = TRUE)
    }
    
  }
}






# In vivo

## MTX

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/24h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_24h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D3", "24h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/24h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_24h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/24h vs D3/dotplot_gsea_in_vivo_MTX_24h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/24h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_24h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D10", "24h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/24h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_24h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/24h vs D10/dotplot_gsea_in_vivo_MTX_24h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 18h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/18h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_18h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D3", "18h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/18h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_18h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/18h vs D3/dotplot_gsea_in_vivo_MTX_18h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/18h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_18h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D10", "18h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/18h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_18h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/18h vs D10/dotplot_gsea_in_vivo_MTX_18h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}




# 12h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/12h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_12h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D3", "12h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/12h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_12h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/12h vs D3/dotplot_gsea_in_vivo_MTX_12h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/12h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_12h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D10", "12h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/12h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_12h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/12h vs D10/dotplot_gsea_in_vivo_MTX_12h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 6h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/6h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_6h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D3", "6h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/6h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_6h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/6h vs D3/dotplot_gsea_in_vivo_MTX_6h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/6h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_6h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D10", "6h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/6h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_6h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/6h vs D10/dotplot_gsea_in_vivo_MTX_6h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}




# 3h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/3h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_3h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D3", "3h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/3h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_3h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/3h vs D3/dotplot_gsea_in_vivo_MTX_3h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/MTX/3h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_3h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "MTX", "D10", "3h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/3h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_MTX_3h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "MTX_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/MTX/3h vs D10/dotplot_gsea_in_vivo_MTX_3h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}





## OXA


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/24h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_24h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D3", "24h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/24h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_24h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/24h vs D3/dotplot_gsea_in_vivo_OXA_24h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/24h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_24h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D10", "24h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/24h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_24h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/24h vs D10/dotplot_gsea_in_vivo_OXA_24h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 18h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/18h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_18h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D3", "18h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/18h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_18h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/18h vs D3/dotplot_gsea_in_vivo_OXA_18h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/18h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_18h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D10", "18h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/18h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_18h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/18h vs D10/dotplot_gsea_in_vivo_OXA_18h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}




# 12h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/12h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_12h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D3", "12h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/12h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_12h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/12h vs D3/dotplot_gsea_in_vivo_OXA_12h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/12h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_12h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D10", "12h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/12h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_12h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/12h vs D10/dotplot_gsea_in_vivo_OXA_12h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 6h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/6h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_6h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D3", "6h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/6h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_6h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/6h vs D3/dotplot_gsea_in_vivo_OXA_6h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/6h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_6h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D10", "6h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/6h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_6h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/6h vs D10/dotplot_gsea_in_vivo_OXA_6h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}




# 3h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/3h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_3h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D3", "3h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/3h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_3h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/3h vs D3/dotplot_gsea_in_vivo_OXA_3h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Oxa/3h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_3h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "OXA", "D10", "3h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/3h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_OXA_3h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "OXA_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/OXA/3h vs D10/dotplot_gsea_in_vivo_OXA_3h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}






## CIS


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/24h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_24h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D3", "24h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/24h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_24h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/24h vs D3/dotplot_gsea_in_vivo_CIS_24h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/24h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_24h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D10", "24h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/24h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_24h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/24h vs D10/dotplot_gsea_in_vivo_CIS_24h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 18h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/18h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_18h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D3", "18h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/18h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_18h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/18h vs D3/dotplot_gsea_in_vivo_CIS_18h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/18h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_18h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D10", "18h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/18h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_18h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/18h vs D10/dotplot_gsea_in_vivo_CIS_18h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}




# 12h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/12h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_12h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D3", "12h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/12h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_12h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/12h vs D3/dotplot_gsea_in_vivo_CIS_12h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/12h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_12h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D10", "12h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/12h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_12h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/12h vs D10/dotplot_gsea_in_vivo_CIS_12h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}



# 6h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/6h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_6h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D3", "6h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/6h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_6h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/6h vs D3/dotplot_gsea_in_vivo_CIS_6h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/6h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_6h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D10", "6h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/6h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_6h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/6h vs D10/dotplot_gsea_in_vivo_CIS_6h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}




# 3h

filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/3h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_3h_D3_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D3", "3h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/3h vs D3/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_3h_D3_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D3_vs_PBS_D3", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/3h vs D3/dotplot_gsea_in_vivo_CIS_3h_vs_D3.pdf", nmax=20)
} else{
  
  print("coucou")
}


filename <- "~/MODICED_link/Enrichment_overlap_in_vitro_in_vivo/GO/Cis/3h_vs_D3+D10+D22/enriched_terms_all_in_vivo_overlap_3h_D10_df.csv"

if (file.exists(filename)){
  
  enrich_plot_go_for_publi(filename, "CIS", "D10", "3h", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/3h vs D10/dotplot_in_vivo_overlap_in_vitro_in_vivo_CIS_3h_D10_all.pdf", "~/MODICED_link/top_table_in_vivo_for_Andrea.csv", "CIS_D10_vs_PBS_D10", "In vivo", "~/MODICED_link/figures_for_publication/In_vivo/all_diff_genes_ORA_GSEA/CIS/3h vs D10/dotplot_gsea_in_vivo_CIS_3h_vs_D10.pdf", nmax=20)
} else{
  
  print("coucou")
}


