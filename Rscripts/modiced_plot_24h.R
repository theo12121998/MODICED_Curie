# Get DEA figures in vitro

library(tidyverse)
library(RColorBrewer)

tt <- read.csv("/home/theo/Documents/MODICED/top_table_drug_control_in_vitro.csv")

cytokines <-c("Ccl5", "Ccl4", "Ccl20", "Il6")


tt_MTX_24 <- tt %>% filter(contrast == "MTX_24h_vs_Control_24h")


volcano_mtx_24 <- ggplot(tt_MTX_24, aes(x=logFC, y=-log10(adj.P.Val), label=genes)) +
                        geom_point(data=filter(tt_MTX_24, logFC > 2 & adj.P.Val <= 0.05), color="dodgerblue", size=1,shape=21, fill="dodgerblue") +
                        geom_point(data=filter(tt_MTX_24, logFC < -2 &  adj.P.Val <= 0.05), color="lightskyblue", size=1, shape=21, fill="lightskyblue") +
                        geom_point(data=filter(tt_MTX_24, abs(logFC) < 2 | adj.P.Val > 0.05), color="lightgray", size=1,shape=21, fill="lightgray") + 
                        ggrepel::geom_text_repel(data=filter(tt_MTX_24, (abs(logFC) > 2 & adj.P.Val <= 10e-20)| genes %in% cytokines),aes(label = genes) , max.overlaps = Inf, fontface="bold",family="Arial") +
                        geom_hline(yintercept = -log10(0.05), colour = "black", linetype
             ="dashed") + geom_vline(xintercept = c(-2,2), colour="black", linetype
             ="dashed") + ggtitle("MTX vs CTRL - 24h") + theme_bw() +theme(plot.title = element_text(hjust=0.5,face="bold")) + theme(axis.title = element_text(face="bold")) + theme(axis.text = element_text(face="bold")) + labs(x = "log2FoldChange", y = "-log10(FDR)")

ggsave("/home/theo/Documents/MODICED/volcano_mtx_24.pdf",volcano_mtx_24,width=9,height=6,device = cairo_pdf)


tt_OXA_24 <- tt %>% filter(contrast == "OXA_24h_vs_Control_24h")



volcano_oxa_24 <- ggplot(tt_OXA_24, aes(x=logFC, y=-log10(adj.P.Val), label=genes)) +
                        geom_point(data=filter(tt_OXA_24, logFC > 2 & adj.P.Val <= 0.05), color="tomato2", size=1,shape=21, fill="tomato2") +
                        geom_point(data=filter(tt_OXA_24, logFC < -2 &  adj.P.Val <= 0.05), color="tan1", size=1, shape=21, fill="tan1") +
                        geom_point(data=filter(tt_OXA_24, abs(logFC) < 2 | adj.P.Val > 0.05), color="lightgray", size=1,shape=21, fill="lightgray") + 
                        ggrepel::geom_text_repel(data=filter(tt_OXA_24, (abs(logFC) > 3.5 & adj.P.Val <= 10e-25) | genes %in% cytokines),aes(label = genes) ,  max.overlaps = Inf, fontface="bold",family="Arial") +
                        geom_hline(yintercept = -log10(0.05), colour = "black", linetype
             ="dashed") + geom_vline(xintercept = c(-2,2), colour="black", linetype
             ="dashed") + ggtitle("OXA vs CTRL - 24h") + theme_bw() +theme(plot.title = element_text(hjust=0.5,face="bold")) + theme(axis.title = element_text(face="bold")) + theme(axis.text = element_text(face="bold")) + labs(x = "log2FoldChange", y = "-log10(FDR)")

ggsave("/home/theo/Documents/MODICED/volcano_oxa_24.pdf",volcano_oxa_24,width=9,height=6,device = cairo_pdf)





tt_CIS_24 <- tt %>% filter(contrast == "CIS_24h_vs_Control_24h")


volcano_cis_24 <- ggplot(tt_CIS_24, aes(x=logFC, y=-log10(adj.P.Val), label=genes)) +
                        geom_point(data=filter(tt_CIS_24, logFC > 2 & adj.P.Val <= 0.05), color="seagreen4", size=1,shape=21, fill="seagreen4") +
                        geom_point(data=filter(tt_CIS_24, logFC < -2 &  adj.P.Val <= 0.05), color="chartreuse3", size=1, shape=21, fill="chartreuse3") +
                        geom_point(data=filter(tt_CIS_24, abs(logFC) < 2 | adj.P.Val > 0.05), color="lightgray", size=1,shape=21, fill="lightgray") + 
                        ggrepel::geom_text_repel(data=filter(tt_CIS_24, (abs(logFC) > 4 & adj.P.Val <= 10e-30) | genes %in% cytokines ),aes(label = genes) , max.overlaps = Inf, fontface="bold",family="Arial") +
                        geom_hline(yintercept = -log10(0.05), colour = "black", linetype
             ="dashed") + geom_vline(xintercept = c(-2,2), colour="black", linetype
             ="dashed") + ggtitle("CIS vs CTRL - 24h") + theme_bw() +theme(plot.title = element_text(hjust=0.5,face="bold")) + theme(axis.title = element_text(face="bold")) + theme(axis.text = element_text(face="bold")) + labs(x = "log2FoldChange", y = "-log10(FDR)")



ggsave("/home/theo/Documents/MODICED/volcano_cis_24.pdf",volcano_cis_24,width=9,height=6,device = cairo_pdf)







