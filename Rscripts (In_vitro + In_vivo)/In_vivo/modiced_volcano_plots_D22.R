# Get DEA figures in vitro

library(tidyverse)
library(RColorBrewer)

tt <- read.csv("top_table_in_vivo_for_Andrea.csv")

cytokines <-c("Ccl5", "Ccl4", "Ccl20", "Il6")


tt_OXA_D22 <- tt %>% filter(contrast == "OXA_D22_vs_PBS_D22")


volcano_OXA_D22 <- ggplot(tt_OXA_D22, aes(x=logFC, y=-log10(adj.P.Val), label=genes)) +
  geom_point(data=filter(tt_OXA_D22, logFC > 2 & adj.P.Val <= 0.05), color="tomato2", size=1,shape=21, fill="tomato2") +
  geom_point(data=filter(tt_OXA_D22, logFC < -2 &  adj.P.Val <= 0.05), color="tan1", size=1, shape=21, fill="tan1") +
  geom_point(data=filter(tt_OXA_D22, abs(logFC) < 2 | adj.P.Val > 0.05), color="lightgray", size=1,shape=21, fill="lightgray") + 
  ggrepel::geom_text_repel(data=filter(tt_OXA_D22, (abs(logFC) > 2 & adj.P.Val <= 0.05)| genes %in% cytokines),aes(label = genes) , max.overlaps = Inf, fontface="bold",family="Arial") +
  geom_hline(yintercept = -log10(0.05), colour = "black", linetype
             ="dashed") + geom_vline(xintercept = c(-2,2), colour="black", linetype
                                     ="dashed") + ggtitle("OXA vs CTRL - D22") + theme_bw() +theme(plot.title = element_text(hjust=0.5,face="bold")) + theme(axis.title = element_text(face="bold")) + theme(axis.text = element_text(face="bold")) + labs(x = "log2FoldChange", y = "-log10(FDR)")

ggsave("volcano_in_vivo_publi/volcano_OXA_D22.pdf",volcano_OXA_D22,width=9,height=6,device = cairo_pdf, create.dir = TRUE)









