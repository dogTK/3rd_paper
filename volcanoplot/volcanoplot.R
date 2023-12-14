# Load necessary libraries
library(DESeq2)
library(tidyverse)

setwd("~/Desktop/研究/3rd_paper")

res_df = read.csv("deseq2_results.csv")

# ボルケーノプロットの作成
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = ifelse(log2FoldChange > 2 & pvalue < 0.05, "red", 
                                ifelse(log2FoldChange < -2 & pvalue < 0.05, "blue", "black"))), alpha=0.5) +
  scale_color_identity(guide = FALSE) +
  theme_minimal() +
  labs(x="log2(FoldChange)", y="-log10(Padj)") +
  geom_hline(yintercept=-log10(0.05), color="blue", linetype="dashed") +
  geom_vline(xintercept=c(-2,2), color="blue", linetype="dashed")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA)) +
  geom_text(aes(label=ifelse((log2FoldChange > 2 | log2FoldChange < -2) & pvalue < 0.05, as.character(SYMBOL), ""),
                color=ifelse((log2FoldChange > 2 | log2FoldChange < -2) & pvalue < 0.05, "black", NA)),
            hjust=1.5, vjust=1.5, check_overlap = TRUE, size=1.5)

setwd("~/Desktop/研究/3rd_paper/volcanoplot")

ggsave("volcanoplot.png", dpi=300)
