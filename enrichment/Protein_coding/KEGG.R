# Necessary libraries
library(readr)
library(patchwork)
library(ggplot2)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("~/Desktop/研究/3rd_paper")

gene_list <- read_csv('deseq2_results_protein_coding.csv')

# Convert gene symbols to Entrez IDs
convertToEntrez <- function(genes) {
  mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
}

# Perform enrichment for KEGG pathway
performKeggEnrichment <- function(entrez_ids) {
  enrichKEGG(gene         = entrez_ids,
             organism     = 'hsa',
             pAdjustMethod = "BH",
             qvalueCutoff  = 0.05)
}

# Assuming you have DESeq2 results stored in the variable "res_df"
differentially_expressed_genes <- subset(gene_list, abs(log2FoldChange) > 2 & pvalue < 0.05)
upregulated_genes <- subset(gene_list, log2FoldChange > 2 & pvalue < 0.05)
downregulated_genes <- subset(gene_list, log2FoldChange < -2 & pvalue < 0.05)

# Convert gene symbols to Entrez IDs
entrez_ids_deg <- convertToEntrez(differentially_expressed_genes$SYMBOL)
entrez_ids_up <- convertToEntrez(upregulated_genes$SYMBOL)
entrez_ids_down <- convertToEntrez(downregulated_genes$SYMBOL)

# Enrichment analysis for KEGG pathway
ekegg_deg <- performKeggEnrichment(entrez_ids_deg)
ekegg_up <- performKeggEnrichment(entrez_ids_up)
ekegg_down <- performKeggEnrichment(entrez_ids_down)

# Convert to data frame
ekegg_deg_df <- as.data.frame(ekegg_deg)
ekegg_up_df <- as.data.frame(ekegg_up)
ekegg_down_df <- as.data.frame(ekegg_down)

# Add a column to indicate whether the enrichment is for upregulated or downregulated genes
ekegg_deg_df$Regulation <- "Deg"
ekegg_up_df$Regulation <- "Up"
ekegg_down_df$Regulation <- "Down"

visualizeKeggEnrichment <- function(ekegg_deg_df, ekegg_up_df, ekegg_down_df) {
  # Combine the results
  combined_results <- rbind(ekegg_deg_df, ekegg_up_df, ekegg_down_df)
  
  # Convert GeneRatio to decimal
  parts <- strsplit(as.character(combined_results$GeneRatio), "/")
  combined_results$GeneRatio <- as.numeric(sapply(parts, `[`, 1)) / as.numeric(sapply(parts, `[`, 2))
  
  # Reorder the levels of 'Description' factor based on 'GeneRatio' in descending order
  combined_results$Description <- reorder(combined_results$Description, combined_results$GeneRatio)
  
  # Plot
  p <- ggplot(combined_results, aes(x=GeneRatio, y=Description)) +
    geom_point(aes(size=Count, color=-log10(p.adjust))) + 
    scale_color_continuous(low="blue", high="red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(Regulation ~ ., scales = "free")  +
    theme_bw()
  
  return(p)
}

# Plot
p <- visualizeKeggEnrichment(ekegg_deg_df, ekegg_up_df, ekegg_down_df)
print(p)
