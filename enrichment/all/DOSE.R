# Necessary libraries
library(readr)
library(patchwork)
library(ggplot2)
library(clusterProfiler)
library(DOSE)

setwd("~/Desktop/研究/3rd_paper")

gene_list <- read_csv('deseq2_results.csv')

# Convert gene symbols to Entrez IDs
convertToEntrez <- function(genes) {
  mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
}

# Perform enrichment and get unique top terms
performDOEnrichment <- function(entrez_ids) {
  enrichDO(gene          = entrez_ids,
            ont           = "DO",
            pvalueCutoff  = 0.05,
            pAdjustMethod = "BH",
            minGSSize     = 5,
            maxGSSize     = 1500,
            qvalueCutoff  = 0.05,
            readable      = FALSE)
}

# Assuming you have DESeq2 results stored in the variable "res_df"
differentially_expressed_genes <- subset(gene_list, abs(log2FoldChange) > 1 & pvalue < 0.05)
upregulated_genes <- subset(gene_list, log2FoldChange > 1 & pvalue < 0.05)
downregulated_genes <- subset(gene_list, log2FoldChange < -1 & pvalue < 0.05)

# Convert gene symbols to Entrez IDs
entrez_ids_deg <- convertToEntrez(differentially_expressed_genes$SYMBOL)
entrez_ids_up <- convertToEntrez(upregulated_genes$SYMBOL)
entrez_ids_down <- convertToEntrez(downregulated_genes$SYMBOL)

# Enrichment analysis
ego_deg <- performDOEnrichment(entrez_ids_deg)
ego_up <- performDOEnrichment(entrez_ids_up)
ego_down <- performDOEnrichment(entrez_ids_down)

visualizeDOEnrichment <- function(enrichment_result) {
  
  # Convert enrichResult object to a dataframe
  df <- as.data.frame(enrichment_result)
  
  # Convert GeneRatio to decimal
  parts <- strsplit(as.character(df$GeneRatio), "/")
  df$GeneRatio <- as.numeric(sapply(parts, `[`, 1)) / as.numeric(sapply(parts, `[`, 2))
  
  # Order by p.adjust in ascending order and select top 30 entries
  df <- df[order(df$p.adjust)[1:30], ]
  
  # Reorder the levels of 'Description' factor based on 'GeneRatio' in descending order
  df$Description <- reorder(df$Description, df$GeneRatio)
  
  # Plot
  p <- ggplot(df, aes(x=GeneRatio, y=Description)) +
    geom_point(aes(size=Count, color=-log10(p.adjust))) + 
    scale_color_continuous(low="blue", high="red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw()
  
  return(p)
}

# Plot
p <- visualizeDOEnrichment(ego_up)
print(p)
