# Necessary libraries
library(readr)
library(patchwork)
library(ggplot2)
library(clusterProfiler)

setwd("~/Desktop/研究/3rd_paper")

gene_list <- read_csv('deseq2_results_protein_coding.csv')

# Convert gene symbols to Entrez IDs
convertToEntrez <- function(genes) {
  mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
}

# Perform enrichment and get unique top terms
performGoEnrichment <- function(entrez_ids) {
  list(
    BP = enrichGO(entrez_ids, OrgDb="org.Hs.eg.db", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05),
    MF = enrichGO(entrez_ids, OrgDb="org.Hs.eg.db", ont="MF", pAdjustMethod="BH", qvalueCutoff=0.05),
    CC = enrichGO(entrez_ids, OrgDb="org.Hs.eg.db", ont="CC", pAdjustMethod="BH", qvalueCutoff=0.05)  )
}

# Assuming you have DESeq2 results stored in the variable "res_df"
differentially_expressed_genes <- subset(gene_list, abs(log2FoldChange) > 2 & pvalue < 0.05)
upregulated_genes <- subset(gene_list, log2FoldChange > 2 & pvalue < 0.05)
downregulated_genes <- subset(gene_list, log2FoldChange < -2 & pvalue < 0.05)

# Convert gene symbols to Entrez IDs
entrez_ids_deg <- convertToEntrez(differentially_expressed_genes$SYMBOL)
entrez_ids_up <- convertToEntrez(upregulated_genes$SYMBOL)
entrez_ids_down <- convertToEntrez(downregulated_genes$SYMBOL)

# Enrichment analysis
ego_deg <- performGoEnrichment(entrez_ids_deg)
ego_up <- performGoEnrichment(entrez_ids_up)
ego_down <- performGoEnrichment(entrez_ids_down)

visualizeGoEnrichment <- function(enrichment_results) {
  combined_plot <- lapply(enrichment_results, function(ego) {
    df <- as.data.frame(ego@result)
    
    # Convert GeneRatio to decimal
    parts <- strsplit(as.character(df$GeneRatio), "/")
    df$GeneRatio <- as.numeric(sapply(parts, `[`, 1)) / as.numeric(sapply(parts, `[`, 2))
    
    df$ont <- names(enrichment_results)[sapply(enrichment_results, function(x) identical(ego, x))]
    
    # Select top 10 rows with smallest p.adjust for each ont
    df <- df[order(df$p.adjust)[1:10], ]
    
    return(df)
  })
  combined_plot <- do.call(rbind, combined_plot)
  
  # Convert 'ont' column to factor with explicit levels
  combined_plot$ont <- factor(combined_plot$ont, levels = c("BP", "MF", "CC"))
  
  # Reorder the levels of 'Description' factor based on 'GeneRatio' in descending order
  combined_plot$Description <- reorder(combined_plot$Description, combined_plot$GeneRatio)
  
  ggplot(combined_plot, aes(x=GeneRatio, y=Description)) +
    geom_point(aes(size=Count, color=-log10(p.adjust))) + 
    scale_color_continuous(low="blue", high="red") +
    facet_grid(ont ~ ., scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw()
}

# Plot
p <- visualizeGoEnrichment(ego_up)
print(p)
