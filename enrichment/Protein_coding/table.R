library(readr)
library(clusterProfiler)

setwd("~/Desktop/研究/3rd_paper")

gene_list <- read_csv('deseq2_results_protein_coding.csv')

# Convert gene symbols to Entrez IDs
convertToEntrez <- function(gene_list) {
  mapIds(org.Hs.eg.db, keys=gene_list, column="ENTREZID", keytype="SYMBOL", multiVals="first")
}

# Assuming you have DESeq2 results stored in the variable "res_df"
differentially_expressed_genes <- subset(gene_list, abs(log2FoldChange) > 2 & pvalue < 0.05)
upregulated_genes <- subset(gene_list, log2FoldChange > 2 & pvalue < 0.05)
downregulated_genes <- subset(gene_list, log2FoldChange < -2 & pvalue < 0.05)

# Convert gene symbols to Entrez IDs
entrez_ids_deg <- convertToEntrez(differentially_expressed_genes$SYMBOL)
entrez_ids_up <- convertToEntrez(upregulated_genes$SYMBOL)
entrez_ids_down <- convertToEntrez(downregulated_genes$SYMBOL)

# Perform enrichment and get unique top terms
performGoEnrichment <- function(entrez_ids) {
  list(
    BP = enrichGO(entrez_ids, OrgDb="org.Hs.eg.db", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05),
    MF = enrichGO(entrez_ids, OrgDb="org.Hs.eg.db", ont="MF", pAdjustMethod="BH", qvalueCutoff=0.05),
    CC = enrichGO(entrez_ids, OrgDb="org.Hs.eg.db", ont="CC", pAdjustMethod="BH", qvalueCutoff=0.05)
  )
}

ego_deg <- performGoEnrichment(entrez_ids_deg)
ego_up <- performGoEnrichment(entrez_ids_up)
ego_down <- performGoEnrichment(entrez_ids_down)

# ... [previous code up to ego_down enrichment analysis]

# Convert enrichment results to a dataframe
convertToDataFrame <- function(enrichment_results) {
  combined_df <- lapply(names(enrichment_results), function(category) {
    ego = enrichment_results[[category]]
    if (is.null(ego)) return(NULL)  # Skip if the enrichment result is NULL
    
    df <- as.data.frame(ego)
    
    # Convert GeneRatio to decimal
    parts <- strsplit(as.character(df$GeneRatio), "/")
    df$Percentage <- as.numeric(sapply(parts, `[`, 1)) / as.numeric(sapply(parts, `[`, 2)) * 100
    
    df$Category <- category
    return(df)
  })
  
  combined_df <- do.call(rbind, combined_df)
  return(combined_df)
}

# Convert the enrichment results to dataframes
df_deg <- convertToDataFrame(ego_deg)
df_up <- convertToDataFrame(ego_up)
df_down <- convertToDataFrame(ego_down)

# Print the table
print(df_up[, c("Category", "ID", "Description", "Count", "Percentage", "p.adjust")])

setwd("~/Desktop/研究/3rd_paper/enrichment/Protein_coding/result/")

write_csv(df_up[, c("Category", "ID", "Description", "Count", "Percentage", "p.adjust")],file = "table_up.csv")
