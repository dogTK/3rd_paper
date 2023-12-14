BiocManager::install("STRINGdb")
library(STRINGdb)
library(tidyverse)

setwd("~/Desktop/研究/3rd_paper")

data<-read_csv('deseq2_results.csv')

deg_genes <- data$SYMBOL[abs(data$log2FoldChange) > 2]
deg_genes_cleaned <- deg_genes[!is.na(deg_genes)]

deg_genes_df <- data.frame(
  Gene = deg_genes_cleaned,
  fold_change = data$log2FoldChange[match(deg_genes_cleaned, data$SYMBOL)]
)

string_db <- STRINGdb$new(version="11", species=9606) # 9606 is the taxon ID for Homo sapiens, adjust if needed
mapped_genes <- string_db$map(deg_genes_df, "Gene", removeUnmappedRows = TRUE)

hits <- mapped_genes$STRING_id
string_db$plot_network(hits)

interactions <- string_db$get_interactions(mapped_genes$STRING_id)

edge_list <- data.frame(
  Edge = 1:nrow(interactions),
  From = interactions$from,
  To = interactions$to,
  Relationship = "Interaction",  
  Score = interactions$combined_score
)

convert_to_symbol <- function(id, mapping_df) {
  symbol <- mapping_df$Gene[which(mapping_df$STRING_id == id)]
  if (length(symbol) > 0) {
    return(symbol[1])
  } else {
    return(NA)
  }
}

edge_list$From <- sapply(edge_list$From, convert_to_symbol, mapping_df = mapped_genes)
edge_list$To <- sapply(edge_list$To, convert_to_symbol, mapping_df = mapped_genes)

filtered_edge_list <- edge_list[edge_list$Score >= 700, ]

# CSVファイルとして保存
write.csv(filtered_edge_list, "interactions_for_cytoscape.csv", row.names=FALSE)
