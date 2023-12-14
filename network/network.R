# 必要なパッケージをインストールします
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("STRINGdb", "igraph"))
library(STRINGdb)
library(igraph)

test_genes <- data.frame(
  Gene = c("UBL5", "NDUFB8", "CHMP2B", "PRPF8", "NOSTRIN", "MFAP1", "CWC22", "PLCH2", "PRPF31", "ATP6AP1", "DSC3"),
  fold_change = c(2.6, 3.5, 4.9, 2.9, 4.2, 3.7, 2.1, 3.1, 3.9, 3, 3.6)
)


example1_mapped <- string_db$map(test_genes, "Gene", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )


# STRINGdbからネットワークのエッジリストを取得
network <- string_db$get_interactions(hits)

# エッジリストからigraphオブジェクトを生成
g <- graph_from_data_frame(network)

# 1. 遺伝子名とSTRING IDの対応関係を作成
gene_to_stringid <- setNames(example1_mapped$Gene, example1_mapped$STRING_id)

# 2. igraphオブジェクトの頂点名を変換
V(g)$name <- gene_to_stringid[V(g)$name]

# ネットワークのプロット
plot(g, layout = layout_with_gem, vertex.size = 10, vertex.label.cex = 0.5, edge.arrow.size = 0.5, main="Gene Interaction Network")
