# 必要なライブラリを読み込む
library(DESeq2)
library(tidyverse)
library(readxl)

setwd("~/Desktop/研究/3rd_paper/gigafile-1008-22137d64530c2b72325cf76f3a119b7b")
# データの読み込み
data <- read_xlsx("PR3462_gene_expression.xlsx")

# TPMデータとサンプル情報を分離
count_data <- round(data[,c("PG5103_01_d.count", "PG5103_03_d.count", "PG5103_04_d.count", "PG5103_05_d.count")])
colData <- data.frame(condition = c("patient","normal", "normal", "normal"))

# DESeqデータセットを作成
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~ condition)

# DESeq2を実行
dds <- DESeq(dds)

# 結果の取得
res <- results(dds)

# Gene Typeが"protein_coding"の遺伝子のインデックスを取得
protein_coding_genes_indices <- which(data$`Gene Type` == "protein_coding")

# 上記のインデックスを使用してresから遺伝子をフィルタリング
res_filtered <- res[protein_coding_genes_indices,]

# `Gene Name` カラムからタンパク質のシンボルを取得して、フィルタリングされた結果セットに追加
res_filtered$SYMBOL <- data$`Gene Name`[protein_coding_genes_indices]

setwd("~/Desktop/研究/3rd_paper/DEseq2/")

# 結果を保存
write.csv(res_filtered, "deseq2_results_protein_coding.csv")
