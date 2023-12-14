# 必要なライブラリを読み込む
library(DESeq2)
library(tidyverse)
library(readxl)

setwd("~/Desktop/research/3rd_paper/rawdata/")
# データの読み込み
data <- read_xlsx("PR3462_gene_expression.xlsx")

# TPMデータとサンプル情報を分離
count_data <- round(data[,c("PG5103_01_d.count", "PG5103_03_d.count", "PG5103_04_d.count", "PG5103_05_d.count")])
colData <- data.frame(condition = c("patient","normal", "normal", "normal"))

# DESeqデータセットを作成
rownames(count_data) <- data$`Transcript ID`
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~ condition)

# count>10でプレフィルターする
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# DESeq2を実行
dds <- DESeq(dds)

# 結果の取得
res <- results(dds)

# `Gene Name` カラムからタンパク質のシンボルを取得して、フィルタリングされた結果セットに追加
res$SYMBOL <- data$`Gene Name`

setwd("~/Desktop/research/3rd_paper/DEseq2/")

# 結果を保存
write.csv(res, "deseq2_results.csv")
