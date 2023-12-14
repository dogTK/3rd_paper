# 必要なライブラリを読み込む
library(DESeq2)
library(tidyverse)
library(readxl)
library(biomaRt)

setwd("~/Desktop/research/3rd_paper/rawdata/")
# データの読み込み
data <- read_xlsx("rna-seq/PR3462_gene_expression.xlsx")

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

# フィルタリングされたデータに基づいて `Transcript ID` を追加
res$TranscriptId <- rownames(res)

df_long <- as.data.frame(res) %>% 
  tidyr::separate_rows(TranscriptId, sep = ";")


# Ensemblのデータベースに接続
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# DESeq2の結果からTranscript IDを取得
transcript_ids <- df_long$TranscriptId
transcript_ids <- gsub("\\..*", "", transcript_ids)

# BioMartを使用してRefSeq IDを取得
transcript_refseq <- getBM(attributes=c('ensembl_transcript_id', 'refseq_mrna'),
                           filters='ensembl_transcript_id',
                           values=transcript_ids,
                           mart=ensembl)

# DESeq2の結果にRefSeq IDを追加
df_long <- left_join(df_long, transcript_refseq, by = c("TranscriptId" = "ensembl_transcript_id"))


setwd("~/Desktop/research/3rd_paper/DEseq2/results/")

# 結果を保存
write.csv(res, "deseq2_results_transcriptId.csv")






