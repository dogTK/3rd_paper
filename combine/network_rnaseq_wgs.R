library(readr)
library(dplyr)

# データを読み込む
string <- read_csv('STRING network - 2 default node.csv')
rnaseq_wgs <- read_csv('../wgs_filtering/result/rna_seq_wgs.csv')

# `string`から必要な列だけを選択する
string_selected <- select(string, `degree.layout`, `display name`)

# 選択した列と`rnaseq_wgs`を結合する
combined_data <- inner_join(string_selected, rnaseq_wgs, by = c("display name" = "GENE"))

# 結合したデータをCSVに書き出す
setwd("~/Desktop/研究/3rd_paper/network")

write_csv(combined_data, 'combined_string_rnaseq_wgs.csv')
