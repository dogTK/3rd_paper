library(dplyr)

# データを読み込む
data <- read.csv("combined_string_rnaseq_wgs.csv")

# IMPACTを数値に変換
impact_scores <- case_when(
  data$IMPACT == "HIGH" ~ 2,
  data$IMPACT == "MODERATE" ~ 1,
  TRUE ~ NA_real_
)

# Zスコアで正規化する関数
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# データフレームにIMPACTスコアを追加
data <- data %>%
  mutate(
    IMPACT_SCORE = impact_scores,
    Z_degree.layout = ifelse(!is.na(degree.layout), z_score(degree.layout), NA),
    Z_ABS_log2FoldChange = z_score(ABS_log2FoldChange),
    Z_padj = -z_score(padj), # p値は低いほうが良いため、符号を反転
    Z_IMPACT_SCORE = ifelse(!is.na(IMPACT_SCORE), z_score(IMPACT_SCORE), NA)
  )

# 各列の重みを定義
weights <- c(
  "Z_degree.layout" = 1,
  "Z_ABS_log2FoldChange" = 1,
  "Z_padj" = 1,
  "Z_IMPACT_SCORE" = 1
)

# スコアを合成
data <- data %>%
  rowwise() %>%
  mutate(
    Combined_Score = sum(c_across(all_of(names(weights))) * weights, na.rm = TRUE)
  ) %>%
  ungroup()

# スコアでランキング
data <- data %>%
  arrange(desc(Combined_Score))

# 上位の遺伝子を見る
head(data)

write_csv(data, 'z-score_string_rnaseq_wgs.csv')
