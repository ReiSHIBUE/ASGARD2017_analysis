### 02_ternary_plots.R
### ASGARD 2017 Processing Site Analysis — Ternary Diagrams
### 三角図（ternary plot）スクリプト
###
### REQUIRES (from 01_data_prep.R):
###   asgard_processing   - 81-sample merged df
###   asgard_processing2  - 78-sample merged df
###   meta_asgard_p2      - metadata for 78 samples
###
### PRODUCES (used by later scripts):
###   ternary_prop         - ternary proportions for 81-sample set (223*3 matrix)
###   ternary_prop_color   - df with dominant filter column (81-sample set)
###   ternary_prop2        - ternary proportions for 78-sample set (221*3 matrix)
###   ternary_prop_color2  - df with dominant filter column (78-sample set)
###   asv_rgb2             - RGB color vector for ASV columns (length 221)
###   sample_rgb2          - RGB color vector for sample rows (length 78)
###   asv_names_0.2        - ASVs dominant in 0.2 µm fraction
###   asv_names_3          - ASVs dominant in 3 µm fraction
###   asv_names_20         - ASVs dominant in 20 µm fraction
###   zero_cols            - ASV names with zero abundance in 0.2 µm samples

# ggtern のインストールと読み込み / Install and load ggtern
if (!requireNamespace("ggtern", quietly = TRUE)) install.packages("ggtern")
library(ggtern)
library(tidyverse)

# ==============================================================================
# Section 1: 81-sample ternary plot (全サンプル)
# ==============================================================================

# フィルターサイズ列とESV列だけ抽出 / Keep filter and ESV columns only
asgard_processing_frac <- asgard_processing %>%
  select(filter, contains("ESV")) # 81*224

# フィルターサイズごとの平均を計算 / Compute per-filter-size mean per ASV
df_for_ternary_ave <- asgard_processing_frac %>%
  group_by(filter) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  select(-filter)
# across(): 全列に同じ操作を適用; na.rm = TRUE で欠損値を除外して平均を計算

rownames(df_for_ternary_ave) <- c("0.2 µm", "20 µm", "3 µm")
df_for_ternary_ave <- t(df_for_ternary_ave) # 転置 / transposition
ternary_prop <- prop.table(df_for_ternary_ave, margin = 1) # margin=1 で行ごとのproportion

# 各ASVを支配的なフィルターサイズで色付け / Colour each ASV by dominant filter size
ternary_prop_color <- as.data.frame(ternary_prop)
ternary_prop_color$dominant <- apply(ternary_prop_color, 1, function(x) {
  colnames(ternary_prop_color)[which.max(x)]
})

ternary_prop_color$dominant <- factor(
  ternary_prop_color$dominant,
  levels = c("0.2 µm", "3 µm", "20 µm")
)

# 三角図を描画 / Draw ternary diagram — 81 samples
p1 <- ggtern(data = ternary_prop_color,
       aes(x = `0.2 µm`, y = `3 µm`, z = `20 µm`, color = dominant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(`0.2 µm` = "red", `3 µm` = "green", `20 µm` = "blue")) +
  labs(x = "0.2 µm", y = "3 µm", z = "20 µm") +
  theme_bw() +
  theme_showarrows() +
  theme_rotate(-90)

# 支配的フィルターサイズで閾値抽出 / Extract ASVs dominant in each size fraction
ext_0.2 <- ternary_prop_color %>% filter(`0.2 µm` >= 0.7)
asv_names_0.2 <- rownames(ext_0.2) # 72 ASVs

ext_3 <- ternary_prop_color %>% filter(`3 µm` >= 0.7)
asv_names_3 <- rownames(ext_3) # 9 ASVs

ext_20 <- ternary_prop_color %>% filter(`20 µm` >= 0.7)
asv_names_20 <- rownames(ext_20) # 36 ASVs

# ==============================================================================
# Section 2: 78-sample ternary plot (≥5000 reads filter applied)
# 5000 reads以上のサンプルのみを使った三角図
# ==============================================================================

asgard_processing_frac2 <- asgard_processing2 %>%
  select(filter, contains("ESV")) # 78*222

df_for_ternary_ave2 <- asgard_processing_frac2 %>%
  group_by(filter) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  select(-filter)

rownames(df_for_ternary_ave2) <- c("0.2 µm", "20 µm", "3 µm")
df_for_ternary_ave2 <- t(df_for_ternary_ave2) # 転置 / transposition
ternary_prop2 <- prop.table(df_for_ternary_ave2, margin = 1)

ternary_prop_color2 <- as.data.frame(ternary_prop2)
ternary_prop_color2$dominant <- apply(ternary_prop_color2, 1, function(x) {
  colnames(ternary_prop_color2)[which.max(x)]
})

ternary_prop_color2$dominant <- factor(
  ternary_prop_color2$dominant,
  levels = c("0.2 µm", "3 µm", "20 µm")
)

p2 <- ggtern(data = ternary_prop_color2,
       aes(x = `0.2 µm`, y = `3 µm`, z = `20 µm`, color = dominant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(`0.2 µm` = "red", `3 µm` = "green", `20 µm` = "blue")) +
  labs(x = "0.2 µm", y = "3 µm", z = "20 µm") +
  theme_bw() +
  theme_showarrows() +
  theme_rotate(-90)

# PDFに保存 / Save both ternary plots to PDF
pdf(file = here::here("output", "ternary", "ternary_plots.pdf"), width = 8, height = 7)
print(p1)
print(p2)
dev.off()

# ==============================================================================
# Section 3: ヒートマップ用カラーベクトル / Colour vectors for heatmaps
# ==============================================================================

# ASV列の色: ternary plotのRGB比率を使う / ASV column colours from ternary RGB ratios
asv_rgb2 <- rgb(
  ternary_prop_color2$`0.2 µm`,
  ternary_prop_color2$`3 µm`,
  ternary_prop_color2$`20 µm`
) # length 221

# サンプル行の色: フィルターサイズ別 / Sample row colours by filter size
filter_colors <- c("0.2 µm" = "red", "3 µm" = "green", "20 µm" = "blue")
sample_rgb2 <- filter_colors[as.character(meta_asgard_p2$filter)] # length 78

# ==============================================================================
# Section 4: zero_cols — 0.2 µm サンプルで全て0のASV列を特定
# Identify ASV columns that are all-zero in the 0.2 µm fraction
# ==============================================================================

# 地図用dfを先に作っておく (mapは04_maps.Rで描画) / Prep the map df (plotting done in 04)
asgard_processing_ggmap <- asgard_processing2 %>%
  select(filter, lat, lon, contains("ESV")) # 78*226

esv <- asgard_processing_ggmap %>%
  select(lat, lon, filter, contains("ESV_")) # 78*224

# 0.2 µm サンプルだけ抽出 / Keep only 0.2 µm samples
esv0.2 <- esv %>%
  filter(filter == "0.2 µm") # 23*224

# ESV列のみ / ESV columns only
esv_only <- esv0.2 %>%
  select(contains("ESV_")) # 23*221

# 全て0の列名を取得 / Names of columns that are entirely zero in 0.2 µm
zero_cols <- names(esv_only)[colSums(esv_only != 0) == 0]
# 3と20 µmで豊富だが0.2 µmにはいないASVたち / ASVs absent in free-living fraction

message("02_ternary_plots.R: done. zero_cols (length=", length(zero_cols), "), asv_rgb2/sample_rgb2 ready.")
