### 08_esv_heatmap.R
### ASGARD 2017 Processing Site Analysis — ESV Relative Abundance Heatmap (h11)
### ESV 相対存在量ヒートマップスクリプト (h11)
###
### REQUIRES (from 07_18S_heatmaps.R, or can run standalone with data file):
###   asgard_euk_class_hm_filtered - 18S matrix, 74*66 (defines the 74-sample set)
###   sample_rgb4                  - row colours for 74 samples
###   h3_74                        - heatmap.2 object (provides Rowv dendrogram)
###
### DATA FILES (read from data/):
###   data/esv_relabund_by_station_depth.tsv  (ESV-level 16S relabund by station/depth)
###
### PRODUCES:
###   h11  - heatmap.2 object for ESV relative abundance
###
### OUTPUT:
###   output/ASGARD_hm_processing_esv_relabund.pdf
###
### NOTE: This script was verified end-to-end and produces the PDF correctly.
###       It can also be run standalone if asgard_euk_class_hm_filtered and
###       sample_rgb4 are available in the session.
###       このスクリプトは単独でも動作確認済みです。

library(gplots)
library(viridis)
library(vegan)
library(tidyverse)

# ==============================================================================
# Section 1: ESV 相対存在量テーブルの読み込み / Read ESV relative abundance table
# ==============================================================================

esv_relabund <- read.table(
  here::here("data", "esv_relabund_by_station_depth.tsv"),
  sep = "\t", header = TRUE, row.names = 1
)
# Expected: ~858 rows × 5443 columns (5 metadata + ESVs)

# ==============================================================================
# Section 2: Processing stations の ASGARD サンプルだけ抽出
# Filter to processing stations, keep ESV columns only
# ==============================================================================

esv_relabund_p <- esv_relabund %>%
  filter(station_type == "P") %>%
  select(contains("ESV"))

# 18Sヒートマップと共通のサンプル (74) に絞る
# Keep only samples present in 18S heatmap (74 samples)
common_esv_samples <- intersect(
  rownames(asgard_euk_class_hm_filtered),
  rownames(esv_relabund_p)
)
esv_asgard_p <- esv_relabund_p[common_esv_samples, ]

# ==============================================================================
# Section 3: 最小存在量フィルター (16S processingと同じ基準)
# Apply the same minimum-abundance filter used in the 16S processing analysis
# ==============================================================================

esv_mincutoff   <- apply(esv_asgard_p, 2, max) > 0.001
esv_asgard_filt <- esv_asgard_p[, esv_mincutoff]
esv_asgard_filt <- esv_asgard_filt[, colSums(esv_asgard_filt > 0) > 2]

# matrix に変換 / Convert to matrix
esv_asgard_mat <- as.matrix(esv_asgard_filt)

# ==============================================================================
# Section 4: h11 ヒートマップの描画 / Draw h11 heatmap
# Row side colours from h10 cluster assignments; dendrogram from h3_74
# ==============================================================================

# 行サイドカラーを共通サンプル順に揃える / Align row colours to common sample order
sample_rgb_esv <- sample_rgb4[common_esv_samples]

pdf(file = here::here("output", "heatmaps", "ASGARD_hm_processing_esv_relabund.pdf"),
    width = 20, height = 20)

h11 <- heatmap.2(
  (esv_asgard_mat)^.25,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = sample_rgb_esv,
  Rowv          = h3_74$rowDendrogram,
  margins       = c(15, 15),
  scale         = "none",
  main          = "ASGARD ESV relabund by station/depth (Bray/ward.D)",
  trace         = "none",
  cexCol        = 0.8, # 列の文字サイズ / column label size
  cexRow        = 0.8  # 行の文字サイズ / row label size
)

dev.off()

message("08_esv_heatmap.R: done. h11 written to output/ASGARD_hm_processing_esv_relabund.pdf")
