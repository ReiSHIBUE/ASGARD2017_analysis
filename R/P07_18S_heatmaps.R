### 07_18S_heatmaps.R
### ASGARD 2017 Processing Site Analysis — 18S Eukaryote Heatmaps
### 18S真核生物ヒートマップスクリプト
###
### REQUIRES (from 01–03):
###   meta_asgard_p2        - metadata for 78 samples
###   asgard_filtered_p_hm2 - 16S ASV matrix, 78*221
###   clusnum_p             - cluster assignments (length 78)
###   sample_rgb3           - row colours by cluster (length 78)
###   h3                    - heatmap.2 object from 03_heatmaps_16S.R
###
### PRODUCES (used by 08_esv_heatmap.R):
###   asgard_euk_class_hm_filtered - 18S class-level heatmap matrix (74*66)
###   sample_rgb4                  - row colours for 74 samples (NAs removed)
###   h3_74                        - heatmap.2 object re-run on 74 samples
###
### DATA FILES (read from data/):
###   data/class_relabund_by_station_depth.tsv   (18S class-level relative abundance)
###   data/phylum_relabund_by_station_depth.tsv  (18S phylum-level relative abundance)
###
### OUTPUT:
###   output/18Sclass_boxplots_processing.pdf
###   output/18Sphylum_boxplots_processing.pdf
###   output/ASGARD_hm_processing_18S.pdf

library(tidyverse)
library(gplots)
library(viridis)
library(vegan)

# ==============================================================================
# Section 1: 18Sデータの読み込み / Read 18S relative abundance data
# ==============================================================================

euk_class  <- read.table(
  here::here("data", "class_relabund_by_station_depth.tsv"),
  sep = "\t", header = TRUE, row.names = 1
) # 851*122

euk_phylum <- read.table(
  here::here("data", "phylum_relabund_by_station_depth.tsv"),
  sep = "\t", header = TRUE, row.names = 1
) # 851*61

# ==============================================================================
# Section 2: Processing station の ASGARD サンプルだけ抽出
# Filter to processing stations
# ==============================================================================

euk_class_p  <- euk_class  %>% filter(station_type == "P") # 291*122
euk_phylum_p <- euk_phylum %>% filter(station_type == "P") # 291*61

# メタデータ列を除去 / Remove metadata columns
euk_class_p  <- euk_class_p  %>%
  select(-c("station", "depth_m", "station_type", "filter")) # 291*118

euk_phylum_p <- euk_phylum_p %>%
  select(-c("station", "depth_m", "station_type", "filter")) # 291*57

# Sample列を追加してleft_joinの準備 / Add Sample column for join
euk_class_p$Sample  <- rownames(euk_class_p)
euk_phylum_p$Sample <- rownames(euk_phylum_p)

# ASGARDサンプル(78)のみに絞る / Keep only ASGARD samples (78 samples from meta_asgard_p2)
meta_asgard_p_euk_class <- left_join(meta_asgard_p2, euk_class_p,  by = "Sample") # 78*164
rownames(meta_asgard_p_euk_class) <- meta_asgard_p_euk_class$Sample

meta_asgard_p_euk_phylum <- left_join(meta_asgard_p2, euk_phylum_p, by = "Sample") # 78*103
rownames(meta_asgard_p_euk_phylum) <- meta_asgard_p_euk_phylum$Sample

# ==============================================================================
# Section 3: クラスター別 boxplots / Cluster-wise boxplots
# ==============================================================================

euk_clusnum_p <- factor(clusnum_p, levels = c("1", "2", "3", "4"))

# 数値列のみ抽出 / Extract numeric columns for plotting
num_cols_class  <- sapply(meta_asgard_p_euk_class,  is.numeric)
num_df_class    <- meta_asgard_p_euk_class[, num_cols_class]

num_cols_phylum <- sapply(meta_asgard_p_euk_phylum, is.numeric)
num_df_phylum   <- meta_asgard_p_euk_phylum[, num_cols_phylum]

# 18S class のboxplots / Boxplots for 18S class-level data
pdf(file = here::here("output", "boxplots", "18Sclass_boxplots_processing.pdf"))

for (var in colnames(num_df_class)) {
  gg <- ggplot(num_df_class, aes(x = euk_clusnum_p, y = .data[[var]])) +
    geom_boxplot(aes(fill = euk_clusnum_p), outlier.shape = NA) +
    geom_jitter(width = .4, height = 0) +
    scale_y_log10() +
    theme(text = element_text(size = 24))
  print(gg)
}

dev.off()

# 18S phylum のboxplots / Boxplots for 18S phylum-level data
pdf(file = here::here("output", "boxplots", "18Sphylum_boxplots_processing.pdf"))

for (var in colnames(num_df_phylum)) {
  gg <- ggplot(num_df_phylum, aes(x = euk_clusnum_p, y = .data[[var]])) +
    geom_boxplot(aes(fill = euk_clusnum_p), outlier.shape = NA) +
    geom_jitter(width = .4, height = 0) +
    scale_y_log10() +
    theme(text = element_text(size = 24))
  print(gg)
}

dev.off()

# ==============================================================================
# Section 4: 18S class ヒートマップ / 18S class-level heatmap
# ==============================================================================

# メタデータ列26列を除いてASVのみのmatrixを作成 / Keep only taxon columns (exclude first 26 metadata cols)
asgard_euk_class_hm <- num_df_class[, -(1:26)]
asgard_euk_class_hm <- as.matrix(asgard_euk_class_hm) # 78*118

# 合計 > 0 の列のみ残す / Keep columns with positive sum
col_sums <- colSums(asgard_euk_class_hm, na.rm = TRUE)
asgard_euk_class_hm_filtered <- asgard_euk_class_hm[, col_sums > 0] # 78*66

# NAを含む行を削除 (行36, 38, 41, 73) / Remove rows with NAs
asgard_euk_class_hm_filtered <- asgard_euk_class_hm_filtered[
  complete.cases(asgard_euk_class_hm_filtered), ] # 74*66

sample_rgb4 <- sample_rgb3[-c(36, 38, 41, 73)] # 74 colours

# 16S heatmap を74サンプルで再描画 / Replot 16S heatmap for 74-sample subset
asgard_filtered_p_hm3 <- asgard_filtered_p_hm2[rownames(asgard_euk_class_hm_filtered), ]

pdf(file = here::here("output", "heatmaps", "ASGARD_hm_processing_18S.pdf"), width = 20, height = 20)

# h3_74: 16S heatmap for 74 samples aligned to 18S data
h3_74 <- heatmap.2(
  (asgard_filtered_p_hm3)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  margins   = c(15, 15),
  scale     = "none",
  main      = "ASGARD_bray/ward.D2 (74 samples)",
  trace     = "none",
  cexCol    = 0.2, # 列の文字サイズ / column label size
  cexRow    = 0.2, # 行の文字サイズ / row label size
  labRow    = meta_asgard_p2$side[match(rownames(asgard_filtered_p_hm3),
                                        meta_asgard_p2$Sample)]
)

# クラスター再抽出 (74サンプル用) / Re-extract clusters for 74-sample set
nclus_p_euk   <- 4
oldclus_p_euk  <- cutree(as.hclust(h3_74$rowDendrogram), k = nclus_p_euk)
oldorder_p_euk <- unname(rle(oldclus_p_euk[as.hclust(h3_74$rowDendrogram)$order])$values)
neworder_p_euk <- (1:nclus_p_euk)
names(neworder_p_euk) <- oldorder_p_euk
clusnum_p_euk  <- unname(neworder_p_euk[as.character(oldclus_p_euk)])
names(clusnum_p_euk) <- names(oldclus_p_euk)

rsc_p_euk   <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb4 <- rsc_p_euk[clusnum_p_euk] # 74 colours

# h10: 18S class ヒートマップ (h3_74 デンドログラム固定)
# h10: 18S class heatmap using dendrogram from h3_74
h10 <- heatmap.2(
  (asgard_euk_class_hm_filtered)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  RowSideColors = sample_rgb4,
  Rowv      = h3_74$rowDendrogram,
  margins   = c(15, 15),
  scale     = "none",
  main      = "ASGARD 18S class heatmap (Bray/ward.D)",
  trace     = "none",
  cexCol    = 0.8, # 列の文字サイズ / column label size
  cexRow    = 0.8  # 行の文字サイズ / row label size
)

dev.off()

message("07_18S_heatmaps.R: done. asgard_euk_class_hm_filtered (74x", ncol(asgard_euk_class_hm_filtered), "), sample_rgb4 (length=", length(sample_rgb4), ") ready.")
