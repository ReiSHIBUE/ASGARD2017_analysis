### 05_beta_diversity_pcoa.R
### ASGARD 2017 Processing Site Analysis — Beta Diversity & PCoA
### ベータ多様性とPCoA（主座標分析）スクリプト
###
### REQUIRES (from 01–03):
###   asgard_filtered_p2      - filtered ASV proportion df, 78*222 (with Sample col)
###   asgard_filtered_p_hm2   - matrix version, 78*221
###   meta_asgard_p2          - metadata for 78 samples
###   clusnum_p               - cluster assignments (length 78)
###   rsc_p                   - 4-colour palette
###   sample_rgb3             - row colours by cluster
###
### PRODUCES:
###   asgard_pcoa_df_p      - PCoA + metadata df for 78 processing samples (78*52)
###   clusnum_p_bp          - clusnum_p as factor, for boxplots
###
### OUTPUT: output/asgard_boxplots_processing.pdf

library(vegan)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# ==============================================================================
# Section 1: ベータ多様性行列の計算 / Compute beta-diversity distance matrices
# ==============================================================================

asgard_beta_p <- asgard_filtered_p_hm2^.25
asgard_beta_p <- asgard_beta_p[, colSums(asgard_beta_p) > 0] # 78*221

# Euclidean, Bray-Curtis, Jaccard の3距離行列 / Three distance matrices
asgard_eucmat_p  <- vegdist(asgard_beta_p, method = "euclidean")
asgard_braymat_p <- vegdist(asgard_beta_p, method = "bray")
asgard_jacmat_p  <- vegdist(asgard_beta_p, method = "jaccard")

# Bray-Curtis 距離のヒートマップ (pheatmap) / Bray-Curtis distance heatmap
pdf(file = here::here("output", "heatmaps", "bray_distance_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  asgard_braymat_p,
  clustering_method = "ward",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  main  = "Bray-Curtis Distance Heatmap",
  border_color = NA
)
dev.off()

# ==============================================================================
# Section 2: PCoA 計算 / Perform PCoA (principal coordinates analysis)
# ==============================================================================

asgard_braydist_p <- as.dist(asgard_braymat_p)
asgard_jacdist_p  <- as.dist(asgard_jacmat_p)
asgard_eucdist_p  <- as.dist(asgard_eucmat_p)

asgard_pcoa_bray_p      <- cmdscale(asgard_braydist_p,  k = 2, eig = TRUE)
asgard_pcoa_jaccard_p   <- cmdscale(asgard_jacdist_p,   k = 2, eig = TRUE)
asgard_pcoa_euclidean_p <- cmdscale(asgard_eucdist_p,   k = 2, eig = TRUE)

# ==============================================================================
# Section 3: PCoA結果をデータフレームへ / Assemble PCoA results into a data frame
# ==============================================================================

asgard_pcoa_df_p <- data.frame(
  Sample           = rownames(asgard_filtered_p2), # 行名だけ取得 / row names only
  PCoA1_Bray       = asgard_pcoa_bray_p$points[, 1],
  PCoA2_Bray       = asgard_pcoa_bray_p$points[, 2],
  PCoA1_Jaccard    = asgard_pcoa_jaccard_p$points[, 1],
  PCoA2_Jaccard    = asgard_pcoa_jaccard_p$points[, 2],
  PCoA1_Euclidean  = asgard_pcoa_euclidean_p$points[, 1],
  PCoA2_Euclidean  = asgard_pcoa_euclidean_p$points[, 2]
) # 78*7

meta_asgard_p2 <- meta_asgard_p2 %>%
  mutate(Sample = rownames(asgard_filtered_p2)) # 78*46

# メタデータとPCoA結果を結合 / Join metadata to PCoA df
asgard_pcoa_df_p <- left_join(asgard_pcoa_df_p, meta_asgard_p2, by = "Sample") # 78*52
asgard_pcoa_df_p$filter <- factor(
  asgard_pcoa_df_p$filter,
  levels = c("0.2 µm", "3 µm", "20 µm")
)

# ==============================================================================
# Section 4: Boxplots — クラスターごとに各変数を可視化
# Boxplots of all metadata & PCoA variables by cluster
# ==============================================================================

# クラスターをfactorに変換 / Convert cluster assignments to factor
rsc_p       <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb3 <- rsc_p[clusnum_p]
clusnum_p_bp <- factor(clusnum_p, levels = c("1", "2", "3", "4"))

pdf(file = here::here("output", "beta_diversity", "asgard_boxplots_processing.pdf"))

# PCoA散布図 / PCoA scatter plot coloured by cluster
plot(asgard_pcoa_df_p$PCoA1_Bray,
     asgard_pcoa_df_p$PCoA2_Bray,
     col = sample_rgb3,
     pch = 19)

# 全変数のboxplot + PCoA点図 / Boxplot + PCoA scatterplot for every variable
for (var in colnames(asgard_pcoa_df_p)) {
  gg <- ggplot(asgard_pcoa_df_p, aes(x = clusnum_p_bp, y = .data[[var]])) +
    geom_boxplot(aes(fill = clusnum_p_bp), outlier.shape = NA) +
    geom_jitter(width = .4, height = 0) +
    theme(text = element_text(size = 24))
  print(gg)

  gp <- ggplot(asgard_pcoa_df_p) +
    geom_point(aes(x = PCoA1_Bray, y = PCoA2_Bray,
                   col = clusnum_p_bp, size = .data[[var]]))
  print(gp)
}

dev.off()

message("05_beta_diversity_pcoa.R: done. asgard_pcoa_df_p (78*", ncol(asgard_pcoa_df_p), ") ready.")
