### S04_beta_diversity_pcoa.R
### ASGARD 2017 Survey Site Analysis — Beta Diversity & PCoA
### ベータ多様性とPCoAスクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02):
###   asgard_frtprop   - fourth-root proportion matrix 181×258
###   meta_asgard      - metadata 181×45
###   clusnum          - cluster assignments (length 181)
###   rsc              - 5-colour palette vector
###
### PRODUCES (consumed by S05, S10):
###   asgard_beta      - filtered fourth-root matrix used for distances
###   asgard_braymat   - Bray-Curtis distance matrix (dist object, 181×181)
###   asgard_eucmat    - Euclidean distance matrix
###   asgard_jacmat    - Jaccard distance matrix
###   asgard_pcoa_df   - PCoA scores + metadata df (181×52+)
###
### OUTPUT:
###   output/survey/beta_diversity/ASGARD_bray_heatmap_survey.pdf
###   output/survey/beta_diversity/ASGARD_boxplots_survey.pdf

library(vegan)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# ==============================================================================
# Section 1: ベータ多様性行列の計算 / Compute beta-diversity distance matrices
# ==============================================================================

asgard_beta     <- asgard_frtprop[, colSums(asgard_frtprop) > 0]
asgard_eucmat   <- vegdist(asgard_beta, method = "euclidean")
asgard_braymat  <- vegdist(asgard_beta, method = "bray")
asgard_jacmat   <- vegdist(asgard_beta, method = "jaccard")

# Bray-Curtis 距離のヒートマップ / Bray-Curtis distance heatmap
pdf(file = here::here("output", "survey", "beta_diversity", "ASGARD_bray_heatmap_survey.pdf"),
    width = 12, height = 10)
pheatmap(
  asgard_braymat,
  clustering_method = "ward",
  color             = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  main              = "Bray-Curtis Distance Heatmap — Survey Samples",
  border_color      = NA
)
dev.off()

# ==============================================================================
# Section 2: PCoA 計算 / Perform PCoA
# ==============================================================================

asgard_braydist <- as.dist(asgard_braymat)
asgard_jacdist  <- as.dist(asgard_jacmat)
asgard_eucdist  <- as.dist(asgard_eucmat)

asgard_pcoa_bray      <- cmdscale(asgard_braydist, k = 2, eig = TRUE)
asgard_pcoa_jaccard   <- cmdscale(asgard_jacdist,  k = 2, eig = TRUE)
asgard_pcoa_euclidean <- cmdscale(asgard_eucdist,  k = 2, eig = TRUE)

# ==============================================================================
# Section 3: PCoA結果をデータフレームへ / Assemble PCoA results into data frame
# ==============================================================================

asgard_pcoa_df <- data.frame(
  Sample           = rownames(asgard_filtered),
  PCoA1_Bray       = asgard_pcoa_bray$points[, 1],
  PCoA2_Bray       = asgard_pcoa_bray$points[, 2],
  PCoA1_Jaccard    = asgard_pcoa_jaccard$points[, 1],
  PCoA2_Jaccard    = asgard_pcoa_jaccard$points[, 2],
  PCoA1_Euclidean  = asgard_pcoa_euclidean$points[, 1],
  PCoA2_Euclidean  = asgard_pcoa_euclidean$points[, 2]
) # 181×7

meta_asgard <- meta_asgard %>%
  mutate(Sample = rownames(asgard_filtered))

asgard_pcoa_df <- left_join(asgard_pcoa_df, meta_asgard, by = "Sample") # 181×52+
asgard_pcoa_df$cluster <- as.factor(clusnum)

# ==============================================================================
# Section 4: Boxplots — クラスターごとに各変数を可視化
# Boxplot + PCoA scatter for every variable, grouped by cluster
# ==============================================================================

clusnum_bp <- factor(clusnum, levels = as.character(1:5))

pdf(file = here::here("output", "survey", "beta_diversity", "ASGARD_boxplots_survey.pdf"))

plot(asgard_pcoa_df$PCoA1_Bray,
     asgard_pcoa_df$PCoA2_Bray,
     col = rsc, pch = 19,
     main = "PCoA Bray-Curtis — Survey Samples",
     xlab = "PCoA1", ylab = "PCoA2")

for (var in colnames(asgard_pcoa_df)) {
  gg <- ggplot(asgard_pcoa_df, aes(x = clusnum_bp, y = .data[[var]])) +
    geom_boxplot(aes(fill = clusnum_bp), outlier.shape = NA) +
    geom_jitter(width = .4, height = 0) +
    theme(text = element_text(size = 24))
  print(gg)

  gp <- ggplot(asgard_pcoa_df) +
    geom_point(aes(x = PCoA1_Bray, y = PCoA2_Bray,
                   col = clusnum_bp, size = .data[[var]]))
  print(gp)
}

dev.off()

message("S04_beta_diversity_pcoa.R: done. asgard_pcoa_df (", nrow(asgard_pcoa_df), "×", ncol(asgard_pcoa_df), ") ready.")
