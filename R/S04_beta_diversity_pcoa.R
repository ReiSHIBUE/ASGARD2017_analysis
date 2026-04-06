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
asgard_pcoa_df$cluster10 <- factor(clusnum10[asgard_pcoa_df$Sample], levels = as.character(1:10))

# ==============================================================================
# Section 3b: PCoA plots / PCoAプロット（Bray-Curtis, Jaccard, Euclidean）
# ==============================================================================

# Calculate % variance explained per axis
bray_eig <- asgard_pcoa_bray$eig
bray_pct1 <- round(bray_eig[1] / sum(bray_eig[bray_eig > 0]) * 100, 1)
bray_pct2 <- round(bray_eig[2] / sum(bray_eig[bray_eig > 0]) * 100, 1)

jac_eig <- asgard_pcoa_jaccard$eig
jac_pct1 <- round(jac_eig[1] / sum(jac_eig[jac_eig > 0]) * 100, 1)
jac_pct2 <- round(jac_eig[2] / sum(jac_eig[jac_eig > 0]) * 100, 1)

euc_eig <- asgard_pcoa_euclidean$eig
euc_pct1 <- round(euc_eig[1] / sum(euc_eig[euc_eig > 0]) * 100, 1)
euc_pct2 <- round(euc_eig[2] / sum(euc_eig[euc_eig > 0]) * 100, 1)

pdf(file = here::here("output", "survey", "beta_diversity", "ASGARD_pcoa_survey.pdf"),
    width = 10, height = 8)

# Page 1: Bray-Curtis PCoA
print(ggplot(asgard_pcoa_df, aes(x = PCoA1_Bray, y = PCoA2_Bray, color = cluster10)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10) +
  coord_cartesian(xlim = range(asgard_pcoa_df$PCoA1_Bray) * 1.1,
                  ylim = range(asgard_pcoa_df$PCoA2_Bray) * 1.1) +
  labs(x = paste0("PCoA1 (", bray_pct1, "%)"),
       y = paste0("PCoA2 (", bray_pct2, "%)"),
       color = "Cluster",
       title = "PCoA — Bray-Curtis (10 clusters)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 2: Jaccard PCoA
print(ggplot(asgard_pcoa_df, aes(x = PCoA1_Jaccard, y = PCoA2_Jaccard, color = cluster10)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10) +
  coord_cartesian(xlim = range(asgard_pcoa_df$PCoA1_Jaccard) * 1.1,
                  ylim = range(asgard_pcoa_df$PCoA2_Jaccard) * 1.1) +
  labs(x = paste0("PCoA1 (", jac_pct1, "%)"),
       y = paste0("PCoA2 (", jac_pct2, "%)"),
       color = "Cluster",
       title = "PCoA — Jaccard (10 clusters)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 3: Euclidean PCoA
print(ggplot(asgard_pcoa_df, aes(x = PCoA1_Euclidean, y = PCoA2_Euclidean, color = cluster10)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10) +
  coord_cartesian(xlim = range(asgard_pcoa_df$PCoA1_Euclidean) * 1.1,
                  ylim = range(asgard_pcoa_df$PCoA2_Euclidean) * 1.1) +
  labs(x = paste0("PCoA1 (", euc_pct1, "%)"),
       y = paste0("PCoA2 (", euc_pct2, "%)"),
       color = "Cluster",
       title = "PCoA — Euclidean (10 clusters)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

dev.off()

# ==============================================================================
# Section 4: Boxplots — 10クラスターごとに各変数を可視化
# Boxplot + PCoA scatter for every variable, grouped by 10 clusters
# ==============================================================================

# 階層名でクラスター列を追加 / Add hierarchical cluster names
asgard_pcoa_df$hier_name <- factor(
  hier_names[as.character(clusnum10[asgard_pcoa_df$Sample])],
  levels = hier_levels
)
asgard_pcoa_df$group <- sub("^(.).*", "\\1", asgard_pcoa_df$hier_name)

pdf(file = here::here("output", "survey", "beta_diversity", "ASGARD_boxplots_survey.pdf"))

plot(asgard_pcoa_df$PCoA1_Bray,
     asgard_pcoa_df$PCoA2_Bray,
     col = rsc10[asgard_pcoa_df$Sample], pch = 19,
     main = "PCoA Bray-Curtis — Survey Samples (10 clusters)",
     xlab = "PCoA1", ylab = "PCoA2")

numeric_vars <- colnames(asgard_pcoa_df)[sapply(asgard_pcoa_df, is.numeric)]
# 有効なデータが十分ある列のみ
numeric_vars <- numeric_vars[sapply(numeric_vars, function(v) sum(!is.na(asgard_pcoa_df[[v]])) > 10)]
for (var in numeric_vars) {
  gg <- ggplot(asgard_pcoa_df, aes(x = hier_name, y = .data[[var]])) +
    geom_boxplot(aes(fill = hier_name), outlier.shape = NA) +
    geom_jitter(width = .4, height = 0) +
    scale_fill_manual(values = cc10) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    labs(title = var) +
    theme(text = element_text(size = 14), legend.position = "none")
  print(gg)

  gp <- ggplot(asgard_pcoa_df) +
    geom_point(aes(x = PCoA1_Bray, y = PCoA2_Bray,
                   col = hier_name, size = .data[[var]])) +
    scale_color_manual(values = cc10)
  print(gp)
}

dev.off()

message("S04_beta_diversity_pcoa.R: done. asgard_pcoa_df (", nrow(asgard_pcoa_df), "×", ncol(asgard_pcoa_df), ") ready.")
