### 06_dbrda.R
### ASGARD 2017 Processing Site Analysis — db-RDA & Environmental Scatter Plots
### 距離ベースRDA（db-RDA）スクリプト
###
### REQUIRES (from 01–05):
###   asgard_filtered_p2    - filtered ASV df, 78*222
###   asgard_filtered_p_hm2 - matrix version, 78*221
###   asgard_pcoa_df_p      - PCoA + metadata df (78*52)
###   meta_asgard_p2        - metadata for 78 samples
###   clusnum_p             - cluster assignments (length 78)
###   sample_rgb3           - row colours by cluster
###
### PRODUCES:
###   asgard_dbrda_model_p   - db-RDA capscale model object
###   asgard_anova_dbrda_p   - ANOVA test of overall model
###   asgard_anova_env_p     - ANOVA test per environmental variable (by margin)
###   asgard_anova_axes_p    - ANOVA test per db-RDA axis
###   asgard_dbrda_merged_p  - merged scores + metadata df (75 samples with complete env data)

library(vegan)
library(tidyverse)

# ==============================================================================
# Section 1: db-RDA 用ベータ多様性行列
# Beta-diversity matrix for db-RDA (78 processing samples)
# ==============================================================================

asgard_beta_p <- asgard_filtered_p_hm2^.25
asgard_beta_p <- asgard_beta_p[, colSums(asgard_beta_p) > 0] # 78*221

asgard_eucmat_p  <- vegdist(asgard_beta_p, method = "euclidean")
asgard_braymat_p <- vegdist(asgard_beta_p, method = "bray")
asgard_jacmat_p  <- vegdist(asgard_beta_p, method = "jaccard")

# ==============================================================================
# Section 2: 環境変数の選択とスケーリング / Select and scale environmental variables
# ==============================================================================

# NAの数を確認 / Check NA counts per env variable
for (var in c("salinity", "temp", "depth_m", "DO", "POC (ug/L)", "chl (ug/l)",
              "PON (ug/L)", "SPM (ug/L)", "FlECO-AFL(mg/m^3)", "chl depth",
              "phaeo (ug/l)", "PO4(uM)", "Sil(uM)", "NO2(uM)", "NH4(uM)",
              "N+N (umol/L)", "NO3(uM)")) {
  na_count <- sum(is.na(asgard_pcoa_df_p[[var]]))
  cat(var, ": NA =", na_count, "\n")
}
# POC+PON+SPM等: NA=27, NO3: NA=3, DO+temp+sal+FlECO: NA=1, depth: NA=0

# PCoAスコア列 (標準化しない) / PCoA score columns (do not scale)
asgard_pcoa_df_sub1_p <- asgard_pcoa_df_p %>%
  select(c("Sample", "PCoA1_Bray", "PCoA2_Bray", "PCoA1_Jaccard",
           "PCoA2_Jaccard", "PCoA1_Euclidean", "PCoA2_Euclidean"))

# 環境変数列をスケーリング (平均0, SD1) / Scale environmental variables
asgard_pcoa_df_sub2_p <- asgard_pcoa_df_p %>%
  select(c("temp", "salinity", "DO", "NO3(uM)", "FlECO-AFL(mg/m^3)")) %>%
  scale()

asgard_pcoa_df_sub_p <- cbind(asgard_pcoa_df_sub1_p, asgard_pcoa_df_sub2_p) # 78*12

# ==============================================================================
# Section 3: complete cases のみ使用 / Keep only rows with complete env data
# ==============================================================================

asgard_filtered_p_frt <- asgard_filtered_p2^.25

asgard_complete_cases_p <- complete.cases(asgard_pcoa_df_sub_p)
asgard_pcoa_cc_p <- asgard_pcoa_df_sub_p[asgard_complete_cases_p, ] # 75*12

asgard_bray_cc_p <- vegdist(
  asgard_filtered_p_frt[asgard_complete_cases_p, ],
  method = "bray"
)

# 列名の特殊文字を修正 / Fix column name with spaces for model formula
colnames(asgard_pcoa_cc_p)[colnames(asgard_pcoa_cc_p) == "chl (ug/l)"] <- "chl_ug_l"

# ==============================================================================
# Section 4: db-RDA モデル構築 / Fit db-RDA model
# ==============================================================================

asgard_dbrda_model_p <- capscale(
  asgard_bray_cc_p ~ temp + salinity + DO + `NO3(uM)` + `FlECO-AFL(mg/m^3)`,
  data = asgard_pcoa_cc_p
)

summary(asgard_dbrda_model_p)

# 全体モデルの有意性検定 / Test overall model significance
asgard_anova_dbrda_p <- anova.cca(asgard_dbrda_model_p, permutations = 999)
print(asgard_anova_dbrda_p) # expected: p < 0.001

# 環境変数ごとの寄与 (margin検定) / Per-variable marginal contribution
asgard_anova_env_p <- anova.cca(asgard_dbrda_model_p, by = "margin", permutations = 999)
print(asgard_anova_env_p)

# 各軸の有意性 / Axis-by-axis significance
asgard_anova_axes_p <- anova.cca(asgard_dbrda_model_p, by = "axis", permutations = 999)
print(asgard_anova_axes_p)

# 調整R² / Adjusted R²
RsquareAdj(asgard_dbrda_model_p)

# ==============================================================================
# Section 5: db-RDA プロット / Plot db-RDA ordination
# ==============================================================================

# サイトスコアの抽出 / Extract site scores
asgard_dbrda_scores_p <- as.data.frame(scores(asgard_dbrda_model_p, display = "sites"))
asgard_dbrda_scores_p$MDS1 <- asgard_dbrda_model_p$CA$u[, 1]
asgard_dbrda_scores_p$MDS2 <- asgard_dbrda_model_p$CA$u[, 2]
asgard_dbrda_scores_p$Sample <- rownames(asgard_dbrda_scores_p)

# バイプロット矢印の抽出 / Extract biplot arrows for environmental vectors
asgard_dbrda_vectors_p <- as.data.frame(scores(asgard_dbrda_model_p, display = "bp"))
asgard_dbrda_vectors_p$Variable <- rownames(asgard_dbrda_vectors_p)

asgard_dbrda_df_p <- merge(asgard_dbrda_scores_p, asgard_pcoa_cc_p,
                           by = "Sample", sort = FALSE)

# フィルターサイズ列を追加 (NAs 3行除く75サンプル) / Add filter size column (75 samples)
df_75_for_dbrda <- meta_asgard_p2 %>%
  filter(!Sample %in% c("BOX_6_6", "BOX_6_7", "BOX_6_26"))

filter_75 <- df_75_for_dbrda$filter
names(filter_75) <- rownames(df_75_for_dbrda)
filter_75 <- factor(filter_75, levels = c("0.2 µm", "3 µm", "20 µm"))

asgard_dbrda_merged_p <- left_join(
  asgard_dbrda_df_p,
  rownames_to_column(as.data.frame(filter_75)),
  by = c("Sample" = "rowname")
)
asgard_dbrda_merged_p$filter <- factor(
  asgard_dbrda_merged_p$filter,
  levels = c("0.2 µm", "3 µm", "20 µm")
)

# 75サンプル用クラスターベクトル / Cluster vector for 75 samples (excl. 3 with NA env data)
clusnum_p_db <- clusnum_p[!names(clusnum_p) %in% c("BOX_6_6", "BOX_6_7", "BOX_6_26")]
clusnum_p_db <- factor(clusnum_p_db, levels = c("1", "2", "3", "4")) # length 75

# ggplot2 で db-RDA を描画 / Plot db-RDA in ggplot2
ggplot() +
  geom_point(data = asgard_dbrda_merged_p,
             aes(x = MDS1, y = MDS2, color = clusnum_p_db, size = `NO3(uM)`)) +
  geom_segment(data = asgard_dbrda_vectors_p,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = asgard_dbrda_vectors_p,
            aes(x = CAP1, y = CAP2, label = Variable),
            vjust = -0.5, hjust = 0.5, size = 5) +
  labs(x = "MDS1", y = "MDS2", color = "cluster") +
  theme_minimal()

# ==============================================================================
# Section 6: 環境変数の散布図 (塩分 vs 温度 / NO3)
# Environmental scatter plots coloured by cluster
# ==============================================================================

rsc_p       <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb3 <- rsc_p[clusnum_p]

a_env <- meta_asgard_p2
a_env$cluster <- clusnum_p

plot(x   = a_env$salinity,
     y   = a_env$temp,
     col = sample_rgb3,
     pch = 19,
     cex = 1.5,
     xlab = "Salinity",
     ylab = "Temperature (°C)")

plot(x   = a_env$salinity,
     y   = a_env$`NO3(uM)`,
     col = sample_rgb3,
     pch = 19,
     cex = 1.5,
     xlab = "Salinity",
     ylab = "NO3 (µM)")

message("06_dbrda.R: done. asgard_dbrda_merged_p (75 samples) and ANOVA results ready.")
