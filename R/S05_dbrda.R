### S05_dbrda.R
### ASGARD 2017 Survey Site Analysis — dbRDA & Environmental Scatter Plots (11 clusters)
### 距離ベースRDAスクリプト（サーベイサイト、11クラスター）
###
### REQUIRES (from S01, S02, S04):
###   asgard_frtprop     - fourth-root proportion matrix 181×258
###   asgard_pcoa_df     - PCoA + metadata df (181×52+)
###   meta_asgard        - metadata 181×45
###   clusnum11          - 11-cluster assignments (factor, from S02)
###   hier_levels_11     - ordered cluster names (from S02)
###   cc11               - 11-colour palette (from S02)
###
### PRODUCES:
###   asgard_dbrda_model - capscale model object
###   asgard_anova_dbrda - ANOVA test of overall model
###   asgard_anova_env   - ANOVA test per environmental variable (by margin)
###   asgard_anova_axes  - ANOVA test per dbRDA axis
###   asgard_pcoa_cc     - complete-case PCoA + env df (used in capscale)
###   asgard_dbrda_merged - merged scores + metadata df
###
### OUTPUT:
###   output/survey/dbrda/ASGARD_dbrda_survey_11clusters.pdf

library(vegan)
library(tidyverse)

# ==============================================================================
# Section 1: NAの確認 / Check NA counts per environmental variable
# ==============================================================================

for (var in c("salinity", "temp", "depth_m", "DO", "POC (ug/L)", "chl (ug/l)",
              "PON (ug/L)", "SPM (ug/L)", "FlECO-AFL(mg/m^3)", "chl depth",
              "phaeo (ug/l)", "PO4(uM)", "Sil(uM)", "NO2(uM)", "NH4(uM)",
              "N+N (umol/L)", "NO3(uM)")) {
  cat(var, ": NA =", sum(is.na(asgard_pcoa_df[[var]])), "\n")
}

# ==============================================================================
# Section 2: VIF確認 / Check VIF for variable selection
# ==============================================================================

# All 9 candidate variables
all_env_vars <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)", "Sil(uM)",
                   "NH4(uM)", "FlECO-AFL(mg/m^3)", "depth_m")
env_all <- asgard_pcoa_df[, all_env_vars]
cc_all <- complete.cases(env_all)
env_all_cc <- scale(env_all[cc_all, ])
bray_all_cc <- vegdist(asgard_frtprop[cc_all, ], method = "bray")

# VIF with all 9 variables
model_all9 <- capscale(bray_all_cc ~ ., data = as.data.frame(env_all_cc))
vif_all9 <- vif.cca(model_all9)
message("\n--- VIF with all 9 variables ---")
print(round(vif_all9, 2))

# VIF with final 5 variables (PO4, Sil, NH4, depth removed due to VIF > 10 or redundancy)
env_5 <- scale(asgard_pcoa_df[cc_all, c("temp", "salinity", "DO", "NO3(uM)", "FlECO-AFL(mg/m^3)")])
model_5 <- capscale(bray_all_cc ~ ., data = as.data.frame(env_5))
vif_5 <- vif.cca(model_5)
message("\n--- VIF with final 5 variables ---")
print(round(vif_5, 2))

# ==============================================================================
# Section 3: 環境変数の選択とスケーリング / Select and scale environmental variables
# temp, salinity, DO, NO3, FlECO-AFL (all VIF < 3)
# ==============================================================================

asgard_pcoa_df_sub1 <- asgard_pcoa_df %>%
  select(c("Sample", "PCoA1_Bray", "PCoA2_Bray", "PCoA1_Jaccard",
           "PCoA2_Jaccard", "PCoA1_Euclidean", "PCoA2_Euclidean"))

asgard_pcoa_df_sub2 <- asgard_pcoa_df %>%
  select(c("temp", "salinity", "DO", "NO3(uM)", "FlECO-AFL(mg/m^3)")) %>%
  scale()

asgard_pcoa_df_sub <- cbind(asgard_pcoa_df_sub1, asgard_pcoa_df_sub2) # 181×12

# ==============================================================================
# Section 4: Complete cases のみ使用 / Keep rows with complete environmental data
# ==============================================================================

asgard_complete_cases <- complete.cases(asgard_pcoa_df_sub) # logical length 181
asgard_pcoa_cc        <- asgard_pcoa_df_sub[asgard_complete_cases, ]

asgard_bray_cc <- vegdist(
  asgard_frtprop[asgard_complete_cases, ],
  method = "bray"
)

# 列名の特殊文字を修正 / Fix column names with spaces for model formula
colnames(asgard_pcoa_cc)[colnames(asgard_pcoa_cc) == "chl (ug/l)"] <- "chl_ug_l"

# ==============================================================================
# Section 5: dbRDA モデル構築 / Fit dbRDA model
# ==============================================================================

asgard_dbrda_model <- capscale(
  asgard_bray_cc ~ temp + salinity + DO + `NO3(uM)` + `FlECO-AFL(mg/m^3)`,
  data = asgard_pcoa_cc
)

summary(asgard_dbrda_model)

asgard_anova_dbrda <- anova.cca(asgard_dbrda_model, permutations = 999)
print(asgard_anova_dbrda)

asgard_anova_env <- anova.cca(asgard_dbrda_model, by = "margin", permutations = 999)
print(asgard_anova_env)

asgard_anova_axes <- anova.cca(asgard_dbrda_model, by = "axis", permutations = 999)
print(asgard_anova_axes)

RsquareAdj(asgard_dbrda_model)

# ==============================================================================
# Section 6: dbRDA スコア抽出とマージ / Extract scores and merge with metadata
# ==============================================================================

asgard_dbrda_scores <- as.data.frame(scores(asgard_dbrda_model, display = "sites"))
asgard_dbrda_scores$MDS1   <- asgard_dbrda_model$CA$u[, 1]
asgard_dbrda_scores$MDS2   <- asgard_dbrda_model$CA$u[, 2]
asgard_dbrda_scores$Sample <- rownames(asgard_dbrda_scores)

asgard_dbrda_vectors <- as.data.frame(scores(asgard_dbrda_model, display = "bp"))
asgard_dbrda_vectors$Variable <- rownames(asgard_dbrda_vectors)
# Shorter labels for plot
asgard_dbrda_vectors$label <- c("Temp", "Salinity", "DO", "NO3", "FlECO")

asgard_dbrda_df <- merge(asgard_dbrda_scores, asgard_pcoa_cc, by = "Sample", sort = FALSE)

asgard_dbrda_merged <- left_join(
  asgard_dbrda_df,
  data.frame(Sample = names(clusnum11),
             cluster11 = factor(as.character(clusnum11), levels = hier_levels_11)),
  by = "Sample"
)

# Add division for faceting
asgard_dbrda_merged$division <- factor(
  sub("^(.).*", "\\1", as.character(asgard_dbrda_merged$cluster11)),
  levels = c("A", "B", "C")
)

# ==============================================================================
# Section 7: dbRDA プロット / Plot dbRDA (11 clusters)
# ==============================================================================

dir.create(here::here("output", "survey", "dbrda"), showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "dbrda", "ASGARD_dbrda_survey_11clusters.pdf"),
    width = 10, height = 8)

# Calculate % variance explained per axis from summary
sm <- summary(asgard_dbrda_model)
prop_expl <- sm$cont$importance["Proportion Explained", ]
cap1_pct <- round(prop_expl["CAP1"] * 100, 1)
cap2_pct <- round(prop_expl["CAP2"] * 100, 1)
mds1_pct <- round(prop_expl["MDS1"] * 100, 1)
mds2_pct <- round(prop_expl["MDS2"] * 100, 1)

# Page 1: CAP axes (constrained) — all clusters
cap1_range <- range(asgard_dbrda_merged$CAP1) * 1.1
cap2_range <- range(asgard_dbrda_merged$CAP2) * 1.1

print(ggplot() +
  geom_point(data = asgard_dbrda_merged,
             aes(x = CAP1, y = CAP2, color = cluster11, size = `NO3(uM)`),
             alpha = 0.7) +
  scale_color_manual(values = cc11, name = "Cluster") +
  scale_size_continuous(range = c(1, 6), name = "NO3 (uM)") +
  geom_segment(data = asgard_dbrda_vectors,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = asgard_dbrda_vectors,
            aes(x = CAP1 * 1.15, y = CAP2 * 1.15, label = label),
            size = 5, fontface = "bold") +
  coord_cartesian(xlim = cap1_range, ylim = cap2_range) +
  labs(x = paste0("CAP1 (", cap1_pct, "%)"),
       y = paste0("CAP2 (", cap2_pct, "%)"),
       title = "dbRDA CAP axes (11 clusters)") +
  theme_grey() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 2: CAP axes faceted by Division
print(ggplot() +
  geom_point(data = asgard_dbrda_merged,
             aes(x = CAP1, y = CAP2, color = cluster11, size = `NO3(uM)`),
             alpha = 0.7) +
  scale_color_manual(values = cc11, name = "Cluster") +
  scale_size_continuous(range = c(1, 6), name = "NO3 (uM)") +
  facet_grid(~ division) +
  coord_cartesian(xlim = cap1_range, ylim = cap2_range) +
  labs(x = paste0("CAP1 (", cap1_pct, "%)"),
       y = paste0("CAP2 (", cap2_pct, "%)"),
       title = "dbRDA CAP axes (by Division A/B/C)") +
  theme_grey() +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 14)))

# Page 3: MDS axes (unconstrained)
mds1_range <- range(asgard_dbrda_merged$MDS1) * 1.1
mds2_range <- range(asgard_dbrda_merged$MDS2) * 1.1

print(ggplot() +
  geom_point(data = asgard_dbrda_merged,
             aes(x = MDS1, y = MDS2, color = cluster11, size = `NO3(uM)`),
             alpha = 0.7) +
  scale_color_manual(values = cc11, name = "Cluster") +
  scale_size_continuous(range = c(1, 6), name = "NO3 (uM)") +
  coord_cartesian(xlim = mds1_range, ylim = mds2_range) +
  labs(x = paste0("MDS1 (", mds1_pct, "%)"),
       y = paste0("MDS2 (", mds2_pct, "%)"),
       title = "dbRDA MDS axes (11 clusters)") +
  theme_grey() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 4: MDS axes faceted by Division
print(ggplot() +
  geom_point(data = asgard_dbrda_merged,
             aes(x = MDS1, y = MDS2, color = cluster11, size = `NO3(uM)`),
             alpha = 0.7) +
  scale_color_manual(values = cc11, name = "Cluster") +
  scale_size_continuous(range = c(1, 6), name = "NO3 (uM)") +
  facet_grid(~ division) +
  coord_cartesian(xlim = mds1_range, ylim = mds2_range) +
  labs(x = paste0("MDS1 (", mds1_pct, "%)"),
       y = paste0("MDS2 (", mds2_pct, "%)"),
       title = "dbRDA MDS axes (by Division A/B/C)") +
  theme_grey() +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 14)))

dev.off()

message("S05_dbrda.R: done. asgard_dbrda_merged (", nrow(asgard_dbrda_merged), " samples) and ANOVA results ready.")
message("  PDF: output/survey/dbrda/ASGARD_dbrda_survey_11clusters.pdf")
