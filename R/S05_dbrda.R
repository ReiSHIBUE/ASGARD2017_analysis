### S05_dbrda.R
### ASGARD 2017 Survey Site Analysis — db-RDA & Environmental Scatter Plots
### 距離ベースRDAスクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02, S04):
###   asgard_frtprop     - fourth-root proportion matrix 181×258
###   asgard_pcoa_df     - PCoA + metadata df (181×52+)
###   meta_asgard        - metadata 181×45
###   clusnum            - cluster assignments (length 181)
###   rsc                - 5-colour palette vector
###
### PRODUCES:
###   asgard_dbrda_model - capscale model object
###   asgard_anova_dbrda - ANOVA test of overall model
###   asgard_anova_env   - ANOVA test per environmental variable (by margin)
###   asgard_anova_axes  - ANOVA test per db-RDA axis
###   asgard_pcoa_cc     - complete-case PCoA + env df (used in capscale)
###   asgard_dbrda_merged - merged scores + metadata df
###
### OUTPUT:
###   output/survey/dbrda/ASGARD_dbrda_survey.pdf
###
### NOTE: envfit code from original script (line 1344) referenced asgard_pcoa_bray_cc
###   which was never defined — that block is omitted here. Use capscale biplot
###   arrows for equivalent visualisation.

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
# Section 2: 環境変数の選択とスケーリング / Select and scale environmental variables
# temp, salinity, DO, NO3, FlECO-AFL を選択 (VIF < 10 を確認済み)
# ==============================================================================

asgard_pcoa_df_sub1 <- asgard_pcoa_df %>%
  select(c("Sample", "PCoA1_Bray", "PCoA2_Bray", "PCoA1_Jaccard",
           "PCoA2_Jaccard", "PCoA1_Euclidean", "PCoA2_Euclidean"))

asgard_pcoa_df_sub2 <- asgard_pcoa_df %>%
  select(c("temp", "salinity", "DO", "NO3(uM)", "FlECO-AFL(mg/m^3)")) %>%
  scale()

asgard_pcoa_df_sub <- cbind(asgard_pcoa_df_sub1, asgard_pcoa_df_sub2) # 181×12

# ==============================================================================
# Section 3: Complete cases のみ使用 / Keep rows with complete environmental data
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
# Section 4: db-RDA モデル構築 / Fit db-RDA model
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
# Section 5: db-RDA スコア抽出とマージ / Extract scores and merge with metadata
# ==============================================================================

asgard_dbrda_scores <- as.data.frame(scores(asgard_dbrda_model, display = "sites"))
asgard_dbrda_scores$MDS1   <- asgard_dbrda_model$CA$u[, 1]
asgard_dbrda_scores$MDS2   <- asgard_dbrda_model$CA$u[, 2]
asgard_dbrda_scores$Sample <- rownames(asgard_dbrda_scores)

asgard_dbrda_vectors <- as.data.frame(scores(asgard_dbrda_model, display = "bp"))
asgard_dbrda_vectors$Variable <- rownames(asgard_dbrda_vectors)

asgard_dbrda_df <- merge(asgard_dbrda_scores, asgard_pcoa_cc, by = "Sample", sort = FALSE)

asgard_dbrda_merged <- left_join(
  asgard_dbrda_df,
  rownames_to_column(as.data.frame(clusnum)),
  by = c("Sample" = "rowname")
)
asgard_dbrda_merged$clusnum <- as.factor(asgard_dbrda_merged$clusnum)

# ==============================================================================
# Section 6: db-RDA プロットと環境変数散布図 / Plot db-RDA and environmental scatter
# ==============================================================================

pdf(file = here::here("output", "survey", "dbrda", "ASGARD_dbrda_survey.pdf"),
    width = 10, height = 8)

print(ggplot() +
  geom_point(data = asgard_dbrda_merged,
             aes(x = MDS1, y = MDS2, color = clusnum, size = `NO3(uM)`)) +
  geom_segment(data = asgard_dbrda_vectors,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 1) +
  geom_text(data = asgard_dbrda_vectors,
            aes(x = CAP1, y = CAP2, label = Variable),
            vjust = -0.5, hjust = 0.5, size = 5) +
  labs(x = "MDS1", y = "MDS2", color = "cluster") +
  theme_minimal())

# 環境変数散布図 / Environmental scatter plots coloured by cluster
clusnum_cc <- factor(
  clusnum[names(clusnum) %in% asgard_pcoa_cc$Sample],
  levels = as.character(1:5)
)
col_cc <- hue_pal()(5)[as.integer(clusnum_cc)]

plot(x    = asgard_pcoa_cc$salinity,
     y    = asgard_pcoa_cc$temp,
     col  = col_cc,
     pch  = 19, cex = 1.5,
     xlab = "Salinity", ylab = "Temperature (°C)",
     main = "Salinity vs. Temperature — Survey Samples")

plot(x    = asgard_pcoa_cc$salinity,
     y    = asgard_pcoa_cc$`NO3(uM)`,
     col  = col_cc,
     pch  = 19, cex = 1.5,
     xlab = "Salinity", ylab = "NO3 (µM)",
     main = "Salinity vs. NO3 — Survey Samples")

dev.off()

message("S05_dbrda.R: done. asgard_dbrda_merged (", nrow(asgard_dbrda_merged), " samples) and ANOVA results ready.")
