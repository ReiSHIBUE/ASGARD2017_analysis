### S11_crosstable.R
### ASGARD 2017 Survey — Cluster x Water Mass Association (11 clusters)
### クラスターと水塊の関連解析（11クラスター）
###
### REQUIRES (from S01, S02):
###   meta_asgard      - metadata 181 samples
###   clusnum11        - 11-cluster assignments (factor, from S02)
###   hier_levels_11   - ordered cluster names (from S02)
###
### OUTPUT:
###   output/survey/crosstable/cluster_watermass_crosstab_11clusters.csv
###   output/survey/crosstable/cluster_watermass_proportions_11clusters.csv
###   output/survey/crosstable/watermass_sample_counts.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Water mass classification / 水塊分類
# Danielson et al. (2020), Table 2 — mutually exclusive
# ==============================================================================

classify_watermass <- function(temp, sal) {
  wm <- rep(NA_character_, length(temp))
  for (i in seq_along(temp)) {
    s <- sal[i]; t <- temp[i]
    if (is.na(s) | is.na(t))                              { wm[i] <- "NA"; next }
    if (s < 30.8  & t >= 3)                               { wm[i] <- "wCW"
    } else if (s < 30.8  & t < 3)                         { wm[i] <- "IMW & cCW"
    } else if (s >= 30.8 & s < 32.5 & t >= 0 & t < 3)    { wm[i] <- "cSW"
    } else if (s >= 30.8 & s < 33.4 & t >= 3)             { wm[i] <- "wSW"
    } else if (s >= 32.5 & s < 33.8 & t >= 0 & t < 3)    { wm[i] <- "AnW"
    } else if (s >= 30.8 & s < 33.8 & t >= -1 & t < 0)   { wm[i] <- "MWW"
    } else if (s >= 30.8 & s < 33.8 & t < -1)             { wm[i] <- "WW"
    } else if ((s >= 33.4 & t >= 3) | (s >= 33.8 & t < 3)) { wm[i] <- "AtlW & BBW"
    } else { wm[i] <- "Unclassified" }
  }
  wm
}

# ==============================================================================
# Section 2: Assign water mass to each sample / 各サンプルに水塊を割り当て
# ==============================================================================

df <- meta_asgard
df$cluster11 <- factor(as.character(clusnum11[rownames(df)]), levels = hier_levels_11)
df$water_mass <- classify_watermass(df$temp, df$salinity)

# ==============================================================================
# Section 3: Cross-tabulation / クロス集計
# ==============================================================================

ct <- table(df$cluster11, df$water_mass)

cat("=== Cross-tabulation: Cluster x Water Mass (11 clusters) ===\n\n")
print(ct)

# Proportions (% of each cluster)
ct_pct <- round(prop.table(ct, margin = 1) * 100, 1)

cat("\n=== Proportions (% of each cluster in each water mass) ===\n\n")
print(ct_pct)

# Dominant water mass per cluster
cat("\n=== Dominant water mass per cluster ===\n\n")
for (cl in hier_levels_11) {
  row <- ct_pct[cl, ]
  row <- sort(row[row > 0], decreasing = TRUE)
  cat(sprintf("%-5s: %s\n", cl,
    paste(paste0(names(row), " (", row, "%)"), collapse = ", ")))
}

# ==============================================================================
# Section 4: Statistical tests / 統計検定
# ==============================================================================

# Chi-squared test (simulated p-value for sparse table)
cat("\n=== Chi-squared test ===\n")
chi <- chisq.test(ct, simulate.p.value = TRUE, B = 9999)
print(chi)

# Fisher exact test (simulated)
cat("\n=== Fisher exact test ===\n")
fisher <- fisher.test(ct, simulate.p.value = TRUE, B = 9999)
print(fisher)

# ==============================================================================
# Section 5: Save results / 結果を保存
# ==============================================================================

dir.create(here("output", "survey", "crosstable"), showWarnings = FALSE, recursive = TRUE)

# Cross-tabulation (counts)
ct_df <- as.data.frame.matrix(ct)
ct_df$cluster <- rownames(ct_df)
ct_df <- ct_df %>% select(cluster, everything())
write.csv(ct_df,
  here("output", "survey", "crosstable", "cluster_watermass_crosstab_11clusters.csv"),
  row.names = FALSE)

# Proportions
ct_pct_df <- as.data.frame.matrix(ct_pct)
ct_pct_df$cluster <- rownames(ct_pct_df)
ct_pct_df <- ct_pct_df %>% select(cluster, everything())
write.csv(ct_pct_df,
  here("output", "survey", "crosstable", "cluster_watermass_proportions_11clusters.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 6: Water mass sample counts (total, regardless of cluster)
# ==============================================================================

wm_total <- table(df$water_mass)
wm_total_df <- data.frame(
  water_mass = names(wm_total),
  n_samples = as.integer(wm_total),
  pct = round(as.numeric(wm_total) / sum(wm_total) * 100, 1)
) %>% arrange(desc(n_samples))

cat("\n=== Water mass sample counts (all 181 samples) ===\n\n")
print(wm_total_df)

# Statistical test results
test_results <- data.frame(
  test = c("Chi-squared (Monte Carlo, B=9999)", "Fisher exact (Monte Carlo, B=9999)"),
  statistic = c(round(chi$statistic, 2), NA),
  p_value = c(chi$p.value, fisher$p.value),
  sig = c(
    ifelse(chi$p.value <= 0.001, "***", ifelse(chi$p.value <= 0.01, "**", ifelse(chi$p.value <= 0.05, "*", "ns"))),
    ifelse(fisher$p.value <= 0.001, "***", ifelse(fisher$p.value <= 0.01, "**", ifelse(fisher$p.value <= 0.05, "*", "ns")))
  )
)

write.csv(test_results,
  here("output", "survey", "crosstable", "cluster_watermass_test_results.csv"),
  row.names = FALSE)

write.csv(wm_total_df,
  here("output", "survey", "crosstable", "watermass_sample_counts.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 7: Water mass position plot / 水塊ポジション図
# Coastal vs Shelf (Y), Warm vs Cool (X)
# ==============================================================================

library(ggrepel)

cluster_wm <- df %>%
  filter(!is.na(water_mass)) %>%
  group_by(cluster11) %>%
  summarise(
    n = n(),
    pct_wSW = sum(water_mass == "wSW") / n() * 100,
    pct_cSW = sum(water_mass == "cSW") / n() * 100,
    pct_AnW = sum(water_mass == "AnW") / n() * 100,
    pct_wCW = sum(water_mass == "wCW") / n() * 100,
    pct_IMW = sum(water_mass == "IMW & cCW") / n() * 100,
    pct_MWW = sum(water_mass == "MWW") / n() * 100,
    pct_WW  = sum(water_mass == "WW") / n() * 100,
    .groups = "drop"
  ) %>%
  mutate(
    coastal_index = pct_wCW + pct_IMW,
    cool_index    = pct_cSW + pct_AnW + pct_MWW + pct_WW
  )

pdf(here("output", "survey", "crosstable", "cluster_watermass_position.pdf"),
    width = 10, height = 8)

print(
  ggplot(cluster_wm, aes(x = cool_index, y = coastal_index)) +
    geom_point(aes(color = cluster11, size = n), alpha = 0.8) +
    geom_text_repel(aes(label = cluster11, color = cluster11),
                    size = 5, fontface = "bold",
                    max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = cc11, guide = "none") +
    scale_size_continuous(name = "# of samples", range = c(3, 10)) +
    scale_x_continuous(limits = c(-5, 105), breaks = seq(0, 100, 25)) +
    scale_y_continuous(limits = c(-5, 75),  breaks = seq(0, 100, 25)) +
    labs(x = "Cool water mass proportion (%)\n(cSW + AnW + MWW + WW)",
         y = "Coastal water mass proportion (%)\n(wCW + IMW & cCW)") +
    annotate("text", x = 0,  y = -3, label = "Warm Shelf",    hjust = 0,
             fontface = "italic", size = 4, color = "gray40") +
    annotate("text", x = 95, y = -3, label = "Cool Shelf",    hjust = 1,
             fontface = "italic", size = 4, color = "gray40") +
    annotate("text", x = 0,  y = 72, label = "Warm Coastal",  hjust = 0,
             fontface = "italic", size = 4, color = "gray40") +
    annotate("text", x = 95, y = 72, label = "Cool Coastal",  hjust = 1,
             fontface = "italic", size = 4, color = "gray40") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      panel.grid.minor = element_blank()
    )
)

dev.off()

message("\nS11_crosstable.R: done.")
message("  Chi-squared p-value: ", chi$p.value)
message("  Fisher p-value: ", fisher$p.value)
message("  CSV: cluster_watermass_crosstab_11clusters.csv")
message("  CSV: cluster_watermass_proportions_11clusters.csv")
message("  CSV: cluster_watermass_test_results.csv")
message("  CSV: watermass_sample_counts.csv")
message("  PDF: cluster_watermass_position.pdf")
