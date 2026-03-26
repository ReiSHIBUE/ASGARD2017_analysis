### cluster_watermass_association.R
### ASGARD 2017 Survey — Cluster x Water Mass Association Analysis
### クラスターと水塊の関連解析
###
### REQUIRES (from S01, S02):
###   meta_asgard   - metadata 181 samples
###   clusnum10     - 10-cluster assignments (from S02)
###
### OUTPUT:
###   output/survey/cluster_watermass_crosstab.csv
###   output/survey/cluster_watermass_proportions.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Water mass classification / Classification des masses d'eau
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
df$cluster <- factor(clusnum10[rownames(df)], levels = as.character(1:10))
df$water_mass <- classify_watermass(df$temp, df$salinity)

# ==============================================================================
# Section 3: Cross-tabulation / クロス集計
# ==============================================================================

ct <- table(df$cluster, df$water_mass)

cat("=== Cross-tabulation: Cluster x Water Mass ===\n\n")
print(ct)

# Proportions (% of each cluster)
ct_pct <- round(prop.table(ct, margin = 1) * 100, 1)

cat("\n=== Proportions (% of each cluster in each water mass) ===\n\n")
print(ct_pct)

# Dominant water mass per cluster
cat("\n=== Dominant water mass per cluster ===\n\n")
for (cl in as.character(1:10)) {
  row <- ct_pct[cl, ]
  row <- sort(row[row > 0], decreasing = TRUE)
  cat(sprintf("Cluster %2s: %s\n", cl,
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

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)

# Cross-tabulation (counts)
ct_df <- as.data.frame.matrix(ct)
ct_df$cluster <- rownames(ct_df)
ct_df <- ct_df %>% select(cluster, everything())
write.csv(ct_df, here("output", "survey", "cluster_watermass_crosstab.csv"), row.names = FALSE)

# Proportions
ct_pct_df <- as.data.frame.matrix(ct_pct)
ct_pct_df$cluster <- rownames(ct_pct_df)
ct_pct_df <- ct_pct_df %>% select(cluster, everything())
write.csv(ct_pct_df, here("output", "survey", "cluster_watermass_proportions.csv"), row.names = FALSE)

message("cluster_watermass_association.R: done.")
message("  Chi-squared p-value: ", chi$p.value)
message("  Fisher p-value: ", fisher$p.value)
message("  CSV: cluster_watermass_crosstab.csv, cluster_watermass_proportions.csv")
