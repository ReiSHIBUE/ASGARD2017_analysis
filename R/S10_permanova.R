### S10_permanova.R
### ASGARD 2017 Survey Site Analysis — PERMANOVA, PERMDISP, Mantel Test (10 clusters)
### PERMANOVA・PERMDISP・マンテル検定スクリプト（サーベイサイト、10クラスター）
###
### REQUIRES (from S01, S02, S04):
###   asgard_braymat   - Bray-Curtis distance matrix (dist object)
###   asgard_eucmat    - Euclidean distance matrix
###   asgard_jacmat    - Jaccard distance matrix
###   asgard_frtprop   - fourth-root proportion matrix 181x258
###   asgard_pcoa_df   - PCoA + metadata df (181x53+)
###   meta_asgard      - metadata 181x45
###   clusnum10        - 10-cluster assignments (from S02)
###
### PRODUCES:
###   asgard_permanova - adonis2() result for 10-cluster model
###   asgard_permdisp  - betadisper() result by 10 clusters
###   asgard_mantel_*  - mantel() results per environmental variable
###
### OUTPUT:
###   output/survey/beta_diversity/ASGARD_permanova_results.txt

library(vegan)
library(tidyverse)

# ==============================================================================
# Section 1: PERMANOVA — 環境変数ごとの寄与 / Per environmental variable
# ==============================================================================

message("\n--- PERMANOVA (adonis2) per environmental variable ---")

for (gvar in c("temp", "salinity", "DO")) {
  res <- adonis2(asgard_braymat ~ asgard_pcoa_df[[gvar]], permutations = 999)
  message(gvar, ": R2 = ", round(res$R2[1], 4), ", p = ", res$`Pr(>F)`[1])
}

# NO3: exclude rows with NA
no3_ok   <- !is.na(asgard_pcoa_df$`NO3(uM)`)
bray_no3 <- vegdist(asgard_frtprop[no3_ok, ], method = "bray")

res_no3 <- adonis2(bray_no3 ~ asgard_pcoa_df$`NO3(uM)`[no3_ok], permutations = 999)
message("NO3(uM): R2 = ", round(res_no3$R2[1], 4), ", p = ", res_no3$`Pr(>F)`[1])

# ==============================================================================
# Section 2: PERMANOVA — 10クラスター / By 10-cluster assignment
# ==============================================================================

message("\n--- PERMANOVA (adonis2) by 10 clusters ---")

cluster10_factor <- factor(clusnum10[rownames(as.matrix(asgard_braymat))],
                           levels = as.character(1:10))

asgard_permanova <- adonis2(
  asgard_braymat ~ cluster10_factor,
  permutations = 999
)
print(asgard_permanova)

# ==============================================================================
# Section 3: Pairwise PERMANOVA — 全ペアワイズ比較 / All pairwise comparisons
# ==============================================================================

message("\n--- Pairwise PERMANOVA (10 clusters, 45 pairs) ---")

clusters <- as.character(1:10)
pairs <- combn(clusters, 2)
pairwise_results <- data.frame(
  pair = character(),
  F_value = numeric(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(pairs))) {
  cl_a <- pairs[1, i]
  cl_b <- pairs[2, i]

  samps_a <- names(clusnum10)[clusnum10 == cl_a]
  samps_b <- names(clusnum10)[clusnum10 == cl_b]
  samps_ab <- c(samps_a, samps_b)

  # Subset distance matrix
  mat_sub <- as.matrix(asgard_braymat)[samps_ab, samps_ab]
  dist_sub <- as.dist(mat_sub)
  group_sub <- factor(c(rep(cl_a, length(samps_a)), rep(cl_b, length(samps_b))))

  res_pair <- adonis2(dist_sub ~ group_sub, permutations = 999)

  pairwise_results <- rbind(pairwise_results, data.frame(
    pair = paste0(cl_a, " vs ", cl_b),
    F_value = round(res_pair$F[1], 2),
    R2 = round(res_pair$R2[1], 4),
    p_value = res_pair$`Pr(>F)`[1]
  ))
}

# BH correction for multiple comparisons
pairwise_results$p_adj <- round(p.adjust(pairwise_results$p_value, method = "BH"), 4)
pairwise_results$sig <- ifelse(pairwise_results$p_adj <= 0.001, "***",
                        ifelse(pairwise_results$p_adj <= 0.01, "**",
                        ifelse(pairwise_results$p_adj <= 0.05, "*", "ns")))

print(pairwise_results, row.names = FALSE)

# ==============================================================================
# Section 4: PERMDISP — グループ内分散の均一性検定 / Homogeneity of dispersion
# ==============================================================================

message("\n--- PERMDISP results (10 clusters) ---")

asgard_permdisp       <- betadisper(asgard_braymat, group = cluster10_factor)
asgard_permdisp_anova <- anova(asgard_permdisp)
print(asgard_permdisp_anova)

# ==============================================================================
# Section 4b: Pairwise PERMDISP — ペアワイズ分散均一性検定
# Pairwise comparison of within-group dispersion (TukeyHSD on betadisper)
# ==============================================================================

message("\n--- Pairwise PERMDISP (TukeyHSD) ---")

asgard_permdisp_tukey <- TukeyHSD(asgard_permdisp)
permdisp_pairwise <- as.data.frame(asgard_permdisp_tukey$group)
permdisp_pairwise$pair <- rownames(permdisp_pairwise)
permdisp_pairwise$sig <- ifelse(permdisp_pairwise$`p adj` <= 0.001, "***",
                          ifelse(permdisp_pairwise$`p adj` <= 0.01, "**",
                          ifelse(permdisp_pairwise$`p adj` <= 0.05, "*", "ns")))

# Reorder columns
permdisp_pairwise <- permdisp_pairwise[, c("pair", "diff", "lwr", "upr", "p adj", "sig")]
colnames(permdisp_pairwise) <- c("pair", "diff_dispersion", "CI_lower", "CI_upper", "p_adj", "sig")
permdisp_pairwise$diff_dispersion <- round(permdisp_pairwise$diff_dispersion, 4)
permdisp_pairwise$CI_lower <- round(permdisp_pairwise$CI_lower, 4)
permdisp_pairwise$CI_upper <- round(permdisp_pairwise$CI_upper, 4)
permdisp_pairwise$p_adj <- round(permdisp_pairwise$p_adj, 4)

# Count significant pairs
n_sig <- sum(permdisp_pairwise$sig != "ns")
message("Significant pairwise PERMDISP: ", n_sig, " / ", nrow(permdisp_pairwise))

print(permdisp_pairwise[permdisp_pairwise$sig != "ns", ], row.names = FALSE)

# ==============================================================================
# Section 5: Mantel test — 距離行列間の相関 / Correlation between distance matrices
# ==============================================================================

message("\n--- Mantel test results ---")

asgard_salinity_dist <- dist(asgard_pcoa_df$salinity)
asgard_temp_dist     <- dist(asgard_pcoa_df$temp)
asgard_do_dist       <- dist(asgard_pcoa_df$DO)
asgard_no3_dist      <- dist(asgard_pcoa_df$`NO3(uM)`[no3_ok])

asgard_mantel_sal  <- mantel(asgard_braymat, asgard_salinity_dist, permutations = 999)
asgard_mantel_temp <- mantel(asgard_braymat, asgard_temp_dist,     permutations = 999)
asgard_mantel_do   <- mantel(asgard_braymat, asgard_do_dist,       permutations = 999)
asgard_mantel_no3  <- mantel(bray_no3,       asgard_no3_dist,      permutations = 999)

for (nm in c("sal", "temp", "do", "no3")) {
  obj <- get(paste0("asgard_mantel_", nm))
  message(nm, ": r = ", round(obj$statistic, 4), ", p = ", round(obj$signif, 4))
}

# ==============================================================================
# Section 6: Save results / 結果をファイルに書き出す
# ==============================================================================

dir.create(here::here("output", "survey", "beta_diversity"), showWarnings = FALSE, recursive = TRUE)

sink(here::here("output", "survey", "beta_diversity", "ASGARD_permanova_results.txt"))

cat("=== PERMANOVA by 10 clusters (adonis2, 999 permutations) ===\n")
print(asgard_permanova)

cat("\n=== Pairwise PERMANOVA (45 pairs, BH-adjusted p-values) ===\n")
print(pairwise_results, row.names = FALSE)

cat("\n=== PERMDISP (betadisper ANOVA, 10 clusters) ===\n")
print(asgard_permdisp_anova)

cat("\n=== Pairwise PERMDISP (TukeyHSD on betadisper) ===\n")
cat("Significant pairs only (p_adj <= 0.05):\n")
print(permdisp_pairwise[permdisp_pairwise$sig != "ns", ], row.names = FALSE)
cat("\nNon-significant pairs: ", sum(permdisp_pairwise$sig == "ns"), " / ", nrow(permdisp_pairwise), "\n")

cat("\n=== Mantel tests (Bray-Curtis vs environmental distance) ===\n")
for (nm in c("sal", "temp", "do", "no3")) {
  obj <- get(paste0("asgard_mantel_", nm))
  cat(nm, ": r =", round(obj$statistic, 4), ", p =", round(obj$signif, 4), "\n")
}

sink()

# Save pairwise results as CSV
write.csv(pairwise_results,
  here::here("output", "survey", "beta_diversity", "pairwise_permanova_10clusters.csv"),
  row.names = FALSE)

# Save pairwise PERMDISP as CSV
write.csv(permdisp_pairwise,
  here::here("output", "survey", "beta_diversity", "pairwise_permdisp_10clusters.csv"),
  row.names = FALSE)

message("\nS10_permanova.R: done.")
message("  TXT: output/survey/beta_diversity/ASGARD_permanova_results.txt")
message("  CSV: output/survey/beta_diversity/pairwise_permanova_10clusters.csv")
message("  CSV: output/survey/beta_diversity/pairwise_permdisp_10clusters.csv")

