### S10_permanova.R
### ASGARD 2017 Survey Site Analysis — PERMANOVA, PERMDISP, Mantel Test
### PERMANOVA・PERMDISP・マンテル検定スクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02, S04):
###   asgard_braymat   - Bray-Curtis distance matrix (dist object)
###   asgard_eucmat    - Euclidean distance matrix
###   asgard_jacmat    - Jaccard distance matrix
###   asgard_pcoa_df   - PCoA + metadata df (181×52+)
###   meta_asgard      - metadata 181×45
###   clusnum          - cluster assignments (length 181)
###
### PRODUCES:
###   asgard_permanova - adonis2() result for overall model
###   asgard_permdisp  - betadisper() result by cluster
###   asgard_mantel_*  - mantel() results per environmental variable
###
### OUTPUT:
###   output/survey/beta_diversity/ASGARD_permanova_results.txt

library(vegan)
library(tidyverse)

# ==============================================================================
# Section 1: PERMANOVA — 環境変数ごとの寄与 / PERMANOVA per environmental variable
# adonis2()で距離行列~変数を検定 / Test each env variable as predictor
# ==============================================================================

message("\n--- PERMANOVA (adonis2) results ---")

for (gvar in c("temp", "salinity", "DO")) {
  res <- adonis2(asgard_braymat ~ asgard_pcoa_df[[gvar]], permutations = 999)
  message(gvar, ": R2 = ", round(res$R2[1], 4), ", p = ", res$`Pr(>F)`[1])
}

# NO3 はNAを含む行を除外 / NO3: must exclude rows with NA
no3_ok  <- !is.na(asgard_pcoa_df$`NO3(uM)`)
bray_no3 <- vegdist(asgard_frtprop[no3_ok, ], method = "bray")

res_no3 <- adonis2(bray_no3 ~ asgard_pcoa_df$`NO3(uM)`[no3_ok], permutations = 999)
message("NO3(uM): R2 = ", round(res_no3$R2[1], 4), ", p = ", res_no3$`Pr(>F)`[1])

# クラスター毎のPERMANOVA / PERMANOVA by cluster assignment
asgard_permanova <- adonis2(
  asgard_braymat ~ factor(clusnum),
  permutations = 999
)
print(asgard_permanova)

# ==============================================================================
# Section 2: PERMDISP — グループ内分散の均一性検定
# Test homogeneity of dispersion within clusters
# ==============================================================================

message("\n--- PERMDISP results ---")

asgard_permdisp      <- betadisper(asgard_braymat, group = factor(clusnum))
asgard_permdisp_anova <- anova(asgard_permdisp)
print(asgard_permdisp_anova)

# ==============================================================================
# Section 3: Mantel test — 距離行列間の相関
# Test correlation between Bray-Curtis and environmental distance matrices
# ==============================================================================

message("\n--- Mantel test results ---")

asgard_salinity_dist <- dist(asgard_pcoa_df$salinity)
asgard_temp_dist     <- dist(asgard_pcoa_df$temp)
asgard_do_dist       <- dist(asgard_pcoa_df$DO)
asgard_no3_dist      <- dist(asgard_pcoa_df$`NO3(uM)`[no3_ok])

asgard_mantel_sal  <- mantel(asgard_braymat,   asgard_salinity_dist, permutations = 999)
asgard_mantel_temp <- mantel(asgard_braymat,   asgard_temp_dist,     permutations = 999)
asgard_mantel_do   <- mantel(asgard_braymat,   asgard_do_dist,       permutations = 999)
asgard_mantel_no3  <- mantel(bray_no3,         asgard_no3_dist,      permutations = 999)

for (nm in c("sal", "temp", "do", "no3")) {
  obj <- get(paste0("asgard_mantel_", nm))
  message(nm, ": r = ", round(obj$statistic, 4), ", p = ", round(obj$signif, 4))
}

# ==============================================================================
# Section 4: 結果をテキストファイルに書き出す / Write results to text file
# ==============================================================================

sink(here::here("output", "survey", "beta_diversity", "ASGARD_permanova_results.txt"))
cat("=== PERMANOVA by cluster (adonis2) ===\n")
print(asgard_permanova)
cat("\n=== PERMDISP (betadisper ANOVA) ===\n")
print(asgard_permdisp_anova)
cat("\n=== Mantel tests ===\n")
for (nm in c("sal", "temp", "do", "no3")) {
  obj <- get(paste0("asgard_mantel_", nm))
  cat(nm, ": r =", round(obj$statistic, 4), ", p =", round(obj$signif, 4), "\n")
}
sink()

message("S10_permanova.R: done. Results in output/survey/beta_diversity/ASGARD_permanova_results.txt")
