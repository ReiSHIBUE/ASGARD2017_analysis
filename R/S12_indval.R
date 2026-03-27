### S12_indval.R
### ASGARD 2017 Survey — Indicator Species Analysis (IndVal)
### 指標種解析（IndVal）
###
### REQUIRES (from S01, S02):
###   asgard_filtered  - 181 x 258 ASV proportion matrix
###   clusnum10        - 10-cluster assignments
###
### OUTPUT:
###   output/survey/indval_results.csv

library(tidyverse)
library(indicspecies)
library(here)

# ==============================================================================
# Section 1: IndVal analysis / IndVal解析
# IndVal.g: group-equalized indicator value
# Accounts for unequal group sizes
# ==============================================================================

cluster_vec <- clusnum10[rownames(asgard_filtered)]

set.seed(42)
indval_result <- multipatt(
  asgard_filtered,
  cluster_vec,
  func = "IndVal.g",
  control = how(nperm = 999)
)

# ==============================================================================
# Section 2: Extract significant results / 有意な結果を抽出
# ==============================================================================

sig <- indval_result$sign
sig$ASV <- rownames(sig)
sig_only <- sig %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

# Identify which cluster(s) each ESV is indicator for
s_cols <- grep("^s[.]", colnames(sig_only), value = TRUE)
sig_only$indicator_cluster <- apply(sig_only[, s_cols], 1, function(x) {
  paste(which(x == 1), collapse = "+")
})

# Single-cluster indicators only
single <- sig_only %>% filter(!grepl("[+]", indicator_cluster))

message("Total significant (p < 0.05): ", nrow(sig_only))
message("Single-cluster indicators: ", nrow(single))

# ==============================================================================
# Section 3: Print top indicators per cluster / 各クラスターの上位指標種を表示
# ==============================================================================

for (cl in as.character(1:10)) {
  cl_ind <- single %>% filter(indicator_cluster == cl) %>% head(10)
  n_total <- sum(single$indicator_cluster == cl)
  message(sprintf("\nCluster %s (%d indicators):", cl, n_total))
  if (nrow(cl_ind) > 0) {
    for (j in seq_len(nrow(cl_ind))) {
      message(sprintf("  %-65s stat=%.3f p=%.3f",
        cl_ind$ASV[j], cl_ind$stat[j], cl_ind$p.value[j]))
    }
  } else {
    message("  (none)")
  }
}

# ==============================================================================
# Section 4: Save results / 結果を保存
# ==============================================================================

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)
write.csv(sig_only, here("output", "survey", "indval_results.csv"), row.names = FALSE)

message("\nS12_indval.R: done.")
message("  CSV: output/survey/indval_results.csv")
