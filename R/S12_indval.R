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
# Section 4: 各クラスターのTop10 ASV (平均相対存在量順)
# Top 10 ASVs per cluster by mean relative abundance
# ==============================================================================

mat <- as.matrix(asgard_filtered)
cl_names <- hier_names[as.character(clusnum10[rownames(mat)])]
names(cl_names) <- rownames(mat)

top10_list <- list()
for (cl_name in hier_levels) {
  samples <- names(cl_names[cl_names == cl_name])
  mean_ra <- sort(colMeans(mat[samples, , drop = FALSE]), decreasing = TRUE)
  top10 <- data.frame(
    cluster = cl_name,
    rank = 1:10,
    ASV = names(mean_ra)[1:10],
    mean_RA = unname(mean_ra[1:10]),
    stringsAsFactors = FALSE
  )
  top10_list[[cl_name]] <- top10

  message(sprintf("\n%s (n=%d) — Top 10 ASVs by mean RA:", cl_name, length(samples)))
  for (i in 1:10) {
    parts <- strsplit(top10$ASV[i], "; ")[[1]]
    short <- paste(tail(parts, 2), collapse = "; ")
    message(sprintf("  %2d. mean RA=%.4f  %s", i, top10$mean_RA[i], short))
  }
}

top10_df <- bind_rows(top10_list)

# ==============================================================================
# Section 5: Save results / 結果を保存
# ==============================================================================

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)
write.csv(sig_only, here("output", "survey", "indval_results.csv"), row.names = FALSE)
write.csv(top10_df, here("output", "survey", "cluster_top10_ASVs.csv"), row.names = FALSE)

message("\nS12_indval.R: done.")
message("  CSV: output/survey/indval_results.csv")
message("  CSV: output/survey/cluster_top10_ASVs.csv")
