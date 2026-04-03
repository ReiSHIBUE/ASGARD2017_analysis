### S18_cluster_specific_ASVs.R
### ASGARD 2017 Survey — Cluster-specific ASV identification (Wilcoxon)
### クラスター特異的ASV抽出（Wilcoxon片側検定）
###
### REQUIRES (from S01, S02):
###   asgard_filtered  - 181 x 258 ASV proportion matrix
###   clusnum10        - 10-cluster assignments
###   hier_names       - hierarchical cluster names
###   hier_levels      - cluster name order
###
### METHOD:
###   For each ASV and each cluster, compare within-cluster RA vs outside-cluster RA
###   using one-sided Wilcoxon rank-sum test (alternative = "greater").
###   BH-corrected p-values. Filter: padj < 0.05 and enrichment ratio > 2.
###
### OUTPUT:
###   output/survey/IndVal/cluster_specific_ASVs_wilcoxon.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Wilcoxon検定 / Wilcoxon rank-sum test per ASV per cluster
# ==============================================================================

mat <- as.matrix(asgard_filtered)
cl <- hier_names[as.character(clusnum10[rownames(mat)])]
names(cl) <- rownames(mat)

results <- list()

for (cl_name in hier_levels) {
  samples_in <- names(cl[cl == cl_name])
  samples_out <- names(cl[cl != cl_name])

  for (asv in colnames(mat)) {
    vals_in <- mat[samples_in, asv]
    vals_out <- mat[samples_out, asv]
    mean_in <- mean(vals_in)
    mean_out <- mean(vals_out)

    if (mean_in < 1e-5 & mean_out < 1e-5) next  # 両方ほぼゼロはスキップ

    wt <- wilcox.test(vals_in, vals_out, alternative = "greater", exact = FALSE)

    results[[length(results) + 1]] <- data.frame(
      cluster = cl_name,
      ASV = asv,
      mean_in = mean_in,
      mean_out = mean_out,
      ratio = (mean_in + 1e-6) / (mean_out + 1e-6),
      pval = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}

res_df <- bind_rows(results)
res_df$padj <- p.adjust(res_df$pval, method = "BH")

# ==============================================================================
# Section 2: 有意なASVを抽出 / Filter significant cluster-specific ASVs
# ==============================================================================

specific <- res_df %>%
  filter(padj < 0.05, ratio > 2) %>%
  arrange(cluster, desc(ratio))

message("=== Cluster-specific ASVs (padj < 0.05, ratio > 2) ===")
message("Total: ", nrow(specific))

for (cl_name in hier_levels) {
  n <- sum(specific$cluster == cl_name)
  message(sprintf("  %s: %d ASVs", cl_name, n))
}

# 短縮名を追加 / Add short name
specific$short_name <- sapply(specific$ASV, function(a) {
  parts <- strsplit(a, "; ")[[1]]
  paste(tail(parts, 2), collapse = "; ")
})

# ==============================================================================
# Section 3: CSV保存 / Save results
# ==============================================================================

specific_out <- specific %>%
  select(cluster, ASV, short_name, mean_in, mean_out, ratio, pval, padj) %>%
  mutate(mean_in = round(mean_in, 6),
         mean_out = round(mean_out, 6),
         ratio = round(ratio, 1),
         pval = signif(pval, 3),
         padj = signif(padj, 3))

dir.create(here("output", "survey", "IndVal"), showWarnings = FALSE, recursive = TRUE)
write.csv(specific_out,
  here("output", "survey", "IndVal", "cluster_specific_ASVs_wilcoxon.csv"),
  row.names = FALSE)

message("\nS18_cluster_specific_ASVs.R: done.")
message("  CSV: output/survey/IndVal/cluster_specific_ASVs_wilcoxon.csv")
