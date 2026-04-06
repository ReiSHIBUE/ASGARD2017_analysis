### S18_cluster_specific_ASVs.R
### ASGARD 2017 Survey — Cluster-specific ASV identification (Wilcoxon)
### クラスター特異的ASV抽出（Wilcoxon片側検定）
###
### REQUIRES (from S01, S02):
###   asgard_filtered  - 181 x 258 ASV proportion matrix
###   clusnum10        - 10-cluster assignments
###   hier_names       - hierarchical cluster names
###   hier_levels      - cluster name order
###   colclusnum       - column cluster assignments (6 clusters)
###
### METHOD:
###   For each ASV and each cluster, compare within-cluster RA vs outside-cluster RA
###   using one-sided Wilcoxon rank-sum test (alternative = "greater").
###   BH-corrected p-values. Filter: padj < 0.05 and enrichment ratio > 2.
###
### OUTPUT:
###   output/survey/IndVal/cluster_specific_ASVs_wilcoxon_with_colcluster.csv
###   output/survey/IndVal/intersection_top5_ASVs.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Wilcoxon検定 / Wilcoxon rank-sum test per ASV per cluster
# ==============================================================================

mat <- as.matrix(asgard_filtered)
cl <- hier_names[as.character(clusnum10[rownames(mat)])]
names(cl) <- rownames(mat)

cc <- colclusnum
names(cc) <- names(colclusnum)

results <- list()

for (cl_name in hier_levels) {
  samples_in <- names(cl[cl == cl_name])
  samples_out <- names(cl[cl != cl_name])

  for (asv in colnames(mat)) {
    vals_in <- mat[samples_in, asv]
    vals_out <- mat[samples_out, asv]
    mean_in <- mean(vals_in)
    mean_out <- mean(vals_out)

    if (mean_in < 1e-5 & mean_out < 1e-5) next

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
# Section 3: Wilcoxon結果にcolumn cluster情報を追加して保存
# Add column cluster (assemblage) to Wilcoxon results and save
# ==============================================================================

specific$col_cluster <- paste0("CC", cc[specific$ASV])

# Sample cluster × Column cluster のクロス集計
cross_tab <- specific %>%
  count(cluster, col_cluster) %>%
  pivot_wider(names_from = col_cluster, values_from = n, values_fill = 0)

message("\n=== Sample cluster x Column cluster (number of specific ASVs) ===")
print(as.data.frame(cross_tab))

dir.create(here("output", "survey", "IndVal"), showWarnings = FALSE, recursive = TRUE)

specific_with_cc <- specific %>%
  select(cluster, col_cluster, ASV, short_name, mean_in, mean_out, ratio, pval, padj) %>%
  arrange(cluster, col_cluster, desc(ratio)) %>%
  mutate(mean_in = round(mean_in, 6),
         mean_out = round(mean_out, 6),
         ratio = round(ratio, 1),
         pval = signif(pval, 3),
         padj = signif(padj, 3))

write.csv(specific_with_cc,
  here("output", "survey", "IndVal", "cluster_specific_ASVs_wilcoxon_with_colcluster.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 4: 交差部分のtop5 ASVs (sample cluster × column cluster)
# Top 5 ASVs per intersection, ranked by mean_in (within-cluster mean RA)
# ==============================================================================

int_results <- list()

for (scl in hier_levels) {
  samples_in <- names(cl[cl == scl])
  samples_out <- names(cl[cl != scl])

  for (ccl in 1:6) {
    asvs_in_cc <- names(cc[cc == ccl])

    for (asv in asvs_in_cc) {
      mean_in <- mean(mat[samples_in, asv])
      mean_out <- mean(mat[samples_out, asv])

      if (mean_in < 1e-5) next

      int_results[[length(int_results) + 1]] <- data.frame(
        sample_cluster = scl,
        col_cluster = paste0("CC", ccl),
        ASV = asv,
        mean_in = mean_in,
        mean_out = mean_out,
        stringsAsFactors = FALSE
      )
    }
  }
}

int_df <- bind_rows(int_results)

int_df$short_name <- sapply(int_df$ASV, function(a) {
  parts <- strsplit(a, "; ")[[1]]
  paste(tail(parts, 2), collapse = "; ")
})

top_intersection <- int_df %>%
  group_by(sample_cluster, col_cluster) %>%
  slice_max(order_by = mean_in, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(sample_cluster, col_cluster, desc(mean_in)) %>%
  select(sample_cluster, col_cluster, ASV, short_name, mean_in, mean_out) %>%
  mutate(mean_in = round(mean_in, 6),
         mean_out = round(mean_out, 6))

write.csv(top_intersection,
  here("output", "survey", "IndVal", "intersection_top5_ASVs.csv"),
  row.names = FALSE)

message("\nS18_cluster_specific_ASVs.R: done.")
message("  CSV: output/survey/IndVal/cluster_specific_ASVs_wilcoxon_with_colcluster.csv")
message("  CSV: output/survey/IndVal/intersection_top5_ASVs.csv")
