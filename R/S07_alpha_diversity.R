### S07_alpha_diversity.R
### ASGARD 2017 Survey Site Analysis — Alpha Diversity (11 clusters)
### アルファ多様性スクリプト（サーベイサイト、11クラスター）
###
### REQUIRES (from S01, S02):
###   asgard_seqcount  - integer count matrix 181×3076 (raw counts, required by estimateR)
###   meta_asgard      - metadata 181×45
###   clusnum11        - 11-cluster assignments (factor, from S02)
###   hier_levels_11   - ordered cluster names (from S02)
###   cc11             - 11-colour palette (from S02)
###
### PRODUCES:
###   asgard_alpha_df  - df with Shannon, Simpson, Chao1 + metadata per sample
###
### OUTPUT:
###   output/survey/alpha_diversity/ASGARD_alpha_diversity_11clusters.pdf
###   output/survey/alpha_diversity/alpha_diversity_summary_11clusters.csv
###   output/survey/alpha_diversity/pairwise_wilcoxon_alpha_11clusters.csv
###
### NOTE: estimateR() requires raw integer counts (not proportions).
###   asgard_seqcount from S01 is correct; asgard_filtered (proportions) must NOT be used.

library(vegan)
library(tidyverse)

# ==============================================================================
# Section 1: アルファ多様性の計算 / Compute alpha diversity indices
# ==============================================================================

shannon_vals <- vegan::diversity(asgard_seqcount, index = "shannon")
simpson_vals <- vegan::diversity(asgard_seqcount, index = "simpson")
chao1_vals   <- estimateR(asgard_seqcount)[1, ]  # row 1 = Chao1

asgard_alpha_df <- data.frame(
  Sample        = rownames(asgard_seqcount),
  Observed_ASVs = rowSums(asgard_seqcount > 0),
  Shannon       = shannon_vals,
  Simpson       = simpson_vals,
  Chao1         = chao1_vals,
  cluster11     = factor(as.character(clusnum11[rownames(asgard_seqcount)]),
                         levels = hier_levels_11)
)

asgard_alpha_df <- left_join(
  asgard_alpha_df,
  rownames_to_column(meta_asgard, var = "Sample"),
  by = "Sample"
)

# ==============================================================================
# Section 2: Kruskal-Wallis 検定 / Kruskal-Wallis tests by 11 clusters
# ==============================================================================

message("\n--- Kruskal-Wallis tests (11 clusters) ---")

for (idx in c("Observed_ASVs", "Shannon", "Simpson", "Chao1")) {
  kt <- kruskal.test(asgard_alpha_df[[idx]] ~ asgard_alpha_df$cluster11)
  message(idx, ": Kruskal-Wallis chi-sq = ", round(kt$statistic, 3),
          ", df = ", kt$parameter, ", p = ", format(kt$p.value, digits = 4))
}

# ==============================================================================
# Section 3: Pairwise Wilcoxon 検定 (BH補正) / Pairwise Wilcoxon with BH correction
# ==============================================================================

message("\n--- Pairwise Wilcoxon tests (BH correction, 11 clusters) ---")

pw_all <- list()
for (idx in c("Observed_ASVs", "Shannon", "Simpson", "Chao1")) {
  pw <- pairwise.wilcox.test(asgard_alpha_df[[idx]], asgard_alpha_df$cluster11,
                              p.adjust.method = "BH")
  pw_mat <- as.data.frame(as.table(pw$p.value))
  pw_mat <- pw_mat[!is.na(pw_mat$Freq), ]
  colnames(pw_mat) <- c("cluster_a", "cluster_b", "p_adj")
  pw_mat$index <- idx
  pw_mat$sig <- ifelse(pw_mat$p_adj <= 0.001, "***",
                ifelse(pw_mat$p_adj <= 0.01,  "**",
                ifelse(pw_mat$p_adj <= 0.05,  "*", "ns")))
  pw_all[[idx]] <- pw_mat

  n_sig <- sum(pw_mat$sig != "ns")
  message(idx, ": significant pairs = ", n_sig, " / ", nrow(pw_mat))
}

pw_combined <- bind_rows(pw_all)
pw_combined$p_adj <- round(pw_combined$p_adj, 4)

# ==============================================================================
# Section 4: Boxplots by 11 clusters / 11クラスター別 boxplot
# ==============================================================================

dir.create(here::here("output", "survey", "alpha_diversity"), showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "alpha_diversity", "ASGARD_alpha_diversity_11clusters.pdf"),
    width = 12, height = 6)

for (idx in c("Observed_ASVs", "Shannon", "Simpson", "Chao1")) {
  gg <- ggplot(asgard_alpha_df, aes(x = cluster11, y = .data[[idx]])) +
    geom_boxplot(aes(fill = cluster11), outlier.shape = NA) +
    geom_jitter(width = 0.3, size = 1.2, alpha = 0.6) +
    scale_fill_manual(values = cc11) +
    labs(title = paste(idx, "by 11 Clusters"),
         x = "Cluster", y = idx) +
    theme_minimal() +
    theme(text = element_text(size = 14), legend.position = "none")
  print(gg)
}

dev.off()

# ==============================================================================
# Section 5: Save CSV / 結果をCSVに保存
# ==============================================================================

# Summary statistics per cluster
summ_all <- list()
for (idx in c("Observed_ASVs", "Shannon", "Simpson", "Chao1")) {
  s <- asgard_alpha_df %>%
    group_by(cluster11) %>%
    summarise(
      n      = n(),
      median = round(median(.data[[idx]]), 3),
      mean   = round(mean(.data[[idx]]), 3),
      sd     = round(sd(.data[[idx]]), 3),
      min    = round(min(.data[[idx]]), 3),
      max    = round(max(.data[[idx]]), 3),
      .groups = "drop"
    )
  s$index <- idx
  summ_all[[idx]] <- s
}
summ_combined <- bind_rows(summ_all)

write.csv(summ_combined,
  here::here("output", "survey", "alpha_diversity", "alpha_diversity_summary_11clusters.csv"),
  row.names = FALSE)

write.csv(pw_combined,
  here::here("output", "survey", "alpha_diversity", "pairwise_wilcoxon_alpha_11clusters.csv"),
  row.names = FALSE)

message("\nS07_alpha_diversity.R: done. asgard_alpha_df (", nrow(asgard_alpha_df), " rows) ready.")
message("  PDF: output/survey/alpha_diversity/ASGARD_alpha_diversity_11clusters.pdf")
message("  CSV: output/survey/alpha_diversity/alpha_diversity_summary_11clusters.csv")
message("  CSV: output/survey/alpha_diversity/pairwise_wilcoxon_alpha_11clusters.csv")
