### S07_alpha_diversity.R
### ASGARD 2017 Survey Site Analysis — Alpha Diversity
### アルファ多様性スクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02):
###   asgard_seqcount  - integer count matrix 181×3076 (raw counts, required by estimateR)
###   meta_asgard      - metadata 181×45
###   clusnum          - cluster assignments (length 181)
###
### PRODUCES:
###   asgard_alpha_df  - df with Shannon, Simpson, Chao1 + metadata per sample
###
### OUTPUT:
###   output/survey/alpha_diversity/ASGARD_alpha_diversity_survey.pdf
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
  Sample          = rownames(asgard_seqcount),
  Observed_ASVs   = rowSums(asgard_seqcount > 0),
  Shannon         = shannon_vals,
  Simpson         = simpson_vals,
  Chao1           = chao1_vals,
  cluster         = as.factor(clusnum[rownames(asgard_seqcount)])
)

asgard_alpha_df <- left_join(
  asgard_alpha_df,
  rownames_to_column(meta_asgard, var = "Sample"),
  by = "Sample"
)

# ==============================================================================
# Section 2: Kruskal-Wallis 検定 / Kruskal-Wallis tests by cluster
# ==============================================================================

for (idx in c("Shannon", "Simpson", "Chao1", "Observed_ASVs")) {
  kt <- kruskal.test(asgard_alpha_df[[idx]] ~ asgard_alpha_df$cluster)
  message(idx, ": Kruskal-Wallis chi-sq = ", round(kt$statistic, 3),
          ", df = ", kt$parameter, ", p = ", round(kt$p.value, 4))
}

# ==============================================================================
# Section 3: Boxplots by cluster / クラスター別 boxplot
# ==============================================================================

pdf(file = here::here("output", "survey", "alpha_diversity", "ASGARD_alpha_diversity_survey.pdf"),
    width = 8, height = 6)

for (idx in c("Observed_ASVs", "Shannon", "Simpson", "Chao1")) {
  gg <- ggplot(asgard_alpha_df, aes(x = cluster, y = .data[[idx]])) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA) +
    geom_jitter(width = .4, height = 0) +
    labs(title = paste(idx, "by Cluster — Survey Samples"),
         x = "Cluster", y = idx) +
    theme_minimal() +
    theme(text = element_text(size = 18))
  print(gg)
}

dev.off()

message("S07_alpha_diversity.R: done. asgard_alpha_df (", nrow(asgard_alpha_df), " rows) ready.")
message("  PDF: output/survey/alpha_diversity/ASGARD_alpha_diversity_survey.pdf")
