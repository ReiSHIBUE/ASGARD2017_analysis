### P15_alpha_diversity.R
### ASGARD 2017 Processing — Shannon diversity by cluster and by filter fraction
### シャノン多様度：4クラスター別・3フィルター画分別のボックスプロット
###
### REQUIRES (from 00_setup.R, P01, P03):
###   asgard_seqcount2  - raw count matrix, 78 samples x 3076 ASVs (P01)
###   meta_asgard_p2    - metadata for the same 78 samples (P01)
###   clusnum_p         - 4-cluster assignment (P03)
###   rsc_p             - 4-colour cluster palette (P03)
###
### Shannon is computed on raw counts, as in S07_alpha_diversity.R.
### シャノン指数は S07 と同様に生カウントから計算する。
###
### OUTPUT:
###   output_p/alpha_diversity/ASGARD_shannon_processing.pdf
###   output_p/alpha_diversity/shannon_by_cluster.csv
###   output_p/alpha_diversity/shannon_by_filter.csv

library(vegan)
library(tidyverse)
library(here)

dir.create(here("output_p", "alpha_diversity"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: Shannon per sample / サンプルごとのシャノン指数
# ==============================================================================

shannon_p <- vegan::diversity(asgard_seqcount2, index = "shannon")

alpha_p <- data.frame(
  Sample  = rownames(asgard_seqcount2),
  Shannon = shannon_p,
  cluster = factor(as.character(clusnum_p[rownames(asgard_seqcount2)]),
                   levels = c("1", "2", "3", "4"),
                   labels = clus_labels_p),
  filter  = factor(meta_asgard_p2[rownames(asgard_seqcount2), "filter"],
                   levels = c("0.2 µm", "3 µm", "20 µm")),
  row.names = NULL
)

cc_p_alpha     <- setNames(rsc_p, clus_levels_p)
filter_colors_p <- c("0.2 µm" = "red", "3 µm" = "green", "20 µm" = "blue")

# ==============================================================================
# Section 2: Kruskal-Wallis tests / クラスカル・ウォリス検定
# ==============================================================================

kt_cluster <- kruskal.test(Shannon ~ cluster, data = alpha_p)
kt_filter  <- kruskal.test(Shannon ~ filter,  data = alpha_p)

message("\n--- Shannon diversity, processing samples (n = ", nrow(alpha_p), ") ---")
message(sprintf("by cluster (k=4):        chi-sq = %.3f, df = %d, p = %s",
                kt_cluster$statistic, kt_cluster$parameter,
                format(kt_cluster$p.value, digits = 4)))
message(sprintf("by filter fraction (3):  chi-sq = %.3f, df = %d, p = %s",
                kt_filter$statistic, kt_filter$parameter,
                format(kt_filter$p.value, digits = 4)))

pw_cluster <- pairwise.wilcox.test(alpha_p$Shannon, alpha_p$cluster,
                                   p.adjust.method = "BH")
pw_filter  <- pairwise.wilcox.test(alpha_p$Shannon, alpha_p$filter,
                                   p.adjust.method = "BH")

# ==============================================================================
# Section 3: Boxplots / ボックスプロット
# ==============================================================================

box_theme <- theme_bw(base_size = 15) +
  theme(axis.title   = element_text(size = 18, face = "bold"),
        axis.text    = element_text(size = 14),
        plot.subtitle = element_text(size = 13),
        legend.position = "none")

pdf(file = here("output_p", "alpha_diversity", "ASGARD_shannon_processing.pdf"),
    width = 9, height = 7)

# Page 1: by cluster
print(
  ggplot(alpha_p, aes(x = cluster, y = Shannon)) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.75, colour = "black") +
    scale_fill_manual(values = cc_p_alpha) +
    labs(x = "Cluster", y = "Shannon diversity",
         subtitle = sprintf("Kruskal-Wallis chi-sq = %.2f, df = %d, p = %s",
                            kt_cluster$statistic, kt_cluster$parameter,
                            format.pval(kt_cluster$p.value, digits = 3))) +
    box_theme
)

# Page 2: by filter fraction
print(
  ggplot(alpha_p, aes(x = filter, y = Shannon)) +
    geom_boxplot(aes(fill = filter), outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.75, colour = "black") +
    scale_fill_manual(values = filter_colors_p) +
    labs(x = "Filter fraction", y = "Shannon diversity",
         subtitle = sprintf("Kruskal-Wallis chi-sq = %.2f, df = %d, p = %s",
                            kt_filter$statistic, kt_filter$parameter,
                            format.pval(kt_filter$p.value, digits = 3))) +
    box_theme
)

dev.off()

# ==============================================================================
# Section 4: CSV summaries / 集計結果の書き出し
# ==============================================================================

summarise_by <- function(df, grp) {
  df %>%
    group_by(.data[[grp]]) %>%
    summarise(n      = n(),
              median = round(median(Shannon), 3),
              mean   = round(mean(Shannon), 3),
              sd     = round(sd(Shannon), 3),
              min    = round(min(Shannon), 3),
              max    = round(max(Shannon), 3),
              .groups = "drop")
}

write.csv(summarise_by(alpha_p, "cluster"),
          here("output_p", "alpha_diversity", "shannon_by_cluster.csv"),
          row.names = FALSE)
write.csv(summarise_by(alpha_p, "filter"),
          here("output_p", "alpha_diversity", "shannon_by_filter.csv"),
          row.names = FALSE)

message("\nP15_alpha_diversity.R: done.")
message("  PDF: output_p/alpha_diversity/ASGARD_shannon_processing.pdf")
message("  CSV: output_p/alpha_diversity/shannon_by_cluster.csv")
message("  CSV: output_p/alpha_diversity/shannon_by_filter.csv")
