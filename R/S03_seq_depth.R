### S03_seq_depth.R
### ASGARD 2017 Survey Site Analysis — Sequencing Depth & Rarefaction (11 clusters)
### シーケンス深度・レアファクションスクリプト（11クラスター）
###
### REQUIRES (from S01, S02):
###   asgard_seqcount  - integer count matrix 181×3076
###   meta_asgard      - metadata 181×45
###   clusnum11        - 11-cluster assignments (factor, from S02)
###   hier_levels_11   - ordered cluster names (from S02)
###   cc11             - 11-colour palette (from S02)
###   rsc11            - row side colours indexed by clusnum11 (from S02)
###
### PRODUCES:
###   asgard_seqdepth_df - df with total reads + cluster + metadata per sample (181 rows)
###
### OUTPUT:
###   output/survey/alpha_diversity/ASGARD_seq_depth_survey.pdf
###   output/survey/alpha_diversity/seq_depth_summary_11clusters.csv

library(tidyverse)
library(vegan)

# ==============================================================================
# Section 1: シーケンス深度の計算 / Compute sequencing depth per sample
# ==============================================================================

seq_depth <- rowSums(asgard_seqcount) # total 16S reads per sample, length 181

asgard_seqdepth_df <- data.frame(
  Sample = rownames(asgard_seqcount),
  Reads  = seq_depth
)
asgard_seqdepth_df <- left_join(
  asgard_seqdepth_df,
  rownames_to_column(meta_asgard, var = "Sample"),
  by = "Sample"
)
asgard_seqdepth_df$cluster11 <- factor(
  as.character(clusnum11[asgard_seqdepth_df$Sample]),
  levels = hier_levels_11
)

# ==============================================================================
# Section 2: クラスター別サマリー / Summary statistics per cluster
# ==============================================================================

depth_summary <- asgard_seqdepth_df %>%
  group_by(cluster11) %>%
  summarise(
    n            = n(),
    median_reads = median(Reads),
    mean_reads   = round(mean(Reads)),
    sd_reads     = round(sd(Reads)),
    min_reads    = min(Reads),
    max_reads    = max(Reads),
    .groups = "drop"
  )

message("\n--- Sequencing depth summary by 11 clusters ---")
print(as.data.frame(depth_summary))

# Kruskal-Wallis test on sequencing depth
kt_depth <- kruskal.test(Reads ~ cluster11, data = asgard_seqdepth_df)
message("\nKruskal-Wallis (seq depth ~ 11 clusters): chi-sq = ",
        round(kt_depth$statistic, 3), ", df = ", kt_depth$parameter,
        ", p = ", format(kt_depth$p.value, digits = 4))

# ==============================================================================
# Section 3: PDFに保存 / Save plots to PDF
# ==============================================================================

dir.create(here::here("output", "survey", "alpha_diversity"), showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "alpha_diversity", "ASGARD_seq_depth_survey.pdf"),
    width = 12, height = 6)

# Page 1: シーケンス深度のdot plot (11クラスター色)
p_depth <- ggplot(
  asgard_seqdepth_df,
  aes(x = reorder(Sample, Reads), y = Reads, colour = cluster11)
) +
  geom_point(size = 2) +
  geom_hline(yintercept = min(seq_depth), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = cc11) +
  labs(title = "Sequencing Depth per Survey Sample (11 clusters)",
       x     = "Sample (ordered by depth)",
       y     = "Total 16S Reads") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

print(p_depth)

# Page 2: クラスター別 boxplot
p_box <- ggplot(asgard_seqdepth_df, aes(x = cluster11, y = Reads)) +
  geom_boxplot(aes(fill = cluster11), outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 1.2, alpha = 0.6) +
  scale_fill_manual(values = cc11) +
  labs(title = "Sequencing Depth by 11 Clusters",
       x = "Cluster", y = "Total 16S Reads") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")

print(p_box)

# Page 3: シーケンス深度のヒストグラム / Histogram
hist(seq_depth, breaks = 20,
     main = "Distribution of 16S Read Depths — Survey Samples",
     xlab = "Total Reads per Sample")
abline(v = min(seq_depth), col = "red", lty = 2)

# Page 4: レアファクション曲線 / Rarefaction curves
rarecurve(asgard_seqcount, step = 500,
          col   = rsc11[rownames(asgard_seqcount)],
          cex   = 0.5,
          label = FALSE,
          main  = "Rarefaction Curves — Survey Samples (11 clusters)")
abline(v = min(seq_depth), col = "red", lty = 2)

dev.off()

# ==============================================================================
# Section 4: CSV保存 / Save summary CSV
# ==============================================================================

write.csv(depth_summary,
  here::here("output", "survey", "alpha_diversity", "seq_depth_summary_11clusters.csv"),
  row.names = FALSE)

message("\nS03_seq_depth.R: done.")
message("  min reads: ", min(seq_depth), ", max reads: ", max(seq_depth))
message("  PDF: output/survey/alpha_diversity/ASGARD_seq_depth_survey.pdf")
message("  CSV: output/survey/alpha_diversity/seq_depth_summary_11clusters.csv")
