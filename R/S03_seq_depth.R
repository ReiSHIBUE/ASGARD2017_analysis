### S03_seq_depth.R
### ASGARD 2017 Survey Site Analysis — Sequencing Depth & Rarefaction
### シーケンス深度・レアファクションスクリプト
###
### REQUIRES (from S01, S02):
###   asgard_seqcount  - integer count matrix 181×3076
###   meta_asgard      - metadata 181×45
###   clusnum          - cluster assignments (length 181)
###   rsc              - 5-colour palette vector
###
### PRODUCES:
###   asgard_seqdepth_df - df with total reads + metadata per sample (181 rows)
###
### OUTPUT:
###   output/survey/alpha_diversity/ASGARD_seq_depth_survey.pdf

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
asgard_seqdepth_df$cluster <- as.factor(clusnum[asgard_seqdepth_df$Sample])

# ==============================================================================
# Section 2: PDFに保存 / Save plots to PDF
# ==============================================================================

pdf(file = here::here("output", "survey", "alpha_diversity", "ASGARD_seq_depth_survey.pdf"),
    width = 12, height = 6)

# シーケンス深度のdot plot / Sequencing depth dot plot coloured by cluster
p_depth <- ggplot(
  asgard_seqdepth_df,
  aes(x = reorder(Sample, Reads), y = Reads, colour = cluster)
) +
  geom_point(size = 2) +
  geom_hline(yintercept = min(seq_depth), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = hue_pal()(5)) +
  labs(title = "Sequencing Depth per Survey Sample",
       x     = "Sample (ordered by depth)",
       y     = "Total 16S Reads") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

print(p_depth)

# シーケンス深度のヒストグラム / Histogram
hist(seq_depth, breaks = 20,
     main = "Distribution of 16S Read Depths — Survey Samples",
     xlab = "Total Reads per Sample")
abline(v = min(seq_depth), col = "red", lty = 2)

# レアファクション曲線 / Rarefaction curves
# rarecurve() can be slow; step=500 balances resolution and speed
rarecurve(asgard_seqcount, step = 500,
          col   = rsc,
          cex   = 0.5,
          label = FALSE,
          main  = "Rarefaction Curves — Survey Samples")
abline(v = min(seq_depth), col = "red", lty = 2)

dev.off()

message("S03_seq_depth.R: done.")
message("  min reads: ", min(seq_depth), ", max reads: ", max(seq_depth))
message("  PDF: output/survey/alpha_diversity/ASGARD_seq_depth_survey.pdf")
