### S19_env_boxplots.R
### ASGARD 2017 Survey Site Analysis — Environmental Boxplots & Bloom/Autotrophy Index (11 clusters)
### 環境変数boxplot・ブルームインデックス（11クラスター）
###
### REQUIRES (from S01, S02):
###   meta_asgard      - metadata 181×45
###   clusnum11        - 11-cluster assignments (factor, from S02)
###   hier_levels_11   - ordered cluster names (from S02)
###   cc11             - 11-colour palette (from S02)
###
### OUTPUT:
###   output/survey/beta_diversity/ASGARD_boxplots_11clusters.pdf
###   output/survey/beta_diversity/ASGARD_bloom_autotrophy_index_11clusters.pdf

library(tidyverse)
library(vegan)
library(gridExtra)

# ==============================================================================
# Section 1: 環境変数boxplot / Environmental variable boxplots by 11 clusters
# ==============================================================================

df <- meta_asgard
df$cluster11 <- factor(as.character(clusnum11[rownames(df)]), levels = hier_levels_11)

env_plot_vars <- list(
  list(var = "temp",                 label = "Temperature (\u00B0C)"),
  list(var = "salinity",             label = "Salinity"),
  list(var = "DO",                   label = "DO"),
  list(var = "NO3(uM)",             label = "NO3 (\u00B5M)"),
  list(var = "PO4(uM)",             label = "PO4 (\u00B5M)"),
  list(var = "Sil(uM)",             label = "Sil (\u00B5M)"),
  list(var = "NH4(uM)",             label = "NH4 (\u00B5M)"),
  list(var = "N+N (umol/L)",        label = "N+N (\u00B5mol/L)"),
  list(var = "FlECO-AFL(mg/m^3)",   label = "FlECO-AFL (mg/m\u00B3)"),
  list(var = "chl (ug/l)",          label = "Chl (ug/l)"),
  list(var = "depth_m",             label = "Depth (m)"),
  list(var = "lat",                  label = "Latitude")
)

df$division <- factor(
  sub("^(.).*", "\\1", as.character(df$cluster11)),
  levels = c("A", "B", "C")
)

dir.create(here::here("output", "survey", "beta_diversity"), showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "beta_diversity", "ASGARD_boxplots_11clusters.pdf"),
    width = 14, height = 6)

for (ev in env_plot_vars) {
  print(
    ggplot(df, aes(x = cluster11, y = .data[[ev$var]])) +
      geom_boxplot(aes(fill = cluster11), outlier.shape = NA) +
      geom_jitter(width = 0.3, size = 1.2, alpha = 0.6) +
      scale_fill_manual(values = cc11, guide = "none") +
      facet_grid(~ division, scales = "free_x", space = "free_x") +
      labs(title = ev$label, x = "Cluster", y = ev$label) +
      theme_bw() +
      theme(text = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 16),
            strip.text = element_text(face = "bold", size = 14))
  )
}

# Combined page: 10 variables in one page
combined_vars <- list(
  list(var = "temp",                 label = "Temp (\u00B0C)"),
  list(var = "salinity",             label = "Salinity"),
  list(var = "DO",                   label = "DO"),
  list(var = "NO3(uM)",             label = "NO3 (\u00B5M)"),
  list(var = "PO4(uM)",             label = "PO4 (\u00B5M)"),
  list(var = "NH4(uM)",             label = "NH4 (\u00B5M)"),
  list(var = "FlECO-AFL(mg/m^3)",   label = "FlECO-AFL"),
  list(var = "chl (ug/l)",          label = "Chl (\u00B5g/l)"),
  list(var = "depth_m",             label = "Depth (m)"),
  list(var = "lat",                  label = "Latitude")
)

plot_list <- list()
for (i in seq_along(combined_vars)) {
  ev <- combined_vars[[i]]
  plot_list[[i]] <- ggplot(df, aes(x = cluster11, y = .data[[ev$var]])) +
    geom_boxplot(aes(fill = cluster11), outlier.shape = NA) +
    geom_jitter(width = 0.3, size = 0.6, alpha = 0.5) +
    scale_fill_manual(values = cc11, guide = "none") +
    facet_grid(~ division, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = ev$label) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(face = "bold", size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}

print(gridExtra::grid.arrange(grobs = plot_list, ncol = 2))

dev.off()

message("  PDF: output/survey/beta_diversity/ASGARD_boxplots_11clusters.pdf")

# ==============================================================================
# Section 2: PCA → Bloom Index / Autotrophy Index
# ==============================================================================

df$sample_id <- rownames(df)
env_vars <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)", "Sil(uM)",
              "NH4(uM)", "FlECO-AFL(mg/m^3)", "depth_m")
env_df <- df[, c("sample_id", "cluster11", env_vars)]
env_complete <- env_df[complete.cases(env_df), ]

env_mat <- env_complete[, env_vars]
log_vars <- c("NO3(uM)", "PO4(uM)", "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)")
for (v in log_vars) env_mat[[v]] <- log1p(env_mat[[v]])
env_scaled <- scale(env_mat)

pca_result <- prcomp(env_scaled, center = FALSE, scale. = FALSE)

pca_scores <- as.data.frame(pca_result$x)
pca_scores$cluster11 <- env_complete$cluster11
pca_scores$sample_id <- env_complete$sample_id

# Rescale PC1 → Bloom Index (0=pre-bloom, 1=post-bloom)
# Rescale PC2 → Autotrophy Index (0=low, 1=high)
rescale01 <- function(x) (x - min(x)) / (max(x) - min(x))
pca_scores$bloom_index      <- rescale01(pca_scores$PC1)
pca_scores$autotrophy_index <- rescale01(pca_scores$PC2)

# Division for faceting
pca_scores$division <- factor(
  sub("^(.).*", "\\1", as.character(pca_scores$cluster11)),
  levels = c("A", "B", "C")
)

# Add raw environmental variables (all 12)
all_env_cols <- c(env_vars, "N+N (umol/L)", "chl (ug/l)", "lat")
pca_scores <- left_join(pca_scores, df[, c("sample_id", all_env_cols)], by = "sample_id")

message("  PCA complete cases: ", nrow(pca_scores))

# ==============================================================================
# Section 3: Bloom/Autotrophy Index plots
# ==============================================================================

pdf(file = here::here("output", "survey", "beta_diversity",
                       "ASGARD_bloom_autotrophy_index_11clusters.pdf"),
    width = 14, height = 6)

# Page 1: Bloom Index boxplot
print(
  ggplot(pca_scores, aes(x = cluster11, y = bloom_index)) +
    geom_boxplot(aes(fill = cluster11), outlier.shape = NA) +
    geom_jitter(aes(color = cluster11), width = 0.3, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = cc11, guide = "none") +
    scale_color_manual(values = cc11, guide = "none") +
    facet_grid(~ division, scales = "free_x", space = "free_x") +
    labs(title = "Bloom Index by Cluster (rescaled PC1)",
         x = "Cluster", y = "Bloom Index (0 = pre-bloom, 1 = post-bloom)") +
    theme_bw() + theme(text = element_text(size = 13),
                        strip.text = element_text(face = "bold", size = 14))
)

# Page 2: Autotrophy Index boxplot
print(
  ggplot(pca_scores, aes(x = cluster11, y = autotrophy_index)) +
    geom_boxplot(aes(fill = cluster11), outlier.shape = NA) +
    geom_jitter(aes(color = cluster11), width = 0.3, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = cc11, guide = "none") +
    scale_color_manual(values = cc11, guide = "none") +
    facet_grid(~ division, scales = "free_x", space = "free_x") +
    labs(title = "Autotrophy Index by Cluster (rescaled PC2)",
         x = "Cluster", y = "Autotrophy Index (0 = low, 1 = high)") +
    theme_bw() + theme(text = element_text(size = 13),
                        strip.text = element_text(face = "bold", size = 14))
)

# Pages 3+: Environmental variables with dot size = Bloom/Autotrophy Index
for (ev in env_plot_vars) {
  p_bloom <- ggplot(pca_scores, aes(x = cluster11, y = .data[[ev$var]])) +
    geom_boxplot(aes(fill = cluster11), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = cluster11, size = bloom_index), width = 0.3, alpha = 0.7) +
    scale_fill_manual(values = cc11, guide = "none") +
    scale_color_manual(values = cc11, guide = "none") +
    scale_size_continuous(name = "Bloom\nIndex", range = c(0.5, 4)) +
    facet_grid(~ division, scales = "free_x", space = "free_x") +
    labs(title = paste(ev$label, "-- dot size: Bloom Index"),
         x = "Cluster", y = ev$label) +
    theme_bw() + theme(text = element_text(size = 12),
                        strip.text = element_text(face = "bold", size = 13))

  p_auto <- ggplot(pca_scores, aes(x = cluster11, y = .data[[ev$var]])) +
    geom_boxplot(aes(fill = cluster11), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = cluster11, size = autotrophy_index), width = 0.3, alpha = 0.7) +
    scale_fill_manual(values = cc11, guide = "none") +
    scale_color_manual(values = cc11, guide = "none") +
    scale_size_continuous(name = "Autotrophy\nIndex", range = c(0.5, 4)) +
    facet_grid(~ division, scales = "free_x", space = "free_x") +
    labs(title = paste(ev$label, "-- dot size: Autotrophy Index"),
         x = "Cluster", y = ev$label) +
    theme_bw() + theme(text = element_text(size = 12),
                        strip.text = element_text(face = "bold", size = 13))

  print(grid.arrange(p_bloom, p_auto, ncol = 1))
}

dev.off()

message("  PDF: output/survey/beta_diversity/ASGARD_bloom_autotrophy_index_11clusters.pdf")

# Clean up temp script output
rm(env_mat, env_scaled, pca_result, pca_scores, env_df, env_complete)

message("\nS19_env_boxplots.R: done.")
