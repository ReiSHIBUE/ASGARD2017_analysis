### P11_env_boxplots.R
### ASGARD 2017 Processing — Environmental boxplots (4 clusters)
### 環境変数 boxplot（処理サイト、k=4 クラスタ + k=2 division）
###
### REQUIRES (from 00_setup.R, P01–P03):
###   meta_asgard_p2  - processing metadata 78 samples
###   clusnum_p       - 4-cluster assignment
###
### OUTPUT:
###   output_p/env_boxplots/ASGARD_boxplots_processing.pdf

library(tidyverse)
library(gridExtra)
library(here)

dir.create(here("output_p", "env_boxplots"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: メタデータ + cluster + k=2 division
# ==============================================================================

df_p <- meta_asgard_p2
df_p$cluster   <- factor(as.character(clusnum_p[rownames(df_p)]),
                         levels = c("1", "2", "3", "4"))
df_p$division2 <- factor(
  ifelse(df_p$cluster == "1", "Free-living", "Particle-associated"),
  levels = c("Free-living", "Particle-associated"))

cc_p <- c("1" = "#E41A1C", "2" = "#377EB8",
          "3" = "#4DAF4A", "4" = "#984EA3")

# ==============================================================================
# Section 2: PDF — single-variable pages + 6-var combined grid
# ==============================================================================

env_plot_vars <- list(
  list(var = "temp",                label = "Temperature (°C)"),
  list(var = "salinity",            label = "Salinity (PSU)"),
  list(var = "DO",                  label = "DO (µmol/kg)"),
  list(var = "FlECO-AFL(mg/m^3)",   label = "FlECO-AFL (mg/m³)"),
  list(var = "NO3(uM)",             label = "NO3 (µM)"),
  list(var = "depth_m",             label = "Sampling depth (m)")
)

pdf(file = here("output_p", "env_boxplots",
                "ASGARD_boxplots_processing.pdf"),
    width = 14, height = 6)

# Single-variable pages (one per variable, with k=2 division facet)
for (ev in env_plot_vars) {
  print(
    ggplot(df_p, aes(x = cluster, y = .data[[ev$var]])) +
      geom_boxplot(aes(fill = cluster), outlier.shape = NA) +
      geom_jitter(width = 0.3, size = 1.2, alpha = 0.6) +
      scale_fill_manual(values = cc_p, guide = "none") +
      facet_grid(~ division2, scales = "free_x", space = "free_x") +
      labs(title = ev$label, x = "Cluster", y = ev$label) +
      theme_bw() +
      theme(text       = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 16),
            strip.text = element_text(face = "bold", size = 14))
  )
}

# 6-variable combined grid (matches survey S19 final page)
six_plot_list <- list()
for (i in seq_along(env_plot_vars)) {
  ev <- env_plot_vars[[i]]
  six_plot_list[[i]] <- ggplot(df_p, aes(x = cluster, y = .data[[ev$var]])) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA) +
    geom_jitter(width = 0.3, size = 0.8, alpha = 0.5) +
    scale_fill_manual(values = cc_p, guide = "none") +
    facet_grid(~ division2, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = ev$label) +
    theme_bw(base_size = 14) +
    theme(strip.text  = element_text(face = "bold", size = 14),
          axis.title  = element_text(size = 14),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11))
}

print(gridExtra::grid.arrange(grobs = six_plot_list, ncol = 3))

dev.off()

message("  PDF: output_p/env_boxplots/ASGARD_boxplots_processing.pdf")

# ==============================================================================
# Section 3: PCA → Bloom Index / Autotrophy Index (k=4 clusters)
# ==============================================================================

df_p$sample_id <- rownames(df_p)
env_vars_p <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)", "Sil(uM)",
                "NH4(uM)", "FlECO-AFL(mg/m^3)", "depth_m")
env_p <- df_p[, c("sample_id", "cluster", env_vars_p)]
env_p_complete <- env_p[complete.cases(env_p), ]

env_mat_p <- env_p_complete[, env_vars_p]
log_vars_p <- c("NO3(uM)", "PO4(uM)", "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)")
for (v in log_vars_p) env_mat_p[[v]] <- log1p(env_mat_p[[v]])
env_scaled_p <- scale(env_mat_p)

pca_p <- prcomp(env_scaled_p, center = FALSE, scale. = FALSE)

# PCA sign is arbitrary; align orientation so the indices have the same meaning
# as in the survey pipeline (S19):
#   PC1: nutrient axis -- want NO3 loading > 0 so 1 - rescale(PC1) = post-bloom
#   PC2: autotrophy axis -- want FlECO loading > 0 so rescale(PC2) = high autotrophy
if (pca_p$rotation["NO3(uM)", "PC1"] < 0) {
  pca_p$rotation[, "PC1"] <- -pca_p$rotation[, "PC1"]
  pca_p$x[, "PC1"]        <- -pca_p$x[, "PC1"]
}
if (pca_p$rotation["FlECO-AFL(mg/m^3)", "PC2"] < 0) {
  pca_p$rotation[, "PC2"] <- -pca_p$rotation[, "PC2"]
  pca_p$x[, "PC2"]        <- -pca_p$x[, "PC2"]
}

pca_scores_p <- as.data.frame(pca_p$x)
pca_scores_p$cluster   <- env_p_complete$cluster
pca_scores_p$sample_id <- env_p_complete$sample_id

# Bloom Index = 1 - rescale(PC1)   (0 = pre-bloom, 1 = post-bloom)
# Autotrophy Index = rescale(PC2)  (0 = low, 1 = high)
rescale01 <- function(x) (x - min(x)) / (max(x) - min(x))
pca_scores_p$bloom_index      <- 1 - rescale01(pca_scores_p$PC1)
pca_scores_p$autotrophy_index <- rescale01(pca_scores_p$PC2)

pca_scores_p$division2 <- factor(
  ifelse(pca_scores_p$cluster == "1", "Free-living", "Particle-associated"),
  levels = c("Free-living", "Particle-associated"))

# Bring raw vars back for the dot-size plots
pca_scores_p <- left_join(pca_scores_p,
                          df_p[, c("sample_id", env_vars_p)],
                          by = "sample_id")

message("  PCA complete cases (processing): ", nrow(pca_scores_p))

# Save PCA outputs
dir.create(here("output_p", "env_boxplots"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(round(as.data.frame(pca_p$rotation), 4),
          here("output_p", "env_boxplots",
               "processing_pca_loadings_for_index.csv"))
write.csv(pca_scores_p %>%
            select(sample_id, cluster, division2,
                   PC1, PC2, bloom_index, autotrophy_index),
          here("output_p", "env_boxplots",
               "processing_pca_indices.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 4: Bloom / Autotrophy Index plots
# ==============================================================================

bloom_theme <- theme_bw() + theme(
  axis.title   = element_text(size = 18, face = "bold"),
  axis.text    = element_text(size = 14),
  strip.text   = element_text(face = "bold", size = 16),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text  = element_text(size = 12))

pdf(file = here("output_p", "env_boxplots",
                "ASGARD_bloom_autotrophy_index_processing.pdf"),
    width = 14, height = 6)

# Page 1: Bloom Index boxplot
print(
  ggplot(pca_scores_p, aes(x = cluster, y = bloom_index)) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA) +
    geom_jitter(aes(color = cluster), width = 0.3, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = cc_p, guide = "none") +
    scale_color_manual(values = cc_p, guide = "none") +
    facet_grid(~ division2, scales = "free_x", space = "free_x") +
    labs(x = "Cluster", y = "Bloom Index (0 = pre-bloom, 1 = post-bloom)") +
    bloom_theme
)

# Page 2: Autotrophy Index boxplot
print(
  ggplot(pca_scores_p, aes(x = cluster, y = autotrophy_index)) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA) +
    geom_jitter(aes(color = cluster), width = 0.3, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = cc_p, guide = "none") +
    scale_color_manual(values = cc_p, guide = "none") +
    facet_grid(~ division2, scales = "free_x", space = "free_x") +
    labs(x = "Cluster", y = "Autotrophy Index (0 = low, 1 = high)") +
    bloom_theme
)

# Pages 3+: Environmental variables with dot size = Bloom / Autotrophy Index
for (ev in env_plot_vars) {
  p_bloom <- ggplot(pca_scores_p, aes(x = cluster, y = .data[[ev$var]])) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = cluster, size = bloom_index),
                width = 0.3, alpha = 0.7) +
    scale_fill_manual(values = cc_p, guide = "none") +
    scale_color_manual(values = cc_p, guide = "none") +
    scale_size_continuous(name = "Bloom\nIndex", range = c(0.5, 4)) +
    facet_grid(~ division2, scales = "free_x", space = "free_x") +
    labs(x = "Cluster", y = ev$label) +
    bloom_theme

  p_auto <- ggplot(pca_scores_p, aes(x = cluster, y = .data[[ev$var]])) +
    geom_boxplot(aes(fill = cluster), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = cluster, size = autotrophy_index),
                width = 0.3, alpha = 0.7) +
    scale_fill_manual(values = cc_p, guide = "none") +
    scale_color_manual(values = cc_p, guide = "none") +
    scale_size_continuous(name = "Autotrophy\nIndex", range = c(0.5, 4)) +
    facet_grid(~ division2, scales = "free_x", space = "free_x") +
    labs(x = "Cluster", y = ev$label) +
    bloom_theme

  print(grid.arrange(p_bloom, p_auto, ncol = 1))
}

dev.off()

message("\nP11_env_boxplots.R: done.")
message("  PDF: output_p/env_boxplots/ASGARD_boxplots_processing.pdf")
message("  PDF: output_p/env_boxplots/ASGARD_bloom_autotrophy_index_processing.pdf")
message("  CSV: output_p/env_boxplots/processing_pca_loadings_for_index.csv")
message("  CSV: output_p/env_boxplots/processing_pca_indices.csv")
