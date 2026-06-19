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
  list(var = "NO3(uM)",             label = "NO3 (µM)"),
  list(var = "FlECO-AFL(mg/m^3)",   label = "FlECO-AFL (mg/m³)"),
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

message("\nP11_env_boxplots.R: done.")
message("  PDF: output_p/env_boxplots/ASGARD_boxplots_processing.pdf")
