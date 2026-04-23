library(here); library(tidyverse); library(vegan)
source(here("R/00_setup.R")); source(here("R/S01_data_prep.R")); source(here("R/S02_heatmaps_16S.R"))

df <- meta_asgard
df$cluster <- factor(as.character(clusnum11[rownames(df)]), levels = hier_levels_11)
df$sample_id <- rownames(df)

# Select environmental variables
env_vars <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)", "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)", "depth_m")

# Keep complete cases only
env_df <- df[, c("sample_id", "cluster", env_vars)]
env_complete <- env_df[complete.cases(env_df), ]
cat("Complete cases:", nrow(env_complete), "/ 181\n")

# Extract numeric matrix
env_mat <- env_complete[, env_vars]

# Check for zero variance
cat("\nVariance per variable:\n")
print(round(apply(env_mat, 2, var), 4))

# Log-transform skewed variables (nutrients), then scale all
# NO3, PO4, Sil, NH4, FlECO are right-skewed -> log(x + 1) transform
log_vars <- c("NO3(uM)", "PO4(uM)", "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)")
for (v in log_vars) {
  env_mat[[v]] <- log1p(env_mat[[v]])
}

# Scale all variables (mean=0, sd=1)
env_scaled <- scale(env_mat)

cat("\nAfter log1p + scale:\n")
print(round(apply(env_scaled, 2, mean), 4))  # should be ~0
print(round(apply(env_scaled, 2, sd), 4))    # should be ~1

# PCA
pca_result <- prcomp(env_scaled, center = FALSE, scale. = FALSE)  # already scaled

cat("\nPCA summary:\n")
print(summary(pca_result))

# Extract scores
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$cluster <- env_complete$cluster
pca_scores$sample_id <- env_complete$sample_id

# Extract loadings (arrows)
pca_loadings <- as.data.frame(pca_result$rotation[, 1:2])
pca_loadings$variable <- rownames(pca_loadings)

# Variance explained
var_expl <- summary(pca_result)$importance["Proportion of Variance", ]
pc1_pct <- round(var_expl[1] * 100, 1)
pc2_pct <- round(var_expl[2] * 100, 1)

cat("\nPC1:", pc1_pct, "%, PC2:", pc2_pct, "%\n")

# Scale loadings for biplot (multiply by sqrt of eigenvalue for better display)
arrow_scale <- 3
pca_loadings$PC1_plot <- pca_loadings$PC1 * arrow_scale
pca_loadings$PC2_plot <- pca_loadings$PC2 * arrow_scale

# Shorter variable names for plot
pca_loadings$label <- c("Temp", "Salinity", "DO", "NO3", "PO4", "Sil", "NH4", "FlECO", "Depth")

# Plot
dir.create(here("output", "survey", "beta_diversity"), showWarnings = FALSE, recursive = TRUE)

pdf(file = here("output", "survey", "beta_diversity", "ASGARD_pca_env_survey.pdf"),
    width = 10, height = 8)

# Page 1: PCA biplot
p1 <- ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = cluster),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = cc11) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_loadings,
            aes(x = PC1_plot * 1.15, y = PC2_plot * 1.15, label = label),
            size = 4, fontface = "bold") +
  coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings$PC1_plot)) * 1.3,
                  ylim = range(c(pca_scores$PC2, pca_loadings$PC2_plot)) * 1.3) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Cluster",
       title = "PCA of Environmental Variables (11 clusters)",
       subtitle = "log1p-transformed nutrients + scaled; 9 variables") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))
print(p1)

# Page 2: Without arrows (cleaner view)
p2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = cc11) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC2) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Cluster",
       title = "PCA of Environmental Variables (points only)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16))
print(p2)

# Pages 3-11: Per-variable PCA plots (dot size = variable value)
env_labels <- c("temp" = "Temperature (\u00B0C)",
                "salinity" = "Salinity",
                "DO" = "Dissolved Oxygen",
                "NO3(uM)" = "NO3 (\u00B5M)",
                "PO4(uM)" = "PO4 (\u00B5M)",
                "Sil(uM)" = "Silicate (\u00B5M)",
                "NH4(uM)" = "NH4 (\u00B5M)",
                "FlECO-AFL(mg/m^3)" = "FlECO-AFL (mg/m\u00B3)",
                "depth_m" = "Depth (m)")

# Add raw env values to pca_scores
for (v in env_vars) {
  pca_scores[[v]] <- env_complete[[v]]
}

for (v in env_vars) {
  p_env <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = cc11) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Cluster",
         title = paste0("PCA — dot size: ", env_labels[v])) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16))
  print(p_env)
}

# Pages 12-20: Per-variable PCA plots faceted by depth_type
# Each page: all + surf + mid + bottom in one view

pca_scores$depth_type <- factor(df$depth_type[match(pca_scores$sample_id, df$sample_id)],
                                levels = c("surf", "mid", "bottom"))

for (v in env_vars) {
  # All depths combined
  p_all <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = cc11) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Cluster",
         title = paste0(env_labels[v], " — all depths")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "right")

  # Faceted by depth_type
  p_facet <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = cc11) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    facet_wrap(~ depth_type) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Cluster",
         title = paste0(env_labels[v], " — by depth type")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold", size = 11),
          legend.position = "right")

  # Combine on one page using patchwork-style grid
  print(gridExtra::grid.arrange(p_all, p_facet, ncol = 1, heights = c(1, 1)))
}

dev.off()

cat("\nDone: output/survey/beta_diversity/ASGARD_pca_env_survey.pdf\n")

# ==============================================================================
# 3グループ (A/B/C) 版PCAプロット
# ==============================================================================

group_colors <- c("A" = "#E31A1C", "B" = "#33A02C", "C" = "#1F78B4")
pca_scores$group <- sub("^(.).*", "\\1", as.character(pca_scores$cluster))

pdf(file = here("output", "survey", "beta_diversity", "ASGARD_pca_env_3groups.pdf"),
    width = 10, height = 8)

# Page 1: PCA biplot with 3 groups
print(ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = group),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_loadings,
            aes(x = PC1_plot * 1.15, y = PC2_plot * 1.15, label = label),
            size = 4, fontface = "bold") +
  coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings$PC1_plot)) * 1.3,
                  ylim = range(c(pca_scores$PC2, pca_loadings$PC2_plot)) * 1.3) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Group",
       title = "PCA of Environmental Variables (3 groups)",
       subtitle = "log1p-transformed nutrients + scaled; 9 variables") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12)))

# Page 2: Points only
print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC2) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Group",
       title = "PCA of Environmental Variables (3 groups, points only)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 3: Faceted by A/B/C
print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ group) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC2) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Group",
       title = "PCA of Environmental Variables (faceted by group)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 13)))

# Pages 4+: Per-variable, dot size = env value
for (v in env_vars) {
  print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = group, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Group",
         title = paste0("PCA -- dot size: ", env_labels[v], " (3 groups)")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16)))
}

dev.off()

cat("Done: output/survey/beta_diversity/ASGARD_pca_env_3groups.pdf\n")

# ==============================================================================
# PC1 vs PC3 版プロット (10クラスター + 3グループ + 4グループ)
# ==============================================================================

# PC3の分散説明率
pc3_pct <- round(summary(pca_result)$importance["Proportion of Variance", 3] * 100, 1)

# PC3 loadings for arrows
pca_loadings_13 <- as.data.frame(pca_result$rotation[, c(1, 3)])
colnames(pca_loadings_13) <- c("PC1", "PC3")
pca_loadings_13$variable <- rownames(pca_loadings_13)
pca_loadings_13$PC1_plot <- pca_loadings_13$PC1 * arrow_scale
pca_loadings_13$PC3_plot <- pca_loadings_13$PC3 * arrow_scale
pca_loadings_13$label <- c("Temp", "Salinity", "DO", "NO3", "PO4", "Sil", "NH4", "FlECO", "Depth")

# PC3 scores
pca_scores$PC3 <- pca_result$x[, 3]

pdf(file = here("output", "survey", "beta_diversity", "ASGARD_pca_env_PC1vsPC3.pdf"),
    width = 10, height = 8)

# Page 1: 11 clusters biplot
print(ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC3, color = cluster),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = cc11) +
  geom_segment(data = pca_loadings_13,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC3_plot),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_loadings_13,
            aes(x = PC1_plot * 1.15, y = PC3_plot * 1.15, label = label),
            size = 4, fontface = "bold") +
  coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings_13$PC1_plot)) * 1.3,
                  ylim = range(c(pca_scores$PC3, pca_loadings_13$PC3_plot)) * 1.3) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC3 (", pc3_pct, "%)"),
       color = "Cluster",
       title = "PCA PC1 vs PC3 (11 clusters)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 2: 3 groups biplot
print(ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC3, color = group),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  geom_segment(data = pca_loadings_13,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC3_plot),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_loadings_13,
            aes(x = PC1_plot * 1.15, y = PC3_plot * 1.15, label = label),
            size = 4, fontface = "bold") +
  coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings_13$PC1_plot)) * 1.3,
                  ylim = range(c(pca_scores$PC3, pca_loadings_13$PC3_plot)) * 1.3) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC3 (", pc3_pct, "%)"),
       color = "Group",
       title = "PCA PC1 vs PC3 (3 groups)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 3: 3 groups faceted
print(ggplot(pca_scores, aes(x = PC1, y = PC3, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ group) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC3) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC3 (", pc3_pct, "%)"),
       color = "Group",
       title = "PCA PC1 vs PC3 (faceted by A/B/C)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 13)))

# Pages 4+: Per-variable, dot size
for (v in env_vars) {
  print(ggplot(pca_scores, aes(x = PC1, y = PC3, color = group, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC3) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC3 (", pc3_pct, "%)"),
         color = "Group",
         title = paste0("PC1 vs PC3 -- dot size: ", env_labels[v])) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16)))
}

dev.off()

cat("Done: output/survey/beta_diversity/ASGARD_pca_env_PC1vsPC3.pdf\n")

# ==============================================================================
# PC1 vs PC3 版プロット (4グループ A/B/C1/C2)
# ==============================================================================

# group4はまだ定義されていないので先に定義
group4_colors_13 <- c("A" = "#E31A1C", "B" = "#33A02C", "C1" = "#1F78B4", "C2" = "#6A3D9A")
group4_map_13 <- c("A1"="A", "A2"="A", "B1"="B", "B2a"="B", "B2b"="B",
                   "C1a"="C1", "C1b1"="C1", "C1b2"="C1", "C2a"="C2", "C2b1"="C2", "C2b2"="C2")
pca_scores$group4 <- group4_map_13[as.character(pca_scores$cluster)]

pdf(file = here("output", "survey", "beta_diversity", "ASGARD_pca_env_4groups_PC1vsPC3.pdf"),
    width = 10, height = 8)

# Page 1: 4 groups biplot
print(ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC3, color = group4),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = group4_colors_13) +
  geom_segment(data = pca_loadings_13,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC3_plot),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_loadings_13,
            aes(x = PC1_plot * 1.15, y = PC3_plot * 1.15, label = label),
            size = 4, fontface = "bold") +
  coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings_13$PC1_plot)) * 1.3,
                  ylim = range(c(pca_scores$PC3, pca_loadings_13$PC3_plot)) * 1.3) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC3 (", pc3_pct, "%)"),
       color = "Group",
       title = "PCA PC1 vs PC3 (4 groups)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 2: Points only
print(ggplot(pca_scores, aes(x = PC1, y = PC3, color = group4)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group4_colors_13) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC3) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC3 (", pc3_pct, "%)"),
       color = "Group",
       title = "PCA PC1 vs PC3 (4 groups, points only)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 3: Faceted
print(ggplot(pca_scores, aes(x = PC1, y = PC3, color = group4)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group4_colors_13) +
  facet_wrap(~ group4) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC3) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC3 (", pc3_pct, "%)"),
       color = "Group",
       title = "PCA PC1 vs PC3 (faceted by 4 groups)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 13)))

# Pages 4+: Per-variable, dot size
for (v in env_vars) {
  print(ggplot(pca_scores, aes(x = PC1, y = PC3, color = group4, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = group4_colors_13) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC3) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC3 (", pc3_pct, "%)"),
         color = "Group",
         title = paste0("PC1 vs PC3 -- dot size: ", env_labels[v], " (4 groups)")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16)))
}

dev.off()

cat("Done: output/survey/beta_diversity/ASGARD_pca_env_4groups_PC1vsPC3.pdf\n")

# ==============================================================================
# 4グループ (A/B/C1/C2) 版PCAプロット
# ==============================================================================

group4_colors <- c("A" = "#E31A1C", "B" = "#33A02C", "C1" = "#1F78B4", "C2" = "#6A3D9A")
group4_map_plot <- c("A1"="A", "A2"="A", "B1"="B", "B2a"="B", "B2b"="B",
                     "C1a"="C1", "C1b1"="C1", "C1b2"="C1", "C2a"="C2", "C2b1"="C2", "C2b2"="C2")
pca_scores$group4 <- group4_map_plot[as.character(pca_scores$cluster)]

pdf(file = here("output", "survey", "beta_diversity", "ASGARD_pca_env_4groups.pdf"),
    width = 10, height = 8)

# Page 1: PCA biplot with 4 groups
print(ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = group4),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = group4_colors) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = pca_loadings,
            aes(x = PC1_plot * 1.15, y = PC2_plot * 1.15, label = label),
            size = 4, fontface = "bold") +
  coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings$PC1_plot)) * 1.3,
                  ylim = range(c(pca_scores$PC2, pca_loadings$PC2_plot)) * 1.3) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Group",
       title = "PCA of Environmental Variables (4 groups)",
       subtitle = "log1p-transformed nutrients + scaled; 9 variables") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12)))

# Page 2: Points only
print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = group4)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group4_colors) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC2) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Group",
       title = "PCA of Environmental Variables (4 groups, points only)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16)))

# Page 3: Faceted by A/B/C1/C2
print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = group4)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = group4_colors) +
  facet_wrap(~ group4) +
  coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                  ylim = range(pca_scores$PC2) * 1.1) +
  labs(x = paste0("PC1 (", pc1_pct, "%)"),
       y = paste0("PC2 (", pc2_pct, "%)"),
       color = "Group",
       title = "PCA of Environmental Variables (faceted by 4 groups)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 13)))

# Pages 4+: Per-variable, dot size = env value
for (v in env_vars) {
  print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = group4, size = .data[[v]])) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = group4_colors) +
    scale_size_continuous(range = c(1, 8), name = env_labels[v]) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Group",
         title = paste0("PCA -- dot size: ", env_labels[v], " (4 groups)")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16)))
}

dev.off()

cat("Done: output/survey/beta_diversity/ASGARD_pca_env_4groups.pdf\n")

# ==============================================================================
# Section: 環境変数による PERMANOVA / PERMDISP (3グループ A/B/C)
# Environmental PERMANOVA / PERMDISP by 3 groups using Euclidean distance
# on log1p-transformed and scaled environmental variables
# ==============================================================================

cat("\n--- Environmental PERMANOVA / PERMDISP by 3 groups (A/B/C) ---\n")

# PERMANOVA用: 全PC軸スコアのdf (174 x 9) + cluster + sample_id
pca_all <- as.data.frame(pca_result$x)  # PC1〜PC9
pca_all$cluster <- env_complete$cluster
pca_all$sample_id <- env_complete$sample_id

# pca_allのPC1〜PC9からEuclidean距離を計算
env_eucdist <- vegdist(pca_all[, paste0("PC", 1:9)], method = "euclidean")

# 3グループ割り当て
env_group3 <- factor(
  sub("^(.).*", "\\1", as.character(env_complete$cluster)),
  levels = c("A", "B", "C")
)

cat("Group sizes: A =", sum(env_group3 == "A"),
    " B =", sum(env_group3 == "B"),
    " C =", sum(env_group3 == "C"), "\n")

# PERMANOVA
env_permanova_3 <- adonis2(env_eucdist ~ env_group3, permutations = 999)
cat("\nPERMANOVA (env, 3 groups):\n")
print(env_permanova_3)

# Pairwise PERMANOVA
groups3 <- c("A", "B", "C")
pairs3 <- combn(groups3, 2)
env_pairwise_3 <- data.frame(
  pair = character(), F_value = numeric(), R2 = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

for (i in seq_len(ncol(pairs3))) {
  cl_a <- pairs3[1, i]
  cl_b <- pairs3[2, i]
  idx_a <- which(env_group3 == cl_a)
  idx_b <- which(env_group3 == cl_b)
  idx_ab <- c(idx_a, idx_b)
  dist_sub <- as.dist(as.matrix(env_eucdist)[idx_ab, idx_ab])
  group_sub <- factor(c(rep(cl_a, length(idx_a)), rep(cl_b, length(idx_b))))
  res_pair <- adonis2(dist_sub ~ group_sub, permutations = 999)
  env_pairwise_3 <- rbind(env_pairwise_3, data.frame(
    pair = paste0(cl_a, " vs ", cl_b),
    F_value = round(res_pair$F[1], 2),
    R2 = round(res_pair$R2[1], 4),
    p_value = res_pair$`Pr(>F)`[1]
  ))
}
env_pairwise_3$p_adj <- round(p.adjust(env_pairwise_3$p_value, method = "BH"), 4)
env_pairwise_3$sig <- ifelse(env_pairwise_3$p_adj <= 0.001, "***",
                      ifelse(env_pairwise_3$p_adj <= 0.01, "**",
                      ifelse(env_pairwise_3$p_adj <= 0.05, "*", "ns")))

cat("\nPairwise PERMANOVA (env, 3 groups):\n")
print(env_pairwise_3, row.names = FALSE)

# PERMDISP
env_permdisp_3 <- betadisper(env_eucdist, group = env_group3)
env_permdisp_3_anova <- anova(env_permdisp_3)
cat("\nPERMDISP (env, 3 groups):\n")
print(env_permdisp_3_anova)

# Pairwise PERMDISP (TukeyHSD)
env_permdisp_3_tukey <- TukeyHSD(env_permdisp_3)
env_permdisp_pw3 <- as.data.frame(env_permdisp_3_tukey$group)
env_permdisp_pw3$pair <- rownames(env_permdisp_pw3)
env_permdisp_pw3$sig <- ifelse(env_permdisp_pw3$`p adj` <= 0.001, "***",
                        ifelse(env_permdisp_pw3$`p adj` <= 0.01, "**",
                        ifelse(env_permdisp_pw3$`p adj` <= 0.05, "*", "ns")))
env_permdisp_pw3 <- env_permdisp_pw3[, c("pair", "diff", "lwr", "upr", "p adj", "sig")]
colnames(env_permdisp_pw3) <- c("pair", "diff_dispersion", "CI_lower", "CI_upper", "p_adj", "sig")
env_permdisp_pw3$diff_dispersion <- round(env_permdisp_pw3$diff_dispersion, 4)
env_permdisp_pw3$CI_lower <- round(env_permdisp_pw3$CI_lower, 4)
env_permdisp_pw3$CI_upper <- round(env_permdisp_pw3$CI_upper, 4)
env_permdisp_pw3$p_adj <- round(env_permdisp_pw3$p_adj, 4)

cat("\nPairwise PERMDISP (env, 3 groups):\n")
print(env_permdisp_pw3, row.names = FALSE)

# ==============================================================================
# Section: 環境変数による PERMANOVA / PERMDISP (4グループ A/B/C1/C2)
# ==============================================================================

cat("\n--- Environmental PERMANOVA / PERMDISP by 4 groups (A/B/C1/C2) ---\n")

group4_map <- c("A1"="A", "A2"="A", "B1"="B", "B2a"="B", "B2b"="B",
                "C1a"="C1", "C1b1"="C1", "C1b2"="C1", "C2a"="C2", "C2b1"="C2", "C2b2"="C2")
env_group4 <- factor(
  group4_map[as.character(env_complete$cluster)],
  levels = c("A", "B", "C1", "C2")
)

cat("Group sizes: A =", sum(env_group4 == "A"),
    " B =", sum(env_group4 == "B"),
    " C1 =", sum(env_group4 == "C1"),
    " C2 =", sum(env_group4 == "C2"), "\n")

# PERMANOVA
env_permanova_4 <- adonis2(env_eucdist ~ env_group4, permutations = 999)
cat("\nPERMANOVA (env, 4 groups):\n")
print(env_permanova_4)

# Pairwise PERMANOVA
groups4 <- c("A", "B", "C1", "C2")
pairs4 <- combn(groups4, 2)
env_pairwise_4 <- data.frame(
  pair = character(), F_value = numeric(), R2 = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

for (i in seq_len(ncol(pairs4))) {
  cl_a <- pairs4[1, i]
  cl_b <- pairs4[2, i]
  idx_a <- which(env_group4 == cl_a)
  idx_b <- which(env_group4 == cl_b)
  idx_ab <- c(idx_a, idx_b)
  dist_sub <- as.dist(as.matrix(env_eucdist)[idx_ab, idx_ab])
  group_sub <- factor(c(rep(cl_a, length(idx_a)), rep(cl_b, length(idx_b))))
  res_pair <- adonis2(dist_sub ~ group_sub, permutations = 999)
  env_pairwise_4 <- rbind(env_pairwise_4, data.frame(
    pair = paste0(cl_a, " vs ", cl_b),
    F_value = round(res_pair$F[1], 2),
    R2 = round(res_pair$R2[1], 4),
    p_value = res_pair$`Pr(>F)`[1]
  ))
}
env_pairwise_4$p_adj <- round(p.adjust(env_pairwise_4$p_value, method = "BH"), 4)
env_pairwise_4$sig <- ifelse(env_pairwise_4$p_adj <= 0.001, "***",
                      ifelse(env_pairwise_4$p_adj <= 0.01, "**",
                      ifelse(env_pairwise_4$p_adj <= 0.05, "*", "ns")))

cat("\nPairwise PERMANOVA (env, 4 groups):\n")
print(env_pairwise_4, row.names = FALSE)

# PERMDISP
env_permdisp_4 <- betadisper(env_eucdist, group = env_group4)
env_permdisp_4_anova <- anova(env_permdisp_4)
cat("\nPERMDISP (env, 4 groups):\n")
print(env_permdisp_4_anova)

# Pairwise PERMDISP
env_permdisp_4_tukey <- TukeyHSD(env_permdisp_4)
env_permdisp_pw4 <- as.data.frame(env_permdisp_4_tukey$group)
env_permdisp_pw4$pair <- rownames(env_permdisp_pw4)
env_permdisp_pw4$sig <- ifelse(env_permdisp_pw4$`p adj` <= 0.001, "***",
                        ifelse(env_permdisp_pw4$`p adj` <= 0.01, "**",
                        ifelse(env_permdisp_pw4$`p adj` <= 0.05, "*", "ns")))
env_permdisp_pw4 <- env_permdisp_pw4[, c("pair", "diff", "lwr", "upr", "p adj", "sig")]
colnames(env_permdisp_pw4) <- c("pair", "diff_dispersion", "CI_lower", "CI_upper", "p_adj", "sig")
env_permdisp_pw4$diff_dispersion <- round(env_permdisp_pw4$diff_dispersion, 4)
env_permdisp_pw4$CI_lower <- round(env_permdisp_pw4$CI_lower, 4)
env_permdisp_pw4$CI_upper <- round(env_permdisp_pw4$CI_upper, 4)
env_permdisp_pw4$p_adj <- round(env_permdisp_pw4$p_adj, 4)

cat("\nPairwise PERMDISP (env, 4 groups):\n")
print(env_permdisp_pw4, row.names = FALSE)

# ==============================================================================
# Section: 環境変数による PERMANOVA / PERMDISP (11クラスター)
# Environmental PERMANOVA / PERMDISP by 11 clusters (C1b split into C1b1/C1b2)
# ==============================================================================

cat("\n--- Environmental PERMANOVA / PERMDISP by 11 clusters ---\n")

# clusnum11のうちenv_completeにあるサンプルのみ
env_clus11 <- factor(
  as.character(clusnum11[env_complete$sample_id]),
  levels = hier_levels_11
)

cat("Group sizes:\n")
print(table(env_clus11))

# PERMANOVA
env_permanova_11 <- adonis2(env_eucdist ~ env_clus11, permutations = 999)
cat("\nPERMANOVA (env, 11 clusters):\n")
print(env_permanova_11)

# Pairwise PERMANOVA (55 pairs)
pairs11 <- combn(hier_levels_11, 2)
env_pairwise_11 <- data.frame(
  pair = character(), F_value = numeric(), R2 = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

for (i in seq_len(ncol(pairs11))) {
  cl_a <- pairs11[1, i]
  cl_b <- pairs11[2, i]
  idx_a <- which(env_clus11 == cl_a)
  idx_b <- which(env_clus11 == cl_b)
  idx_ab <- c(idx_a, idx_b)
  dist_sub <- as.dist(as.matrix(env_eucdist)[idx_ab, idx_ab])
  group_sub <- factor(c(rep(cl_a, length(idx_a)), rep(cl_b, length(idx_b))))
  res_pair <- adonis2(dist_sub ~ group_sub, permutations = 999)
  env_pairwise_11 <- rbind(env_pairwise_11, data.frame(
    pair    = paste0(cl_a, " vs ", cl_b),
    F_value = round(res_pair$F[1], 2),
    R2      = round(res_pair$R2[1], 4),
    p_value = res_pair$`Pr(>F)`[1]
  ))
}
env_pairwise_11$p_adj <- round(p.adjust(env_pairwise_11$p_value, method = "BH"), 4)
env_pairwise_11$sig <- ifelse(env_pairwise_11$p_adj <= 0.001, "***",
                       ifelse(env_pairwise_11$p_adj <= 0.01,  "**",
                       ifelse(env_pairwise_11$p_adj <= 0.05,  "*", "ns")))

cat("\nPairwise PERMANOVA (env, 11 clusters):\n")
print(env_pairwise_11, row.names = FALSE)

# PERMDISP
env_permdisp_11       <- betadisper(env_eucdist, group = env_clus11)
env_permdisp_11_anova <- anova(env_permdisp_11)
cat("\nPERMDISP (env, 11 clusters):\n")
print(env_permdisp_11_anova)

# Pairwise PERMDISP (TukeyHSD)
env_permdisp_11_tukey <- TukeyHSD(env_permdisp_11)
env_permdisp_pw11 <- as.data.frame(env_permdisp_11_tukey$group)
env_permdisp_pw11$pair <- rownames(env_permdisp_pw11)
env_permdisp_pw11$sig <- ifelse(env_permdisp_pw11$`p adj` <= 0.001, "***",
                          ifelse(env_permdisp_pw11$`p adj` <= 0.01,  "**",
                          ifelse(env_permdisp_pw11$`p adj` <= 0.05,  "*", "ns")))
env_permdisp_pw11 <- env_permdisp_pw11[, c("pair", "diff", "lwr", "upr", "p adj", "sig")]
colnames(env_permdisp_pw11) <- c("pair", "diff_dispersion", "CI_lower", "CI_upper", "p_adj", "sig")
env_permdisp_pw11$diff_dispersion <- round(env_permdisp_pw11$diff_dispersion, 4)
env_permdisp_pw11$CI_lower        <- round(env_permdisp_pw11$CI_lower, 4)
env_permdisp_pw11$CI_upper        <- round(env_permdisp_pw11$CI_upper, 4)
env_permdisp_pw11$p_adj           <- round(env_permdisp_pw11$p_adj, 4)

cat("\nPairwise PERMDISP (env, 11 clusters):\n")
print(env_permdisp_pw11, row.names = FALSE)

# ==============================================================================
# CSV保存: 3グループ + 4グループの結果を1つのCSVにまとめる
# ==============================================================================

# 3グループ
summary_3 <- rbind(
  data.frame(groups = "3 (A/B/C)", test = "PERMANOVA (overall)", pair = "A/B/C",
             R2 = round(env_permanova_3$R2[1], 4), F_value = round(env_permanova_3$F[1], 2),
             p_value = env_permanova_3$`Pr(>F)`[1], sig = ifelse(env_permanova_3$`Pr(>F)`[1] <= 0.001, "***", ifelse(env_permanova_3$`Pr(>F)`[1] <= 0.01, "**", ifelse(env_permanova_3$`Pr(>F)`[1] <= 0.05, "*", "ns")))),
  data.frame(groups = "3 (A/B/C)", test = "PERMDISP (overall)", pair = "A/B/C",
             R2 = NA, F_value = round(env_permdisp_3_anova$`F value`[1], 2),
             p_value = round(env_permdisp_3_anova$`Pr(>F)`[1], 4), sig = ifelse(env_permdisp_3_anova$`Pr(>F)`[1] <= 0.001, "***", ifelse(env_permdisp_3_anova$`Pr(>F)`[1] <= 0.01, "**", ifelse(env_permdisp_3_anova$`Pr(>F)`[1] <= 0.05, "*", "ns")))),
  data.frame(groups = "3 (A/B/C)", test = "Pairwise PERMANOVA", pair = env_pairwise_3$pair,
             R2 = env_pairwise_3$R2, F_value = env_pairwise_3$F_value,
             p_value = env_pairwise_3$p_adj, sig = env_pairwise_3$sig),
  data.frame(groups = "3 (A/B/C)", test = "Pairwise PERMDISP", pair = c("A vs B", "A vs C", "B vs C"),
             R2 = NA, F_value = NA,
             p_value = env_permdisp_pw3$p_adj, sig = env_permdisp_pw3$sig)
)

# 4グループ
summary_4 <- rbind(
  data.frame(groups = "4 (A/B/C1/C2)", test = "PERMANOVA (overall)", pair = "A/B/C1/C2",
             R2 = round(env_permanova_4$R2[1], 4), F_value = round(env_permanova_4$F[1], 2),
             p_value = env_permanova_4$`Pr(>F)`[1], sig = ifelse(env_permanova_4$`Pr(>F)`[1] <= 0.001, "***", ifelse(env_permanova_4$`Pr(>F)`[1] <= 0.01, "**", ifelse(env_permanova_4$`Pr(>F)`[1] <= 0.05, "*", "ns")))),
  data.frame(groups = "4 (A/B/C1/C2)", test = "PERMDISP (overall)", pair = "A/B/C1/C2",
             R2 = NA, F_value = round(env_permdisp_4_anova$`F value`[1], 2),
             p_value = round(env_permdisp_4_anova$`Pr(>F)`[1], 4), sig = ifelse(env_permdisp_4_anova$`Pr(>F)`[1] <= 0.001, "***", ifelse(env_permdisp_4_anova$`Pr(>F)`[1] <= 0.01, "**", ifelse(env_permdisp_4_anova$`Pr(>F)`[1] <= 0.05, "*", "ns")))),
  data.frame(groups = "4 (A/B/C1/C2)", test = "Pairwise PERMANOVA", pair = env_pairwise_4$pair,
             R2 = env_pairwise_4$R2, F_value = env_pairwise_4$F_value,
             p_value = env_pairwise_4$p_adj, sig = env_pairwise_4$sig),
  data.frame(groups = "4 (A/B/C1/C2)", test = "Pairwise PERMDISP",
             pair = c("A vs B", "A vs C1", "A vs C2", "B vs C1", "B vs C2", "C1 vs C2"),
             R2 = NA, F_value = NA,
             p_value = env_permdisp_pw4$p_adj, sig = env_permdisp_pw4$sig)
)

# 11クラスター
pairs11_labels <- apply(combn(hier_levels_11, 2), 2, function(x) paste(x, collapse = " vs "))
summary_11 <- rbind(
  data.frame(groups = "11 clusters", test = "PERMANOVA (overall)",
             pair = paste(hier_levels_11, collapse = "/"),
             R2 = round(env_permanova_11$R2[1], 4),
             F_value = round(env_permanova_11$F[1], 2),
             p_value = env_permanova_11$`Pr(>F)`[1],
             sig = ifelse(env_permanova_11$`Pr(>F)`[1] <= 0.001, "***",
                   ifelse(env_permanova_11$`Pr(>F)`[1] <= 0.01, "**",
                   ifelse(env_permanova_11$`Pr(>F)`[1] <= 0.05, "*", "ns")))),
  data.frame(groups = "11 clusters", test = "PERMDISP (overall)",
             pair = paste(hier_levels_11, collapse = "/"),
             R2 = NA,
             F_value = round(env_permdisp_11_anova$`F value`[1], 2),
             p_value = round(env_permdisp_11_anova$`Pr(>F)`[1], 4),
             sig = ifelse(env_permdisp_11_anova$`Pr(>F)`[1] <= 0.001, "***",
                   ifelse(env_permdisp_11_anova$`Pr(>F)`[1] <= 0.01, "**",
                   ifelse(env_permdisp_11_anova$`Pr(>F)`[1] <= 0.05, "*", "ns")))),
  data.frame(groups = "11 clusters", test = "Pairwise PERMANOVA",
             pair = env_pairwise_11$pair,
             R2 = env_pairwise_11$R2, F_value = env_pairwise_11$F_value,
             p_value = env_pairwise_11$p_adj, sig = env_pairwise_11$sig),
  data.frame(groups = "11 clusters", test = "Pairwise PERMDISP",
             pair = env_permdisp_pw11$pair,
             R2 = NA, F_value = NA,
             p_value = env_permdisp_pw11$p_adj, sig = env_permdisp_pw11$sig)
)

env_summary <- rbind(summary_3, summary_4, summary_11)

write.csv(env_summary,
  here("output", "survey", "beta_diversity", "env_permanova_permdisp_summary.csv"),
  row.names = FALSE)

cat("\nDone: env PERMANOVA/PERMDISP results saved.\n")
cat("  CSV: output/survey/beta_diversity/env_permanova_permdisp_summary.csv\n")

# ==============================================================================
# Section: C階層的分割 PERMANOVA/PERMDISP (k=2,4,5,6 on each PC axis)
# Hierarchical subdivision of C (n=97) with per-PC-axis tests
# ==============================================================================

cat("\n--- C hierarchical PERMANOVA/PERMDISP ---\n")

hc_full <- as.hclust(h1$rowDendrogram)
c_samples <- names(clusnum11)[clusnum11 %in% c("C1a", "C1b1", "C1b2", "C2a", "C2b1", "C2b2")]
d_c <- as.dist(as.matrix(cophenetic(hc_full))[c_samples, c_samples])
hc_c <- hclust(d_c, method = "ward.D")

common_c <- intersect(c_samples, pca_all$sample_id)
pca_c <- pca_all[match(common_c, pca_all$sample_id), ]

all_results_c <- list()

for (k in c(2, 4, 5, 6)) {
  cuts <- cutree(hc_c, k = k)
  c_order <- hc_c$order
  ord <- unname(rle(cuts[hc_c$labels[c_order]])$values)
  remap <- setNames(1:k, ord)
  cuts <- remap[as.character(cuts)]
  names(cuts) <- names(cutree(hc_c, k = k))

  c_names <- paste0("C", 1:k)
  group_c <- factor(c_names[cuts[common_c]], levels = c_names)

  for (pc in paste0("PC", 1:9)) {
    pc_dist <- vegdist(as.matrix(pca_c[, pc, drop = FALSE]), method = "euclidean")

    res_perm <- adonis2(pc_dist ~ group_c, permutations = 999)
    all_results_c[[length(all_results_c) + 1]] <- data.frame(
      k = k, PC = pc, test = "PERMANOVA (overall)", pair = paste(c_names, collapse = "/"),
      R2 = round(res_perm$R2[1], 4), F_value = round(res_perm$F[1], 2),
      p_value = res_perm$`Pr(>F)`[1])

    disp <- betadisper(pc_dist, group = group_c)
    disp_anova <- anova(disp)
    all_results_c[[length(all_results_c) + 1]] <- data.frame(
      k = k, PC = pc, test = "PERMDISP (overall)", pair = paste(c_names, collapse = "/"),
      R2 = NA, F_value = round(disp_anova$`F value`[1], 2),
      p_value = round(disp_anova$`Pr(>F)`[1], 4))

    pairs_k <- combn(c_names, 2)
    pw_pvals <- c(); pw_rows <- list()
    for (i in seq_len(ncol(pairs_k))) {
      cl_a <- pairs_k[1, i]; cl_b <- pairs_k[2, i]
      idx_a <- which(group_c == cl_a); idx_b <- which(group_c == cl_b)
      idx_ab <- c(idx_a, idx_b)
      dist_sub <- as.dist(as.matrix(pc_dist)[idx_ab, idx_ab])
      group_sub <- factor(c(rep(cl_a, length(idx_a)), rep(cl_b, length(idx_b))))
      res_pair <- adonis2(dist_sub ~ group_sub, permutations = 999)
      pw_pvals <- c(pw_pvals, res_pair$`Pr(>F)`[1])
      pw_rows[[i]] <- data.frame(
        k = k, PC = pc, test = "Pairwise PERMANOVA",
        pair = paste0(cl_a, " vs ", cl_b),
        R2 = round(res_pair$R2[1], 4), F_value = round(res_pair$F[1], 2),
        p_value = res_pair$`Pr(>F)`[1])
    }
    pw_padj <- p.adjust(pw_pvals, method = "BH")
    for (i in seq_along(pw_rows)) {
      pw_rows[[i]]$p_value <- round(pw_padj[i], 4)
      all_results_c[[length(all_results_c) + 1]] <- pw_rows[[i]]
    }

    tukey <- TukeyHSD(disp)
    pw_disp <- as.data.frame(tukey$group)
    pair_labels <- apply(combn(c_names, 2), 2, function(x) paste(x, collapse = " vs "))
    for (j in seq_len(nrow(pw_disp))) {
      all_results_c[[length(all_results_c) + 1]] <- data.frame(
        k = k, PC = pc, test = "Pairwise PERMDISP", pair = pair_labels[j],
        R2 = NA, F_value = NA, p_value = round(pw_disp$`p adj`[j], 4))
    }
  }
}

result_c <- bind_rows(all_results_c)
result_c$sig <- ifelse(result_c$p_value <= 0.001, "***",
                ifelse(result_c$p_value <= 0.01, "**",
                ifelse(result_c$p_value <= 0.05, "*", "ns")))

write.csv(result_c,
  here("output", "survey", "beta_diversity", "env_permanova_permdisp_C_hierarchical.csv"),
  row.names = FALSE)

# デンドログラムノードごとのpairwise結果
node_pairs <- list(
  list(k = 2, pair = "C1 vs C2", node = "C1 vs C2"),
  list(k = 4, pair = "C1 vs C2", node = "C1a vs C1b"),
  list(k = 4, pair = "C3 vs C4", node = "C2a vs C2b"),
  list(k = 5, pair = "C4 vs C5", node = "C2b1 vs C2b2"),
  list(k = 6, pair = "C2 vs C3", node = "C1b1 vs C1b2")
)

node_results <- list()
for (np in node_pairs) {
  for (pc in paste0("PC", 1:9)) {
    perm <- result_c[result_c$k == np$k & result_c$PC == pc & result_c$test == "Pairwise PERMANOVA" & result_c$pair == np$pair, ]
    disp <- result_c[result_c$k == np$k & result_c$PC == pc & result_c$test == "Pairwise PERMDISP" & result_c$pair == np$pair, ]
    if (nrow(perm) == 0 | nrow(disp) == 0) next
    node_results[[length(node_results) + 1]] <- data.frame(
      node = np$node, k = np$k, PC = pc,
      permanova_R2 = perm$R2, permanova_p = perm$p_value, permanova_sig = perm$sig,
      permdisp_p = disp$p_value, permdisp_sig = disp$sig, stringsAsFactors = FALSE)
  }
}

node_df <- bind_rows(node_results)

# BH補正: 各ノード内でPC1〜PC9の9回検定に対して補正
node_df <- node_df %>%
  group_by(node) %>%
  mutate(
    permanova_p = round(p.adjust(permanova_p, method = "BH"), 4),
    permanova_sig = ifelse(permanova_p <= 0.001, "***",
                    ifelse(permanova_p <= 0.01,  "**",
                    ifelse(permanova_p <= 0.05,  "*", "ns")))
  ) %>%
  ungroup()

# 各ノードの最重要PC (BH補正後PERMANOVA sig & PERMDISP ns, max R²)
top_pcs <- node_df %>%
  filter(permanova_p <= 0.05, permdisp_sig == "ns") %>%
  group_by(node) %>%
  slice_max(order_by = permanova_R2, n = 1) %>%
  ungroup()

node_df$top_PC <- ""
for (i in seq_len(nrow(node_df))) {
  tp <- top_pcs[top_pcs$node == node_df$node[i] & top_pcs$PC == node_df$PC[i], ]
  if (nrow(tp) > 0) node_df$top_PC[i] <- "***TOP***"
}

write.csv(node_df,
  here("output", "survey", "beta_diversity", "env_C_dendrogram_node_pairwise.csv"),
  row.names = FALSE)

cat("Done: C hierarchical PERMANOVA/PERMDISP saved.\n")
cat("  CSV: output/survey/beta_diversity/env_permanova_permdisp_C_hierarchical.csv\n")
cat("  CSV: output/survey/beta_diversity/env_C_dendrogram_node_pairwise.csv\n")

# ==============================================================================
# Section: ノード×PC軸 ヒートマップ / Node x PC axis heatmap
# ==============================================================================

node_order <- c("C1 vs C2", "C1a vs C1b", "C2a vs C2b", "C2b1 vs C2b2", "C1b1 vs C1b2")

plot_df_node <- node_df %>%
  mutate(
    node = factor(node, levels = rev(node_order)),
    PC = factor(PC, levels = paste0("PC", 1:9)),
    sig_color = ifelse(permanova_sig == "ns", "ns", "sig"),
    label = ifelse(permanova_sig == "ns", "",
            paste0("R\u00B2=", sprintf("%.2f", permanova_R2), "\n", permanova_sig)),
    is_top = ifelse(top_PC == "***TOP***", TRUE, FALSE)
  )

pdf(here("output", "survey", "beta_diversity", "env_C_node_pairwise_heatmap.pdf"),
    width = 12, height = 5)

print(
  ggplot(plot_df_node, aes(x = PC, y = node)) +
    geom_tile(aes(fill = sig_color), color = "white", linewidth = 1) +
    scale_fill_manual(values = c("sig" = "#4292C6", "ns" = "grey90"), guide = "none") +
    geom_text(aes(label = label), size = 3, lineheight = 0.8) +
    geom_tile(data = plot_df_node %>% filter(is_top),
              aes(x = PC, y = node), fill = NA, color = "red", linewidth = 1.5) +
    labs(x = "PC axis", y = "Dendrogram node",
         subtitle = "Blue = significant (BH-adjusted p < 0.05), Red border = top PC") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank())
)

dev.off()

cat("  PDF: output/survey/beta_diversity/env_C_node_pairwise_heatmap.pdf\n")

# ==============================================================================
# Section: PCスコア地図 (PC1-PC9, depth_typeでfacet)
# Map of PC scores by depth type
# size = |PC score|, shape = sign (+/-), color = 11 clusters
# ==============================================================================

library(ggmap)
ggmap::register_stadiamaps(Sys.getenv("STADIA_MAPS_KEY"))

pca_map <- pca_all
pca_map$hier_name <- factor(
  as.character(clusnum11[pca_map$sample_id]),
  levels = hier_levels_11
)
pca_map$lat <- meta_asgard[pca_map$sample_id, "lat"]
pca_map$lon <- meta_asgard[pca_map$sample_id, "lon"]
pca_map$depth_type <- factor(
  meta_asgard[pca_map$sample_id, "depth_type"],
  levels = c("surf", "mid", "bottom")
)

pca_map <- pca_map %>% filter(!is.na(lat), !is.na(lon))

set.seed(42)
j <- 0.02
pca_map$lon_j <- pca_map$lon + runif(nrow(pca_map), -j, j)
pca_map$lat_j <- pca_map$lat + runif(nrow(pca_map), -j, j)

bbox_pc <- make_bbox(lon = pca_map$lon, lat = pca_map$lat, f = 0.1)
mapz_pc <- get_stadiamap(bbox_pc, maptype = "stamen_terrain", zoom = 4)

loadings_rot <- pca_result$rotation
ve_pct <- round(summary(pca_result)$importance["Proportion of Variance", ] * 100, 1)

pdf(here("output", "survey", "maps", "map_PC_scores.pdf"), width = 16, height = 10)

for (pc in paste0("PC", 1:9)) {
  pc_num <- as.integer(sub("PC", "", pc))

  pca_map$pc_val <- pca_map[[pc]]
  pca_map$pc_abs <- abs(pca_map$pc_val)
  pca_map$pc_sign <- ifelse(pca_map$pc_val >= 0, "positive", "negative")

  ld <- loadings_rot[, pc_num]
  top3 <- names(sort(abs(ld), decreasing = TRUE))[1:3]
  top3_str <- paste(sapply(top3, function(v) {
    sprintf("%s (%+.2f)", v, ld[v])
  }), collapse = ", ")

  print(ggmap(mapz_pc) +
    geom_point(data = pca_map,
               aes(x = lon_j, y = lat_j, color = hier_name,
                   size = pc_abs, shape = pc_sign),
               alpha = 0.8) +
    scale_color_manual(values = cc11, name = "Cluster") +
    scale_size_continuous(range = c(0.5, 7), name = paste0("|", pc, "|")) +
    scale_shape_manual(values = c("positive" = 16, "negative" = 17),
                       name = "Sign") +
    facet_grid(~ depth_type) +
    labs(title = paste0(pc, " scores (", ve_pct[pc_num], "% variance)"),
         subtitle = paste0("Top loadings: ", top3_str)) +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text = element_text(face = "bold", size = 13),
          panel.background = element_rect(fill = "grey85"),
          legend.position = "right"))
}

dev.off()

cat("Done: output/survey/maps/map_PC_scores.pdf\n")

# ==============================================================================
# Section: 11クラスター版 PCA biplot + 8変数まとめページ
# PCA biplot (11 clusters) + combined 8-variable page
# ==============================================================================

# Add Chl and NH4 to pca_scores for combined plot
pca_scores$`chl (ug/l)` <- meta_asgard[pca_scores$sample_id, "chl (ug/l)"]
pca_scores$`NH4(uM)`    <- meta_asgard[pca_scores$sample_id, "NH4(uM)"]

pdf(file = here("output", "survey", "beta_diversity", "ASGARD_pca_env_11clusters.pdf"),
    width = 10, height = 8)

# 8変数のみ表示（PO4, Sil除外、Chl追加だが矢印はPCA変数のみ）
biplot_vars <- c("Temp", "Salinity", "DO", "NO3", "NH4", "FlECO", "Depth")
pca_loadings_8 <- pca_loadings[pca_loadings$label %in% biplot_vars, ]

# Page 1: PCA biplot with arrows (8 variables)
print(
  ggplot() +
    geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = cluster),
               size = 3, alpha = 0.7) +
    scale_color_manual(values = cc11) +
    geom_segment(data = pca_loadings_8,
                 aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "black", linewidth = 0.8) +
    geom_text(data = pca_loadings_8,
              aes(x = PC1_plot * 1.15, y = PC2_plot * 1.15, label = label),
              size = 4, fontface = "bold") +
    coord_cartesian(xlim = range(c(pca_scores$PC1, pca_loadings_8$PC1_plot)) * 1.3,
                    ylim = range(c(pca_scores$PC2, pca_loadings_8$PC2_plot)) * 1.3) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Cluster",
         title = "PCA of Environmental Variables (11 clusters)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16))
)

# Page 2: Without arrows
print(
  ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = cc11) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         color = "Cluster",
         title = "PCA of Environmental Variables (points only)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16))
)

# Page 3: 8 variables combined in one page
pca_8vars <- list(
  list(var = "temp",                 label = "Temperature (\u00B0C)"),
  list(var = "salinity",             label = "Salinity (PSU)"),
  list(var = "DO",                   label = "DO (\u00B5mol/kg)"),
  list(var = "NO3(uM)",             label = "NO3 (\u00B5M)"),
  list(var = "FlECO-AFL(mg/m^3)",   label = "FlECO-AFL (mg/m\u00B3)"),
  list(var = "chl (ug/l)",          label = "Chl a (\u00B5g/l)"),
  list(var = "NH4(uM)",             label = "NH4 (\u00B5M)"),
  list(var = "depth_m",             label = "Sampling depth (m)")
)

pca_8_plots <- list()
for (i in seq_along(pca_8vars)) {
  ev <- pca_8vars[[i]]
  pca_8_plots[[i]] <- ggplot(pca_scores, aes(x = PC1, y = PC2,
                                              color = cluster, size = .data[[ev$var]])) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = cc11, guide = "none") +
    scale_size_continuous(range = c(0.5, 5), name = ev$label) +
    coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                    ylim = range(pca_scores$PC2) * 1.1) +
    labs(x = paste0("PC1 (", pc1_pct, "%)"),
         y = paste0("PC2 (", pc2_pct, "%)"),
         title = ev$label) +
    theme_minimal(base_size = 8) +
    theme(plot.title = element_text(face = "bold", size = 9),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
}

print(gridExtra::grid.arrange(grobs = pca_8_plots, ncol = 2))

# Pages 4-11: Individual 8 variables (full size)
for (ev in pca_8vars) {
  print(
    ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster, size = .data[[ev$var]])) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = cc11) +
      scale_size_continuous(range = c(1, 8), name = ev$label) +
      coord_cartesian(xlim = range(pca_scores$PC1) * 1.1,
                      ylim = range(pca_scores$PC2) * 1.1) +
      labs(x = paste0("PC1 (", pc1_pct, "%)"),
           y = paste0("PC2 (", pc2_pct, "%)"),
           color = "Cluster",
           title = paste0("PCA - dot size: ", ev$label)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 16))
  )
}

dev.off()

cat("Done: output/survey/beta_diversity/ASGARD_pca_env_11clusters.pdf\n")
