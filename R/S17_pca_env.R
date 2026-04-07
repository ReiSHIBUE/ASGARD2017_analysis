library(here); library(tidyverse); library(vegan)
source(here("R/00_setup.R")); source(here("R/S01_data_prep.R")); source(here("R/S02_heatmaps_16S.R"))

df <- meta_asgard
df$cluster <- factor(hier_names[as.character(clusnum10[rownames(df)])], levels = hier_levels)
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
  scale_color_manual(values = cc10) +
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
       title = "PCA of Environmental Variables (10 clusters)",
       subtitle = "log1p-transformed nutrients + scaled; 9 variables") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))
print(p1)

# Page 2: Without arrows (cleaner view)
p2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10) +
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
    scale_color_manual(values = cc10) +
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
    scale_color_manual(values = cc10) +
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
    scale_color_manual(values = cc10) +
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

cat("\nDone: output/survey/beta_diversity/ASGARD_pca_env_survey.pdf (20 pages)\n")

# ==============================================================================
# Section: 環境変数による PERMANOVA / PERMDISP (3グループ A/B/C)
# Environmental PERMANOVA / PERMDISP by 3 groups using Euclidean distance
# on log1p-transformed and scaled environmental variables
# ==============================================================================

cat("\n--- Environmental PERMANOVA / PERMDISP by 3 groups (A/B/C) ---\n")

# env_scaled は既に log1p + scale 済み (174 samples x 9 vars)
env_eucdist <- vegdist(env_scaled, method = "euclidean")

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

# CSV保存
write.csv(env_pairwise_3,
  here("output", "survey", "beta_diversity", "env_pairwise_permanova_3groups.csv"),
  row.names = FALSE)
write.csv(env_permdisp_pw3,
  here("output", "survey", "beta_diversity", "env_pairwise_permdisp_3groups.csv"),
  row.names = FALSE)

cat("\nDone: env PERMANOVA/PERMDISP results saved.\n")
