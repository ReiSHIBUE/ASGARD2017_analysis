### P09_cluster_summary.R
### ASGARD 2017 Processing — Per-cluster geographic + environmental summary
### k=4クラスタごとの地理・環境特徴まとめ
###
### REQUIRES (from P01–P03):
###   meta_asgard_p2  - metadata for 78 samples
###   clusnum_p       - k=4 cluster assignment
###
### OUTPUT:
###   output_p/cluster_summary/processing_cluster_summary.csv
###   output_p/cluster_summary/processing_cluster_stations_long.csv
###   output_p/cluster_summary/processing_cluster_x_station_counts.csv
###   output_p/cluster_summary/processing_cluster_env_summary.csv

library(dplyr)
library(tidyr)
library(here)

dir.create(here("output_p", "cluster_summary"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: メタデータ + クラスタ + 海域 (Bering vs Chukchi, lat >= 66°N)
# ==============================================================================

meta_p <- meta_asgard_p2[names(clusnum_p), ]
meta_p$cluster <- as.integer(clusnum_p)
meta_p$sea <- ifelse(meta_p$lat >= 66, "Chukchi",
                     ifelse(is.na(meta_p$lat), NA, "Bering"))

# ==============================================================================
# Section 2: per (cluster × station) long format
# ==============================================================================

long_df <- meta_p %>%
  group_by(cluster, station) %>%
  summarise(sea       = first(sea),
            lat       = round(mean(lat, na.rm = TRUE), 3),
            lon       = round(mean(lon, na.rm = TRUE), 3),
            n_samples = n(),
            n_surf    = sum(depth_type == "surf",   na.rm = TRUE),
            n_mid     = sum(depth_type == "mid",    na.rm = TRUE),
            n_bottom  = sum(depth_type == "bottom", na.rm = TRUE),
            n_0.2um   = sum(filter == "0.2 µm", na.rm = TRUE),
            n_3um     = sum(filter == "3 µm",   na.rm = TRUE),
            n_20um    = sum(filter == "20 µm",  na.rm = TRUE),
            .groups   = "drop") %>%
  arrange(cluster, station)

write.csv(long_df,
          here("output_p", "cluster_summary",
               "processing_cluster_stations_long.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 3: cluster × station wide cross-tab
# ==============================================================================

wide_df <- meta_p %>%
  count(cluster, station) %>%
  pivot_wider(names_from = station, values_from = n, values_fill = 0) %>%
  arrange(cluster)

write.csv(wide_df,
          here("output_p", "cluster_summary",
               "processing_cluster_x_station_counts.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 4: per-cluster summary (one row per cluster)
# ==============================================================================

summary_df <- meta_p %>%
  group_by(cluster) %>%
  summarise(
    n_total      = n(),
    n_stations   = n_distinct(station),
    stations     = paste(sort(unique(station)), collapse = ", "),
    n_Bering     = sum(sea == "Bering",  na.rm = TRUE),
    n_Chukchi    = sum(sea == "Chukchi", na.rm = TRUE),
    n_surf       = sum(depth_type == "surf",   na.rm = TRUE),
    n_mid        = sum(depth_type == "mid",    na.rm = TRUE),
    n_bottom     = sum(depth_type == "bottom", na.rm = TRUE),
    n_0.2um      = sum(filter == "0.2 µm", na.rm = TRUE),
    n_3um        = sum(filter == "3 µm",   na.rm = TRUE),
    n_20um       = sum(filter == "20 µm",  na.rm = TRUE),
    mean_lat     = round(mean(lat,     na.rm = TRUE), 2),
    mean_lon     = round(mean(lon,     na.rm = TRUE), 2),
    mean_depth_m = round(mean(depth_m, na.rm = TRUE), 1),
    .groups      = "drop"
  )

write.csv(summary_df,
          here("output_p", "cluster_summary",
               "processing_cluster_summary.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 5: per-cluster environmental summary (mean + median + IQR)
# ==============================================================================

env_vars_p <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)",
                "Sil(uM)", "NH4(uM)", "chl (ug/l)", "FlECO-AFL(mg/m^3)",
                "depth_m", "lat", "lon")

env_long <- meta_p %>%
  select(cluster, all_of(env_vars_p)) %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "value")

env_summary <- env_long %>%
  group_by(cluster, variable) %>%
  summarise(n      = sum(!is.na(value)),
            mean   = round(mean(value,   na.rm = TRUE), 3),
            median = round(median(value, na.rm = TRUE), 3),
            sd     = round(sd(value,     na.rm = TRUE), 3),
            q1     = round(quantile(value, 0.25, na.rm = TRUE), 3),
            q3     = round(quantile(value, 0.75, na.rm = TRUE), 3),
            min    = round(min(value,    na.rm = TRUE), 3),
            max    = round(max(value,    na.rm = TRUE), 3),
            .groups = "drop") %>%
  arrange(cluster, variable)

write.csv(env_summary,
          here("output_p", "cluster_summary",
               "processing_cluster_env_summary.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 6: PERMANOVA + PERMDISP (overall + pairwise) on k=4 clusters
# Bray–Curtis on fourth-root proportions, 999 permutations
# ==============================================================================

library(vegan)

cluster_p <- factor(clusnum_p[rownames(asgard_filtered_p_hm2)],
                    levels = c("1", "2", "3", "4"))
bray_p    <- vegdist((asgard_filtered_p_hm2) ^ 0.25, method = "bray")

set.seed(42)

# --- Overall PERMANOVA ---
overall_perm <- adonis2(bray_p ~ cluster_p, permutations = 999)
overall_perm_df <- data.frame(
  term     = rownames(overall_perm),
  Df       = overall_perm$Df,
  SumOfSqs = round(overall_perm$SumOfSqs, 4),
  R2       = round(overall_perm$R2, 4),
  F        = round(overall_perm$F, 3),
  p_value  = overall_perm$`Pr(>F)`
)
write.csv(overall_perm_df,
          here("output_p", "cluster_summary",
               "processing_permanova_overall.csv"),
          row.names = FALSE)

# --- Overall PERMDISP ---
overall_disp <- betadisper(bray_p, group = cluster_p)
disp_anova   <- anova(overall_disp)
overall_disp_df <- data.frame(
  term     = rownames(disp_anova),
  Df       = disp_anova$Df,
  SumSq    = round(disp_anova$`Sum Sq`, 4),
  MeanSq   = round(disp_anova$`Mean Sq`, 4),
  F        = round(disp_anova$`F value`, 3),
  p_value  = disp_anova$`Pr(>F)`
)
write.csv(overall_disp_df,
          here("output_p", "cluster_summary",
               "processing_permdisp_overall.csv"),
          row.names = FALSE)

# --- Pairwise PERMANOVA (6 pairs) ---
cl_levels  <- levels(cluster_p)
pairs_p    <- combn(cl_levels, 2)
pairwise_perm <- data.frame(pair = character(), F = numeric(),
                            R2 = numeric(), p_value = numeric())
mat_bray_p <- as.matrix(bray_p)

for (i in seq_len(ncol(pairs_p))) {
  a <- pairs_p[1, i]; b <- pairs_p[2, i]
  smp <- names(cluster_p)[cluster_p %in% c(a, b)]
  d   <- as.dist(mat_bray_p[smp, smp])
  g   <- droplevels(cluster_p[smp])
  res <- adonis2(d ~ g, permutations = 999)
  pairwise_perm <- rbind(pairwise_perm, data.frame(
    pair    = paste0(a, " vs ", b),
    F       = round(res$F[1], 3),
    R2      = round(res$R2[1], 4),
    p_value = res$`Pr(>F)`[1]
  ))
}
pairwise_perm$p_adj_BH <- round(p.adjust(pairwise_perm$p_value, "BH"), 4)
pairwise_perm$sig <- ifelse(pairwise_perm$p_adj_BH <= 0.001, "***",
                     ifelse(pairwise_perm$p_adj_BH <= 0.01,  "**",
                     ifelse(pairwise_perm$p_adj_BH <= 0.05,  "*", "ns")))
write.csv(pairwise_perm,
          here("output_p", "cluster_summary",
               "processing_permanova_pairwise.csv"),
          row.names = FALSE)

# --- Pairwise PERMDISP (TukeyHSD on betadisper) ---
tk_disp <- TukeyHSD(overall_disp)
pairwise_disp <- as.data.frame(tk_disp$group)
pairwise_disp$pair <- rownames(pairwise_disp)
pairwise_disp <- pairwise_disp[, c("pair", "diff", "lwr", "upr", "p adj")]
colnames(pairwise_disp) <- c("pair", "diff_dispersion",
                             "CI_lower", "CI_upper", "p_adj")
pairwise_disp$diff_dispersion <- round(pairwise_disp$diff_dispersion, 4)
pairwise_disp$CI_lower        <- round(pairwise_disp$CI_lower, 4)
pairwise_disp$CI_upper        <- round(pairwise_disp$CI_upper, 4)
pairwise_disp$p_adj           <- round(pairwise_disp$p_adj, 4)
pairwise_disp$sig <- ifelse(pairwise_disp$p_adj <= 0.001, "***",
                     ifelse(pairwise_disp$p_adj <= 0.01,  "**",
                     ifelse(pairwise_disp$p_adj <= 0.05,  "*", "ns")))
write.csv(pairwise_disp,
          here("output_p", "cluster_summary",
               "processing_permdisp_pairwise.csv"),
          row.names = FALSE)

message("\n--- PERMANOVA / PERMDISP ---")
message(sprintf("Overall PERMANOVA: F = %.2f, R2 = %.3f, p = %.4f",
                overall_perm_df$F[1], overall_perm_df$R2[1],
                overall_perm_df$p_value[1]))
message(sprintf("Overall PERMDISP:  F = %.2f, p = %.4f",
                overall_disp_df$F[1], overall_disp_df$p_value[1]))
message(sprintf("Pairwise PERMANOVA sig (BH<0.05): %d / %d",
                sum(pairwise_perm$sig != "ns"), nrow(pairwise_perm)))
message(sprintf("Pairwise PERMDISP sig  (Tukey<0.05): %d / %d",
                sum(pairwise_disp$sig != "ns"), nrow(pairwise_disp)))

message("\nP09_cluster_summary.R: done.")
message("  CSV: output_p/cluster_summary/processing_cluster_summary.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_stations_long.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_x_station_counts.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_env_summary.csv")
message("  CSV: output_p/cluster_summary/processing_permanova_overall.csv")
message("  CSV: output_p/cluster_summary/processing_permdisp_overall.csv")
message("  CSV: output_p/cluster_summary/processing_permanova_pairwise.csv")
message("  CSV: output_p/cluster_summary/processing_permdisp_pairwise.csv")
