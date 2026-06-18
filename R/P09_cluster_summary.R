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

message("\nP09_cluster_summary.R: done.")
message("  CSV: output_p/cluster_summary/processing_cluster_summary.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_stations_long.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_x_station_counts.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_env_summary.csv")
