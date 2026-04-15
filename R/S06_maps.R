### S06_maps.R
### ASGARD 2017 Survey Site Analysis — Maps (11 clusters)
### 地図スクリプト（サーベイサイト、11クラスター）
###
### REQUIRES (from S01, S02):
###   meta_asgard      - metadata 181×45 (must have lat, lon, depth_type, station)
###   clusnum11        - 11-cluster assignments (factor, from S02)
###   hier_levels_11   - ordered cluster names (from S02)
###   cc11             - 11-colour palette (from S02)
###
### PRODUCES:
###   mapz_survey      - ggmap basemap object for survey stations
###
### OUTPUT:
###   output/survey/maps/ASGARD_survey_map_11clusters.pdf
###   output/survey/maps/ASGARD_survey_map_11clusters_detail.pdf
###
### NOTE: requires ggmap (Stadia Maps). API key loaded from .Renviron:
###   STADIA_MAPS_KEY=your_key_here

library(tidyverse)
library(ggmap)
library(ggrepel)

ggmap::register_stadiamaps(Sys.getenv("STADIA_MAPS_KEY"))

# ==============================================================================
# Section 1: ベースマップの定義 / Define base map
# ==============================================================================

a_map <- meta_asgard %>% filter(!is.na(lat), !is.na(lon))
a_map$cluster11  <- factor(as.character(clusnum11[rownames(a_map)]),
                           levels = hier_levels_11)
a_map$depth_type <- factor(a_map$depth_type, levels = c("surf", "mid", "bottom"))
a_map$ID         <- seq_len(nrow(a_map))

n_per <- table(a_map$cluster11)
clbl  <- paste0(names(n_per), " (n=", n_per, ")")
names(clbl) <- names(n_per)

bbox        <- make_bbox(lon = a_map$lon, lat = a_map$lat, f = 0.1)
mapz_survey <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom = 4)

# ==============================================================================
# Section 2: 基本マップ / Basic maps
# ==============================================================================

dir.create(here::here("output", "survey", "maps"), showWarnings = FALSE, recursive = TRUE)

pdf(here::here("output", "survey", "maps", "ASGARD_survey_map_11clusters.pdf"),
    width = 16, height = 12)

# 全体マップ / Full map
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 3, alpha = 0.7) +
    scale_color_manual(values = cc11, labels = clbl) +
    labs(title = "ASGARD 2017 Survey - 11 clusters", color = NULL) +
    theme(plot.title = element_text(face = "bold", size = 16),
          legend.text = element_text(size = 10))
)

# クラスターでfacet / Faceted by cluster
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = cc11, guide = "none") +
    facet_wrap(~ cluster11, ncol = 4) +
    labs(title = "ASGARD 2017 Survey - 11 clusters (faceted)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 12))
)

dev.off()

# ==============================================================================
# Section 3: 詳細マップ / Detailed maps
# ==============================================================================

pdf(here::here("output", "survey", "maps", "ASGARD_survey_map_11clusters_detail.pdf"),
    width = 20, height = 20)

# station名ラベル付き / With station labels
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 3, alpha = 0.7) +
    scale_color_manual(values = cc11, labels = clbl) +
    geom_text_repel(data = a_map, aes(x = lon, y = lat, label = station),
                    size = 2, max.overlaps = 20, alpha = 0.7) +
    labs(title = "ASGARD 2017 Survey - 11 clusters", color = NULL) +
    theme(plot.title = element_text(face = "bold", size = 16))
)

# depth_typeでfacet / Faceted by depth type
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 3, alpha = 0.7) +
    scale_color_manual(values = cc11, labels = clbl) +
    geom_text_repel(data = a_map, aes(x = lon, y = lat, label = station),
                    size = 1.8, max.overlaps = 15, alpha = 0.6) +
    facet_wrap(~ depth_type) +
    labs(title = "ASGARD 2017 Survey - 11 clusters by depth type", color = NULL) +
    theme(plot.title = element_text(face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 13))
)

# クラスターでfacet / Faceted by cluster
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = cc11, guide = "none") +
    facet_wrap(~ cluster11, ncol = 4) +
    labs(title = "ASGARD 2017 Survey - 11 clusters (faceted)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 12))
)

# depth_type x cluster グリッド / Grid
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 2, alpha = 0.8) +
    scale_color_manual(values = cc11, guide = "none") +
    facet_grid(depth_type ~ cluster11) +
    labs(title = "ASGARD 2017 Survey - depth type x cluster") +
    theme(plot.title = element_text(face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 10))
)

dev.off()

# ==============================================================================
# Section 4: 地理サマリーCSV / Geographic summary CSV per cluster
# ==============================================================================

dt <- as.data.frame.matrix(table(a_map$cluster11, a_map$depth_type))
dt$cluster11 <- rownames(dt)

geo_summary <- a_map %>%
  group_by(cluster11) %>%
  summarise(
    n          = n(),
    n_stations = n_distinct(station),
    stations   = paste(sort(unique(station)), collapse = ", "),
    lat_min    = round(min(lat, na.rm = TRUE), 2),
    lat_max    = round(max(lat, na.rm = TRUE), 2),
    lon_min    = round(min(lon, na.rm = TRUE), 2),
    lon_max    = round(max(lon, na.rm = TRUE), 2),
    depth_m_min = round(min(depth_m, na.rm = TRUE), 1),
    depth_m_max = round(max(depth_m, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  left_join(dt, by = "cluster11") %>%
  mutate(
    surf_pct   = round(surf / n * 100, 1),
    mid_pct    = round(mid / n * 100, 1),
    bottom_pct = round(bottom / n * 100, 1)
  )

write.csv(geo_summary,
  here::here("output", "survey", "maps", "cluster11_geographic_summary.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 5: Depth type 検定 / Statistical test for depth type × cluster
# ==============================================================================

message("\n--- Depth type analysis (11 clusters) ---")

# Contingency table: cluster × depth_type
ct <- table(a_map$cluster11, a_map$depth_type)
message("\nContingency table (cluster x depth_type):")
print(ct)

# 全体検定: Fisher exact test (expected counts small for some cells)
fisher_all <- fisher.test(ct, simulate.p.value = TRUE, B = 9999)
message("\nFisher exact test (11 clusters x 3 depth types):")
message("  p = ", round(fisher_all$p.value, 4))

# Chi-squared test (参考)
chisq_all <- chisq.test(ct, simulate.p.value = TRUE, B = 9999)
message("Chi-squared test (simulated p):")
message("  chi-sq = ", round(chisq_all$statistic, 3), ", p = ", round(chisq_all$p.value, 4))

# 各クラスターの depth type 分布が全体分布と異なるかを検定
# (各クラスター vs 残り全部のFisher exact test)
message("\n--- Per-cluster Fisher exact test (cluster vs rest) ---")

depth_fisher <- data.frame(
  cluster = character(), p_value = numeric(), stringsAsFactors = FALSE
)

overall_dist <- colSums(ct)

for (cl in hier_levels_11) {
  cl_counts  <- ct[cl, ]
  rest_counts <- overall_dist - cl_counts
  mat2x3 <- rbind(cl_counts, rest_counts)
  ft <- fisher.test(mat2x3)
  depth_fisher <- rbind(depth_fisher, data.frame(
    cluster = cl, p_value = ft$p.value
  ))
}

depth_fisher$p_adj <- round(p.adjust(depth_fisher$p_value, method = "BH"), 4)
depth_fisher$p_value <- round(depth_fisher$p_value, 4)
depth_fisher$sig <- ifelse(depth_fisher$p_adj <= 0.001, "***",
                    ifelse(depth_fisher$p_adj <= 0.01,  "**",
                    ifelse(depth_fisher$p_adj <= 0.05,  "*", "ns")))

message("\nPer-cluster Fisher test (BH-adjusted):")
print(depth_fisher, row.names = FALSE)

# Depth type の割合テーブル
depth_pct <- as.data.frame.matrix(round(prop.table(ct, margin = 1) * 100, 1))
depth_pct$cluster11 <- rownames(depth_pct)
depth_pct <- depth_pct[, c("cluster11", "surf", "mid", "bottom")]
colnames(depth_pct) <- c("cluster11", "surf_pct", "mid_pct", "bottom_pct")

depth_result <- merge(depth_pct, depth_fisher, by.x = "cluster11", by.y = "cluster")
depth_result <- depth_result[match(hier_levels_11, depth_result$cluster11), ]

write.csv(depth_result,
  here::here("output", "survey", "maps", "cluster11_depth_type_test.csv"),
  row.names = FALSE)

message("\nDepth type results saved.")

message("\nS06_maps.R: done.")
message("  PDF: ASGARD_survey_map_11clusters.pdf, ASGARD_survey_map_11clusters_detail.pdf")
message("  CSV: output/survey/maps/cluster11_geographic_summary.csv")
message("  CSV: output/survey/maps/cluster11_depth_type_test.csv")
