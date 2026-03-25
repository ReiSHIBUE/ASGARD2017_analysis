### S06_maps.R
### ASGARD 2017 Survey Site Analysis — Maps (10 clusters)
### 地図スクリプト（サーベイサイト、10クラスター）
###
### REQUIRES (from S01, S02):
###   meta_asgard      - metadata 181×45 (must have lat, lon, depth_type, station)
###   clusnum10        - 10-cluster assignments (length 181, from S02)
###   cc10             - 10-colour palette (from S02)
###
### PRODUCES:
###   mapz_survey      - ggmap basemap object for survey stations
###
### OUTPUT:
###   output/survey/maps/ASGARD_survey_map_10clusters.pdf
###   output/survey/maps/ASGARD_survey_map_10clusters_detail.pdf
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

lvls <- as.character(1:10)

a_map <- meta_asgard %>% filter(!is.na(lat), !is.na(lon))
a_map$cluster    <- factor(clusnum10[rownames(a_map)], levels = lvls)
a_map$depth_type <- factor(a_map$depth_type, levels = c("surf", "mid", "bottom"))
a_map$ID         <- seq_len(nrow(a_map))

n_per <- table(a_map$cluster)
clbl  <- paste0("Cluster ", names(n_per), " (n=", n_per, ")")
names(clbl) <- names(n_per)

bbox        <- make_bbox(lon = a_map$lon, lat = a_map$lat, f = 0.1)
mapz_survey <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom = 4)

# ==============================================================================
# Section 2: 基本マップ / Basic maps
# ==============================================================================

pdf(here::here("output", "survey", "maps", "ASGARD_survey_map_10clusters.pdf"),
    width = 16, height = 12)

# 全体マップ / Full map
ggmap(mapz_survey) +
  geom_point(data = a_map, aes(x = lon, y = lat, color = cluster), size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10, labels = clbl) +
  labs(title = "ASGARD 2017 Survey - 10 clusters", color = NULL) +
  theme(plot.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 10))

# クラスターでfacet / Faceted by cluster
ggmap(mapz_survey) +
  geom_point(data = a_map, aes(x = lon, y = lat, color = cluster), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = cc10, guide = "none") +
  facet_wrap(~ cluster, ncol = 5) +
  labs(title = "ASGARD 2017 Survey - 10 clusters (faceted)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 12))

dev.off()

# ==============================================================================
# Section 3: 詳細マップ / Detailed maps
# ==============================================================================

pdf(here::here("output", "survey", "maps", "ASGARD_survey_map_10clusters_detail.pdf"),
    width = 20, height = 20)

# station名ラベル付き / With station labels
ggmap(mapz_survey) +
  geom_point(data = a_map, aes(x = lon, y = lat, color = cluster), size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10, labels = clbl) +
  geom_text_repel(data = a_map, aes(x = lon, y = lat, label = station),
                  size = 2, max.overlaps = 20, alpha = 0.7) +
  labs(title = "ASGARD 2017 Survey - 10 clusters", color = NULL) +
  theme(plot.title = element_text(face = "bold", size = 16))

# depth_typeでfacet / Faceted by depth type
ggmap(mapz_survey) +
  geom_point(data = a_map, aes(x = lon, y = lat, color = cluster), size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10, labels = clbl) +
  geom_text_repel(data = a_map, aes(x = lon, y = lat, label = station),
                  size = 1.8, max.overlaps = 15, alpha = 0.6) +
  facet_wrap(~ depth_type) +
  labs(title = "ASGARD 2017 Survey - 10 clusters by depth type", color = NULL) +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 13))

# クラスターでfacet / Faceted by cluster
ggmap(mapz_survey) +
  geom_point(data = a_map, aes(x = lon, y = lat, color = cluster), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = cc10, guide = "none") +
  facet_wrap(~ cluster, ncol = 5) +
  labs(title = "ASGARD 2017 Survey - 10 clusters (faceted)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 12))

# depth_type × cluster グリッド / Grid
ggmap(mapz_survey) +
  geom_point(data = a_map, aes(x = lon, y = lat, color = cluster), size = 2, alpha = 0.8) +
  scale_color_manual(values = cc10, guide = "none") +
  facet_grid(depth_type ~ cluster) +
  labs(title = "ASGARD 2017 Survey - depth type x cluster") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 10))

dev.off()

message("S06_maps.R: done.")
message("  PDF: ASGARD_survey_map_10clusters.pdf, ASGARD_survey_map_10clusters_detail.pdf")
