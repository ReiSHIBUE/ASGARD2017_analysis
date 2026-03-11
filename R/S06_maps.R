### S06_maps.R
### ASGARD 2017 Survey Site Analysis — Maps
### 地図スクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02):
###   meta_asgard      - metadata 181×45 (must have lat, lon, depth_type, station)
###   clusnum          - cluster assignments (length 181)
###   asgard_pcoa_df   - PCoA + metadata df (from S04, for lat/lon join)
###
### PRODUCES:
###   mapz_survey      - ggmap basemap object for survey stations
###
### OUTPUT:
###   output/survey/maps/ASGARD_survey_map.pdf
###
### NOTE: requires ggmap (Stadia Maps). API key loaded from .Renviron:
###   STADIA_MAPS_KEY=your_key_here
###   Bug fix: hardcoded API key in original script (line 1267) replaced with
###   Sys.getenv("STADIA_MAPS_KEY").

library(tidyverse)
library(ggmap)
library(ggrepel)

ggmap::register_stadiamaps(Sys.getenv("STADIA_MAPS_KEY"))

# ==============================================================================
# Section 1: ベースマップの定義 / Define base map
# ==============================================================================

a_map <- meta_asgard
a_map$cluster    <- as.factor(clusnum)
a_map$depth_type <- factor(a_map$depth_type, levels = c("surf", "mid", "bottom"))
a_map$ID         <- seq_len(nrow(a_map))

bbox     <- make_bbox(lon = a_map$lon, lat = a_map$lat, f = 0.1)
mapz_survey <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom = 4)

# ==============================================================================
# Section 2: サンプルの地図 / Sample maps coloured by cluster
# ==============================================================================

pdf(here::here("output", "survey", "maps", "ASGARD_survey_map.pdf"),
    width = 20, height = 10)

# クラスターと水深でfacet / Faceted by depth type and cluster
map_plot_1 <- ggmap(mapz_survey) +
  geom_point(data = a_map,
             aes(x = lon, y = lat, color = cluster),
             size = 3, alpha = 1) +
  geom_text_repel(data = a_map,
                  aes(x = lon, y = lat, label = ID),
                  size = 3, max.overlaps = Inf) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ depth_type) +
  labs(color = "cluster", title = "ASGARD 2017 Survey Stations by Depth Type")

print(map_plot_1)

# クラスターのみでfacet / Faceted by cluster only
map_plot_2 <- ggmap(mapz_survey) +
  geom_point(data = a_map,
             aes(x = lon, y = lat, color = cluster),
             size = 3, alpha = 1) +
  geom_text_repel(data = a_map,
                  aes(x = lon, y = lat, label = ID),
                  size = 3, max.overlaps = Inf) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(depth_type ~ cluster) +
  labs(color = "cluster", title = "ASGARD 2017 Survey Stations by Cluster × Depth")

print(map_plot_2)

dev.off()

message("S06_maps.R: done. PDF: output/survey/maps/ASGARD_survey_map.pdf")
