### 04_maps.R
### ASGARD 2017 Processing Site Analysis — Maps
### 地図スクリプト
###
### REQUIRES (from 01–03):
###   asgard_processing2       - 78-sample merged df
###   meta_asgard_p2           - metadata for 78 samples
###   asgard_filtered_p_hm2    - ASV proportion matrix 78*221
###   clusnum_p                - cluster assignment vector (length 78)
###   zero_cols                - ASV names absent in 0.2 µm fraction
###
### PRODUCES:
###   mapz                     - ggmap basemap object
###   asgard_processing_ggmap  - map-ready df (78*226)
###   esv_zero_only            - 0.2 µm-absent ASVs with lat/lon/filter (78*52)
###
### OUTPUT:
###   output_p/processing_map.pdf
###   output_p/maps_pa.pdf
###
### NOTE: requires ggmap (uses Stadia Maps). Set API key if needed:
###   ggmap::register_stadiamaps("YOUR_API_KEY")

library(tidyverse)
library(ggmap)
library(ggrepel)

ggmap::register_stadiamaps(Sys.getenv("STADIA_MAPS_KEY"))

# ==============================================================================
# Section 1: ベースマップの定義 / Define base map
# ==============================================================================

asgard_processing_ggmap <- asgard_processing2 %>%
  select(filter, lat, lon, contains("ESV")) # 78*226

bbox <- make_bbox(
  lon = asgard_processing_ggmap$lon,
  lat = asgard_processing_ggmap$lat,
  f   = 0.1
)

mapz <- get_stadiamap(bbox, maptype = "stamen_terrain", zoom = 4)

asgard_processing_ggmap$filter <- factor(
  asgard_processing_ggmap$filter,
  levels = c("0.2 µm", "3 µm", "20 µm")
)

# ==============================================================================
# Section 2: サンプルの地図 (クラスター色付き) / Processing station sample map
# ==============================================================================

# サンプル番号を付与 / Assign sequential IDs
a <- meta_asgard_p2
a$ID <- 1:nrow(a)
a$cluster <- as.factor(clusnum_p)
a$depth_type <- factor(a$depth_type, levels = c("surf", "mid", "bottom"))
a <- a %>% filter(Sample != "BOX_6_26") # 77*48

pdf(here::here("output_p", "maps", "processing_map.pdf"), width = 20, height = 20)

# クラスター × 水深で分割 / Faceted by depth and cluster
map_plot_3 <- ggmap(mapz) +
  geom_point(data = a,
             aes(x = lon, y = lat, color = cluster),
             alpha = 1, size = 3) +
  geom_text_repel(data = a,
                  aes(x = lon, y = lat, label = ID),
                  size = 3, max.overlaps = Inf) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(depth_type ~ cluster) +
  labs(color = "cluster")

print(map_plot_3)

# クラスターのみで分割 / Faceted by cluster only
map_plot_4 <- ggmap(mapz) +
  geom_point(data = a,
             aes(x = lon, y = lat, color = cluster),
             alpha = 1, size = 3) +
  geom_text_repel(data = a,
                  aes(x = lon, y = lat, label = ID),
                  size = 3, max.overlaps = Inf) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ cluster) +
  labs(color = "cluster")

print(map_plot_4)

dev.off()

# ==============================================================================
# Section 3: Particle-associated ASVの地図 (全て0のASVのみ)
# Per-ESV presence/absence maps for ASVs absent in 0.2 µm
# ==============================================================================

esv <- asgard_processing_ggmap %>%
  select(lat, lon, filter, contains("ESV_")) # 78*224

# 0.2 µm で全て0だったASVのみ抽出 / Subset to ASVs absent in free-living fraction
esv_zero_only <- asgard_processing_ggmap %>%
  select(lat, lon, filter, any_of(zero_cols)) # 78*52

pdf(here::here("output_p", "maps", "maps_pa.pdf"), width = 8, height = 6)

asv_cols <- colnames(esv_zero_only)[4:ncol(esv_zero_only)]

for (asv in asv_cols) {
  map_plot <- ggmap(mapz) +
    geom_point(data = esv_zero_only,
               aes(x = lon, y = lat, size = .data[[asv]]),
               alpha = 1) +
    theme(
      text                = element_text(size = 14),
      axis.text.x         = element_text(angle = 45, hjust = 1),
      legend.position     = "right",
      legend.key.size     = unit(0.5, "cm"),
      legend.text         = element_text(size = 2),
      legend.title        = element_text(size = 2)
    ) +
    facet_grid(. ~ filter)
  print(map_plot)
}

dev.off()

# ==============================================================================
# Section 4: S06-style hero maps (overall + faceted) — 4-cluster main figure
# ==============================================================================

a_map_p <- meta_asgard_p2
a_map_p$cluster    <- factor(as.character(clusnum_p[rownames(a_map_p)]),
                             levels = c("1", "2", "3", "4"))
a_map_p$depth_type <- factor(a_map_p$depth_type,
                             levels = c("surf", "mid", "bottom"))
a_map_p$division2  <- factor(ifelse(a_map_p$cluster == "1",
                                    "Free-living", "Particle-associated"),
                             levels = c("Free-living", "Particle-associated"))
a_map_p <- a_map_p[!is.na(a_map_p$lat) & !is.na(a_map_p$lon), ]

n_per_p     <- table(a_map_p$cluster)
clbl_p      <- paste0(names(n_per_p), " (n=", n_per_p, ")")
names(clbl_p) <- names(n_per_p)
n_per_d2    <- table(a_map_p$division2)
clbl_d2     <- paste0(names(n_per_d2), " (n=", n_per_d2, ")")
names(clbl_d2) <- names(n_per_d2)

cc_p   <- c("1" = "#E41A1C", "2" = "#377EB8",
            "3" = "#4DAF4A", "4" = "#984EA3")
cc_d2  <- c("Free-living" = "#E41A1C",
            "Particle-associated" = "#377EB8")

# ==============================================================================
# Section 5: Detail maps (station labels + depth facets) — main paper figure
# ==============================================================================

pdf(here::here("output_p", "maps",
               "ASGARD_processing_map_4clusters_detail.pdf"),
    width = 20, height = 20)

detail_theme_p <- theme(
  plot.title   = element_text(face = "bold", size = 20),
  axis.title   = element_text(size = 18, face = "bold"),
  axis.text    = element_text(size = 14),
  legend.title = element_text(size = 16, face = "bold"),
  legend.text  = element_text(size = 14),
  strip.text   = element_text(face = "bold", size = 18))

# With station labels
print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = cluster),
               size = 3.5, alpha = 0.7) +
    scale_color_manual(values = cc_p, labels = clbl_p) +
    geom_text_repel(data = a_map_p, aes(x = lon, y = lat, label = station),
                    size = 3, max.overlaps = 20, alpha = 0.8) +
    labs(title = "ASGARD 2017 Processing - 4 clusters", color = NULL) +
    detail_theme_p
)

# By depth_type
print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = cluster),
               size = 3, alpha = 0.7) +
    scale_color_manual(values = cc_p, labels = clbl_p) +
    geom_text_repel(data = a_map_p, aes(x = lon, y = lat, label = station),
                    size = 2.5, max.overlaps = 15, alpha = 0.7) +
    facet_wrap(~ depth_type) +
    labs(title = "ASGARD 2017 Processing - 4 clusters by depth type",
         color = NULL) +
    detail_theme_p
)

# depth_type x cluster grid
print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = cluster),
               size = 2, alpha = 0.8) +
    scale_color_manual(values = cc_p, guide = "none") +
    facet_grid(depth_type ~ cluster) +
    detail_theme_p +
    theme(strip.text = element_text(face = "bold", size = 13))
)

dev.off()

message("04_maps.R: done. PDFs written:")
message("  output_p/maps/processing_map.pdf (legacy)")
message("  output_p/maps/maps_pa.pdf (per-ASV PA maps)")
message("  output_p/maps/ASGARD_processing_map_4clusters_detail.pdf (S06-style)")
