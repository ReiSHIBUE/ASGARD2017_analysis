### S23_dbo3_envmap.R
### ASGARD 2017 Survey — Chl-a & NO3 maps at DBO3 stations
### DBO3 ステーションの Chl-a と NO3 濃度マップ
###
### REQUIRES (from 00_setup.R, S01, S02, S06):
###   meta_asgard      - metadata 181×45
###   clusnum11        - 11-cluster assignments
###   hier_levels_11   - ordered cluster names
###   cc11             - 11-colour palette
###   mapz_survey      - Stadia basemap (from S06)
###
### OUTPUT:
###   output/survey/maps/ASGARD_dbo3_envmap.pdf

library(tidyverse)
library(ggmap)
library(ggrepel)

# ==============================================================================
# Section 1: Filter to DBO3 samples
# ==============================================================================

m <- meta_asgard %>% rownames_to_column("Sample")
m$cluster11  <- factor(as.character(clusnum11[m$Sample]), levels = hier_levels_11)
m$depth_type <- factor(m$depth_type, levels = c("surf", "mid", "bottom"))

dbo3 <- m %>%
  filter(grepl("^DBO3", station)) %>%
  rename(NO3 = `NO3(uM)`, chl = `chl (ug/l)`)

cat("=== DBO3 samples ===\n")
cat("Total: ", nrow(dbo3), "\n")
print(table(dbo3$station, dbo3$depth_type, useNA = "always"))

# ==============================================================================
# Section 2: Per-station mean (across depths)
# ステーション単位で平均（depth_type の3層を統合）
# ==============================================================================

dbo3_station <- dbo3 %>%
  group_by(station) %>%
  summarise(
    lat        = mean(lat, na.rm = TRUE),
    lon        = mean(lon, na.rm = TRUE),
    n_samples  = n(),
    mean_NO3   = round(mean(NO3, na.rm = TRUE), 2),
    mean_chl   = round(mean(chl, na.rm = TRUE), 2),
    clusters   = paste(sort(unique(as.character(cluster11))), collapse = ", "),
    .groups    = "drop"
  )

cat("\n=== DBO3 station-level means ===\n")
print(as.data.frame(dbo3_station), row.names = FALSE)

# ==============================================================================
# Section 3: Per (station × depth_type)
# ==============================================================================

dbo3_depth <- dbo3 %>%
  group_by(station, depth_type) %>%
  summarise(lat       = mean(lat, na.rm = TRUE),
            lon       = mean(lon, na.rm = TRUE),
            mean_NO3  = round(mean(NO3, na.rm = TRUE), 2),
            mean_chl  = round(mean(chl, na.rm = TRUE), 2),
            n_samples = n(),
            .groups   = "drop") %>%
  filter(!is.na(depth_type))

# ==============================================================================
# Section 4: DBO3-focused basemap
# DBO3 ライン周辺だけにズーム
# ==============================================================================

bbox_dbo3 <- make_bbox(lon = dbo3_station$lon, lat = dbo3_station$lat, f = 0.3)
ggmap::register_stadiamaps(Sys.getenv("STADIA_MAPS_KEY"))
mapz_dbo3 <- get_stadiamap(bbox_dbo3, maptype = "stamen_terrain", zoom = 7)

# ==============================================================================
# Section 5: Map plots
# ==============================================================================

dir.create(here::here("output", "survey", "maps"), showWarnings = FALSE, recursive = TRUE)

pdf(here::here("output", "survey", "maps", "ASGARD_dbo3_envmap.pdf"),
    width = 16, height = 10)

# Page 1: per-station NO3 (depth-averaged)
print(
  ggmap(mapz_dbo3) +
    geom_point(data = dbo3_station,
               aes(x = lon, y = lat, size = mean_NO3),
               color = "#1f77b4", alpha = 0.7) +
    geom_text_repel(data = dbo3_station,
                    aes(x = lon, y = lat,
                        label = paste0(station, "\n", mean_NO3, " µM")),
                    size = 3, max.overlaps = 20) +
    scale_size_continuous(range = c(2, 14), name = "NO3 (µM)") +
    labs(title = "DBO3 station mean NO3 (depth-averaged)",
         x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(face = "bold", size = 14))
)

# Page 2: per-station Chl-a (depth-averaged)
print(
  ggmap(mapz_dbo3) +
    geom_point(data = dbo3_station,
               aes(x = lon, y = lat, size = mean_chl),
               color = "#2ca02c", alpha = 0.7) +
    geom_text_repel(data = dbo3_station,
                    aes(x = lon, y = lat,
                        label = paste0(station, "\n", mean_chl, " µg/L")),
                    size = 3, max.overlaps = 20) +
    scale_size_continuous(range = c(2, 14), name = "Chl-a (µg/L)") +
    labs(title = "DBO3 station mean Chl-a (depth-averaged)",
         x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(face = "bold", size = 14))
)

# Page 2: per (station × depth_type) NO3 — faceted by depth
print(
  ggmap(mapz_dbo3) +
    geom_point(data = dbo3_depth,
               aes(x = lon, y = lat, size = mean_NO3),
               color = "#1f77b4", alpha = 0.7) +
    geom_text_repel(data = dbo3_depth,
                    aes(x = lon, y = lat, label = paste0(station, ": ", mean_NO3)),
                    size = 2.5, max.overlaps = 25) +
    scale_size_continuous(range = c(2, 10), name = "NO3 (µM)") +
    facet_wrap(~ depth_type, nrow = 1) +
    labs(title = "DBO3 NO3 by depth type",
         x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12))
)

# Page 3: per (station × depth_type) Chl-a
print(
  ggmap(mapz_dbo3) +
    geom_point(data = dbo3_depth,
               aes(x = lon, y = lat, size = mean_chl),
               color = "#2ca02c", alpha = 0.7) +
    geom_text_repel(data = dbo3_depth,
                    aes(x = lon, y = lat, label = paste0(station, ": ", mean_chl)),
                    size = 2.5, max.overlaps = 25) +
    scale_size_continuous(range = c(2, 10), name = "Chl-a (µg/L)") +
    facet_wrap(~ depth_type, nrow = 1) +
    labs(title = "DBO3 Chl-a by depth type",
         x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12))
)

dev.off()

write.csv(dbo3_station,
  here::here("output", "survey", "maps", "ASGARD_dbo3_envmap_station.csv"),
  row.names = FALSE)
write.csv(dbo3_depth,
  here::here("output", "survey", "maps", "ASGARD_dbo3_envmap_depth.csv"),
  row.names = FALSE)

message("\nS23_dbo3_envmap.R: done.")
message("  PDF: output/survey/maps/ASGARD_dbo3_envmap.pdf")
message("  CSV: output/survey/maps/ASGARD_dbo3_envmap_station.csv")
message("  CSV: output/survey/maps/ASGARD_dbo3_envmap_depth.csv")
