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
# Section 4b: Hero overview maps (no station labels) — supplementary figure
# Page 2 title intentionally removed for slide / manuscript use
# ==============================================================================

pdf(here::here("output_p", "maps",
               "ASGARD_processing_map_4clusters.pdf"),
    width = 16, height = 12)

# Page 1: Full map, all clusters
print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = cluster),
               size = 3, alpha = 0.7) +
    scale_color_manual(values = cc_p, labels = clbl_p) +
    labs(title = "ASGARD 2017 Processing - 4 clusters", color = NULL) +
    theme(plot.title  = element_text(face = "bold", size = 16),
          legend.text = element_text(size = 11))
)

# Page 2: Faceted by cluster (no title)
print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = cluster),
               size = 2.5, alpha = 0.8) +
    scale_color_manual(values = cc_p, guide = "none") +
    facet_wrap(~ cluster, ncol = 4) +
    theme(strip.text = element_text(face = "bold", size = 13))
)

# Page 3: k=2 division (Free-living vs Particle-associated)
print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = division2),
               size = 3, alpha = 0.7) +
    scale_color_manual(values = cc_d2, labels = clbl_d2) +
    facet_wrap(~ division2) +
    labs(color = NULL) +
    theme(strip.text = element_text(face = "bold", size = 13))
)

dev.off()

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

# ==============================================================================
# Section 6: division x depth_type sample maps and ASV community contribution
# Equivalent of S06 Page 5 & 6 for processing (k=2: Free-living vs Particle-associated)
# ==============================================================================

# Page A: division2 x depth_type grid of sample locations
div2_colors <- c("Free-living" = "#E41A1C",
                 "Particle-associated" = "#377EB8")

pdf(here::here("output_p", "maps",
               "ASGARD_processing_division_depth_map.pdf"),
    width = 16, height = 14)

print(
  ggmap(mapz) +
    geom_point(data = a_map_p,
               aes(x = lon, y = lat, color = division2),
               size = 2.5, alpha = 0.8) +
    scale_color_manual(values = div2_colors, name = "Division") +
    facet_grid(depth_type ~ division2) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text  = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12))
)

# Page B: Per-division ASV community contribution x depth_type
# Assign each ASV to the division (1 vs 2-4) in which its mean RA is higher
mat_div_p   <- sweep(asgard_filtered_p_hm2, 1,
                     rowSums(asgard_filtered_p_hm2), "/")
fl_samples_p <- names(clusnum_p)[clusnum_p == 1]
pa_samples_p <- names(clusnum_p)[clusnum_p != 1]
asv_mean_fl <- colMeans(mat_div_p[fl_samples_p, , drop = FALSE])
asv_mean_pa <- colMeans(mat_div_p[pa_samples_p, , drop = FALSE])
asv_div_assign <- ifelse(asv_mean_fl >= asv_mean_pa,
                         "Free-living", "Particle-associated")
asv_div_f <- factor(asv_div_assign,
                    levels = c("Free-living", "Particle-associated"))
n_asv_div_p <- table(asv_div_f)

comp_RA_mat_p <- sapply(c("Free-living", "Particle-associated"), function(d) {
  rowSums(mat_div_p[, asv_div_f == d, drop = FALSE])
})

comp_df_p <- as.data.frame(comp_RA_mat_p) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(cluster    = factor(as.character(clusnum_p[Sample]),
                                    levels = c("1", "2", "3", "4")),
                lat        = meta_asgard_p2[Sample, "lat"],
                lon        = meta_asgard_p2[Sample, "lon"],
                depth_type = meta_asgard_p2[Sample, "depth_type"]) %>%
  dplyr::filter(!is.na(lat), !is.na(lon), !is.na(depth_type)) %>%
  tidyr::pivot_longer(c("Free-living", "Particle-associated"),
                      names_to = "Division", values_to = "RA") %>%
  dplyr::mutate(RA_pct     = RA * 100,
                Division   = factor(Division,
                                    levels = c("Free-living",
                                               "Particle-associated")),
                depth_type = factor(depth_type,
                                    levels = c("surf", "mid", "bottom")))

comp_df_p$Division_lab <- factor(
  paste0(as.character(comp_df_p$Division),
         " (", n_asv_div_p[as.character(comp_df_p$Division)], " ASVs)"),
  levels = paste0(c("Free-living", "Particle-associated"),
                  " (", n_asv_div_p[c("Free-living",
                                      "Particle-associated")], " ASVs)")
)

print(
  ggmap(mapz) +
    geom_point(data = comp_df_p,
               aes(x = lon, y = lat,
                   size = RA_pct, color = cluster),
               alpha = 0.75) +
    scale_color_manual(values = cc_p, name = "Cluster") +
    scale_size_continuous(range = c(1, 12),
                          name = "Relative abundance (%)",
                          breaks = c(0, 25, 50, 75)) +
    facet_grid(depth_type ~ Division_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text   = element_text(face = "bold", size = 10),
          axis.title   = element_text(size = 16, face = "bold"),
          axis.text    = element_text(size = 12),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text  = element_text(size = 14))
)

# Page C: Per-cluster (k=4) ASV community contribution x depth_type
asv_cluster_assign <- character(ncol(mat_div_p))
names(asv_cluster_assign) <- colnames(mat_div_p)
asv_max_mean_c <- rep(-Inf, ncol(mat_div_p))
for (k in c("1", "2", "3", "4")) {
  smp_k <- names(clusnum_p)[as.character(clusnum_p) == k]
  asv_mean_k <- colMeans(mat_div_p[smp_k, , drop = FALSE])
  upd <- asv_mean_k > asv_max_mean_c
  asv_cluster_assign[upd] <- k
  asv_max_mean_c[upd] <- asv_mean_k[upd]
}
asv_cluster_f <- factor(asv_cluster_assign, levels = c("1","2","3","4"))
n_asv_cluster_p <- table(asv_cluster_f)

comp_RA_mat_c <- sapply(c("1","2","3","4"), function(k) {
  rowSums(mat_div_p[, asv_cluster_f == k, drop = FALSE])
})

comp_df_c <- as.data.frame(comp_RA_mat_c) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(cluster    = factor(as.character(clusnum_p[Sample]),
                                    levels=c("1","2","3","4")),
                lat        = meta_asgard_p2[Sample,"lat"],
                lon        = meta_asgard_p2[Sample,"lon"],
                depth_type = meta_asgard_p2[Sample,"depth_type"]) %>%
  dplyr::filter(!is.na(lat), !is.na(lon), !is.na(depth_type)) %>%
  tidyr::pivot_longer(c("1","2","3","4"),
                      names_to="CC", values_to="RA") %>%
  dplyr::mutate(RA_pct=RA*100,
                CC=factor(CC, levels=c("1","2","3","4")),
                depth_type=factor(depth_type,
                                  levels=c("surf","mid","bottom")))

comp_df_c$CC_lab <- factor(
  paste0("Cluster ", as.character(comp_df_c$CC),
         " (", n_asv_cluster_p[as.character(comp_df_c$CC)], " ASVs)"),
  levels = paste0("Cluster ", c("1","2","3","4"),
                  " (", n_asv_cluster_p[c("1","2","3","4")], " ASVs)")
)

print(
  ggmap(mapz) +
    geom_point(data = comp_df_c,
               aes(x = lon, y = lat, size = RA_pct, color = cluster),
               alpha = 0.75) +
    scale_color_manual(values = cc_p, name = "Cluster") +
    scale_size_continuous(range = c(1, 12),
                          name = "Relative abundance (%)",
                          breaks = c(0, 25, 50, 75)) +
    facet_grid(depth_type ~ CC_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text   = element_text(face = "bold", size = 14),
          axis.title   = element_text(size = 16, face = "bold"),
          axis.text    = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12))
)

dev.off()

# ==============================================================================
# Section 7: S06 Page 6-style community-contribution map (per-cluster, k=4)
# Equivalent of S06 map_colclus_abundance Page 6, adapted to processing.
# Each ASV is assigned to the cluster (1-4) where its mean RA is highest;
# per-sample component RA is then mapped by depth_type x cluster-component,
# coloured by the sample's own cluster, sized by component RA.
# Overlap mitigation: large bubbles drawn first + alpha 0.5 + small jitter.
# NOTE: Page 1 = per-cluster component RA (as above). Page 2 = per-sample ASV
#       richness (Observed ASVs) faceted by depth x cluster, no ASV assignment.
# ==============================================================================

# Assign each ASV to the cluster (1-4) with highest mean RA
asv_assign_p7 <- character(ncol(mat_div_p))
names(asv_assign_p7) <- colnames(mat_div_p)
asv_best_p7 <- rep(-Inf, ncol(mat_div_p))
for (k in c("1", "2", "3", "4")) {
  smp_k <- names(clusnum_p)[as.character(clusnum_p) == k]
  mk    <- colMeans(mat_div_p[smp_k, , drop = FALSE])
  upd   <- mk > asv_best_p7
  asv_assign_p7[upd] <- k
  asv_best_p7[upd]   <- mk[upd]
}
asv_f_p7 <- factor(asv_assign_p7, levels = c("1", "2", "3", "4"))
n_asv_p7 <- table(asv_f_p7)

comp_RA_p7 <- sapply(c("1", "2", "3", "4"),
                     function(k) rowSums(mat_div_p[, asv_f_p7 == k, drop = FALSE]))

comp_df_p7 <- as.data.frame(comp_RA_p7) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(cluster    = factor(as.character(clusnum_p[Sample]),
                                    levels = c("1", "2", "3", "4")),
                lat        = meta_asgard_p2[Sample, "lat"],
                lon        = meta_asgard_p2[Sample, "lon"],
                depth_type = meta_asgard_p2[Sample, "depth_type"]) %>%
  dplyr::filter(!is.na(lat), !is.na(lon), !is.na(depth_type)) %>%
  tidyr::pivot_longer(c("1", "2", "3", "4"),
                      names_to = "Component", values_to = "RA") %>%
  dplyr::mutate(RA_pct     = RA * 100,
                Component  = factor(Component, levels = c("1", "2", "3", "4")),
                depth_type = factor(depth_type,
                                    levels = c("surf", "mid", "bottom"))) %>%
  dplyr::arrange(desc(RA_pct))   # large bubbles first

comp_df_p7$Component_lab <- factor(
  paste0("Cluster ", as.character(comp_df_p7$Component),
         " (", n_asv_p7[as.character(comp_df_p7$Component)], " ASVs)"),
  levels = paste0("Cluster ", c("1", "2", "3", "4"),
                  " (", n_asv_p7[c("1", "2", "3", "4")], " ASVs)")
)

pdf(here::here("output_p", "maps",
               "processing_colclus_abundance_4clusters.pdf"),
    width = 16, height = 14)

print(
  ggmap(mapz) +
    geom_point(data = comp_df_p7,
               aes(x = lon, y = lat, size = RA_pct, color = cluster),
               alpha = 0.5,
               position = position_jitter(width = 0.05, height = 0.05,
                                          seed = 1)) +
    scale_color_manual(values = cc_p, name = "Cluster") +
    scale_size_continuous(range = c(0.5, 11),
                          name = "Relative abundance (%)",
                          breaks = c(0, 25, 50, 75)) +
    facet_grid(depth_type ~ Component_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text   = element_text(face = "bold", size = 14),
          axis.title   = element_text(size = 16, face = "bold"),
          axis.text    = element_text(size = 11),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12))
)

# Page 2: per-sample ASV richness (Observed ASVs) — faceted by depth x cluster
# 割り当て無し。各サンプルを自クラスター列に1回だけ配置し、richness(Observed ASVs)を
# ドットサイズで表示。richness = rowSums(asgard_filtered_p_hm2 > 0)（221 ASVセット内の存在種数）。
rich_p_smp <- rowSums(asgard_filtered_p_hm2 > 0)
rich_df_p <- data.frame(
  Sample        = rownames(asgard_filtered_p_hm2),
  Observed_ASVs = as.integer(rich_p_smp),
  cluster       = factor(as.character(clusnum_p[rownames(asgard_filtered_p_hm2)]),
                         levels = c("1", "2", "3", "4")),
  lon           = meta_asgard_p2[rownames(asgard_filtered_p_hm2), "lon"],
  lat           = meta_asgard_p2[rownames(asgard_filtered_p_hm2), "lat"],
  depth_type    = factor(meta_asgard_p2[rownames(asgard_filtered_p_hm2), "depth_type"],
                         levels = c("surf", "mid", "bottom")),
  stringsAsFactors = FALSE
)
rich_df_p <- rich_df_p[!is.na(rich_df_p$lon) & !is.na(rich_df_p$lat) &
                       !is.na(rich_df_p$depth_type) & !is.na(rich_df_p$cluster), ]
n_cl_p <- table(rich_df_p$cluster)
rich_df_p$cluster_lab <- factor(
  paste0("Cluster ", as.character(rich_df_p$cluster),
         " (n=", n_cl_p[as.character(rich_df_p$cluster)], ")"),
  levels = paste0("Cluster ", c("1", "2", "3", "4"),
                  " (n=", n_cl_p[c("1", "2", "3", "4")], ")"))

print(
  ggmap(mapz) +
    geom_point(data = rich_df_p,
               aes(x = lon, y = lat, size = Observed_ASVs, color = cluster),
               alpha = 0.85) +
    scale_color_manual(values = cc_p, name = "Cluster") +
    scale_size_continuous(range = c(1, 11), name = "Observed ASVs") +
    facet_grid(depth_type ~ cluster_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text      = element_text(face = "bold", size = 18),
          axis.title      = element_text(size = 18, face = "bold"),
          axis.text       = element_text(size = 13),
          axis.text.x     = element_text(size = 13, angle = 45, hjust = 1),
          legend.title    = element_text(size = 17, face = "bold"),
          legend.text     = element_text(size = 15),
          legend.key.size = unit(0.9, "cm"))
)

dev.off()

# ==============================================================================
# Section 8: Free-living vs Particle-associated (k=2) ASV richness map
# 各サンプルの richness (Observed ASVs) を depth × FL/PA division グリッドで地図化。
# FL = cluster 1, PA = clusters 2-4（サンプル所属）。割り当て(winner-take-all)は不要。
# richness = rowSums(asgard_filtered_p_hm2 > 0)。rich_df_p (Section 7) を再利用。
# ==============================================================================

rich_df_flpa <- rich_df_p
rich_df_flpa$division2 <- factor(
  ifelse(rich_df_flpa$cluster == "1", "Free-living", "Particle-associated"),
  levels = c("Free-living", "Particle-associated"))
n_d2_p <- table(rich_df_flpa$division2)
rich_df_flpa$division2_lab <- factor(
  paste0(as.character(rich_df_flpa$division2),
         " (n=", n_d2_p[as.character(rich_df_flpa$division2)], ")"),
  levels = paste0(c("Free-living", "Particle-associated"),
                  " (n=", n_d2_p[c("Free-living", "Particle-associated")], ")"))

pdf(here::here("output_p", "maps",
               "processing_division_contribution_FLPA.pdf"),
    width = 14, height = 14)

print(
  ggmap(mapz) +
    geom_point(data = rich_df_flpa,
               aes(x = lon, y = lat, size = Observed_ASVs, color = cluster),
               alpha = 0.85) +
    scale_color_manual(values = cc_p, name = "Cluster") +
    scale_size_continuous(range = c(1, 11), name = "Observed ASVs") +
    facet_grid(depth_type ~ division2_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text.x    = element_text(face = "bold", size = 11),
          strip.text.y    = element_text(face = "bold", size = 18),
          axis.title      = element_text(size = 18, face = "bold"),
          axis.text       = element_text(size = 13),
          axis.text.x     = element_text(size = 13, angle = 45, hjust = 1),
          legend.title    = element_text(size = 17, face = "bold"),
          legend.text     = element_text(size = 15),
          legend.key.size = unit(0.9, "cm"))
)

dev.off()

# ==============================================================================
# Section 9: Kruskal-Wallis test: FL/PA component RA ~ depth_type
# Tests whether the free-living / particle-associated component RA differs
# among surf / mid / bottom (cf. S06 Section 8b). Bonferroni alpha = 0.05/2.
# ==============================================================================

flpa_kw_results <- comp_df_p %>%
  dplyr::select(Sample, Division, RA_pct, depth_type) %>%
  dplyr::group_by(Division) %>%
  dplyr::summarise(
    n_surf      = sum(depth_type == "surf"),
    n_mid       = sum(depth_type == "mid"),
    n_bottom    = sum(depth_type == "bottom"),
    mean_surf   = round(mean(RA_pct[depth_type == "surf"]),   2),
    mean_mid    = round(mean(RA_pct[depth_type == "mid"]),    2),
    mean_bottom = round(mean(RA_pct[depth_type == "bottom"]), 2),
    chi_sq      = round(kruskal.test(RA_pct ~ depth_type)$statistic, 3),
    df          = kruskal.test(RA_pct ~ depth_type)$parameter,
    p_value     = round(kruskal.test(RA_pct ~ depth_type)$p.value, 4),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    bonferroni_alpha = 0.05 / 2,
    sig_bonferroni   = ifelse(p_value < bonferroni_alpha, "*", "ns")
  )

cat("\n--- Kruskal-Wallis: FL/PA component RA ~ depth_type ---\n")
cat("(Bonferroni alpha = 0.05/2 = 0.025)\n")
print(as.data.frame(flpa_kw_results), row.names = FALSE)

write.csv(flpa_kw_results,
  here::here("output_p", "maps", "processing_FLPA_componentRA_depth_kruskal.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 10: Kruskal-Wallis test: per-cluster (k=4) component RA ~ depth_type
# Tests whether each cluster's component RA differs among surf / mid / bottom
# (processing analogue of S06's A/B/C division test). Reuses comp_df_p7.
# Bonferroni alpha = 0.05/4 = 0.0125.
# ==============================================================================

dir.create(here::here("output_p", "cluster_summary"),
           showWarnings = FALSE, recursive = TRUE)

cluster_kw_results <- comp_df_p7 %>%
  dplyr::filter(depth_type %in% c("surf", "mid", "bottom")) %>%
  dplyr::mutate(depth_type = factor(depth_type,
                                    levels = c("surf", "mid", "bottom"))) %>%
  dplyr::group_by(Component) %>%
  dplyr::summarise(
    n_surf      = sum(depth_type == "surf"),
    n_mid       = sum(depth_type == "mid"),
    n_bottom    = sum(depth_type == "bottom"),
    mean_surf   = round(mean(RA_pct[depth_type == "surf"]),   2),
    mean_mid    = round(mean(RA_pct[depth_type == "mid"]),    2),
    mean_bottom = round(mean(RA_pct[depth_type == "bottom"]), 2),
    chi_sq      = round(kruskal.test(RA_pct ~ depth_type)$statistic, 3),
    df          = kruskal.test(RA_pct ~ depth_type)$parameter,
    p_value     = round(kruskal.test(RA_pct ~ depth_type)$p.value, 4),
    .groups     = "drop"
  ) %>%
  dplyr::rename(Cluster = Component) %>%
  dplyr::mutate(
    bonferroni_alpha = round(0.05 / 4, 4),
    sig_bonferroni   = ifelse(p_value < 0.05 / 4, "*", "ns")
  )

cat("\n--- Kruskal-Wallis: per-cluster component RA ~ depth_type ---\n")
cat("(Bonferroni alpha = 0.05/4 = 0.0125)\n")
print(as.data.frame(cluster_kw_results), row.names = FALSE)

write.csv(cluster_kw_results,
  here::here("output_p", "cluster_summary",
             "processing_cluster_componentRA_depth_kruskal.csv"),
  row.names = FALSE)

message("04_maps.R: done. PDFs written:")
message("  output_p/maps/processing_map.pdf (legacy)")
message("  output_p/maps/maps_pa.pdf (per-ASV PA maps)")
message("  output_p/maps/ASGARD_processing_map_4clusters.pdf (hero)")
message("  output_p/maps/ASGARD_processing_map_4clusters_detail.pdf (detail)")
message("  output_p/maps/ASGARD_processing_division_depth_map.pdf (division x depth)")
message("  output_p/maps/processing_colclus_abundance_4clusters.pdf (p1 component RA; p2 richness)")
message("  output_p/maps/processing_division_contribution_FLPA.pdf (FL/PA richness)")
