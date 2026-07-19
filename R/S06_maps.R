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
library(sf)          # coastline/land polygons for filled-contour maps (Section 13)

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

detail_theme <- theme(
  plot.title   = element_text(face = "bold", size = 20),
  axis.title   = element_text(size = 18, face = "bold"),
  axis.text    = element_text(size = 14),
  legend.title = element_text(size = 16, face = "bold"),
  legend.text  = element_text(size = 14),
  strip.text   = element_text(face = "bold", size = 18))

# station名ラベル付き / With station labels
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 3.5, alpha = 0.7) +
    scale_color_manual(values = cc11, labels = clbl) +
    geom_text_repel(data = a_map, aes(x = lon, y = lat, label = station),
                    size = 3, max.overlaps = 20, alpha = 0.8) +
    labs(color = NULL) +
    detail_theme
)

# depth_typeでfacet / Faceted by depth type
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 3, alpha = 0.7) +
    scale_color_manual(values = cc11, labels = clbl) +
    geom_text_repel(data = a_map, aes(x = lon, y = lat, label = station),
                    size = 2.5, max.overlaps = 15, alpha = 0.7) +
    facet_wrap(~ depth_type) +
    labs(color = NULL) +
    detail_theme
)

# クラスターでfacet / Faceted by cluster
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = cc11, guide = "none") +
    facet_wrap(~ cluster11, ncol = 4) +
    detail_theme
)

# depth_type x cluster グリッド / Grid
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map, aes(x = lon, y = lat, color = cluster11), size = 2, alpha = 0.8) +
    scale_color_manual(values = cc11, guide = "none") +
    facet_grid(depth_type ~ cluster11) +
    detail_theme + theme(strip.text = element_text(face = "bold", size = 13))
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

# ==============================================================================
# Section 6: 全ステーションマップ（トランセクト色分け）
# All stations map colored by transect prefix
# ==============================================================================

stations_all <- meta_asgard %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct(station, .keep_all = TRUE) %>%
  select(station, lat, lon)

# Assign transect group (DBO split into DBO1, DBO2, DBO3)
stations_all$prefix <- sub("[0-9.]+$", "", stations_all$station)
stations_all$prefix <- ifelse(grepl("^DBO1", stations_all$station), "DBO1",
                       ifelse(grepl("^DBO2", stations_all$station), "DBO2",
                       ifelse(grepl("^DBO3", stations_all$station), "DBO3",
                       stations_all$prefix)))
stations_all$prefix <- factor(stations_all$prefix,
  levels = c("SLY", "CBE", "CB", "CBW", "CPL", "DBO1", "DBO2", "DBO3", "CNL", "IL", "KL", "CL")
)

prefix_colors <- c(
  "SLY" = "#D32F2F", "CBE" = "#F57C00", "CB" = "#000000",
  "CBW" = "#1565C0", "CPL" = "#2E7D32", "DBO1" = "#C2185B",
  "DBO2" = "#E91E63", "DBO3" = "#880E4F",
  "CNL" = "#7B1FA2", "IL" = "#00838F", "KL" = "#FFEB3B",
  "CL" = "#4E342E"
)

pdf(here::here("output", "survey", "maps", "ASGARD_survey_all_stations.pdf"),
    width = 14, height = 10)

print(
  ggmap(mapz_survey) +
    geom_point(data = stations_all, aes(x = lon, y = lat, color = prefix),
               size = 3.5, alpha = 0.8) +
    geom_text_repel(data = stations_all, aes(x = lon, y = lat, label = station, color = prefix),
                    size = 3, max.overlaps = 40, alpha = 0.8, show.legend = FALSE) +
    scale_color_manual(values = prefix_colors, name = "Transect") +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.text = element_text(size = 11),
          legend.title = element_text(face = "bold", size = 12))
)

dev.off()

# ==============================================================================
# Section 7: Assemblage abundance × depth-type grid map
# 6 ASV assemblage × 3 depth type (surf/mid/bottom) のグリッドマップ
# REQUIRES: colclusnum (from S02), asgard_filtered (from S01)
# ==============================================================================

# Per-sample total RA in each ASV assemblage (1..6), in %
assem_ra <- sapply(1:6, function(k) {
  asv_k <- names(colclusnum)[colclusnum == k]
  rowSums(asgard_filtered[, asv_k, drop = FALSE])
}) * 100
colnames(assem_ra) <- paste0("Assemblage_", 1:6)

# Long table with lat/lon, cluster11, depth_type
asm_long <- as.data.frame(assem_ra) %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::pivot_longer(-Sample, names_to = "Assemblage", values_to = "RA_pct") %>%
  dplyr::inner_join(
    a_map %>% dplyr::select(Sample = ID, lat, lon, cluster11, depth_type) %>%
      dplyr::mutate(Sample = rownames(a_map)[match(Sample, a_map$ID)]),
    by = "Sample"
  )

# Assemblage labels with ASV count
assem_n_tab <- table(colclusnum)
asm_long$Assemblage <- factor(
  asm_long$Assemblage,
  levels = paste0("Assemblage_", 1:6),
  labels = paste0("Assemblage ", 1:6,
                  " (n=", assem_n_tab[as.character(1:6)], ")")
)

pdf(here::here("output", "survey", "maps", "map_colclus_abundance_11clusters.pdf"),
    width = 24, height = 12)

# Page 1: 6 assemblage panels (all depths combined)
print(
  ggmap(mapz_survey) +
    geom_point(data = asm_long,
               aes(x = lon, y = lat, color = cluster11, size = RA_pct),
               alpha = 0.7) +
    scale_color_manual(values = cc11, name = "Cluster") +
    scale_size_continuous(range = c(0.5, 6), name = "Assemblage RA (%)") +
    facet_wrap(~ Assemblage, nrow = 2) +
    labs(title = "Per-assemblage abundance per sample",
         x = "Longitude", y = "Latitude") +
    theme(strip.text  = element_text(face = "bold", size = 12),
          plot.title  = element_text(face = "bold", size = 16))
)

# Page 2: Assemblage × depth_type grid
print(
  ggmap(mapz_survey) +
    geom_point(data = asm_long,
               aes(x = lon, y = lat, color = cluster11, size = RA_pct),
               alpha = 0.7) +
    scale_color_manual(values = cc11, name = "Cluster") +
    scale_size_continuous(range = c(0.5, 6), name = "Assemblage RA (%)") +
    facet_grid(depth_type ~ Assemblage) +
    labs(title = "Per-assemblage abundance × depth type",
         x = "Longitude", y = "Latitude") +
    theme(strip.text  = element_text(face = "bold", size = 10),
          plot.title  = element_text(face = "bold", size = 16))
)

# Page 3: 11-cluster × depth_type grid (sample locations)
# Sample counts per cluster for panel labels
cluster_n_tab <- table(a_map$cluster11)
a_map$cluster11_lab <- factor(
  paste0(as.character(a_map$cluster11),
         " (n=", cluster_n_tab[as.character(a_map$cluster11)], ")"),
  levels = paste0(hier_levels_11,
                  " (n=", cluster_n_tab[hier_levels_11], ")")
)

print(
  ggmap(mapz_survey) +
    geom_point(data = a_map,
               aes(x = lon, y = lat, color = cluster11),
               size = 2, alpha = 0.8) +
    scale_color_manual(values = cc11, guide = "none") +
    facet_grid(depth_type ~ cluster11_lab) +
    labs(title = "Sample locations by cluster × depth type (11 clusters)",
         x = "Longitude", y = "Latitude") +
    theme(strip.text  = element_text(face = "bold", size = 9),
          plot.title  = element_text(face = "bold", size = 16))
)

# Page 4: k=3 division (A / B / C) sample locations
a_map$division3 <- factor(substr(as.character(a_map$cluster11), 1, 1),
                          levels = c("A", "B", "C"))
div3_n <- table(a_map$division3)
a_map$division3_lab <- factor(
  paste0(as.character(a_map$division3),
         " (n=", div3_n[as.character(a_map$division3)], ")"),
  levels = paste0(c("A", "B", "C"),
                  " (n=", div3_n[c("A", "B", "C")], ")")
)
div3_colors <- c("A" = "#E31A1C", "B" = "#33A02C", "C" = "#1F78B4")

print(
  ggmap(mapz_survey) +
    geom_point(data = a_map,
               aes(x = lon, y = lat, color = division3),
               size = 2.5, alpha = 0.8) +
    scale_color_manual(values = div3_colors, name = "Division") +
    facet_wrap(~ division3_lab, nrow = 1) +
    labs(title = "Sample locations by k=3 division (A / B / C)",
         x = "Longitude", y = "Latitude") +
    theme(strip.text  = element_text(face = "bold", size = 14),
          plot.title  = element_text(face = "bold", size = 16))
)

# Page 5: k=3 division × depth_type grid
print(
  ggmap(mapz_survey) +
    geom_point(data = a_map,
               aes(x = lon, y = lat, color = division3),
               size = 2.5, alpha = 0.8) +
    scale_color_manual(values = div3_colors, name = "Division") +
    facet_grid(depth_type ~ division3_lab) +
    labs(title = "Sample locations by k=3 division × depth type",
         x = "Longitude", y = "Latitude") +
    theme(strip.text  = element_text(face = "bold", size = 12),
          plot.title  = element_text(face = "bold", size = 16))
)

# Page 6: Per-division ASV community contribution × depth
# サンプル別 RA 化 → 各ASVを「mean RAが最大のdivision」に割り当て
# → 各サンプルの A/B/C 成分RA を depth × division でmap化
mat_div <- sweep(asgard_filtered, 1, rowSums(asgard_filtered), "/")
asv_division_assignment <- character(ncol(mat_div))
names(asv_division_assignment) <- colnames(mat_div)
asv_max_mean_d <- rep(-Inf, ncol(mat_div))
for (d in c("A", "B", "C")) {
  smp <- names(clusnum11)[substr(as.character(clusnum11), 1, 1) == d]
  asv_mean_d <- colMeans(mat_div[smp, , drop = FALSE])
  upd <- asv_mean_d > asv_max_mean_d
  asv_division_assignment[upd] <- d
  asv_max_mean_d[upd] <- asv_mean_d[upd]
}
asv_div_f <- factor(asv_division_assignment, levels = c("A", "B", "C"))
n_asv_div <- table(asv_div_f)

comp_RA_mat <- sapply(c("A", "B", "C"), function(d) {
  rowSums(mat_div[, asv_div_f == d, drop = FALSE])
})

comp_df <- as.data.frame(comp_RA_mat) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(cluster11  = factor(as.character(clusnum11[Sample]),
                                    levels = hier_levels_11),
                lat        = meta_asgard[Sample, "lat"],
                lon        = meta_asgard[Sample, "lon"],
                depth_type = meta_asgard[Sample, "depth_type"]) %>%
  dplyr::filter(!is.na(lat), !is.na(lon), !is.na(depth_type)) %>%
  tidyr::pivot_longer(c("A", "B", "C"),
                      names_to = "Division", values_to = "RA") %>%
  dplyr::mutate(RA_pct     = RA * 100,
                Division   = factor(Division, levels = c("A", "B", "C")),
                depth_type = factor(depth_type,
                                    levels = c("surf", "mid", "bottom")))

comp_df$Division_lab <- factor(
  paste0(as.character(comp_df$Division),
         " (", n_asv_div[as.character(comp_df$Division)], " ASVs)"),
  levels = paste0(c("A", "B", "C"),
                  " (", n_asv_div[c("A", "B", "C")], " ASVs)")
)

print(
  ggmap(mapz_survey) +
    geom_point(data = comp_df,
               aes(x = lon, y = lat,
                   size = RA_pct, color = cluster11),
               alpha = 0.75) +
    scale_color_manual(values = cc11, name = "Cluster") +
    scale_size_continuous(range = c(0.5, 10),
                          name = "Relative abundance (%)",
                          breaks = c(0, 25, 50, 75)) +
    facet_grid(depth_type ~ Division_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text   = element_text(face = "bold", size = 20),
          axis.title   = element_text(size = 18),
          axis.text    = element_text(size = 14),
          legend.title = element_text(size = 17, face = "bold"),
          legend.text  = element_text(size = 15))
)

dev.off()

# ==============================================================================
# Section 8: Kruskal-Wallis test: assemblage RA ~ depth_type
# Assemblage 別に「相対存在量が depth_type で異なるか」を検定
# ==============================================================================

kw_long <- asm_long %>%
  dplyr::filter(!is.na(depth_type)) %>%
  dplyr::mutate(Assemblage_clean = sub(" \\(n=.+\\)$", "", as.character(Assemblage)))

kw_results <- kw_long %>%
  dplyr::group_by(Assemblage_clean) %>%
  dplyr::summarise(
    n_surf   = sum(depth_type == "surf"),
    n_mid    = sum(depth_type == "mid"),
    n_bottom = sum(depth_type == "bottom"),
    mean_surf   = round(mean(RA_pct[depth_type == "surf"]),   2),
    mean_mid    = round(mean(RA_pct[depth_type == "mid"]),    2),
    mean_bottom = round(mean(RA_pct[depth_type == "bottom"]), 2),
    chi_sq  = round(kruskal.test(RA_pct ~ depth_type)$statistic, 3),
    df      = kruskal.test(RA_pct ~ depth_type)$parameter,
    p_value = round(kruskal.test(RA_pct ~ depth_type)$p.value, 4),
    .groups = "drop"
  ) %>%
  dplyr::mutate(sig = ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01,  "**",
                       ifelse(p_value < 0.05,  "*", "ns"))))

cat("\n--- Kruskal-Wallis: assemblage RA ~ depth_type ---\n")
print(as.data.frame(kw_results), row.names = FALSE)

write.csv(kw_results,
  here::here("output", "survey", "maps", "assemblage_RA_depth_kruskal.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 8b: Kruskal-Wallis test: ASV-division component RA ~ depth_type
# (Page 6 の数値検定。各サンプルにおける A/B/C 成分RAが depth で異なるか)
# Bonferroni 補正: α = 0.05 / 3 = 0.0167
# ==============================================================================

div_kw_df <- comp_df %>%
  dplyr::select(Sample, Division, RA_pct, depth_type)

div_kw_results <- div_kw_df %>%
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
    bonferroni_alpha = 0.05 / 3,
    sig_bonferroni   = ifelse(p_value < bonferroni_alpha, "*", "ns")
  )

cat("\n--- Kruskal-Wallis: ASV-division component RA ~ depth_type ---\n")
cat("(Bonferroni alpha = 0.05/3 = 0.0167)\n")
print(as.data.frame(div_kw_results), row.names = FALSE)

write.csv(div_kw_results,
  here::here("output", "survey", "maps", "division_componentRA_depth_kruskal.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 9: Per-station / per-transect breakdown for selected clusters
# B2a, B2b, C2b1, C2b2 がどのステーション・トランセクトで多く採取されたか
# ==============================================================================

mb <- meta_asgard %>% tibble::rownames_to_column("Sample")
mb$cluster <- factor(as.character(clusnum11[mb$Sample]), levels = hier_levels_11)

mb$transect <- ifelse(grepl("^DBO3", mb$station), "DBO3",
               ifelse(grepl("^DBO2", mb$station), "DBO2",
               ifelse(grepl("^DBO1", mb$station), "DBO1",
               ifelse(grepl("^CBE",  mb$station), "CBE",
               ifelse(grepl("^CBW",  mb$station), "CBW",
               ifelse(grepl("^CB[0-9]", mb$station), "CB",
               ifelse(grepl("^CNL",  mb$station), "CNL",
               ifelse(grepl("^CPL",  mb$station), "CPL",
               ifelse(grepl("^IL",   mb$station), "IL",
               ifelse(grepl("^KL",   mb$station), "KL",
               ifelse(grepl("^CL",   mb$station), "CL",
               ifelse(grepl("^SLY",  mb$station), "SLY", "Other"))))))))))))

target_clusters <- c("B2a", "B2b", "C2b1", "C2b2")

station_breakdown  <- list()
transect_breakdown <- list()

for (cl in target_clusters) {
  sub <- mb %>% dplyr::filter(cluster == cl)
  if (nrow(sub) == 0) next

  st <- sub %>%
    dplyr::group_by(station) %>%
    dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(pct_of_cluster = round(n_samples / sum(n_samples) * 100, 1),
                  cluster = cl) %>%
    dplyr::arrange(-n_samples)
  station_breakdown[[cl]] <- st

  tr <- sub %>%
    dplyr::group_by(transect) %>%
    dplyr::summarise(n_samples  = dplyr::n(),
                     n_stations = dplyr::n_distinct(station),
                     stations   = paste(sort(unique(station)), collapse = ", "),
                     .groups = "drop") %>%
    dplyr::mutate(pct_of_cluster = round(n_samples / sum(n_samples) * 100, 1),
                  cluster = cl) %>%
    dplyr::arrange(-n_samples)
  transect_breakdown[[cl]] <- tr

  cat("\n--- Cluster", cl, "---\n")
  cat("Total samples:", nrow(sub),
      "; unique stations:", length(unique(sub$station)), "\n")
  cat("\nPer station:\n")
  print(as.data.frame(st), row.names = FALSE)
  cat("\nPer transect:\n")
  print(as.data.frame(tr), row.names = FALSE)
}

station_breakdown_df  <- dplyr::bind_rows(station_breakdown) %>%
  dplyr::select(cluster, station, n_samples, pct_of_cluster)
transect_breakdown_df <- dplyr::bind_rows(transect_breakdown) %>%
  dplyr::select(cluster, transect, n_samples, n_stations, pct_of_cluster, stations)

write.csv(station_breakdown_df,
  here::here("output", "survey", "maps", "B2_C2b_per_station.csv"),
  row.names = FALSE)
write.csv(transect_breakdown_df,
  here::here("output", "survey", "maps", "B2_C2b_per_transect.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 10: ASV richness map (Observed ASVs) by depth x division
# 各サンプルの観測ASV数(richness)を depth × division (A/B/C) グリッドで地図化。
#   richness = rowSums(asgard_filtered > 0)  ... クラスタリング用258 ASVセット内の存在種数
#   ドット色 = 11クラスター(cc11)、サイズ = Observed ASVs。ASVのグループ割り当ては不要。
# NOTE: 全16Sプール(3076)での観測種数・Chao1 は S07 (asgard_alpha_df) を参照。
# ==============================================================================

a_rich <- a_map
a_rich$Observed_ASVs <- rowSums(asgard_filtered[rownames(a_rich), ] > 0)
a_rich$division <- factor(substr(as.character(a_rich$cluster11), 1, 1),
                          levels = c("A", "B", "C"))
a_rich <- a_rich[!is.na(a_rich$division), ]

n_div_r <- table(a_rich$division)
a_rich$division_lab <- factor(
  paste0(as.character(a_rich$division), " (n=", n_div_r[as.character(a_rich$division)], ")"),
  levels = paste0(c("A", "B", "C"), " (n=", n_div_r[c("A", "B", "C")], ")"))

pdf(here::here("output", "survey", "maps",
               "ASGARD_survey_richness_by_depth_division.pdf"),
    width = 16, height = 14)
print(
  ggmap(mapz_survey) +
    geom_point(data = a_rich,
               aes(lon, lat, size = Observed_ASVs, color = cluster11),
               alpha = 0.85) +
    scale_color_manual(values = cc11, name = "Cluster") +
    scale_size_continuous(range = c(1, 11), name = "Observed ASVs") +
    facet_grid(depth_type ~ division_lab) +
    labs(x = "Longitude", y = "Latitude") +
    theme(strip.text      = element_text(face = "bold", size = 18),
          axis.title      = element_text(size = 18, face = "bold"),
          axis.text       = element_text(size = 13),
          legend.title    = element_text(size = 17, face = "bold"),
          legend.text     = element_text(size = 15),
          legend.key.size = unit(0.9, "cm"))
)
dev.off()

# ==============================================================================
# Section 11: Kruskal-Wallis: Observed ASVs (richness) ~ depth_type
# Section 10 の richness が depth (surf/mid/bottom) で異なるかを検定。
# 全体 + division別 (A/B/C, Bonferroni alpha = 0.05/3)。a_rich を再利用。
# ==============================================================================

kw_rich <- a_rich
kw_rich$depth_type <- factor(kw_rich$depth_type, levels = c("surf", "mid", "bottom"))
kw_rich <- kw_rich[!is.na(kw_rich$depth_type), ]

.kw_summary <- function(sub, label, alpha) {
  k  <- kruskal.test(Observed_ASVs ~ depth_type, data = sub)
  md <- tapply(sub$Observed_ASVs, sub$depth_type, median)
  nn <- table(sub$depth_type)
  data.frame(group = label, n = nrow(sub),
             n_surf = as.integer(nn["surf"]), n_mid = as.integer(nn["mid"]),
             n_bottom = as.integer(nn["bottom"]),
             median_surf = md["surf"], median_mid = md["mid"], median_bottom = md["bottom"],
             chi_sq = round(k$statistic, 3), df = k$parameter,
             p_value = round(k$p.value, 4), bonferroni_alpha = round(alpha, 4),
             sig = ifelse(k$p.value < alpha, "*", "ns"), row.names = NULL)
}

richness_kw <- .kw_summary(kw_rich, "Overall", 0.05)
for (d in c("A", "B", "C")) {
  sub <- kw_rich[kw_rich$division == d, ]
  if (nrow(sub) >= 3 && length(unique(as.character(sub$depth_type))) >= 2)
    richness_kw <- rbind(richness_kw, .kw_summary(sub, paste0("Division_", d), 0.05 / 3))
}

cat("\n--- Kruskal-Wallis: Observed ASVs (258-set richness) ~ depth_type ---\n")
print(richness_kw, row.names = FALSE)
write.csv(richness_kw,
  here::here("output", "survey", "maps", "richness_depth_kruskal.csv"),
  row.names = FALSE)

# ==============================================================================
# Section 12: Cruise progression map — stations coloured by sampling date,
# with near-surface temperature / salinity overlaid (as point size) for context.
#   REQUIRES: a_map (date, temp, salinity, lat, lon), mapz_survey (Section 1).
#   date format e.g. "Jun 19 2017 21:22"; one (shallowest) sample per station.
# ==============================================================================

a_prog <- a_map
# parse English month names + build date labels regardless of session locale
.old_lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
a_prog$date_p <- as.POSIXct(a_prog$date, format = "%b %d %Y %H:%M", tz = "UTC")
a_prog <- a_prog[!is.na(a_prog$lat) & !is.na(a_prog$lon) & !is.na(a_prog$date_p), ]
a_prog <- a_prog[order(a_prog$station, a_prog$depth_m), ]   # shallowest first
stn_prog <- a_prog[!duplicated(a_prog$station), ]           # one point per station
stn_prog$date_day <- as.Date(stn_prog$date_p)
stn_prog$date_num <- as.numeric(stn_prog$date_day)
stn_track <- stn_prog[order(stn_prog$date_p), ]            # ship route in time order
date_brks   <- pretty(stn_prog$date_day, n = 5)
date_labels <- format(date_brks, "%B %d")                  # e.g. "June 17"
Sys.setlocale("LC_TIME", .old_lct)

prog_theme <- theme(axis.title   = element_text(size = 16, face = "bold"),
                    axis.text    = element_text(size = 12),
                    legend.title = element_text(size = 13, face = "bold"),
                    legend.text  = element_text(size = 11))

prog_date_scale <- scale_color_viridis_c(name = "Date", option = "C",
                     breaks = as.numeric(date_brks), labels = date_labels)

pdf(here::here("output", "survey", "maps", "ASGARD_survey_cruise_progression.pdf"),
    width = 12, height = 9)

# Page 1: cruise progression (colour = date) with route line
print(
  ggmap(mapz_survey) +
    geom_path(data = stn_track, aes(lon, lat), color = "grey45",
              linewidth = 0.4, alpha = 0.6) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num), size = 3.5, alpha = 0.9) +
    prog_date_scale +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 2: cruise progression with station-name labels
print(
  ggmap(mapz_survey) +
    geom_path(data = stn_track, aes(lon, lat), color = "grey45",
              linewidth = 0.4, alpha = 0.6) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num), size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(data = stn_prog, aes(lon, lat, label = station),
                             size = 2.5, max.overlaps = Inf, alpha = 0.85) +
    prog_date_scale +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 3: colour = date, size = near-surface temperature
print(
  ggmap(mapz_survey) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num, size = temp), alpha = 0.85) +
    prog_date_scale +
    scale_size_continuous(name = "Temp (°C)", range = c(1, 9)) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 4: temperature, with station-name labels
print(
  ggmap(mapz_survey) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num, size = temp), alpha = 0.85) +
    ggrepel::geom_text_repel(data = stn_prog, aes(lon, lat, label = station),
                             size = 2.5, max.overlaps = Inf, alpha = 0.85) +
    prog_date_scale +
    scale_size_continuous(name = "Temp (°C)", range = c(1, 9)) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 5: colour = date, size = near-surface salinity
print(
  ggmap(mapz_survey) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num, size = salinity), alpha = 0.85) +
    prog_date_scale +
    scale_size_continuous(name = "Salinity", range = c(1, 9)) +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 6: salinity, with station-name labels
print(
  ggmap(mapz_survey) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num, size = salinity), alpha = 0.85) +
    ggrepel::geom_text_repel(data = stn_prog, aes(lon, lat, label = station),
                             size = 2.5, max.overlaps = Inf, alpha = 0.85) +
    prog_date_scale +
    scale_size_continuous(name = "Salinity", range = c(1, 9)) +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 7: temperature, with station labels and cruise route
print(
  ggmap(mapz_survey) +
    geom_path(data = stn_track, aes(lon, lat), color = "grey45",
              linewidth = 0.4, alpha = 0.6) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num, size = temp), alpha = 0.85) +
    ggrepel::geom_text_repel(data = stn_prog, aes(lon, lat, label = station),
                             size = 2.5, max.overlaps = Inf, alpha = 0.85) +
    prog_date_scale +
    scale_size_continuous(name = "Temp (°C)", range = c(1, 9)) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# Page 8: salinity, with station labels and cruise route
print(
  ggmap(mapz_survey) +
    geom_path(data = stn_track, aes(lon, lat), color = "grey45",
              linewidth = 0.4, alpha = 0.6) +
    geom_point(data = stn_prog, aes(lon, lat, color = date_num, size = salinity), alpha = 0.85) +
    ggrepel::geom_text_repel(data = stn_prog, aes(lon, lat, label = station),
                             size = 2.5, max.overlaps = Inf, alpha = 0.85) +
    prog_date_scale +
    scale_size_continuous(name = "Salinity", range = c(1, 9)) +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
)

# ==============================================================================
# Section 13: distance-masked (30 km) filled-contour maps of near-surface
# salinity and temperature, added as pages to the cruise progression PDF.
# IDW (power 2, latitude-corrected km distance); land masked white via
# rnaturalearth. NOTE: interpolated from sparse ship stations -> only filled
# within ~30 km of a station (non-synoptic; context only).
# ==============================================================================

fc_stn <- stn_prog[complete.cases(stn_prog[, c("lon", "lat", "temp", "salinity")]), ]
fc_radius_km <- 30
fc_xr <- range(fc_stn$lon) + c(-0.5, 0.5)
fc_yr <- range(fc_stn$lat) + c(-0.3, 0.3)
fc_kx <- 111 * cos(mean(fc_stn$lat) * pi / 180); fc_ky <- 111   # km per degree
fc_g  <- expand.grid(lon = seq(fc_xr[1], fc_xr[2], length.out = 240),
                     lat = seq(fc_yr[1], fc_yr[2], length.out = 240))
fc_D  <- sqrt((outer(fc_g$lon, fc_stn$lon, "-") * fc_kx)^2 +
              (outer(fc_g$lat, fc_stn$lat, "-") * fc_ky)^2)
fc_W  <- 1 / (fc_D^2); fc_W[!is.finite(fc_W)] <- max(fc_W[is.finite(fc_W)], na.rm = TRUE) * 1e6
fc_nn <- apply(fc_D, 1, min)
fc_land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
fc_jet  <- colorRampPalette(c("#00007F", "#0000FF", "#007FFF", "#00FFFF",
                              "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"))

fc_plot <- function(vals, brks, legend_name, lab_every = 1, show_labels = FALSE) {
  gg <- fc_g
  gg$z <- as.vector((fc_W %*% vals) / rowSums(fc_W))
  gg$z[fc_nn > fc_radius_km] <- NA
  p <- ggplot() +
    geom_contour_filled(data = gg[!is.na(gg$z), ],
                        aes(lon, lat, z = z, fill = after_stat(level_mid)), breaks = brks) +
    scale_fill_stepsn(colours = fc_jet(length(brks) - 1), breaks = brks,
                      limits = range(brks), name = legend_name,
                      labels = function(b) ifelse(b %% lab_every == 0, as.character(b), ""),
                      guide = guide_colorsteps(barheight = grid::unit(6, "cm"),
                                               show.limits = TRUE)) +
    geom_sf(data = fc_land, fill = "white", color = "grey30", linewidth = 0.3) +
    geom_point(data = fc_stn, aes(lon, lat), color = "black", size = 0.6, alpha = 0.6)
  if (show_labels)
    p <- p + ggrepel::geom_text_repel(data = fc_stn, aes(lon, lat, label = station),
                                      size = 2.2, max.overlaps = Inf, alpha = 0.85)
  p +
    coord_sf(xlim = fc_xr, ylim = fc_yr, crs = 4326, expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_line(color = "grey85", linewidth = 0.2))
}

# Pages 9-10: salinity / temperature filled contours (clean, final figures)
print(fc_plot(fc_stn$salinity, seq(25, 33, by = 0.5), "Salinity (psu)", lab_every = 1))
print(fc_plot(fc_stn$temp,     seq(-2, 11, by = 1),   "Temp (°C)",  lab_every = 2))
# Pages 11-12: same, with station-name labels (working reference, not for final figs)
print(fc_plot(fc_stn$salinity, seq(25, 33, by = 0.5), "Salinity (psu)", lab_every = 1, show_labels = TRUE))
print(fc_plot(fc_stn$temp,     seq(-2, 11, by = 1),   "Temp (°C)",  lab_every = 2, show_labels = TRUE))

# Pages 13-14: filled contour + labels on the terrain basemap (semi-transparent,
# so the Stadia terrain shows through, like pages 1-8)
fc_plot_gg <- function(vals, brks, legend_name, lab_every = 1) {
  gg <- fc_g
  gg$z <- as.vector((fc_W %*% vals) / rowSums(fc_W))
  gg$z[fc_nn > fc_radius_km] <- NA
  ggmap(mapz_survey) +
    geom_contour_filled(data = gg[!is.na(gg$z), ],
                        aes(lon, lat, z = z, fill = after_stat(level_mid)),
                        breaks = brks, alpha = 0.6) +
    scale_fill_stepsn(colours = fc_jet(length(brks) - 1), breaks = brks,
                      limits = range(brks), name = legend_name,
                      labels = function(b) ifelse(b %% lab_every == 0, as.character(b), ""),
                      guide = guide_colorsteps(barheight = grid::unit(6, "cm"), show.limits = TRUE)) +
    geom_point(data = fc_stn, aes(lon, lat), color = "black", size = 0.6, alpha = 0.7) +
    ggrepel::geom_text_repel(data = fc_stn, aes(lon, lat, label = station),
                             size = 2.2, max.overlaps = Inf, alpha = 0.85) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
}
print(fc_plot_gg(fc_stn$salinity, seq(25, 33, by = 0.5), "Salinity (psu)", lab_every = 1))
print(fc_plot_gg(fc_stn$temp,     seq(-2, 11, by = 1),   "Temp (°C)",  lab_every = 2))

# Pages 15-16: terrain basemap + filled field (fill) + station points coloured by
# sampling date (colour) + cruise route + labels -> field AND cruise timing.
fc_plot_gg_date <- function(vals, brks, legend_name, lab_every = 1) {
  gg <- fc_g
  gg$z <- as.vector((fc_W %*% vals) / rowSums(fc_W))
  gg$z[fc_nn > fc_radius_km] <- NA
  ggmap(mapz_survey) +
    geom_contour_filled(data = gg[!is.na(gg$z), ],
                        aes(lon, lat, z = z, fill = after_stat(level_mid)),
                        breaks = brks, alpha = 0.55) +
    scale_fill_stepsn(colours = fc_jet(length(brks) - 1), breaks = brks,
                      limits = range(brks), name = legend_name,
                      labels = function(b) ifelse(b %% lab_every == 0, as.character(b), ""),
                      guide = guide_colorsteps(barheight = grid::unit(5, "cm"),
                                               show.limits = TRUE, order = 2)) +
    geom_path(data = stn_track, aes(lon, lat), color = "grey30",
              linewidth = 0.4, alpha = 0.6) +
    geom_point(data = fc_stn, aes(lon, lat, color = date_num), size = 2.4, alpha = 0.95) +
    scale_color_viridis_c(name = "Date", option = "C",
                          breaks = as.numeric(date_brks), labels = date_labels,
                          guide = guide_colorbar(order = 1)) +
    ggrepel::geom_text_repel(data = fc_stn, aes(lon, lat, label = station),
                             size = 2, max.overlaps = Inf, alpha = 0.8) +
    labs(x = "Longitude", y = "Latitude") +
    prog_theme
}
print(fc_plot_gg_date(fc_stn$salinity, seq(25, 33, by = 0.5), "Salinity (psu)", lab_every = 1))
print(fc_plot_gg_date(fc_stn$temp,     seq(-2, 11, by = 1),   "Temp (°C)",  lab_every = 2))

dev.off()

message("\nS06_maps.R: done.")
message("  PDF: ASGARD_survey_map_11clusters.pdf, ASGARD_survey_map_11clusters_detail.pdf")
message("  PDF: ASGARD_survey_all_stations.pdf")
message("  PDF: map_colclus_abundance_11clusters.pdf")
message("  PDF: ASGARD_survey_richness_by_depth_division.pdf")
message("  PDF: ASGARD_survey_cruise_progression.pdf")
message("  CSV: output/survey/maps/richness_depth_kruskal.csv")
message("  CSV: output/survey/maps/cluster11_geographic_summary.csv")
message("  CSV: output/survey/maps/cluster11_depth_type_test.csv")
message("  CSV: output/survey/maps/assemblage_RA_depth_kruskal.csv")
message("  CSV: output/survey/maps/B2_C2b_per_station.csv")
message("  CSV: output/survey/maps/B2_C2b_per_transect.csv")
