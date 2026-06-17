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
  paste0("Division ", as.character(comp_df$Division),
         " (", n_asv_div[as.character(comp_df$Division)], " ASVs)"),
  levels = paste0("Division ", c("A", "B", "C"),
                  " (", n_asv_div[c("A", "B", "C")], " ASVs)")
)

print(
  ggmap(mapz_survey) +
    geom_point(data = comp_df,
               aes(x = lon, y = lat,
                   size = RA_pct, color = cluster11),
               alpha = 0.75) +
    scale_color_manual(values = cc11, name = "Cluster") +
    scale_size_continuous(range = c(0.5, 8),
                          name = "Component RA (%)",
                          breaks = c(0, 25, 50, 75)) +
    facet_grid(depth_type ~ Division_lab) +
    labs(title = "Per-division ASV community contribution x depth type",
         x = "Longitude", y = "Latitude") +
    theme(strip.text = element_text(face = "bold", size = 13),
          plot.title = element_text(face = "bold", size = 16))
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

message("\nS06_maps.R: done.")
message("  PDF: ASGARD_survey_map_11clusters.pdf, ASGARD_survey_map_11clusters_detail.pdf")
message("  PDF: ASGARD_survey_all_stations.pdf")
message("  PDF: map_colclus_abundance_11clusters.pdf")
message("  CSV: output/survey/maps/cluster11_geographic_summary.csv")
message("  CSV: output/survey/maps/cluster11_depth_type_test.csv")
message("  CSV: output/survey/maps/assemblage_RA_depth_kruskal.csv")
message("  CSV: output/survey/maps/B2_C2b_per_station.csv")
message("  CSV: output/survey/maps/B2_C2b_per_transect.csv")
