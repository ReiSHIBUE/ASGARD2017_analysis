### S22_dbo3_waffle.R
### ASGARD 2017 Survey — Class composition at DBO3 stations
### DBO3 ステーション別のクラス組成ワッフル
###
### REQUIRES (from 00_setup.R, S01, S02):
###   asgard_filtered  - 181×258 ASV proportion matrix
###   meta_asgard      - metadata 181×45
###   clusnum11        - 11-cluster assignments (factor)
###   hier_levels_11   - ordered cluster names
###   cc11             - 11-colour palette
###   shorternames     - "Order; Family; Genus; ESV_N"
###   fullnameboot     - taxonomy + bootstrap df
###
### OUTPUT:
###   output/survey/waffle/ASGARD_dbo3_class_waffle.pdf

library(tidyverse)
library(RColorBrewer)
library(waffle)

# ==============================================================================
# Section 1: Identify DBO3 samples
# ==============================================================================

m <- meta_asgard %>% rownames_to_column("Sample")
m$cluster11 <- factor(as.character(clusnum11[m$Sample]), levels = hier_levels_11)

# Filter to DBO3 stations only
dbo3_m <- m %>% filter(grepl("^DBO3", station))

cat("=== DBO3 samples per station ===\n")
dbo3_summary <- dbo3_m %>%
  group_by(station) %>%
  summarise(n_samples = n(),
            clusters  = paste(sort(unique(as.character(cluster11))), collapse = ", "),
            .groups = "drop") %>%
  arrange(station)
print(as.data.frame(dbo3_summary), row.names = FALSE)

dbo3_samples <- dbo3_m$Sample

# ==============================================================================
# Section 2: ASV → Class mapping
# ==============================================================================

asv_idx        <- match(colnames(asgard_filtered), shorternames)
strip_boot     <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
asv_class      <- strip_boot(fullnameboot$Class[asv_idx])
names(asv_class) <- colnames(asgard_filtered)

# ==============================================================================
# Section 3: Compute Class composition per DBO3 station
# サンプル単位で Class 合計 → ステーション内のサンプル間で平均
# ==============================================================================

mat_dbo3 <- asgard_filtered[dbo3_samples, , drop = FALSE]

long_dbo3 <- as.data.frame(mat_dbo3) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(class = asv_class[ASV]) %>%
  inner_join(dbo3_m %>% select(Sample, station, cluster11), by = "Sample") %>%
  group_by(Sample) %>%
  mutate(RA = Abundance / sum(Abundance)) %>%
  ungroup()

# Per-sample sum per class
sample_class_dbo3 <- long_dbo3 %>%
  group_by(Sample, station, cluster11, class) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

# Per-station mean
station_class <- sample_class_dbo3 %>%
  group_by(station, class) %>%
  summarise(mean_RA = mean(sample_RA), .groups = "drop") %>%
  group_by(station) %>%
  mutate(pct = mean_RA / sum(mean_RA) * 100,
         units = round(pct)) %>%
  ungroup() %>%
  filter(units > 0)

# Top 12 classes + Other (use Class waffle palette)
class_order <- station_class %>%
  group_by(class) %>%
  summarise(total = sum(mean_RA), .groups = "drop") %>%
  arrange(-total) %>% pull(class)

n_top_class <- 12
top_classes <- class_order[seq_len(min(n_top_class, length(class_order)))]
station_class <- station_class %>%
  mutate(class_grp = ifelse(class %in% top_classes, class, "Other")) %>%
  group_by(station, class_grp) %>%
  summarise(units = sum(units), .groups = "drop") %>%
  rename(class = class_grp)

station_class$class <- factor(station_class$class, levels = c(top_classes, "Other"))

# ==============================================================================
# Section 4: Station labels with cluster + n
# ==============================================================================

station_label <- dbo3_m %>%
  group_by(station) %>%
  summarise(n_samples = n(),
            clusters  = paste(sort(unique(as.character(cluster11))), collapse = ","),
            .groups = "drop") %>%
  mutate(label = paste0(station, "\n(n=", n_samples, "; ", clusters, ")"))

station_class <- station_class %>%
  left_join(station_label %>% select(station, label), by = "station")
station_class$label <- factor(station_class$label,
  levels = station_label$label[order(station_label$station)])

# Same colour scheme as the Class waffle for consistency
class_emphasis <- c(
  "Alphaproteobacteria" = "#E41A1C",
  "Gammaproteobacteria" = "#FF7F00",
  "Flavobacteriia"      = "#377EB8",
  "Bacteroidia"         = "#4DAF4A",
  "Verrucomicrobiae"    = "#984EA3",
  "Planctomycetacia"    = "#A65628",
  "Mollicutes"          = "#F781BF",
  "Nitrososphaeria"     = "#17BECF",
  "Actinobacteria"      = "#FFFF33",
  "Spirochaetia"        = "#666666",
  "Acidimicrobiia"      = "#A6CEE3"
)
remaining <- setdiff(top_classes, names(class_emphasis))
fallback  <- suppressWarnings(
  brewer.pal(max(3, length(remaining)), "Set3"))[seq_along(remaining)]
class_colors <- c(class_emphasis[intersect(top_classes, names(class_emphasis))],
                  setNames(fallback, remaining),
                  "Other" = "gray70")

# ==============================================================================
# Section 5: Waffle plot
# ==============================================================================

dir.create(here::here("output", "survey", "waffle"),
           showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_dbo3_class_waffle.pdf"),
    width = 16, height = 10)

print(
  ggplot(station_class, aes(fill = class, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ label, nrow = 2) +
    coord_equal() +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(strip.text  = element_text(face = "bold", size = 16,
                                     margin = margin(t = 6, b = 6)),
          legend.text  = element_text(size = 14),
          legend.title = element_text(size = 15, face = "bold"))
)

# ==============================================================================
# Section 6: Per-cluster (DBO3 samples only) class composition
# 各クラスター内の DBO3 サンプルのみで Class 組成
# ==============================================================================

cluster_class_dbo3 <- sample_class_dbo3 %>%
  group_by(cluster11, class) %>%
  summarise(mean_RA = mean(sample_RA), .groups = "drop") %>%
  group_by(cluster11) %>%
  mutate(pct   = mean_RA / sum(mean_RA) * 100,
         units = round(pct)) %>%
  ungroup() %>%
  filter(units > 0)

# Re-use the same Top 12 classes selected from station_class so palette matches
cluster_class_dbo3 <- cluster_class_dbo3 %>%
  mutate(class_grp = ifelse(class %in% top_classes, class, "Other")) %>%
  group_by(cluster11, class_grp) %>%
  summarise(units = sum(units), .groups = "drop") %>%
  rename(class = class_grp)
cluster_class_dbo3$class <- factor(cluster_class_dbo3$class,
                                   levels = c(top_classes, "Other"))

# Cluster labels with (n=X) — DBO3 samples count
n_per_cluster_dbo3 <- dbo3_m %>% group_by(cluster11) %>%
  summarise(n = n(), .groups = "drop")
label_map <- setNames(paste0(n_per_cluster_dbo3$cluster11,
                             " (n=", n_per_cluster_dbo3$n, ")"),
                      as.character(n_per_cluster_dbo3$cluster11))
cluster_class_dbo3$cluster_label <- factor(
  label_map[as.character(cluster_class_dbo3$cluster11)],
  levels = label_map[as.character(n_per_cluster_dbo3$cluster11)]
)

print(
  ggplot(cluster_class_dbo3, aes(fill = class, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ cluster_label, nrow = 2) +
    coord_equal() +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    labs(title = "Class composition of DBO3 samples within each cluster",
         subtitle = "Each waffle = 100 units; n = number of DBO3 samples in cluster",
         x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text    = element_text(face = "bold", size = 11,
                                       margin = margin(t = 6, b = 6)),
          legend.text   = element_text(size = 9))
)

dev.off()

# CSV summary
write.csv(station_class %>% select(station, class, units),
  here::here("output", "survey", "waffle", "ASGARD_dbo3_class_waffle.csv"),
  row.names = FALSE)
write.csv(cluster_class_dbo3 %>% select(cluster11, class, units),
  here::here("output", "survey", "waffle", "ASGARD_dbo3_per_cluster_waffle.csv"),
  row.names = FALSE)

# Bacteroidia % per station (highlight)
bact_pct <- station_class %>%
  filter(class == "Bacteroidia") %>%
  select(station, Bacteroidia_pct = units) %>%
  arrange(station)

cat("\n=== Bacteroidia % per DBO3 station ===\n")
print(as.data.frame(bact_pct), row.names = FALSE)

message("\nS22_dbo3_waffle.R: done.")
message("  PDF: output/survey/waffle/ASGARD_dbo3_class_waffle.pdf")
message("  CSV: output/survey/waffle/ASGARD_dbo3_class_waffle.csv")
