### S08_taxonomy.R
### ASGARD 2017 Survey Site Analysis — Taxonomy Composition
### 分類組成スクリプト（サーベイサイト）
###
### REQUIRES (from 00_setup.R, S01, S02):
###   asgard_filtered  - 181×258 ASV proportion matrix (shorternames as colnames)
###   meta_asgard      - metadata 181×45
###   clusnum11        - 11 sample cluster assignments (factor, length 181)
###   hier_levels_11   - ordered 11-cluster names
###   cc11             - 11-colour palette
###   colclusnum       - 6 ASV column cluster assignments (length 258)
###   shorternames     - "Order; Family; Genus; ESV_N" vector (length 18535)
###   fullnameboot     - taxonomy + bootstrap df (length 18535)
###
### OUTPUT:
###   output/survey/taxonomy/ASGARD_taxonomy_barchart_11clusters.pdf
###   output/survey/taxonomy/ASGARD_class_big3_per_11clusters.csv
###   output/survey/waffle/ASGARD_taxonomy_waffle_11clusters.pdf           (Phylum)
###   output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_class.pdf     (Class)
###   output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_order.pdf     (Order)
###   output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_family.pdf    (Family)
###   output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_genus.pdf     (Genus)
###   output/survey/waffle/ASGARD_taxonomy_waffle_6assemblages.pdf         (Phylum)
###   output/survey/waffle/ASGARD_taxonomy_waffle_6assemblages_class.pdf   (Class)
###
### NOTE: waffle package — install via:
###   install.packages("waffle", repos = "https://cinc.rud.is")
###   or: remotes::install_github("hrbrmstr/waffle")

library(tidyverse)
library(scales)
library(RColorBrewer)

dir.create(here::here("output", "survey", "taxonomy"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("output", "survey", "waffle"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: ASV → Phylum / Order マッピング
# Map each filtered ASV column to its Phylum and Order
# ==============================================================================

# shorternames の位置でfullnamebootを引く / Look up taxonomy by shortername index
asv_idx        <- match(colnames(asgard_filtered), shorternames)
asv_phylum_raw <- fullnameboot$Phylum[asv_idx]
asv_class_raw  <- fullnameboot$Class[asv_idx]
asv_order_raw  <- fullnameboot$Order[asv_idx]
asv_family_raw <- fullnameboot$Family[asv_idx]
asv_genus_raw  <- fullnameboot$Genus[asv_idx]

# "Name (bootstrap)" → "Name" にする / Strip bootstrap suffix
strip_boot <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
asv_phylum <- strip_boot(asv_phylum_raw)
asv_class  <- strip_boot(asv_class_raw)
asv_order  <- strip_boot(asv_order_raw)
asv_family <- strip_boot(asv_family_raw)
asv_genus  <- strip_boot(asv_genus_raw)
names(asv_phylum) <- colnames(asgard_filtered)
names(asv_family) <- colnames(asgard_filtered)
names(asv_genus)  <- colnames(asgard_filtered)
names(asv_class)  <- colnames(asgard_filtered)
names(asv_order)  <- colnames(asgard_filtered)

# ==============================================================================
# Section 2: 11クラスター別 Order レベル積み上げ棒グラフ
# Stacked bar chart per sample (ordered by 11-cluster)
# ==============================================================================

asgard_long <- as.data.frame(asgard_filtered) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(
    Order   = asv_order[ASV],
    Phylum  = asv_phylum[ASV],
    cluster = factor(as.character(clusnum11[Sample]), levels = hier_levels_11)
  ) %>%
  group_by(Sample) %>%
  mutate(RA = Abundance / sum(Abundance)) %>%
  ungroup()

# Cluster labels with sample counts (used by 11-cluster waffle facets)
cluster_n_tab        <- table(clusnum11)
make_cluster_label   <- function(k) paste0(k, " (n=", cluster_n_tab[as.character(k)], ")")
cluster_label_levels <- make_cluster_label(hier_levels_11)

pdf(file = here::here("output", "survey", "taxonomy", "ASGARD_taxonomy_barchart_11clusters.pdf"),
    width = 20, height = 10)

asgard_ordered <- asgard_long %>%
  arrange(cluster, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

print(
  ggplot(asgard_ordered, aes(x = Sample, y = RA, fill = Order)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = percent) +
    facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
    labs(title = "Microbial Community Composition — 11 clusters",
         x = "Sample (ordered by cluster)", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text  = element_text(size = 6),
          strip.text   = element_text(face = "bold", size = 11))
)

dev.off()

# ==============================================================================
# Section 3: Waffle chart — 11 sample clusters
# 11 サンプルクラスター毎の Phylum 組成
# ==============================================================================

if (!requireNamespace("waffle", quietly = TRUE)) {
  message("S08_taxonomy.R: waffle package not installed. Skipping waffle charts.")
  message("  Install: remotes::install_github('hrbrmstr/waffle')")
} else {

library(waffle)

# ---- 3a: 11 sample cluster waffle (Phylum) ----

# 正しい集計: per-sample に Phylum 合計を出してから cluster mean を取る
# (mean(RA) を直接取ると ASV 数の少ない分類群が過大評価されるため)
sample_phylum_sum <- asgard_long %>%
  group_by(Sample, cluster, Phylum) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

waffle_clu_df <- sample_phylum_sum %>%
  group_by(cluster, Phylum) %>%
  summarise(mean_RA = mean(sample_RA), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(pct = mean_RA / sum(mean_RA) * 100,
         units = round(pct)) %>%
  ungroup() %>%
  filter(units > 0)

# Phylum を全体存在量順に並べる / Order Phyla by overall abundance
phylum_order <- waffle_clu_df %>%
  group_by(Phylum) %>%
  summarise(total = sum(mean_RA), .groups = "drop") %>%
  arrange(-total) %>%
  pull(Phylum)

# 上位 n_top + "Other" にまとめる / Keep top phyla, lump rest into "Other"
n_top <- 10
top_phyla <- phylum_order[seq_len(min(n_top, length(phylum_order)))]
waffle_clu_df <- waffle_clu_df %>%
  mutate(Phylum_grp = ifelse(Phylum %in% top_phyla, Phylum, "Other")) %>%
  group_by(cluster, Phylum_grp) %>%
  summarise(units = sum(units), .groups = "drop") %>%
  rename(Phylum = Phylum_grp)

waffle_clu_df$Phylum <- factor(waffle_clu_df$Phylum,
                               levels = c(top_phyla, "Other"))
# Add sample counts to cluster panel labels
waffle_clu_df$cluster <- factor(make_cluster_label(as.character(waffle_clu_df$cluster)),
                                levels = cluster_label_levels)

# Color palette: top phyla colored + Other gray
phylum_colors <- setNames(
  c(brewer.pal(min(n_top, 9), "Set1"),
    if (n_top > 9) brewer.pal(n_top - 9, "Set3") else NULL,
    "gray70"),
  c(top_phyla, "Other")
)

pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_11clusters.pdf"),
    width = 16, height = 12)

print(
  ggplot(waffle_clu_df, aes(fill = Phylum, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ cluster, nrow = 3) +
    coord_equal() +
    scale_fill_manual(values = phylum_colors, drop = FALSE) +
    labs(title = "Phylum Composition per Sample Cluster (11 clusters)",
         subtitle = "Each waffle = 100 units (~ 100% relative abundance)",
         x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text    = element_text(face = "bold", size = 12),
          legend.text   = element_text(size = 9))
)

dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_11clusters.pdf")

# ==============================================================================
# Section 4: Waffle chart — 6 ASV column clusters (assemblages)
# 6 列クラスター毎の Phylum 組成（assemblage 内の各 ASV を平均存在量で重み付け）
# ==============================================================================

# 各 ASV の全体平均 RA を重みとして Phylum 別に集計
asv_mean_ra <- colMeans(asgard_filtered)  # length 258

# Assemblage 別 ASV 数 / ASV counts per assemblage for panel labels
assem_n <- table(colclusnum)
assem_label <- function(k) paste0("Assemblage ", k, " (n=", assem_n[as.character(k)], ")")

waffle_col_df <- data.frame(
  ASV    = colnames(asgard_filtered),
  Phylum = asv_phylum,
  assemblage = assem_label(as.character(colclusnum[colnames(asgard_filtered)])),
  weight = asv_mean_ra,
  stringsAsFactors = FALSE
) %>%
  group_by(assemblage, Phylum) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  group_by(assemblage) %>%
  mutate(pct = weight / sum(weight) * 100,
         units = round(pct)) %>%
  ungroup() %>%
  filter(units > 0)

# Phylum を全体存在量で並べる / Order Phyla by overall weight
col_phylum_order <- waffle_col_df %>%
  group_by(Phylum) %>%
  summarise(total = sum(weight), .groups = "drop") %>%
  arrange(-total) %>%
  pull(Phylum)

n_top_col <- 10
top_phyla_col <- col_phylum_order[seq_len(min(n_top_col, length(col_phylum_order)))]
waffle_col_df <- waffle_col_df %>%
  mutate(Phylum_grp = ifelse(Phylum %in% top_phyla_col, Phylum, "Other")) %>%
  group_by(assemblage, Phylum_grp) %>%
  summarise(units = sum(units), .groups = "drop") %>%
  rename(Phylum = Phylum_grp)

waffle_col_df$Phylum <- factor(waffle_col_df$Phylum,
                               levels = c(top_phyla_col, "Other"))

col_phylum_colors <- setNames(
  c(brewer.pal(min(n_top_col, 9), "Set1"),
    if (n_top_col > 9) brewer.pal(n_top_col - 9, "Set3") else NULL,
    "gray70"),
  c(top_phyla_col, "Other")
)

waffle_col_df$assemblage <- factor(waffle_col_df$assemblage,
  levels = assem_label(sort(unique(colclusnum))))

pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_6assemblages.pdf"),
    width = 14, height = 10)

print(
  ggplot(waffle_col_df, aes(fill = Phylum, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ assemblage, nrow = 2) +
    coord_equal() +
    scale_fill_manual(values = col_phylum_colors, drop = FALSE) +
    labs(title = "Phylum Composition per ASV Assemblage (6 column clusters)",
         subtitle = "Each waffle = 100 units (~ 100% within-assemblage abundance)",
         x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text    = element_text(face = "bold", size = 12,
                                       margin = margin(t = 6, b = 6)),
          legend.text   = element_text(size = 9))
)

dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_6assemblages.pdf")

# ==============================================================================
# Section 5: Waffle chart — 11 sample clusters at CLASS level
# 11 サンプルクラスター毎の Class レベル組成
# ==============================================================================

# Class情報をlongに追加 / Add Class to long table
long_class <- as.data.frame(asgard_filtered) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(Class = asv_class[ASV],
         cluster = factor(make_cluster_label(as.character(clusnum11[Sample])),
                          levels = cluster_label_levels)) %>%
  group_by(Sample) %>%
  mutate(RA = Abundance / sum(Abundance)) %>%
  ungroup()

# per-sample に Class 合計 → cluster mean (Phylum版と同じ正しい集計)
sample_class_sum <- long_class %>%
  group_by(Sample, cluster, Class) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

waffle_clu_class <- sample_class_sum %>%
  group_by(cluster, Class) %>%
  summarise(mean_RA = mean(sample_RA), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(pct = mean_RA / sum(mean_RA) * 100,
         units = round(pct)) %>%
  ungroup() %>%
  filter(units > 0)

class_order <- waffle_clu_class %>%
  group_by(Class) %>%
  summarise(total = sum(mean_RA), .groups = "drop") %>%
  arrange(-total) %>%
  pull(Class)

n_top_class <- 12
top_classes <- class_order[seq_len(min(n_top_class, length(class_order)))]
waffle_clu_class <- waffle_clu_class %>%
  mutate(Class_grp = ifelse(Class %in% top_classes, Class, "Other")) %>%
  group_by(cluster, Class_grp) %>%
  summarise(units = sum(units), .groups = "drop") %>%
  rename(Class = Class_grp)

waffle_clu_class$Class <- factor(waffle_clu_class$Class,
                                 levels = c(top_classes, "Other"))

# Class color palette: emphasize bloom-related taxa (Flavobacteriia, Gammaproteobacteria, Alphaproteobacteria)
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
# Fallback for classes not in emphasis
remaining_classes <- setdiff(top_classes, names(class_emphasis))
fallback_palette  <- brewer.pal(max(3, length(remaining_classes)), "Set3")[seq_along(remaining_classes)]
class_colors <- c(class_emphasis[intersect(top_classes, names(class_emphasis))],
                  setNames(fallback_palette, remaining_classes),
                  "Other" = "gray70")

pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_11clusters_class.pdf"),
    width = 16, height = 12)

print(
  ggplot(waffle_clu_class, aes(fill = Class, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ cluster, nrow = 3) +
    coord_equal() +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    labs(title = "Class Composition per Sample Cluster (11 clusters)",
         subtitle = "Each waffle = 100 units (~ 100% relative abundance); top 12 classes shown",
         x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text    = element_text(face = "bold", size = 12),
          legend.text   = element_text(size = 9))
)

dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_class.pdf")

# ==============================================================================
# Section 5b: CSV table — Alpha / Gamma / Bacteroidia per 11 cluster
# 11クラスター毎の Big-3 (Alpha + Gamma + Bacteroidia) % をテーブルとして保存
# ==============================================================================

# Build a clean two-step aggregation directly from clusnum11 (avoids the
# "(n=X)"-labelled factor levels used in the waffle plots)
big3_long <- as.data.frame(asgard_filtered) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(Class   = asv_class[ASV],
         cluster = factor(as.character(clusnum11[Sample]), levels = hier_levels_11)) %>%
  group_by(Sample) %>% mutate(RA = Abundance / sum(Abundance)) %>% ungroup() %>%
  group_by(Sample, cluster, Class) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

big3_table <- big3_long %>%
  filter(Class %in% c("Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia")) %>%
  group_by(cluster, Class) %>%
  summarise(mean_pct = round(mean(sample_RA) * 100, 1), .groups = "drop") %>%
  pivot_wider(names_from = Class, values_from = mean_pct, values_fill = 0) %>%
  rename(Alpha = Alphaproteobacteria,
         Gamma = Gammaproteobacteria,
         Bact  = Bacteroidia)

# Sample counts per cluster
n_per_cluster <- as.data.frame(table(clusnum11)) %>%
  rename(cluster = clusnum11, n = Freq)
n_per_cluster$cluster <- as.character(n_per_cluster$cluster)
big3_table$cluster <- as.character(big3_table$cluster)

big3_table <- big3_table %>%
  left_join(n_per_cluster, by = "cluster") %>%
  mutate(Big3_total  = round(Alpha + Gamma + Bact, 1),
         Other       = round(100 - Big3_total, 1),
         Alpha_share = round(Alpha / (Alpha + Gamma) * 100, 1)) %>%
  select(cluster, n, Alpha, Gamma, Bact, Big3_total, Other, Alpha_share)

big3_table <- big3_table[match(hier_levels_11, big3_table$cluster), ]

cat("\n--- 11 cluster: Alpha / Gamma / Bacteroidia (mean %) ---\n")
print(as.data.frame(big3_table), row.names = FALSE)

write.csv(big3_table,
  here::here("output", "survey", "taxonomy",
             "ASGARD_class_big3_per_11clusters.csv"),
  row.names = FALSE)
message("  CSV: output/survey/taxonomy/ASGARD_class_big3_per_11clusters.csv")

# ==============================================================================
# Section 6: Waffle chart — 6 assemblages at CLASS level
# 6 列クラスター毎の Class レベル組成
# ==============================================================================

waffle_col_class <- data.frame(
  ASV    = colnames(asgard_filtered),
  Class  = asv_class,
  assemblage = assem_label(as.character(colclusnum[colnames(asgard_filtered)])),
  weight = asv_mean_ra,
  stringsAsFactors = FALSE
) %>%
  group_by(assemblage, Class) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  group_by(assemblage) %>%
  mutate(pct = weight / sum(weight) * 100,
         units = round(pct)) %>%
  ungroup() %>%
  filter(units > 0)

col_class_order <- waffle_col_class %>%
  group_by(Class) %>%
  summarise(total = sum(weight), .groups = "drop") %>%
  arrange(-total) %>%
  pull(Class)

n_top_col_class <- 12
top_col_classes <- col_class_order[seq_len(min(n_top_col_class, length(col_class_order)))]
waffle_col_class <- waffle_col_class %>%
  mutate(Class_grp = ifelse(Class %in% top_col_classes, Class, "Other")) %>%
  group_by(assemblage, Class_grp) %>%
  summarise(units = sum(units), .groups = "drop") %>%
  rename(Class = Class_grp)

waffle_col_class$Class <- factor(waffle_col_class$Class,
                                 levels = c(top_col_classes, "Other"))
waffle_col_class$assemblage <- factor(waffle_col_class$assemblage,
  levels = assem_label(sort(unique(colclusnum))))

remaining_col <- setdiff(top_col_classes, names(class_emphasis))
fallback_col  <- brewer.pal(max(3, length(remaining_col)), "Set3")[seq_along(remaining_col)]
col_class_colors <- c(class_emphasis[intersect(top_col_classes, names(class_emphasis))],
                      setNames(fallback_col, remaining_col),
                      "Other" = "gray70")

pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_6assemblages_class.pdf"),
    width = 14, height = 10)

print(
  ggplot(waffle_col_class, aes(fill = Class, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ assemblage, nrow = 2) +
    coord_equal() +
    scale_fill_manual(values = col_class_colors, drop = FALSE) +
    labs(title = "Class Composition per ASV Assemblage (6 column clusters)",
         subtitle = "Each waffle = 100 units; top 12 classes shown",
         x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text    = element_text(face = "bold", size = 12,
                                       margin = margin(t = 6, b = 6)),
          legend.text   = element_text(size = 9))
)

dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_6assemblages_class.pdf")

# ==============================================================================
# Section 7-9: Waffle charts — 11 sample clusters at ORDER / FAMILY / GENUS
# 11 サンプルクラスター × Order / Family / Genus レベル
# ==============================================================================

# Helper: build cluster-level waffle data at any taxonomic rank
build_waffle_taxon <- function(taxon_vec, n_top = 15) {
  # Per-sample sum per taxon → per-cluster mean
  long_t <- as.data.frame(asgard_filtered) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
    mutate(taxon   = taxon_vec[ASV],
           cluster = factor(make_cluster_label(as.character(clusnum11[Sample])),
                            levels = cluster_label_levels)) %>%
    group_by(Sample) %>%
    mutate(RA = Abundance / sum(Abundance)) %>%
    ungroup()

  sample_sum <- long_t %>%
    group_by(Sample, cluster, taxon) %>%
    summarise(sample_RA = sum(RA), .groups = "drop")

  cluster_mean <- sample_sum %>%
    group_by(cluster, taxon) %>%
    summarise(mean_RA = mean(sample_RA), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(pct = mean_RA / sum(mean_RA) * 100,
           units = round(pct)) %>%
    ungroup() %>%
    filter(units > 0)

  # Top-n by overall abundance
  taxon_order <- cluster_mean %>%
    group_by(taxon) %>%
    summarise(total = sum(mean_RA), .groups = "drop") %>%
    arrange(-total) %>%
    pull(taxon)

  top_t <- taxon_order[seq_len(min(n_top, length(taxon_order)))]
  cluster_mean <- cluster_mean %>%
    mutate(taxon_grp = ifelse(taxon %in% top_t, taxon, "Other")) %>%
    group_by(cluster, taxon_grp) %>%
    summarise(units = sum(units), .groups = "drop") %>%
    rename(taxon = taxon_grp)
  cluster_mean$taxon <- factor(cluster_mean$taxon, levels = c(top_t, "Other"))
  list(df = cluster_mean, top = top_t)
}

# Helper: build a colour vector with "Other" as gray
make_palette <- function(top_taxa) {
  n <- length(top_taxa)
  cols <- if (n <= 12) {
    suppressWarnings(c(brewer.pal(min(n, 8), "Set1"),
                       brewer.pal(min(max(n - 8, 3), 8), "Set2")))[seq_len(n)]
  } else {
    colorRampPalette(c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"),
                       brewer.pal(8, "Set3")))(n)
  }
  setNames(c(cols, "gray70"), c(top_taxa, "Other"))
}

# Helper: draw waffle
draw_taxon_waffle <- function(wdat, palette, title) {
  ggplot(wdat$df, aes(fill = taxon, values = units)) +
    geom_waffle(n_rows = 10, size = 0.2, colour = "white", flip = TRUE) +
    facet_wrap(~ cluster, nrow = 3) +
    coord_equal() +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(title = title,
         subtitle = "Each waffle = 100 units (~ 100% relative abundance)",
         x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          strip.text    = element_text(face = "bold", size = 12),
          legend.text   = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
}

# ---- Section 7: Order ----
wd_order <- build_waffle_taxon(asv_order, n_top = 15)
pal_order <- make_palette(wd_order$top)
pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_11clusters_order.pdf"),
    width = 18, height = 12)
print(draw_taxon_waffle(wd_order, pal_order,
                        "Order Composition per Sample Cluster (11 clusters)"))
dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_order.pdf")

# ---- Section 8: Family ----
wd_family <- build_waffle_taxon(asv_family, n_top = 15)
pal_family <- make_palette(wd_family$top)
pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_11clusters_family.pdf"),
    width = 18, height = 12)
print(draw_taxon_waffle(wd_family, pal_family,
                        "Family Composition per Sample Cluster (11 clusters)"))
dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_family.pdf")

# ---- Section 9: Genus ----
wd_genus <- build_waffle_taxon(asv_genus, n_top = 15)
pal_genus <- make_palette(wd_genus$top)
pdf(file = here::here("output", "survey", "waffle",
                     "ASGARD_taxonomy_waffle_11clusters_genus.pdf"),
    width = 18, height = 12)
print(draw_taxon_waffle(wd_genus, pal_genus,
                        "Genus Composition per Sample Cluster (11 clusters)"))
dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_11clusters_genus.pdf")

}  # end if(waffle installed)

message("\nS08_taxonomy.R: done.")
message("  PDF: output/survey/taxonomy/ASGARD_taxonomy_barchart_11clusters.pdf")
