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
###   output/survey/waffle/ASGARD_taxonomy_waffle_11clusters.pdf
###   output/survey/waffle/ASGARD_taxonomy_waffle_6assemblages.pdf
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
asv_idx       <- match(colnames(asgard_filtered), shorternames)
asv_phylum_raw <- fullnameboot$Phylum[asv_idx]
asv_order_raw  <- fullnameboot$Order[asv_idx]

# "Name (bootstrap)" → "Name" にする / Strip bootstrap suffix
strip_boot <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
asv_phylum <- strip_boot(asv_phylum_raw)
asv_order  <- strip_boot(asv_order_raw)
names(asv_phylum) <- colnames(asgard_filtered)
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

# 各クラスターの平均RA per Phylum / Mean RA per Phylum within each sample cluster
waffle_clu_df <- asgard_long %>%
  group_by(cluster, Phylum) %>%
  summarise(mean_RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
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

waffle_col_df <- data.frame(
  ASV    = colnames(asgard_filtered),
  Phylum = asv_phylum,
  assemblage = paste0("Assemblage ", as.character(colclusnum[colnames(asgard_filtered)])),
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
  levels = paste0("Assemblage ", sort(unique(colclusnum))))

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
          strip.text    = element_text(face = "bold", size = 12),
          legend.text   = element_text(size = 9))
)

dev.off()
message("  PDF: output/survey/waffle/ASGARD_taxonomy_waffle_6assemblages.pdf")

}  # end if(waffle installed)

message("\nS08_taxonomy.R: done.")
message("  PDF: output/survey/taxonomy/ASGARD_taxonomy_barchart_11clusters.pdf")
