### S08_taxonomy.R
### ASGARD 2017 Survey Site Analysis — Taxonomy Composition
### 分類組成スクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02):
###   asgard_filtered  - 181×258 ASV proportion matrix (shorternames as colnames)
###   meta_asgard      - metadata 181×45
###   clusnum          - sample cluster assignments (length 181)
###   colclusnum       - ASV column cluster assignments (length 258)
###   fullnameboot     - taxonomy + bootstrap df (built in 00_setup.R)
###
### PRODUCES: (PDFs only — no new session objects required downstream)
###
### OUTPUT:
###   output/survey/taxonomy/ASGARD_taxonomy_barchart_survey.pdf
###   output/survey/taxonomy/ASGARD_taxonomy_waffle_survey.pdf
###
### NOTE: waffle package — install via:
###   install.packages("waffle", repos = "https://cinc.rud.is")
###   or: remotes::install_github("hrbrmstr/waffle")

library(tidyverse)
library(vegan)
library(scales)

# ==============================================================================
# Section 1: クラスレベルの分類別集計 / Aggregate ASVs by Order-level taxonomy
# shorternames形式: "Order (boot); Family (boot); Genus (boot); ESV_N"
# Extract Order label from colnames of asgard_filtered
# ==============================================================================

# 各ASV列のOrderを抽出 / Extract Order from ASV column names
asv_order <- sapply(strsplit(colnames(asgard_filtered), "; "), `[`, 1)
names(asv_order) <- colnames(asgard_filtered)

# サンプルxASVのdfをlongにする / Pivot to long format
asgard_long <- asgard_filtered %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(
    Order   = asv_order[ASV],
    cluster = as.factor(clusnum[Sample])
  )

# サンプル内相対存在量を計算 / Compute relative abundance within each sample
asgard_long <- asgard_long %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# クラスター × Order の平均 / Mean relative abundance per cluster × Order
asgard_class_clus <- asgard_long %>%
  group_by(cluster, Order) %>%
  summarise(mean_RA = mean(Relative_Abundance, na.rm = TRUE), .groups = "drop")

# ==============================================================================
# Section 2: Stacked bar chart / 積み上げ棒グラフ
# ==============================================================================

pdf(file = here::here("output", "survey", "taxonomy", "ASGARD_taxonomy_barchart_survey.pdf"),
    width = 20, height = 10)

# クラスター毎のサンプルの積み上げ棒グラフ / Per-sample stacked bar ordered by cluster
asgard_ordered <- asgard_long %>%
  arrange(cluster, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

gg_bar <- ggplot(asgard_ordered,
                 aes(x = Sample, y = Relative_Abundance, fill = Order)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent) +
  facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
  labs(title = "Microbial Community Composition — ASGARD 2017 Survey Samples",
       x     = "Sample (ordered by cluster)",
       y     = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text  = element_text(size = 6))

print(gg_bar)

dev.off()

# ==============================================================================
# Section 3: Waffle chart — phylum per column cluster
# ワッフルチャート (列クラスター毎の門レベル組成)
# ==============================================================================

# waffle パッケージが利用可能か確認 / Check if waffle is available
if (requireNamespace("waffle", quietly = TRUE)) {
  library(waffle)

  # 各ASVの門 (Phylum) を取得 / Get Phylum for each ASV column
  # fullnameboot はrow名がASV配列、Phylum列を持つ (from 00_setup.R)
  asv_phylum <- fullnameboot[
    match(colnames(asgard_filtered), rownames(fullnameboot)), "Phylum"
  ]
  names(asv_phylum) <- colnames(asgard_filtered)

  # 各列クラスターのphylum組成 / Phylum composition per column cluster
  waffle_df <- data.frame(
    ASV       = names(colclusnum),
    colclus   = colclusnum,
    Phylum    = asv_phylum[names(colclusnum)],
    weight    = colMeans(asgard_filtered)[names(colclusnum)]
  ) %>%
    group_by(colclus, Phylum) %>%
    summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
    group_by(colclus) %>%
    mutate(units = round(weight / sum(weight) * 100)) %>%
    ungroup()

  pdf(file = here::here("output", "survey", "taxonomy", "ASGARD_taxonomy_waffle_survey.pdf"),
      width = 20, height = 12)

  gg_waffle <- ggplot(waffle_df, aes(fill = Phylum, values = units)) +
    geom_waffle(n_rows = 10, size = 0.00, colour = NA) +
    facet_wrap(~ colclus) +
    coord_equal() +
    scale_fill_viridis_d(option = "turbo") +
    labs(title = "Phylum Composition per ASV Column Cluster — Survey Samples") +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1))

  print(gg_waffle)
  dev.off()

  message("S08_taxonomy.R: done. Waffle chart written.")
} else {
  message("S08_taxonomy.R: waffle package not installed. Skipping waffle chart.")
  message("  Install with: remotes::install_github('hrbrmstr/waffle')")
}

message("S08_taxonomy.R: done. PDFs in output/survey/taxonomy/")
