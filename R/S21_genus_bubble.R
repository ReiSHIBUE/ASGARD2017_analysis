### S21_genus_bubble.R
### ASGARD 2017 Survey — Top 100 genus bubble plot
### 上位 100 属の bubble plot（11 sample clusters）
###
### REQUIRES (from 00_setup.R, S01, S02):
###   asgard_filtered  - 181×258 ASV proportion matrix
###   clusnum11        - 11-cluster assignments (factor)
###   hier_levels_11   - ordered cluster names
###   cc11             - 11-colour palette
###   shorternames     - "Order; Family; Genus; ESV_N" vector
###   fullnameboot     - taxonomy + bootstrap df
###
### OUTPUT:
###   output/survey/bubble/ASGARD_top100_genus_bubble.pdf

library(tidyverse)
library(RColorBrewer)

# ==============================================================================
# Section 1: ASV → Genus / Phylum mapping
# ==============================================================================

asv_idx        <- match(colnames(asgard_filtered), shorternames)
strip_boot     <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
asv_genus      <- strip_boot(fullnameboot$Genus[asv_idx])
asv_phylum     <- strip_boot(fullnameboot$Phylum[asv_idx])
asv_class      <- strip_boot(fullnameboot$Class[asv_idx])
asv_order      <- strip_boot(fullnameboot$Order[asv_idx])
asv_family     <- strip_boot(fullnameboot$Family[asv_idx])
names(asv_genus)  <- colnames(asgard_filtered)
names(asv_phylum) <- colnames(asgard_filtered)
names(asv_class)  <- colnames(asgard_filtered)

# ==============================================================================
# Section 2: Compute per-cluster mean RA per genus
# サンプル単位で genus 合計 → cluster 平均（同じ二段集計）
# ==============================================================================

long_g <- as.data.frame(asgard_filtered) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(genus   = asv_genus[ASV],
         phylum  = asv_phylum[ASV],
         class   = asv_class[ASV],
         cluster = factor(as.character(clusnum11[Sample]), levels = hier_levels_11)) %>%
  group_by(Sample) %>%
  mutate(RA = Abundance / sum(Abundance)) %>%
  ungroup()

# Per-sample sum per genus
sample_genus <- long_g %>%
  group_by(Sample, cluster, phylum, class, genus) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

# Per-cluster mean
cluster_genus <- sample_genus %>%
  group_by(cluster, phylum, class, genus) %>%
  summarise(mean_RA_pct = round(mean(sample_RA) * 100, 3),
            .groups = "drop")

# ==============================================================================
# Section 3: Select Top 100 genera by overall mean RA
# ==============================================================================

overall_genus <- long_g %>%
  group_by(Sample, genus) %>%
  summarise(sample_RA = sum(RA), .groups = "drop") %>%
  group_by(genus) %>%
  summarise(overall_mean = mean(sample_RA), .groups = "drop") %>%
  arrange(-overall_mean)

top100 <- overall_genus %>% slice_head(n = 100) %>% pull(genus)

cat("Top 100 genera selected; first 10:\n")
print(head(overall_genus, 10))

# Limit to top 100
plot_df <- cluster_genus %>%
  filter(genus %in% top100) %>%
  filter(mean_RA_pct > 0)  # drop zero entries (no bubble drawn)

# Phylum for each top-100 genus (for color)
genus_phylum <- plot_df %>%
  group_by(genus) %>%
  summarise(phylum = first(phylum),
            class  = first(class),
            overall = sum(mean_RA_pct),
            .groups = "drop")

# Order genera by phylum, then by overall abundance (top first within phylum)
genus_order <- genus_phylum %>%
  arrange(phylum, -overall) %>%
  pull(genus)
plot_df$genus <- factor(plot_df$genus, levels = rev(genus_order))

# Color by phylum (limit to top 10 phyla, others "Other")
top_phyla <- genus_phylum %>%
  group_by(phylum) %>%
  summarise(total = sum(overall), .groups = "drop") %>%
  arrange(-total) %>%
  slice_head(n = 10) %>%
  pull(phylum)

plot_df <- plot_df %>%
  mutate(phylum_grp = ifelse(phylum %in% top_phyla, phylum, "Other"))
phylum_colors <- setNames(
  c(brewer.pal(min(length(top_phyla), 8), "Set1"),
    if (length(top_phyla) > 8) brewer.pal(length(top_phyla) - 8, "Set2") else NULL,
    "gray70"),
  c(top_phyla, "Other")
)

# ==============================================================================
# Section 4: Bubble plot
# ==============================================================================

dir.create(here::here("output", "survey", "bubble"),
           showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "bubble",
                     "ASGARD_top100_genus_bubble.pdf"),
    width = 12, height = 22)

print(
  ggplot(plot_df, aes(x = cluster, y = genus,
                      size = mean_RA_pct, color = phylum_grp)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(range = c(0.5, 8),
                          name = "Mean RA (%)",
                          breaks = c(0.5, 1, 2, 5, 10, 20)) +
    scale_color_manual(values = phylum_colors, name = "Phylum") +
    labs(title = "Top 100 genera across 11 sample clusters",
         subtitle = paste0("Bubble size = mean relative abundance (%) within cluster; ",
                           "genera ordered by phylum then overall abundance"),
         x = "Cluster", y = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          axis.text.y   = element_text(size = 6),
          axis.text.x   = element_text(face = "bold", size = 11),
          panel.grid.major = element_line(color = "gray92"),
          panel.grid.minor = element_blank(),
          legend.position = "right")
)

dev.off()

# Also save CSV
write.csv(plot_df,
  here::here("output", "survey", "bubble", "ASGARD_top100_genus_bubble.csv"),
  row.names = FALSE)

message("\nS21_genus_bubble.R: done.")
message("  Top 100 genera selected from ", length(unique(asv_genus)), " unique genera")
message("  PDF: output/survey/bubble/ASGARD_top100_genus_bubble.pdf")
message("  CSV: output/survey/bubble/ASGARD_top100_genus_bubble.csv")
