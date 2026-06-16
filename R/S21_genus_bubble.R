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
###   output/survey/bubble/ASGARD_top50_genus_bubble.pdf
###   (plus accompanying CSVs)

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

# Combined "Family; Genus" label per ASV (clarifies "uncultured" entries)
asv_fam_gen <- paste(asv_family, asv_genus, sep = "; ")
names(asv_genus)   <- colnames(asgard_filtered)
names(asv_phylum)  <- colnames(asgard_filtered)
names(asv_class)   <- colnames(asgard_filtered)
names(asv_family)  <- colnames(asgard_filtered)
names(asv_fam_gen) <- colnames(asgard_filtered)

# ==============================================================================
# Section 2: Compute per-cluster mean RA per genus
# サンプル単位で genus 合計 → cluster 平均（同じ二段集計）
# ==============================================================================

long_g <- as.data.frame(asgard_filtered) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(fam_gen = asv_fam_gen[ASV],
         phylum  = asv_phylum[ASV],
         class   = asv_class[ASV],
         cluster = factor(as.character(clusnum11[Sample]), levels = hier_levels_11)) %>%
  group_by(Sample) %>%
  mutate(RA = Abundance / sum(Abundance)) %>%
  ungroup()

# Per-sample sum per family-genus combination
sample_genus <- long_g %>%
  group_by(Sample, cluster, phylum, class, fam_gen) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

# Per-cluster mean
cluster_genus <- sample_genus %>%
  group_by(cluster, phylum, class, fam_gen) %>%
  summarise(mean_RA_pct = round(mean(sample_RA) * 100, 3),
            .groups = "drop")

# ==============================================================================
# Section 3: Overall ranking + colour palette for Class
# ==============================================================================

overall_genus <- long_g %>%
  group_by(Sample, fam_gen) %>%
  summarise(sample_RA = sum(RA), .groups = "drop") %>%
  group_by(fam_gen) %>%
  summarise(overall_mean = mean(sample_RA), .groups = "drop") %>%
  arrange(-overall_mean)

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

# ==============================================================================
# Section 4: Bubble-plot helper (parametrized by N)
# ==============================================================================

dir.create(here::here("output", "survey", "bubble"),
           showWarnings = FALSE, recursive = TRUE)

make_bubble <- function(n_top, pdf_h) {
  top_n <- overall_genus %>% slice_head(n = n_top) %>% pull(fam_gen)

  plot_df <- cluster_genus %>%
    filter(fam_gen %in% top_n, mean_RA_pct > 0)

  genus_class <- plot_df %>%
    group_by(fam_gen) %>%
    summarise(phylum = first(phylum), class = first(class),
              overall = sum(mean_RA_pct), .groups = "drop")

  genus_order <- genus_class %>% arrange(class, -overall) %>% pull(fam_gen)
  plot_df$fam_gen <- factor(plot_df$fam_gen, levels = rev(genus_order))

  top_classes <- genus_class %>%
    group_by(class) %>% summarise(total = sum(overall), .groups = "drop") %>%
    arrange(-total) %>% slice_head(n = 12) %>% pull(class)

  plot_df <- plot_df %>%
    mutate(class_grp = ifelse(class %in% top_classes, class, "Other"))

  remaining <- setdiff(top_classes, names(class_emphasis))
  fallback  <- suppressWarnings(
    brewer.pal(max(3, length(remaining)), "Set3"))[seq_along(remaining)]
  class_colors <- c(class_emphasis[intersect(top_classes, names(class_emphasis))],
                    setNames(fallback, remaining),
                    "Other" = "gray70")

  out_pdf <- here::here("output", "survey", "bubble",
                       paste0("ASGARD_top", n_top, "_genus_bubble.pdf"))
  out_csv <- here::here("output", "survey", "bubble",
                       paste0("ASGARD_top", n_top, "_genus_bubble.csv"))

  pdf(file = out_pdf, width = 12, height = pdf_h)
  print(
    ggplot(plot_df, aes(x = cluster, y = fam_gen,
                        size = mean_RA_pct, color = class_grp)) +
      geom_point(alpha = 0.8) +
      scale_size_continuous(range = c(0.5, 8),
                            name = "Mean RA (%)",
                            breaks = c(0.5, 1, 2, 5, 10, 20)) +
      scale_color_manual(values = class_colors, name = "Class") +
      labs(title = paste0("Top ", n_top, " genera across 11 sample clusters"),
           subtitle = paste0("Bubble size = mean relative abundance (%) within cluster; ",
                             "genera ordered by class then overall abundance"),
           x = "Cluster", y = NULL) +
      theme_bw(base_size = 9) +
      theme(plot.title    = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            axis.text.y   = element_text(size = 7),
            axis.text.x   = element_text(face = "bold", size = 11),
            panel.grid.major = element_line(color = "gray92"),
            panel.grid.minor = element_blank(),
            legend.position = "right")
  )
  dev.off()
  write.csv(plot_df, out_csv, row.names = FALSE)
  message("  PDF: ", sub(".*/", "output/survey/bubble/", out_pdf))
}

# ==============================================================================
# Section 5: Generate Top 100 and Top 50 bubble plots
# ==============================================================================

make_bubble(n_top = 100, pdf_h = 22)
make_bubble(n_top = 50,  pdf_h = 14)

# ==============================================================================
# Section 6: 3 bloom-responsive classes with objective filter
#   - Bacteroidia, Gammaproteobacteria, Alphaproteobacteria
#   - Filter: occurrence > 10% AND cumulative within-class RA ≤ 95%
# ==============================================================================

target_classes_3 <- c("Bacteroidia", "Gammaproteobacteria", "Alphaproteobacteria")

keep_asvs_3   <- names(asv_class)[asv_class %in% target_classes_3]
mat3          <- asgard_filtered[, keep_asvs_3, drop = FALSE]
gen_vec_3     <- asv_fam_gen[keep_asvs_3]
class_vec_3   <- asv_class[keep_asvs_3]
gen_levels_3  <- unique(gen_vec_3)
genus_mat_3   <- sapply(gen_levels_3, function(g) {
  asvs <- names(gen_vec_3)[gen_vec_3 == g]
  if (length(asvs) == 1) mat3[, asvs] else rowSums(mat3[, asvs, drop = FALSE])
})
rownames(genus_mat_3) <- rownames(mat3)

gen2class_3       <- tapply(class_vec_3, gen_vec_3,
                            function(x) names(sort(table(x), decreasing = TRUE))[1])
sample_total_3cl  <- rowSums(genus_mat_3)
genus_RA_3cl      <- sweep(genus_mat_3, 1, sample_total_3cl, "/")
genus_RA_3cl[is.na(genus_RA_3cl)] <- 0

genus_summary_3 <- tibble(
  fam_gen     = colnames(genus_mat_3),
  class       = gen2class_3[colnames(genus_mat_3)],
  occurrence  = colMeans(genus_mat_3 > 0),
  mean_RA     = colMeans(genus_RA_3cl)
) %>%
  filter(occurrence > 0.10) %>%
  arrange(-mean_RA) %>%
  mutate(cumRA = cumsum(mean_RA) / sum(mean_RA))

# 累積95%まで（超える直前まで含む）
selected_3 <- genus_summary_3 %>%
  filter(cumRA <= 0.95 | lag(cumRA, default = 0) < 0.95)

selected_genera_3 <- selected_3$fam_gen

plot_df_3 <- as.data.frame(genus_mat_3) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "fam_gen", values_to = "abundance") %>%
  mutate(cluster = factor(as.character(clusnum11[Sample]), levels = hier_levels_11),
         class   = gen2class_3[fam_gen],
         RA      = abundance / sample_total_3cl[Sample]) %>%
  filter(fam_gen %in% selected_genera_3) %>%
  group_by(cluster, fam_gen, class) %>%
  summarise(mean_RA_pct = round(mean(RA, na.rm = TRUE) * 100, 3), .groups = "drop") %>%
  filter(mean_RA_pct > 0)

gen_order_3 <- selected_3 %>% arrange(class, -mean_RA) %>% pull(fam_gen)
plot_df_3$fam_gen <- factor(plot_df_3$fam_gen, levels = rev(gen_order_3))

class_colors_3 <- c(
  "Alphaproteobacteria" = "#E41A1C",
  "Gammaproteobacteria" = "#FF7F00",
  "Bacteroidia"         = "#4DAF4A"
)

n_genus_3 <- length(unique(plot_df_3$fam_gen))
pdf_h_3   <- max(8, n_genus_3 * 0.32)

out_pdf_3 <- here::here("output", "survey", "bubble",
                        "ASGARD_3classes_objective_bubble.pdf")
out_csv_3 <- here::here("output", "survey", "bubble",
                        "ASGARD_3classes_objective_bubble.csv")

pdf(file = out_pdf_3, width = 13, height = pdf_h_3)
print(
  ggplot(plot_df_3, aes(x = cluster, y = fam_gen,
                        size = mean_RA_pct, color = class)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(1, 12),
                          name = "Mean RA within 3 classes (%)",
                          breaks = c(0.5, 1, 2, 5, 10, 20)) +
    scale_color_manual(values = class_colors_3, name = "Class") +
    labs(title = "Major genera in Bacteroidia / Gamma / Alphaproteobacteria",
         subtitle = paste0("Filter: occurrence > 10% AND cumulative RA <= 95% (",
                           n_genus_3, " genera)"),
         x = "Cluster", y = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 18),
          plot.subtitle = element_text(size = 12),
          axis.text.y   = element_text(size = 11),
          axis.text.x   = element_text(face = "bold", size = 15),
          legend.title  = element_text(size = 12),
          legend.text   = element_text(size = 11),
          panel.grid.major = element_line(color = "gray92"),
          panel.grid.minor = element_blank(),
          legend.position = "right")
)
dev.off()
write.csv(plot_df_3, out_csv_3, row.names = FALSE)
message("  PDF: output/survey/bubble/ASGARD_3classes_objective_bubble.pdf (",
        n_genus_3, " genera)")

message("\nS21_genus_bubble.R: done.")
message("  Total unique family;genus combinations: ", length(unique(asv_fam_gen)))
