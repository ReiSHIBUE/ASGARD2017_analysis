### P12_genus_bubble.R
### ASGARD 2017 Processing — Genus bubble plots (4 clusters)
### 上位 genus の bubble plot（処理サイト、4 クラスタ）
###
### REQUIRES (from 00_setup.R, P01–P03):
###   asgard_filtered_p_hm2 - ASV matrix 78 x 221
###   clusnum_p             - 4-cluster assignment
###   shorternames          - ESV short labels
###   fullnameboot          - taxonomy + bootstrap df
###
### OUTPUT:
###   output_p/bubble/ASGARD_top100_genus_bubble.pdf  (+ csv)
###   output_p/bubble/ASGARD_top50_genus_bubble.pdf   (+ csv)
###   output_p/bubble/ASGARD_3classes_objective_bubble.pdf  (+ csv)

library(tidyverse)
library(RColorBrewer)

dir.create(here::here("output_p", "bubble"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: ASV → Genus / Class mapping
# ==============================================================================

asv_idx_p     <- match(colnames(asgard_filtered_p_hm2), shorternames)
strip_boot    <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
asv_genus_p   <- strip_boot(fullnameboot$Genus[asv_idx_p])
asv_phylum_p  <- strip_boot(fullnameboot$Phylum[asv_idx_p])
asv_class_p   <- strip_boot(fullnameboot$Class[asv_idx_p])
asv_family_p  <- strip_boot(fullnameboot$Family[asv_idx_p])
asv_fam_gen_p <- paste(asv_family_p, asv_genus_p, sep = "; ")
names(asv_genus_p)   <- colnames(asgard_filtered_p_hm2)
names(asv_phylum_p)  <- colnames(asgard_filtered_p_hm2)
names(asv_class_p)   <- colnames(asgard_filtered_p_hm2)
names(asv_family_p)  <- colnames(asgard_filtered_p_hm2)
names(asv_fam_gen_p) <- colnames(asgard_filtered_p_hm2)

# ==============================================================================
# Section 2: Per-cluster mean RA per genus (4 clusters)
# ==============================================================================

cc_p <- c("1" = "#E41A1C", "2" = "#377EB8",
          "3" = "#4DAF4A", "4" = "#984EA3")

long_g_p <- as.data.frame(asgard_filtered_p_hm2) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(fam_gen = asv_fam_gen_p[ASV],
         phylum  = asv_phylum_p[ASV],
         class   = asv_class_p[ASV],
         cluster = factor(as.character(clusnum_p[Sample]),
                          levels = c("1", "2", "3", "4"))) %>%
  group_by(Sample) %>%
  mutate(RA = Abundance / sum(Abundance)) %>%
  ungroup()

sample_genus_p <- long_g_p %>%
  group_by(Sample, cluster, phylum, class, fam_gen) %>%
  summarise(sample_RA = sum(RA), .groups = "drop")

cluster_genus_p <- sample_genus_p %>%
  group_by(cluster, phylum, class, fam_gen) %>%
  summarise(mean_RA_pct = round(mean(sample_RA) * 100, 3),
            .groups = "drop")

# ==============================================================================
# Section 3: Overall ranking + Class colour palette
# ==============================================================================

overall_genus_p <- long_g_p %>%
  group_by(Sample, fam_gen) %>%
  summarise(sample_RA = sum(RA), .groups = "drop") %>%
  group_by(fam_gen) %>%
  summarise(overall_mean = mean(sample_RA), .groups = "drop") %>%
  arrange(-overall_mean)

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
# Section 4: Top-N bubble helper
# ==============================================================================

make_bubble_p <- function(n_top, pdf_h) {
  top_n <- overall_genus_p %>% slice_head(n = n_top) %>% pull(fam_gen)
  plot_df <- cluster_genus_p %>%
    filter(fam_gen %in% top_n, mean_RA_pct > 0)

  genus_class <- plot_df %>%
    group_by(fam_gen) %>%
    summarise(phylum = first(phylum), class = first(class),
              overall = sum(mean_RA_pct), .groups = "drop")
  genus_order <- genus_class %>%
    arrange(class, -overall) %>% pull(fam_gen)
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

  out_pdf <- here::here("output_p", "bubble",
                       paste0("ASGARD_top", n_top, "_genus_bubble.pdf"))
  out_csv <- here::here("output_p", "bubble",
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
      labs(title = paste0("Top ", n_top, " genera across 4 processing clusters"),
           subtitle = "Bubble size = mean relative abundance (%) within cluster",
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
  message("  PDF: output_p/bubble/", basename(out_pdf))
}

make_bubble_p(n_top = 100, pdf_h = 22)
make_bubble_p(n_top = 50,  pdf_h = 14)

# ==============================================================================
# Section 5: 3 bloom-responsive classes + objective filter
#   - Bacteroidia, Gammaproteobacteria, Alphaproteobacteria
#   - occurrence > 10% AND cumulative within-class RA <= 95%
# ==============================================================================

target_classes_3 <- c("Bacteroidia", "Gammaproteobacteria", "Alphaproteobacteria")

keep_asvs_p   <- names(asv_class_p)[asv_class_p %in% target_classes_3]
mat3_p        <- asgard_filtered_p_hm2[, keep_asvs_p, drop = FALSE]
gen_vec_p     <- asv_fam_gen_p[keep_asvs_p]
class_vec_p   <- asv_class_p[keep_asvs_p]
gen_levels_p  <- unique(gen_vec_p)
genus_mat_p   <- sapply(gen_levels_p, function(g) {
  asvs <- names(gen_vec_p)[gen_vec_p == g]
  if (length(asvs) == 1) mat3_p[, asvs] else rowSums(mat3_p[, asvs, drop = FALSE])
})
rownames(genus_mat_p) <- rownames(mat3_p)

gen2class_p       <- tapply(class_vec_p, gen_vec_p,
                            function(x) names(sort(table(x), decreasing = TRUE))[1])
sample_total_3cl_p <- rowSums(genus_mat_p)
genus_RA_3cl_p     <- sweep(genus_mat_p, 1, sample_total_3cl_p, "/")
genus_RA_3cl_p[is.na(genus_RA_3cl_p)] <- 0

genus_summary_p <- tibble(
  fam_gen     = colnames(genus_mat_p),
  class       = gen2class_p[colnames(genus_mat_p)],
  occurrence  = colMeans(genus_mat_p > 0),
  mean_RA     = colMeans(genus_RA_3cl_p)
) %>%
  filter(occurrence > 0.10) %>%
  arrange(-mean_RA) %>%
  mutate(cumRA = cumsum(mean_RA) / sum(mean_RA))

selected_p <- genus_summary_p %>%
  filter(cumRA <= 0.95 | lag(cumRA, default = 0) < 0.95)
selected_genera_p <- selected_p$fam_gen

plot_df_3p <- as.data.frame(genus_mat_p) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "fam_gen", values_to = "abundance") %>%
  mutate(cluster = factor(as.character(clusnum_p[Sample]),
                          levels = c("1", "2", "3", "4")),
         class   = gen2class_p[fam_gen],
         RA      = abundance / sample_total_3cl_p[Sample]) %>%
  filter(fam_gen %in% selected_genera_p) %>%
  group_by(cluster, fam_gen, class) %>%
  summarise(mean_RA_pct = round(mean(RA, na.rm = TRUE) * 100, 3),
            .groups = "drop") %>%
  filter(mean_RA_pct > 0)

gen_order_p <- selected_p %>% arrange(class, -mean_RA) %>% pull(fam_gen)
plot_df_3p$fam_gen <- factor(plot_df_3p$fam_gen, levels = rev(gen_order_p))

class_colors_p <- c(
  "Alphaproteobacteria" = "#E41A1C",
  "Gammaproteobacteria" = "#FF7F00",
  "Bacteroidia"         = "#4DAF4A"
)

n_genus_3p <- length(unique(plot_df_3p$fam_gen))
pdf_h_3p   <- max(8, n_genus_3p * 0.32)

out_pdf_3p <- here::here("output_p", "bubble",
                         "ASGARD_3classes_objective_bubble.pdf")
out_csv_3p <- here::here("output_p", "bubble",
                         "ASGARD_3classes_objective_bubble.csv")

pdf(file = out_pdf_3p, width = 13, height = pdf_h_3p)
print(
  ggplot(plot_df_3p, aes(x = cluster, y = fam_gen,
                         size = mean_RA_pct, color = class)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(1, 12),
                          name = "Mean RA within 3 classes (%)",
                          breaks = c(0.5, 1, 2, 5, 10, 20)) +
    scale_color_manual(values = class_colors_p, name = "Class") +
    labs(title = "Major genera in Bacteroidia / Gamma / Alphaproteobacteria",
         subtitle = paste0("Filter: occurrence > 10% AND cumulative RA <= 95% (",
                           n_genus_3p, " genera)"),
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
write.csv(plot_df_3p, out_csv_3p, row.names = FALSE)
message("  PDF: output_p/bubble/ASGARD_3classes_objective_bubble.pdf (",
        n_genus_3p, " genera)")

message("\nP12_genus_bubble.R: done.")
message("  Total unique family;genus combinations: ", length(unique(asv_fam_gen_p)))
