### P13_ternary_v2.R
### ASGARD 2017 Processing — Improved ternary plots (3 pages)
### 改良版 ternary plot（class色／cluster色／RGB+主要属ラベル）
###
### REQUIRES (from 00_setup.R, P01-P03):
###   ternary_prop_color2   - df with 0.2/3/20 µm proportions per ASV (from P02)
###   asgard_filtered_p_hm2 - ASV matrix (from P01)
###   clusnum_p             - 4-cluster assignment (from P03)
###   shorternames          - ESV short labels (from 00_setup)
###   fullnameboot          - taxonomy df (from 00_setup)
###
### OUTPUT:
###   output_p/ternary/ASGARD_ternary_3pages.pdf

library(ggtern)
library(gridExtra)
library(dplyr)
library(tibble)

dir.create(here::here("output_p", "ternary"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: Per-ASV data frame
# ==============================================================================

asv_ids   <- rownames(ternary_prop_color2)
if (is.null(asv_ids) || all(asv_ids == as.character(seq_along(asv_ids)))) {
  asv_ids <- colnames(asgard_filtered_p_hm2)
}

ternary_df <- data.frame(
  p02 = ternary_prop_color2$`0.2 µm`,
  p03 = ternary_prop_color2$`3 µm`,
  p20 = ternary_prop_color2$`20 µm`,
  asv = asv_ids,
  stringsAsFactors = FALSE
)

# Mean RA across all 78 samples (sample-normalized)
mat_p <- asgard_filtered_p_hm2
mat_p <- sweep(mat_p, 1, rowSums(mat_p), "/")
ternary_df$mean_RA_pct <- 100 * colMeans(mat_p[, ternary_df$asv, drop = FALSE])

# Cluster-dominance per ASV
cluster_means_p <- sapply(c("1","2","3","4"), function(k) {
  smp <- names(clusnum_p)[as.character(clusnum_p) == k]
  colMeans(mat_p[smp, ternary_df$asv, drop = FALSE])
})
ternary_df$cluster_dom <- factor(
  c("1","2","3","4")[apply(cluster_means_p, 1, which.max)],
  levels = c("1","2","3","4"), labels = clus_labels_p)

# Taxonomy
asv_idx     <- match(ternary_df$asv, shorternames)
strip_boot  <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
ternary_df$class   <- strip_boot(fullnameboot$Class[asv_idx])
ternary_df$genus   <- strip_boot(fullnameboot$Genus[asv_idx])
ternary_df$family  <- strip_boot(fullnameboot$Family[asv_idx])
ternary_df$fam_gen <- paste(ternary_df$family, ternary_df$genus, sep = "; ")

# ==============================================================================
# Section 2: Palettes
# ==============================================================================

target_classes <- c("Bacteroidia", "Alphaproteobacteria", "Gammaproteobacteria")
ternary_df$class_grp <- ifelse(ternary_df$class %in% target_classes,
                               ternary_df$class, "Other")
class_levels <- c(target_classes, "Other")
ternary_df$class_grp <- factor(ternary_df$class_grp, levels = class_levels)

class_colors <- c("Bacteroidia"         = "#4DAF4A",
                  "Alphaproteobacteria" = "#E41A1C",
                  "Gammaproteobacteria" = "#FF7F00",
                  "Other"               = "gray70")[class_levels]

cc_p <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                 clus_levels_p)   # FL1, PA1, PA2, PA3

# RGB blend per ASV (same as asv_rgb2 from P02)
ternary_df$rgb <- rgb(ternary_df$p02, ternary_df$p03, ternary_df$p20)

# Top 20 genera labels (Page 3)
ternary_df$label <- ""
top20_idx <- order(-ternary_df$mean_RA_pct)[1:20]
for (i in top20_idx) {
  g <- ternary_df$genus[i]
  if (is.na(g) || g %in% c("", "uncultured", "uncultured_marine_group")) {
    ternary_df$label[i] <- ternary_df$family[i]
  } else {
    ternary_df$label[i] <- g
  }
}

# ==============================================================================
# Section 3: 3-page PDF
# ==============================================================================

base_tern_theme <- theme_bw(base_size = 18) +
  theme(plot.title    = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 14),
        legend.title  = element_text(face = "bold", size = 16),
        legend.text   = element_text(size = 14),
        axis.title    = element_text(face = "bold", size = 18))

p1 <- ggtern(data = ternary_df,
       aes(x = p02, y = p03, z = p20,
           color = class_grp, size = mean_RA_pct)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = class_colors, name = "Class") +
  scale_size_continuous(range = c(1, 10),
                        name = "Mean RA (%)",
                        breaks = c(0.1, 0.5, 1, 2, 5)) +
  labs(title = "ASVs by Class",
       subtitle = "Position = proportion of total RA in each filter fraction",
       x = "0.2 µm", y = "3 µm", z = "20 µm") +
  base_tern_theme + theme_showarrows() + theme_rotate(-90)

p2 <- ggtern(data = ternary_df,
       aes(x = p02, y = p03, z = p20,
           color = cluster_dom, size = mean_RA_pct)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = cc_p, name = "Dominant cluster") +
  scale_size_continuous(range = c(1, 10),
                        name = "Mean RA (%)",
                        breaks = c(0.1, 0.5, 1, 2, 5)) +
  labs(title = "ASVs by dominant sample cluster (k = 4)",
       subtitle = "Each ASV assigned to the cluster with highest mean RA",
       x = "0.2 µm", y = "3 µm", z = "20 µm") +
  base_tern_theme + theme_showarrows() + theme_rotate(-90)

p3 <- ggtern(data = ternary_df,
       aes(x = p02, y = p03, z = p20, size = mean_RA_pct)) +
  geom_point(aes(color = I(rgb)), alpha = 0.85) +
  geom_text(aes(label = label),
            size = 4, fontface = "bold",
            vjust = -1.0, check_overlap = TRUE) +
  scale_size_continuous(range = c(1, 10),
                        name = "Mean RA (%)",
                        breaks = c(0.1, 0.5, 1, 2, 5)) +
  labs(title = "ASVs with RGB blend + top 20 genera",
       subtitle = "Red = 0.2 µm, Green = 3 µm, Blue = 20 µm",
       x = "0.2 µm", y = "3 µm", z = "20 µm") +
  base_tern_theme + theme_showarrows() + theme_rotate(-90)

pdf(file = here::here("output_p", "ternary",
                      "ASGARD_ternary_3pages.pdf"),
    width = 12, height = 11)
print(p1)
print(p2)
print(p3)

# Page 4: 3 panels combined
small_theme <- theme_bw(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_blank(),
        legend.title  = element_text(face = "bold", size = 10),
        legend.text   = element_text(size = 9),
        legend.position = "bottom",
        legend.box      = "vertical")

p1s <- p1 + small_theme +
  labs(title = "(b) ASVs by Class") +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)),
         size  = guide_legend(nrow = 1))
p2s <- p2 + small_theme +
  labs(title = "(a) ASVs by dominant sample cluster (k = 4)") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 4)),
         size  = guide_legend(nrow = 1))
p3s <- p3 + small_theme +
  labs(title = "(c) ASVs with RGB blend + top 20 genera") +
  geom_text(aes(label = label), size = 2.5, fontface = "bold",
            vjust = -1.0, check_overlap = TRUE)

grid.arrange(p2s, p1s, p3s, nrow = 1)

# ============================================================================
# Page 5: alternative panel using the 3-class objective filter
#   (c') = filter to Bacteroidia / Alphaproteobacteria / Gammaproteobacteria
#          with occurrence > 10% AND cumulative within-class RA <= 95%
# ============================================================================

# Compute objective filter for the 3 bloom classes
target_3 <- c("Bacteroidia", "Alphaproteobacteria", "Gammaproteobacteria")
keep_3   <- ternary_df$class %in% target_3
sub3_df  <- ternary_df[keep_3, ]

# Per-sample 3-class normalized RA
keep_asvs_3 <- sub3_df$asv
mat3        <- asgard_filtered_p_hm2[, keep_asvs_3, drop = FALSE]
sample_total_3 <- rowSums(mat3)
ra3 <- sweep(mat3, 1, sample_total_3, "/")
ra3[is.na(ra3)] <- 0

sub3_df$occurrence <- colMeans(mat3 > 0)
sub3_df$mean_RA_3cl <- colMeans(ra3)

sub3_filt <- sub3_df %>%
  filter(occurrence > 0.10) %>%
  arrange(-mean_RA_3cl) %>%
  mutate(cumRA = cumsum(mean_RA_3cl) / sum(mean_RA_3cl))
sub3_final <- sub3_filt %>%
  filter(cumRA <= 0.95 | lag(cumRA, default = 0) < 0.95)

n_objective <- nrow(sub3_final)
message("Objective filter (page 5 c'): ", n_objective, " ASVs")

# Page-5 (c'): labelled with family;genus
sub3_final$label_obj <- ifelse(
  is.na(sub3_final$genus) |
    sub3_final$genus %in% c("", "uncultured", "uncultured_marine_group"),
  sub3_final$family, sub3_final$genus
)

p3_obj <- ggtern(data = sub3_final,
       aes(x = p02, y = p03, z = p20, size = 100 * mean_RA_3cl)) +
  geom_point(aes(color = I(rgb)), alpha = 0.85) +
  geom_text(aes(label = label_obj), size = 3, fontface = "bold",
            vjust = -1.0, check_overlap = TRUE) +
  scale_size_continuous(range = c(1, 10),
                        name = "Mean RA in 3 classes (%)",
                        breaks = c(0.5, 1, 2, 5, 10)) +
  labs(title = "ASVs by 3-class objective filter",
       subtitle = paste0("Occurrence > 10% AND cumulative RA <= 95% (n = ",
                         n_objective, " ASVs)"),
       x = "0.2 µm", y = "3 µm", z = "20 µm") +
  base_tern_theme + theme_showarrows() + theme_rotate(-90)

print(p3_obj)

# Page 6: combined (a) cluster + (b) class + (c') objective
p2s2 <- p2 + small_theme +
  labs(title = "(a) ASVs in four clusters") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 4)),
         size  = guide_legend(nrow = 1))
p1s2 <- p1 + small_theme +
  labs(title = "(b) ASVs by Class") +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)),
         size  = guide_legend(nrow = 1))
p3obj_s <- p3_obj + small_theme +
  labs(title = "(c) Objective-filter ASVs (3 classes)") +
  geom_text(aes(label = label_obj), size = 2.5, fontface = "bold",
            vjust = -1.0, check_overlap = TRUE)

grid.arrange(p2s2, p1s2, p3obj_s, nrow = 1)

# Pages 7-9: (A) clusters and (B) class, each as a full-page ternary with
# labelled arrows (like page 5), plus a combined (A)+(B) sheet. ggtern arrows are
# clipped by grid.arrange, so the combined page rasterizes each standalone panel
# (arrows preserved) before placing them side by side.
p2s7 <- p2 + labs(title = "(A) ASVs in four clusters")
p1s7 <- p1 + labs(title = "(B) ASVs by Class")

print(p2s7)   # Page 7: (A) standalone
print(p1s7)   # Page 8: (B) standalone

# Page 9: combined — rasterize each standalone panel, then arrange side by side
.tf_a <- tempfile(fileext = ".png"); .tf_b <- tempfile(fileext = ".png")
png(.tf_a, width = 1100, height = 1000, res = 130); print(p2s7); dev.off()
png(.tf_b, width = 1100, height = 1000, res = 130); print(p1s7); dev.off()
grid.arrange(grid::rasterGrob(png::readPNG(.tf_a), interpolate = TRUE),
             grid::rasterGrob(png::readPNG(.tf_b), interpolate = TRUE), nrow = 1)

dev.off()

message("\nP13_ternary_v2.R: done.")
message("  PDF: output_p/ternary/ASGARD_ternary_3pages.pdf")
