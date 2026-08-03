### P17_fraction_logratio.R
### ASGARD 2017 Processing — Double-volcano plot of the filter-fraction ratios
### フィルター画分間の比を2軸にとった「ダブル・ボルケーノ」図
###
### For every ASV the mean relative abundance within each filter fraction is
### computed, and two log2 ratios are plotted against each other:
###
###   x = log2( 3 µm / 0.2 µm )    small-particle vs free-living step
###   y = log2( 20 µm / 3 µm )     large-particle vs small-particle step
###
### Point size is sqrt(mean relative abundance across all 78 samples), so area
### scales with abundance. The quadrants read as:
###
###   upper right  0.2 < 3 < 20   abundance rises with particle size
###   lower left   0.2 > 3 > 20   abundance falls with particle size
###   lower right  peak at 3 µm
###   upper left   dip at 3 µm (highest in the free-living and largest fractions)
###
### 各ASVについて画分ごとの平均相対存在量から2つのlog2比を計算し、
### 点の大きさは平均相対存在量の平方根で表す。
###
### Zeros: a pseudocount of half the smallest non-zero fraction mean is added to
### every fraction mean before taking the ratios, so ASVs absent from a fraction
### stay on the plot at a finite, clearly extreme position.
### ゼロ対策として最小非ゼロ値の半分を擬似カウントとして加える。
###
### Page 3 keeps the same processing-derived axes but sizes each point by the
### ASV's mean relative abundance in the *survey* samples. Survey samples are
### whole water filtered directly onto 0.2 µm Sterivex cartridges — nothing is
### pre-filtered out — so they should contain the particle-attached cells too,
### and the page shows how much of each fraction behaviour whole-water sampling
### actually recovers. It is drawn only when the survey objects are in session.
### サーベイ試料は前処理なしの全海水を0.2 µmで濾過したもので、粒子付着菌も
### 含むはずである。3ページ目はそれがどの程度回収されているかを示す。
### 3ページ目は軸はそのままに、点の大きさをサーベイ試料での平均相対存在量にする。
###
### REQUIRES (from 00_setup.R, P01):
###   asgard_filtered_p_hm2, meta_asgard_p2, shorternames, fullnameboot
### OPTIONAL (from S01):
###   asgard_filtered  - survey relative abundances, for page 3
###
### OUTPUT:
###   output_p/fraction_order/ASGARD_fraction_logratio_volcano.pdf
###     page 1  volcano, coloured by fraction rank order
###     page 2  volcano, coloured by class
###     page 3  volcano, sized by whole-water (survey) abundance
###     page 4  recovery of these ASVs in the whole-water samples
###   output_p/fraction_order/fraction_logratio_per_ASV.csv

library(tidyverse)
library(here)

dir.create(here("output_p", "fraction_order"), showWarnings = FALSE, recursive = TRUE)

FR <- c("0.2 µm", "3 µm", "20 µm")

# ==============================================================================
# Section 1: fraction means, ratios and sizes
# ==============================================================================

frac_mean <- sapply(FR, function(fr) {
  smp <- rownames(meta_asgard_p2)[meta_asgard_p2$filter == fr]
  colMeans(asgard_filtered_p_hm2[smp, , drop = FALSE])
})

pseudo <- min(frac_mean[frac_mean > 0]) / 2
message(sprintf("P17: pseudocount = %.3g (half the smallest non-zero fraction mean)",
                pseudo))

adj <- frac_mean + pseudo

vol <- data.frame(
  ASV        = rownames(frac_mean),
  mean_0.2um = frac_mean[, "0.2 µm"],
  mean_3um   = frac_mean[, "3 µm"],
  mean_20um  = frac_mean[, "20 µm"],
  x_log2_3_over_0.2  = log2(adj[, "3 µm"]  / adj[, "0.2 µm"]),
  y_log2_20_over_3   = log2(adj[, "20 µm"] / adj[, "3 µm"]),
  RA_pct     = 100 * colMeans(asgard_filtered_p_hm2)[rownames(frac_mean)],
  row.names  = NULL
)

# ASVs missing from a fraction sit at a pseudocount-driven position, so mark
# them: the streak along y = 0 at negative x, for instance, is ASVs absent from
# both the 3 and 20 µm fractions.
# 画分に検出されないASVは擬似カウント由来の座標になるため、記号で区別する。
vol$detected_all <- factor(
  ifelse(rowSums(frac_mean == 0) == 0, "all three fractions", "zero in >=1 fraction"),
  levels = c("all three fractions", "zero in >=1 fraction"))

# taxonomy, as in P12/P13
asv_idx    <- match(vol$ASV, shorternames)
strip_boot <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
vol$class  <- strip_boot(fullnameboot$Class[asv_idx])
vol$genus  <- strip_boot(fullnameboot$Genus[asv_idx])

target_classes <- c("Bacteroidia", "Alphaproteobacteria", "Gammaproteobacteria")
vol$class_grp  <- factor(ifelse(vol$class %in% target_classes, vol$class, "Other"),
                         levels = c(target_classes, "Other"))
class_cols <- c("Bacteroidia" = "#4DAF4A", "Alphaproteobacteria" = "#E41A1C",
                "Gammaproteobacteria" = "#FF7F00", "Other" = "#9E9E9E")

# rank-order group, the same top-three-plus-Other scheme as S26
rank_pattern <- function(x, nms = FR) {
  r <- rank(-x, ties.method = "min")
  paste(vapply(sort(unique(r)),
               function(k) paste(sort(nms[r == k]), collapse = " = "),
               character(1)), collapse = " > ")
}
pattern_all   <- apply(frac_mean, 1, rank_pattern)
top3          <- names(sort(table(pattern_all), decreasing = TRUE))[1:3]
shorten_order <- function(z) paste0(gsub(" µm", "", z), " µm")
vol$rank_grp  <- factor(ifelse(pattern_all %in% top3,
                               shorten_order(pattern_all), "Other"),
                        levels = c(shorten_order(top3), "Other"))
rank_cols <- setNames(
  c(unname(c("0.2" = "#E41A1C", "3" = "#4DAF4A", "20" = "#377EB8")[
      sub(" .*$", "", top3)]), "#808080"),
  levels(vol$rank_grp))

# ==============================================================================
# Section 2: quadrant counts
# ==============================================================================

vol$quadrant <- with(vol, ifelse(
  x_log2_3_over_0.2 >= 0 & y_log2_20_over_3 >= 0, "rises with size (0.2 < 3 < 20)",
  ifelse(x_log2_3_over_0.2 >= 0 & y_log2_20_over_3 < 0, "peak at 3 µm",
  ifelse(x_log2_3_over_0.2 < 0 & y_log2_20_over_3 >= 0, "dip at 3 µm",
         "falls with size (0.2 > 3 > 20)"))))

quad_tab <- vol %>%
  dplyr::group_by(quadrant) %>%
  dplyr::summarise(n_ASVs = dplyr::n(),
                   summed_RA_pct = round(sum(RA_pct), 2),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_ASVs))

message("\n--- ASVs per quadrant ---")
for (i in seq_len(nrow(quad_tab)))
  message(sprintf("  %-32s n = %3d   summed RA = %5.1f %%",
                  quad_tab$quadrant[i], quad_tab$n_ASVs[i],
                  quad_tab$summed_RA_pct[i]))

# ==============================================================================
# Section 2b: how well whole-water (survey) sampling recovers these ASVs
#
# Survey samples are whole water on 0.2 µm Sterivex, so every processing ASV
# could in principle appear there. Look each one up in the RAW survey table —
# before the survey abundance filter — to separate "never seen" from "seen but
# filtered out".
# サーベイ（全海水）でどれだけ回収できているかを、フィルター前の生データで確認する。
# ==============================================================================

have_survey <- exists("asgard_filtered") && exists("seqtab_16Sprop")

if (have_survey) {
  raw_surv  <- seqtab_16Sprop[rownames(asgard_filtered), ]
  raw_cols  <- colnames(seqtab_filt)[match(vol$ASV, shorternames)]
  in_raw    <- raw_cols %in% colnames(raw_surv)
  raw_m     <- raw_surv[, raw_cols[in_raw], drop = FALSE]

  vol$survey_n_detected <- NA_integer_
  vol$survey_mean_RA_pct <- NA_real_
  vol$survey_max_RA_pct  <- NA_real_
  vol$survey_n_detected[in_raw]  <- colSums(raw_m > 0)
  vol$survey_mean_RA_pct[in_raw] <- 100 * colMeans(raw_m)
  vol$survey_max_RA_pct[in_raw]  <- 100 * apply(raw_m, 2, max)
  vol$in_filtered_survey <- vol$ASV %in% colnames(asgard_filtered)

  vol$survey_recovery <- factor(
    ifelse(vol$in_filtered_survey, "passed the survey filter",
           ifelse(!is.na(vol$survey_n_detected) & vol$survey_n_detected > 0,
                  "detected, below filter", "never detected")),
    levels = c("passed the survey filter", "detected, below filter",
               "never detected"))

  survey_depth  <- rowSums(asgard_seqcount)
  det_limit_pct <- 100 / median(survey_depth)

  message(sprintf(paste0("P17: of %d processing ASVs, %d are detected in the raw ",
                         "survey table, %d pass the survey abundance filter, ",
                         "%d are never detected in %d whole-water samples"),
                  nrow(vol), sum(vol$survey_n_detected > 0, na.rm = TRUE),
                  sum(vol$in_filtered_survey),
                  sum(vol$survey_n_detected == 0, na.rm = TRUE),
                  nrow(asgard_filtered)))
  message(sprintf("     median survey depth %.0f reads -> ~%.3g %% RA gives one expected read",
                  median(survey_depth), det_limit_pct))
}

write.csv(vol[order(-vol$RA_pct), ],
          here("output_p", "fraction_order", "fraction_logratio_per_ASV.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 3: the plot
# ==============================================================================

lim_x <- max(abs(vol$x_log2_3_over_0.2)) * 1.05
lim_y <- max(abs(vol$y_log2_20_over_3))  * 1.05

quad_labels <- data.frame(
  x = c( lim_x,  lim_x, -lim_x, -lim_x) * 0.95,
  y = c( lim_y, -lim_y,  lim_y, -lim_y) * 0.95,
  hjust = c(1, 1, 0, 0),
  vjust = c(1, 0, 1, 0),
  label = c("0.2 < 3 < 20\nrises with size", "peak at 3 µm",
            "dip at 3 µm", "0.2 > 3 > 20\nfalls with size")
)

volcano_plot <- function(colour_by, palette, legend_title) {
  ggplot(vol, aes(x = x_log2_3_over_0.2, y = y_log2_20_over_3)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_point(aes(size = RA_pct, fill = .data[[colour_by]],
                   shape = detected_all),
               colour = "grey20", stroke = 0.25, alpha = 0.8) +
    scale_shape_manual(values = c("all three fractions" = 21,
                                  "zero in >=1 fraction" = 24),
                       name = "Detected in") +
    geom_text(data = quad_labels,
              aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 4.5, fontface = "bold",
              colour = "grey35", lineheight = 0.95) +
    scale_fill_manual(values = palette, name = legend_title) +
    # area scales with abundance: radius ~ sqrt(RA)
    scale_size(range = c(1, 14), trans = "sqrt",
               name = "Mean RA (%)\n(area ~ RA)") +
    coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
    labs(x = expression(bold(log[2]~"( 3 µm / 0.2 µm )")),
         y = expression(bold(log[2]~"( 20 µm / 3 µm )")),
         subtitle = sprintf(paste0("221 ASVs; fraction means with a pseudocount of ",
                                   "%.2g added before the ratios"), pseudo)) +
    guides(fill  = guide_legend(override.aes = list(size = 5, shape = 21)),
           shape = guide_legend(override.aes = list(size = 5, fill = "grey70"))) +
    theme_bw(base_size = 15) +
    theme(axis.title    = element_text(size = 17),
          axis.text     = element_text(size = 13),
          plot.subtitle = element_text(size = 12),
          legend.title  = element_text(size = 13, face = "bold"))
}

pdf(here("output_p", "fraction_order", "ASGARD_fraction_logratio_volcano.pdf"),
    width = 12, height = 10)
print(volcano_plot("rank_grp",  rank_cols,  "Fraction rank order"))
print(volcano_plot("class_grp", class_cols, "Class"))

# ------------------------------------------------------------------------------
# Page 3: same axes, sized by survey abundance instead of processing abundance
# 軸は同じで、点の大きさをサーベイ試料での平均相対存在量にした版
# ------------------------------------------------------------------------------

if (have_survey) {
  survey_mean <- 100 * colMeans(asgard_filtered)
  vol$RA_survey_pct <- unname(survey_mean[vol$ASV])
  vol$in_survey <- !is.na(vol$RA_survey_pct)
  vol$RA_survey_pct[!vol$in_survey] <- 0

  n_shared <- sum(vol$in_survey)
  message(sprintf("P17: %d of %d ASVs are also in the survey table (%d are not)",
                  n_shared, nrow(vol), nrow(vol) - n_shared))

  print(
    ggplot(vol[vol$in_survey, ],
           aes(x = x_log2_3_over_0.2, y = y_log2_20_over_3)) +
      geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
      geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
      # ASVs absent from the survey table, marked but not sized
      geom_point(data = vol[!vol$in_survey, ],
                 aes(x = x_log2_3_over_0.2, y = y_log2_20_over_3),
                 inherit.aes = FALSE, shape = 4, size = 2,
                 colour = "grey55", stroke = 0.6) +
      geom_point(aes(size = RA_survey_pct, fill = rank_grp, shape = detected_all),
                 colour = "grey20", stroke = 0.25, alpha = 0.8) +
      scale_shape_manual(values = c("all three fractions" = 21,
                                    "zero in >=1 fraction" = 24),
                         name = "Detected in
(processing)") +
      geom_text(data = quad_labels,
                aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
                inherit.aes = FALSE, size = 4.5, fontface = "bold",
                colour = "grey35", lineheight = 0.95) +
      scale_fill_manual(values = rank_cols, name = "Fraction rank order") +
      scale_size(range = c(1, 14), trans = "sqrt",
                 name = "Mean RA (%)
in survey samples") +
      coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
      labs(x = expression(bold(log[2]~"( 3 µm / 0.2 µm )")),
           y = expression(bold(log[2]~"( 20 µm / 3 µm )")),
           subtitle = sprintf(paste0("Axes from the processing fractions; point ",
                                     "area is the mean relative abundance across ",
                                     "the 181 whole-water survey samples.\n%d of ",
                                     "%d ASVs occur in the survey table; the %d ",
                                     "that do not are drawn as crosses"),
                              n_shared, nrow(vol), nrow(vol) - n_shared)) +
      guides(fill  = guide_legend(override.aes = list(size = 5, shape = 21)),
             shape = guide_legend(override.aes = list(size = 5, fill = "grey70"))) +
      theme_bw(base_size = 15) +
      theme(axis.title    = element_text(size = 17),
            axis.text     = element_text(size = 13),
            plot.subtitle = element_text(size = 12),
            legend.title  = element_text(size = 13, face = "bold"))
  )
  # ----------------------------------------------------------------------------
  # Page 4: how much of each quadrant whole-water sampling recovers
  # 全海水サンプリングで各象限がどれだけ回収されているか
  # ----------------------------------------------------------------------------

  quad_short <- c("falls with size (0.2 > 3 > 20)" = "falls with size\n0.2 > 3 > 20",
                  "peak at 3 µm"                   = "peak at\n3 µm",
                  "rises with size (0.2 < 3 < 20)" = "rises with size\n0.2 < 3 < 20",
                  "dip at 3 µm"                    = "dip at\n3 µm")
  quad_order <- names(quad_short)

  rec_df <- vol %>%
    dplyr::mutate(quad = factor(quad_short[quadrant],
                                levels = unname(quad_short[quad_order]))) %>%
    dplyr::count(quad, survey_recovery, name = "n_ASVs")

  p_rec <- ggplot(rec_df, aes(x = quad, y = n_ASVs, fill = survey_recovery)) +
    geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = ifelse(n_ASVs > 0, n_ASVs, "")),
              position = position_stack(vjust = 0.5), size = 4.2,
              fontface = "bold", colour = "grey15") +
    scale_fill_manual(values = c("passed the survey filter" = "#4DAF4A",
                                 "detected, below filter"   = "#FDBF6F",
                                 "never detected"           = "#E41A1C"),
                      name = NULL) +
    labs(x = NULL, y = "Number of ASVs",
         title = "(A) Recovery of processing ASVs in whole-water samples") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          plot.title  = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 11),
          axis.title  = element_text(size = 14, face = "bold"))

  ra_cmp <- vol %>%
    dplyr::mutate(quad = factor(quad_short[quadrant],
                                levels = unname(quad_short[quad_order]))) %>%
    dplyr::group_by(quad) %>%
    dplyr::summarise(`processing (within fractions)` = sum(RA_pct),
                     `whole water (survey)` = sum(survey_mean_RA_pct, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_longer(-quad, names_to = "dataset", values_to = "RA_pct") %>%
    dplyr::mutate(dataset = factor(dataset,
                                   levels = c("processing (within fractions)",
                                              "whole water (survey)")))

  p_ra <- ggplot(ra_cmp, aes(x = quad, y = RA_pct, fill = dataset)) +
    geom_col(width = 0.72, position = position_dodge(width = 0.78),
             colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.1f", RA_pct)),
              position = position_dodge(width = 0.78), vjust = -0.35,
              size = 4.2, fontface = "bold") +
    scale_fill_manual(values = c("processing (within fractions)" = "#6A3D9A",
                                 "whole water (survey)"          = "#1F78B4"),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Summed mean relative abundance (%)",
         title = "(B) Where that abundance sits",
         caption = paste0(
           "Fraction abundances are computed within each filter, so they\n",
           "overstate the whole-water share. At the median survey depth\n",
           "(", format(round(median(survey_depth)), big.mark = ","),
           " reads) an ASV needs ~", signif(det_limit_pct, 2),
           " % RA for one expected read,\n",
           "so non-detection is also consistent with dilution alone.")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          plot.title   = element_text(face = "bold", size = 14),
          plot.caption = element_text(hjust = 0, size = 9.5, colour = "grey30",
                                      lineheight = 1.05),
          axis.text.x  = element_text(size = 11),
          axis.title   = element_text(size = 14, face = "bold"))

  print(gridExtra::grid.arrange(p_rec, p_ra, ncol = 2))

} else {
  message("P17: asgard_filtered not found — skipping the survey pages ",
          "(run S01_data_prep.R first to include them).")
}

dev.off()

message("\nP17_fraction_logratio.R: done.")
message("  PDF: output_p/fraction_order/ASGARD_fraction_logratio_volcano.pdf")
message("  CSV: output_p/fraction_order/fraction_logratio_per_ASV.csv")
