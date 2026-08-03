### P16_fraction_rank_order.R
### ASGARD 2017 Processing — Rank order of the three filter fractions per ASV
### ASVごとに、3つのフィルター画分を存在量の高い順に並べた順序を数える
###
### For each of the 221 ASVs the relative abundance is summarised within each
### filter fraction, the three fractions are ranked from most to least abundant,
### and the ASVs are counted per rank pattern — the six strict permutations
### (0.2 > 3 > 20, 0.2 > 20 > 3, 3 > 0.2 > 20, 3 > 20 > 0.2, 20 > 0.2 > 3,
### 20 > 3 > 0.2) plus the tie patterns.
###
### Two summaries are reported:
###   median - as requested. Most ASVs are absent from most samples, so the
###            median is 0 in at least one fraction for 202 of 221 ASVs and the
###            ranking is dominated by ties.
###   mean   - the same calculation on fraction means, where nearly every ASV
###            gets a strict ordering. This is the quantity the ternary diagram
###            and S26's fraction grouping use.
###
### 中央値は0が多く同順位になりやすいため、平均値による集計も併記する。
###
### REQUIRES (from P01):
###   asgard_filtered_p_hm2 - 78 x 221 relative abundance matrix
###   meta_asgard_p2        - metadata with the `filter` column
###   shorternames          - taxonomy labels (00_setup.R)
###
### OUTPUT:
###   output_p/fraction_order/ASGARD_fraction_rank_order.pdf
###   output_p/fraction_order/fraction_rank_counts.csv
###   output_p/fraction_order/fraction_rank_per_ASV.csv

library(tidyverse)
library(here)

dir.create(here("output_p", "fraction_order"), showWarnings = FALSE, recursive = TRUE)

FR <- c("0.2 µm", "3 µm", "20 µm")

frac_samples <- lapply(FR, function(fr)
  rownames(meta_asgard_p2)[meta_asgard_p2$filter == fr])
names(frac_samples) <- FR

# ==============================================================================
# Section 1: per-ASV summaries within each fraction
# 各ASVについて画分ごとの中央値・平均値・検出数を計算
# ==============================================================================

summarise_fraction <- function(FUN) {
  out <- sapply(FR, function(fr)
    apply(asgard_filtered_p_hm2[frac_samples[[fr]], , drop = FALSE], 2, FUN))
  colnames(out) <- FR
  out
}

med_ra <- summarise_fraction(median)
mean_ra <- summarise_fraction(mean)
n_det  <- sapply(FR, function(fr)
  colSums(asgard_filtered_p_hm2[frac_samples[[fr]], , drop = FALSE] > 0))

# ==============================================================================
# Section 2: rank pattern per ASV
# 高い順に並べた文字列を作る（同値は "=" でつなぐ）
# ==============================================================================

rank_pattern <- function(x, nms = FR) {
  o <- order(x, decreasing = TRUE)
  s <- nms[o]
  v <- x[o]
  paste0(s[1],
         if (v[2] == v[1]) " = " else " > ", s[2],
         if (v[3] == v[2]) " = " else " > ", s[3])
}

# a tie anywhere makes the pattern non-strict; order ties alphabetically so the
# same tie always produces the same string
canonical <- function(x, nms = FR) {
  r <- rank(-x, ties.method = "min")
  parts <- vapply(sort(unique(r)), function(k) paste(sort(nms[r == k]), collapse = " = "),
                  character(1))
  paste(parts, collapse = " > ")
}

pattern_med  <- apply(med_ra,  1, canonical)
pattern_mean <- apply(mean_ra, 1, canonical)

PERMUTATIONS <- c("0.2 µm > 3 µm > 20 µm", "0.2 µm > 20 µm > 3 µm",
                  "3 µm > 0.2 µm > 20 µm", "3 µm > 20 µm > 0.2 µm",
                  "20 µm > 0.2 µm > 3 µm", "20 µm > 3 µm > 0.2 µm")

# ==============================================================================
# Section 3: counts per pattern
# ==============================================================================

count_patterns <- function(pattern, statistic) {
  tab <- as.data.frame(table(pattern), stringsAsFactors = FALSE)
  colnames(tab) <- c("pattern", "n_ASVs")
  tab$statistic <- statistic
  tab$strict    <- tab$pattern %in% PERMUTATIONS
  tab[order(-tab$strict, match(tab$pattern, PERMUTATIONS), -tab$n_ASVs),
      c("statistic", "pattern", "strict", "n_ASVs")]
}

counts <- rbind(count_patterns(pattern_med,  "median"),
                count_patterns(pattern_mean, "mean"))

# every permutation listed even when no ASV matches it
counts <- counts %>%
  dplyr::right_join(
    expand.grid(statistic = c("median", "mean"), pattern = PERMUTATIONS,
                stringsAsFactors = FALSE) %>%
      dplyr::mutate(strict = TRUE) %>%
      dplyr::full_join(counts, by = c("statistic", "pattern", "strict")),
    by = c("statistic", "pattern", "strict", "n_ASVs")) %>%
  dplyr::mutate(n_ASVs = ifelse(is.na(n_ASVs), 0L, as.integer(n_ASVs))) %>%
  dplyr::distinct()

write.csv(counts[order(counts$statistic, -counts$strict, -counts$n_ASVs), ],
          here("output_p", "fraction_order", "fraction_rank_counts.csv"),
          row.names = FALSE)

message("\n--- ASVs per fraction rank pattern (median RA, as requested) ---")
med_tab <- counts[counts$statistic == "median", ]
for (i in seq_len(nrow(med_tab)))
  message(sprintf("  %-26s %s%d", med_tab$pattern[i],
                  if (med_tab$strict[i]) "" else "(tie) ", med_tab$n_ASVs[i]))
message("  strict permutations total: ",
        sum(med_tab$n_ASVs[med_tab$strict]), " / ", ncol(asgard_filtered_p_hm2))

message("\n--- ASVs per fraction rank pattern (mean RA) ---")
mean_tab <- counts[counts$statistic == "mean", ]
for (i in seq_len(nrow(mean_tab)))
  message(sprintf("  %-26s %s%d", mean_tab$pattern[i],
                  if (mean_tab$strict[i]) "" else "(tie) ", mean_tab$n_ASVs[i]))
message("  strict permutations total: ",
        sum(mean_tab$n_ASVs[mean_tab$strict]), " / ", ncol(asgard_filtered_p_hm2))

# ==============================================================================
# Section 4: per-ASV table
# ==============================================================================

asv_tab <- data.frame(
  ASV            = colnames(asgard_filtered_p_hm2),
  median_0.2um   = round(100 * med_ra[, "0.2 µm"], 4),
  median_3um     = round(100 * med_ra[, "3 µm"],   4),
  median_20um    = round(100 * med_ra[, "20 µm"],  4),
  mean_0.2um     = round(100 * mean_ra[, "0.2 µm"], 4),
  mean_3um       = round(100 * mean_ra[, "3 µm"],   4),
  mean_20um      = round(100 * mean_ra[, "20 µm"],  4),
  n_detected_0.2um = n_det[, "0.2 µm"],
  n_detected_3um   = n_det[, "3 µm"],
  n_detected_20um  = n_det[, "20 µm"],
  pattern_median = pattern_med,
  pattern_mean   = pattern_mean,
  row.names = NULL
)

write.csv(asv_tab,
          here("output_p", "fraction_order", "fraction_rank_per_ASV.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 5: figure
# ==============================================================================

plot_df <- counts
plot_df$pattern <- factor(
  plot_df$pattern,
  levels = c(PERMUTATIONS,
             setdiff(unique(plot_df$pattern[!plot_df$strict]), PERMUTATIONS)))
plot_df$statistic <- factor(plot_df$statistic, levels = c("median", "mean"),
                            labels = c("median relative abundance (as requested)",
                                       "mean relative abundance"))
plot_df$kind <- ifelse(plot_df$strict, "strict order", "tie")

pdf(here("output_p", "fraction_order", "ASGARD_fraction_rank_order.pdf"),
    width = 12, height = 8)

print(
  ggplot(plot_df, aes(x = pattern, y = n_ASVs, fill = kind)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = n_ASVs), vjust = -0.35, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("strict order" = "#377EB8", "tie" = "#BDBDBD"),
                      name = NULL) +
    facet_wrap(~ statistic, ncol = 1, scales = "free_y") +
    labs(x = "Filter fractions ordered from most to least abundant",
         y = "Number of ASVs",
         subtitle = paste0("221 ASVs; ties arise when an ASV has the same summary ",
                           "value in two or more fractions (usually zero)")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 14) +
    theme(plot.margin   = margin(10, 10, 10, 26),
          axis.text.x   = element_text(angle = 30, hjust = 1, size = 12),
          axis.title    = element_text(size = 16, face = "bold"),
          strip.text    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 12),
          legend.position = "top")
)

# ==============================================================================
# Section 6: relative abundance distribution within each mean-rank category
# 平均値ランクのカテゴリごとに、ASVの相対存在量のヒストグラムを描く
#
# Categories come from the mean-based ranking, where 193 of 221 ASVs are
# strictly ordered. Each ASV contributes one value per panel: page 1 its overall
# mean relative abundance across all 78 samples, page 2 its mean within each
# fraction. Abundances span several orders of magnitude, so the x axis is log10.
# ==============================================================================

overall_mean <- colMeans(asgard_filtered_p_hm2)   # per-ASV mean RA (proportion)

cat_levels <- c(PERMUTATIONS,
                setdiff(sort(unique(pattern_mean)), PERMUTATIONS))
cat_n      <- table(factor(pattern_mean, levels = cat_levels))
cat_labels <- sprintf("%s\n(n = %d)", cat_levels, cat_n)

hist_df <- data.frame(
  ASV      = names(overall_mean),
  RA_pct   = 100 * overall_mean,
  category = factor(pattern_mean[names(overall_mean)],
                    levels = cat_levels, labels = cat_labels),
  strict   = pattern_mean[names(overall_mean)] %in% PERMUTATIONS,
  row.names = NULL
)

# same values split by fraction (page 2)
hist_frac <- do.call(rbind, lapply(FR, function(fr) data.frame(
  ASV      = rownames(mean_ra),
  fraction = factor(fr, levels = FR),
  RA_pct   = 100 * mean_ra[, fr],
  category = factor(pattern_mean[rownames(mean_ra)],
                    levels = cat_levels, labels = cat_labels),
  row.names = NULL)))
hist_frac <- hist_frac[hist_frac$RA_pct > 0, ]

frac_cols <- c("0.2 µm" = "#E41A1C", "3 µm" = "#4DAF4A", "20 µm" = "#377EB8")

pdf(here("output_p", "fraction_order", "ASGARD_fraction_rank_abundance_hist.pdf"),
    width = 13, height = 8)

# Page 1: overall mean RA per ASV, one panel per rank category
print(
  ggplot(hist_df, aes(x = RA_pct, fill = strict)) +
    geom_histogram(bins = 22, colour = "grey25", linewidth = 0.2) +
    scale_x_log10(labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#BDBDBD"),
                      guide = "none") +
    facet_wrap(~ category, ncol = 5) +
    labs(x = "Mean relative abundance per ASV across all 78 samples (%, log scale)",
         y = "Number of ASVs",
         subtitle = paste0("ASVs grouped by the order of their fraction means; ",
                           "grey panels are the tie categories")) +
    theme_bw(base_size = 13) +
    theme(strip.text    = element_text(face = "bold", size = 11),
          axis.title    = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 12))
)

# Page 2: the same ASVs, mean RA within each fraction
print(
  ggplot(hist_frac, aes(x = RA_pct, fill = fraction)) +
    geom_histogram(bins = 22, colour = "grey25", linewidth = 0.2,
                   position = "identity", alpha = 0.65) +
    scale_x_log10(labels = function(x) format(x, drop0trailing = TRUE)) +
    scale_fill_manual(values = frac_cols, name = "Filter fraction") +
    facet_wrap(~ category, ncol = 5) +
    labs(x = "Mean relative abundance within each filter fraction (%, log scale)",
         y = "Number of ASVs",
         subtitle = "Zero values omitted so they can be shown on a log axis") +
    theme_bw(base_size = 13) +
    theme(strip.text    = element_text(face = "bold", size = 11),
          axis.title    = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "top")
)

# Page 3: same as page 1 on the fourth-root scale the pipeline uses for
# distances and heatmaps. Zeros are representable here, unlike on a log axis.
# ページ3：パイプラインと同じ4乗根変換（0も表示できる）
print(
  ggplot(hist_df, aes(x = (RA_pct / 100)^0.25, fill = strict)) +
    geom_histogram(bins = 22, colour = "grey25", linewidth = 0.2) +
    scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#BDBDBD"),
                      guide = "none") +
    facet_wrap(~ category, ncol = 5) +
    labs(x = expression(bold("Mean relative abundance per ASV across all 78 samples (proportion"^0.25*")")),
         y = "Number of ASVs",
         subtitle = paste0("Fourth-root transformed, as used for the Bray-Curtis ",
                           "distances and heatmaps")) +
    theme_bw(base_size = 13) +
    theme(strip.text    = element_text(face = "bold", size = 11),
          axis.title    = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 12))
)

# Page 4: per-fraction means, fourth-root transformed (zeros retained)
hist_frac_all <- do.call(rbind, lapply(FR, function(fr) data.frame(
  ASV      = rownames(mean_ra),
  fraction = factor(fr, levels = FR),
  RA_prop  = mean_ra[, fr],
  category = factor(pattern_mean[rownames(mean_ra)],
                    levels = cat_levels, labels = cat_labels),
  row.names = NULL)))

print(
  ggplot(hist_frac_all, aes(x = RA_prop^0.25, fill = fraction)) +
    geom_histogram(bins = 22, colour = "grey25", linewidth = 0.2,
                   position = "identity", alpha = 0.65) +
    scale_fill_manual(values = frac_cols, name = "Filter fraction") +
    facet_wrap(~ category, ncol = 5) +
    labs(x = expression(bold("Mean relative abundance within each filter fraction (proportion"^0.25*")")),
         y = "Number of ASVs",
         subtitle = "Fourth-root transformed; zero values retained") +
    theme_bw(base_size = 13) +
    theme(strip.text    = element_text(face = "bold", size = 11),
          axis.title    = element_text(size = 15, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "top")
)

dev.off()

# per-category abundance summary
ra_by_cat <- hist_df %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(n_ASVs = dplyr::n(),
                   median_RA_pct = round(median(RA_pct), 4),
                   mean_RA_pct   = round(mean(RA_pct), 4),
                   min_RA_pct    = round(min(RA_pct), 4),
                   max_RA_pct    = round(max(RA_pct), 4),
                   total_RA_pct  = round(sum(RA_pct), 2),
                   .groups = "drop")

write.csv(ra_by_cat,
          here("output_p", "fraction_order", "fraction_rank_abundance_summary.csv"),
          row.names = FALSE)

message("\n--- mean RA per rank category (mean-based ranks) ---")
for (i in seq_len(nrow(ra_by_cat)))
  message(sprintf("  %-26s n = %3d   median %6.4f %%   sum %6.2f %%",
                  gsub("\n", " ", ra_by_cat$category[i]), ra_by_cat$n_ASVs[i],
                  ra_by_cat$median_RA_pct[i], ra_by_cat$total_RA_pct[i]))

message("\nP16_fraction_rank_order.R: done.")
message("  PDF: output_p/fraction_order/ASGARD_fraction_rank_order.pdf")
message("  PDF: output_p/fraction_order/ASGARD_fraction_rank_abundance_hist.pdf")
message("  CSV: output_p/fraction_order/fraction_rank_abundance_summary.csv")
message("  CSV: output_p/fraction_order/fraction_rank_counts.csv")
message("  CSV: output_p/fraction_order/fraction_rank_per_ASV.csv")
