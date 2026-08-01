### P18_accumulation_curve.R
### ASGARD 2017 — Concatenated ASV accumulation curve across sample types
### サンプル種別を連結したASV蓄積曲線
###
### Samples are added to the pool in blocks — first every whole-water (survey)
### sample, then the processing 0.2 µm fraction, then 3 µm, then 20 µm — and the
### cumulative number of distinct ASVs is tracked. The curve is monotonic by
### construction; the step at each block boundary is the number of ASVs that
### sample type contributes beyond everything collected before it.
###
### 全海水 → 0.2 µm → 3 µm → 20 µm の順に試料を加え、累積ASV数を追う。
### 境界での立ち上がりが、その画分だけが捉えているASVの数にあたる。
###
### Sample order within each block is random, so every curve is averaged over
### N_PERM permutations (mean +/- sd shown).
### ブロック内の順序はランダムなので、N_PERM 回の平均を描く。
###
### Page 3 repeats the first curve on data rarefied to RAREFY_DEPTH reads per
### sample (samples below that are dropped), so the yield of the fractions
### cannot be explained by their sequencing depth.
### ページ3は各試料を RAREFY_DEPTH リードにレアファイした版（不足試料は除外）。
###
### REQUIRES:
###   seqtab_16Smat    - raw 16S counts (00_setup.R)
###   asgard_filtered  - survey samples (S01)
###   meta_asgard_p2   - processing metadata with `filter` (P01)
###
### OUTPUT:
###   output_p/accumulation/ASGARD_accumulation_concatenated.pdf
###     page 1  whole water first, then the fractions
###     page 2  fractions first, then whole water
###     page 3  whole water first, rarefied to RAREFY_DEPTH reads
###   output_p/accumulation/accumulation_block_summary.csv
###   output_p/accumulation/accumulation_curve.csv

library(tidyverse)
library(vegan)
library(here)

dir.create(here("output_p", "accumulation"), showWarnings = FALSE, recursive = TRUE)

set.seed(42)
N_PERM       <- 200
RAREFY_DEPTH <- 5000

# ==============================================================================
# Section 1: sample blocks and the presence/absence matrix
# ==============================================================================

survey_samples <- rownames(asgard_filtered)
frac_samples   <- lapply(c("0.2 µm", "3 µm", "20 µm"), function(fr)
  rownames(meta_asgard_p2)[meta_asgard_p2$filter == fr])
names(frac_samples) <- c("0.2 µm", "3 µm", "20 µm")

blocks <- c(list("whole water" = survey_samples), frac_samples)
blocks <- lapply(blocks, function(x) intersect(x, rownames(seqtab_16Smat)))

message("P18: samples per block — ",
        paste(sprintf("%s: %d", names(blocks), lengths(blocks)), collapse = ", "),
        " (total ", sum(lengths(blocks)), ")")

pa <- seqtab_16Smat[unlist(blocks, use.names = FALSE), ] > 0
pa <- pa[, colSums(pa) > 0, drop = FALSE]
message("     ASVs detected in at least one of these samples: ", ncol(pa))

# rarefied sample set: drop anything below RAREFY_DEPTH reads
sample_depth <- rowSums(seqtab_16Smat)
blocks_rar   <- lapply(blocks, function(s) s[sample_depth[s] >= RAREFY_DEPTH])

message(sprintf("     rarefying to %d reads keeps %d of %d samples (%s)",
                RAREFY_DEPTH, sum(lengths(blocks_rar)), sum(lengths(blocks)),
                paste(sprintf("%s: %d/%d", names(blocks), lengths(blocks_rar),
                              lengths(blocks)), collapse = ", ")))

# ==============================================================================
# Section 2: accumulation averaged over permutations within each block
# ==============================================================================

# rarefy = TRUE re-rarefies the counts inside every permutation, so the curve
# carries both the ordering and the subsampling uncertainty
# rarefy=TRUE では並び順と抽出の両方の不確実性を含める
accumulate <- function(block_list, n_perm = N_PERM, rarefy = FALSE) {
  samples <- unlist(block_list, use.names = FALSE)
  n_tot   <- length(samples)
  counts  <- if (rarefy) seqtab_16Smat[samples, , drop = FALSE] else NULL
  curves  <- matrix(NA_real_, nrow = n_perm, ncol = n_tot)

  for (p in seq_len(n_perm)) {
    m_pa <- if (rarefy) {
      suppressWarnings(vegan::rrarefy(counts, RAREFY_DEPTH)) > 0
    } else pa

    ord   <- unlist(lapply(block_list, function(s) sample(s)), use.names = FALSE)
    m     <- m_pa[ord, , drop = FALSE]
    # row of first detection for every ASV; ASVs never seen get NA
    first <- max.col(t(m), ties.method = "first")
    first[colSums(m) == 0] <- NA_integer_
    curves[p, ] <- cumsum(tabulate(first[!is.na(first)], nbins = n_tot))
  }

  list(
    summary = data.frame(
      n_samples = seq_len(n_tot),
      mean      = colMeans(curves),
      sd        = apply(curves, 2, sd),
      block     = factor(rep(names(block_list), lengths(block_list)),
                         levels = names(block_list)),
      row.names = NULL),
    curves = curves)
}

acc_fwd   <- accumulate(blocks)
acc_rev   <- accumulate(rev(blocks))
acc_rar   <- accumulate(blocks_rar, rarefy = TRUE)
curve_fwd <- acc_fwd$summary
curve_rev <- acc_rev$summary
curve_rar <- acc_rar$summary

block_summary <- function(cv, label) {
  cv %>%
    dplyr::group_by(block) %>%
    dplyr::summarise(n_samples = dplyr::n(),
                     cumulative_ASVs = round(dplyr::last(mean), 1),
                     .groups = "drop") %>%
    dplyr::mutate(order = label,
                  new_ASVs = round(cumulative_ASVs - dplyr::lag(cumulative_ASVs, default = 0), 1),
                  pct_of_total = round(100 * new_ASVs / max(cumulative_ASVs), 1))
}

RAR_LABEL <- paste0("whole water first, rarefied to ", RAREFY_DEPTH)

summ <- rbind(block_summary(curve_fwd, "whole water first"),
              block_summary(curve_rev, "size fractions first"),
              block_summary(curve_rar, RAR_LABEL))

write.csv(summ, here("output_p", "accumulation", "accumulation_block_summary.csv"),
          row.names = FALSE)
write.csv(rbind(dplyr::mutate(curve_fwd, order = "whole water first"),
                dplyr::mutate(curve_rev, order = "size fractions first"),
                dplyr::mutate(curve_rar, order = RAR_LABEL)),
          here("output_p", "accumulation", "accumulation_curve.csv"),
          row.names = FALSE)

report_blocks <- function(label) {
  message("\n--- ASVs added by each block (", label, ") ---")
  d <- summ[summ$order == label, ]
  for (i in seq_len(nrow(d)))
    message(sprintf("  %-12s %3d samples   +%6.1f ASVs   cumulative %6.1f",
                    d$block[i], d$n_samples[i], d$new_ASVs[i], d$cumulative_ASVs[i]))
}
report_blocks("whole water first")
report_blocks("size fractions first")
report_blocks(RAR_LABEL)

# ==============================================================================
# Section 3: figures
# ==============================================================================

# fractions keep the ternary/filter colours; whole water is near-black so it is
# not confused with the 20 µm blue
# 画分は三角図と同じ配色、全海水は区別のため黒に近い色にする
block_cols <- c("whole water" = "#252525", "0.2 µm" = "#E41A1C",
                "3 µm" = "#4DAF4A", "20 µm" = "#377EB8")

accum_plot <- function(cv, title) {
  bnd  <- cumsum(as.integer(table(cv$block)[levels(cv$block)]))
  labs <- data.frame(
    x     = bnd - as.integer(table(cv$block)[levels(cv$block)]) / 2,
    y     = max(cv$mean) * 0.06,
    block = factor(levels(cv$block), levels = levels(cv$block)),
    n     = paste0("+", round(diff(c(0, cv$mean[bnd]))), " ASVs")
  )

  ggplot(cv, aes(x = n_samples, y = mean)) +
    geom_vline(xintercept = head(bnd, -1), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = block),
                alpha = 0.25) +
    geom_line(aes(colour = block, group = 1), linewidth = 1.1) +
    geom_text(data = labs, aes(x = x, y = y, label = n, colour = block),
              inherit.aes = FALSE, fontface = "bold", size = 4.5) +
    scale_colour_manual(values = block_cols, name = "Sample type") +
    scale_fill_manual(values = block_cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Samples added (in block order)",
         y = "Cumulative number of ASVs",
         title = title) +
    theme_bw(base_size = 15) +
    theme(plot.title   = element_text(face = "bold", size = 14),
          axis.title   = element_text(size = 16, face = "bold"),
          axis.text    = element_text(size = 13),
          legend.position = "top")
}

pdf(here("output_p", "accumulation", "ASGARD_accumulation_concatenated.pdf"),
    width = 12, height = 8)

print(accum_plot(curve_fwd,
                 paste0("Whole water first, then the size fractions (mean of ",
                        N_PERM, " within-block orderings, ribbon = 1 sd)")))
print(accum_plot(curve_rev,
                 paste0("Size fractions first, then whole water (mean of ",
                        N_PERM, " within-block orderings, ribbon = 1 sd)")))

print(accum_plot(
  curve_rar,
  sprintf(paste0("Whole water first, rarefied to %d reads per sample (mean of %d ",
                 "within-block orderings, ribbon = 1 sd)\n%d of %d samples reach ",
                 "that depth and are kept: %s"),
          RAREFY_DEPTH, N_PERM, sum(lengths(blocks_rar)), sum(lengths(blocks)),
          paste(sprintf("%s %d/%d", names(blocks), lengths(blocks_rar),
                        lengths(blocks)), collapse = ", "))))

dev.off()

message("\nP18_accumulation_curve.R: done.")
message("  PDF: output_p/accumulation/ASGARD_accumulation_concatenated.pdf")
message("  CSV: output_p/accumulation/accumulation_block_summary.csv")
message("  CSV: output_p/accumulation/accumulation_curve.csv")
