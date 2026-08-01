### S26_processing_cluster_asv_map.R
### ASGARD 2017 — Where the processing-dataset ASVs occur at the survey stations
### プロセス側のASVが、サーベイ測点でどれだけ存在するかの地図
###
### ASVs from the processing dataset are put into four groups, the summed
### relative abundance of each group is computed for every survey sample, and
### that sum is shown as point size on a map (depth bin x group) and as a
### boxplot (group x depth bin).
###
### Three groupings are produced:
###
###   0. By fraction rank order (the main figure). The fractions are ranked by
###      each ASV's fraction means (see P16_fraction_rank_order.R); the three
###      most common orders get their own group and everything else — the three
###      rare orders and the ties — is pooled as "Other".
###   1. By filter fraction. Each ASV's relative abundance is
###      averaged within each filter fraction and normalised across the three
###      fractions — the same quantity the ternary diagram plots. An ASV is
###      assigned to a fraction when more than 50 % of that normalised abundance
###      sits there, and to "Mixed" when no fraction reaches 50 %.
###   2. By processing sample cluster (kept for comparison). Each ASV goes to the
###      cluster in which it reaches its highest mean relative abundance.
###
### 0. 画分平均の順位パターン別（上位3パターン＋その他）
### 1. フィルター画分別（>50%で1画分に偏るASV、どれも50%未満なら Mixed）
### 2. プロセスクラスター別（平均相対存在量が最大のクラスター）
###
### REQUIRES:
###   processing: asgard_filtered_p_hm2, meta_asgard_p2, clusnum_p (P01, P03)
###   survey:     asgard_filtered, meta_asgard                     (S01)
###   optional:   mapz_survey                                      (S06)
###
### OUTPUT:
###   output/survey/maps/ASGARD_survey_ASVgroup_RA.pdf
###   output/survey/maps/ASGARD_survey_ASVgroup_RA_boxplot.pdf
###   output/survey/maps/survey_ASVgroup_RA.csv
###   output/survey/maps/survey_ASVgroup_RA_summary.csv

library(tidyverse)
library(here)

dir.create(here("output", "survey", "maps"), showWarnings = FALSE, recursive = TRUE)

FRACTIONS <- c("0.2 µm", "3 µm", "20 µm")
CLUSTERS  <- c("1", "2", "3", "4")

# ==============================================================================
# Section 1: Three ways of grouping the processing ASVs
# ==============================================================================

# --- 1a. by filter fraction --------------------------------------------------
# mean relative abundance per fraction, then normalised per ASV across fractions
frac_mean <- sapply(FRACTIONS, function(fr) {
  smp <- rownames(meta_asgard_p2)[meta_asgard_p2$filter == fr]
  colMeans(asgard_filtered_p_hm2[smp, , drop = FALSE])
})
frac_prop <- prop.table(frac_mean, margin = 1)          # rows sum to 1

frac_max   <- apply(frac_prop, 1, max)
frac_which <- FRACTIONS[apply(frac_prop, 1, which.max)]

assign_frac <- ifelse(frac_max > 0.5, frac_which, "Mixed")
names(assign_frac) <- rownames(frac_prop)

FRAC_LEVELS <- c(FRACTIONS, "Mixed")
cc_frac <- c("0.2 µm" = "#E41A1C", "3 µm" = "#4DAF4A",
             "20 µm"  = "#377EB8", "Mixed" = "#808080")

message("S26: ASVs by filter fraction (>50 % in one fraction): ",
        paste(sprintf("%s=%d", FRAC_LEVELS,
                      table(factor(assign_frac, levels = FRAC_LEVELS))),
              collapse = ", "))

# --- 1a2. by the rank order of the three fraction means -----------------------
# top three orders kept, everything else pooled as "Other"
# 上位3つの順序をそのまま使い、残りは "Other" にまとめる
rank_pattern <- function(x, nms = FRACTIONS) {
  r <- rank(-x, ties.method = "min")
  paste(vapply(sort(unique(r)),
               function(k) paste(sort(nms[r == k]), collapse = " = "),
               character(1)),
        collapse = " > ")
}

pattern_all <- apply(frac_mean, 1, rank_pattern)
top3        <- names(sort(table(pattern_all), decreasing = TRUE))[1:3]

# "0.2 µm > 3 µm > 20 µm" -> "0.2 > 3 > 20 µm" so the facet strips fit
# パネル見出しに収まるようラベルを短縮する
shorten_order <- function(x) paste0(gsub(" µm", "", x), " µm")

assign_rank <- ifelse(pattern_all %in% top3, shorten_order(pattern_all), "Other")
names(assign_rank) <- rownames(frac_mean)

RANK_LEVELS <- c(shorten_order(top3), "Other")

# coloured by whichever fraction ranks first, "Other" grey
cc_rank     <- setNames(
  c(unname(c("0.2" = "#E41A1C", "3" = "#4DAF4A", "20" = "#377EB8")[
      sub(" .*$", "", top3)]), "#808080"),
  RANK_LEVELS)

message("S26: ASVs by fraction rank order (top 3 + Other): ",
        paste(sprintf("%s=%d", RANK_LEVELS,
                      table(factor(assign_rank, levels = RANK_LEVELS))),
              collapse = ", "))

# --- 1b. by processing sample cluster ----------------------------------------
assign_clus  <- character(ncol(asgard_filtered_p_hm2))
names(assign_clus) <- colnames(asgard_filtered_p_hm2)
asv_max_mean <- rep(-Inf, ncol(asgard_filtered_p_hm2))

for (k in CLUSTERS) {
  smp_k      <- names(clusnum_p)[as.character(clusnum_p) == k]
  asv_mean_k <- colMeans(asgard_filtered_p_hm2[smp_k, , drop = FALSE])
  upd        <- asv_mean_k > asv_max_mean
  assign_clus[upd]  <- k
  asv_max_mean[upd] <- asv_mean_k[upd]
}
assign_clus <- unname(clus_labels_p[assign_clus])
names(assign_clus) <- colnames(asgard_filtered_p_hm2)

cc_clus <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), clus_levels_p)

message("S26: ASVs by processing cluster: ",
        paste(sprintf("%s=%d", clus_levels_p,
                      table(factor(assign_clus, levels = clus_levels_p))),
              collapse = ", "))

# ==============================================================================
# Section 2: Summed relative abundance of each group in every survey sample
# サーベイ試料ごとに、各グループのASVの相対存在量を合計する
# ==============================================================================

group_ra <- function(assign, levels_in_order) {
  shared <- intersect(names(assign), colnames(asgard_filtered))
  n_asv  <- sapply(levels_in_order, function(g)
    length(intersect(names(assign)[assign == g], shared)))

  ra <- sapply(levels_in_order, function(g) {
    cols <- intersect(names(assign)[assign == g], shared)
    if (length(cols) == 0) return(rep(0, nrow(asgard_filtered)))
    rowSums(asgard_filtered[, cols, drop = FALSE])
  })
  rownames(ra) <- rownames(asgard_filtered)

  df <- as.data.frame(ra) %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(all_of(levels_in_order),
                        names_to = "group", values_to = "RA") %>%
    dplyr::mutate(RA_pct = 100 * RA)

  df$lon        <- meta_asgard[df$Sample, "lon"]
  df$lat        <- meta_asgard[df$Sample, "lat"]
  df$depth_type <- factor(meta_asgard[df$Sample, "depth_type"],
                          levels = c("surf", "mid", "bottom"))
  df <- df[!is.na(df$lon) & !is.na(df$lat) & !is.na(df$depth_type), ]

  # count on its own line so long rank-order labels fit the facet strips
  # 長いラベルが見出しに収まるよう、ASV数は改行して表示する
  labs_g <- paste0(levels_in_order, "\n(", n_asv[levels_in_order], " ASVs)")
  df$group_lab <- factor(paste0(df$group, "\n(", n_asv[df$group], " ASVs)"),
                         levels = labs_g)
  list(df = df, n_asv = n_asv, labels = labs_g, levels = levels_in_order)
}

ra_rank <- group_ra(assign_rank, RANK_LEVELS)
ra_frac <- group_ra(assign_frac, FRAC_LEVELS)
ra_clus <- group_ra(assign_clus, clus_levels_p)

message("     of which also present in the survey ASV table — by fraction: ",
        paste(sprintf("%s=%d", FRAC_LEVELS, ra_frac$n_asv[FRAC_LEVELS]),
              collapse = ", "))

write.csv(dplyr::bind_rows(
            dplyr::mutate(ra_rank$df, grouping = "fraction_rank_order"),
            dplyr::mutate(ra_frac$df, grouping = "filter_fraction"),
            dplyr::mutate(ra_clus$df, grouping = "processing_cluster"))[
            , c("grouping", "Sample", "group", "RA_pct", "lon", "lat", "depth_type")],
          here("output", "survey", "maps", "survey_ASVgroup_RA.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 3: Map — depth bin (rows) x ASV group (columns)
# ==============================================================================

base_map <- function(df) {
  if (exists("mapz_survey")) {
    list(ggmap::ggmap(mapz_survey))
  } else {
    message("     mapz_survey not found (no Stadia Maps key) — ",
            "drawing a coastline basemap instead.")
    land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    list(ggplot(),
         geom_sf(data = land, fill = "grey92", colour = "grey60", linewidth = 0.2),
         coord_sf(xlim = range(df$lon) + c(-1, 1),
                  ylim = range(df$lat) + c(-0.5, 0.5), expand = FALSE))
  }
}

group_palette <- function(res, palette) {
  setNames(unname(palette[res$levels]), res$labels)
}

ra_map <- function(res, palette, subtitle) {
  Reduce(`+`, base_map(res$df)) +
    geom_point(data = res$df,
               aes(x = lon, y = lat, size = RA_pct, colour = group_lab),
               alpha = 0.65) +
    scale_colour_manual(values = group_palette(res, palette), guide = "none") +
    scale_size_continuous(range = c(0.3, 7),
                          name = "Summed relative\nabundance (%)") +
    facet_grid(depth_type ~ group_lab) +
    labs(x = "Longitude", y = "Latitude", subtitle = subtitle) +
    theme_bw(base_size = 14) +
    theme(strip.text    = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 13),
          axis.title    = element_text(size = 16, face = "bold"),
          axis.text     = element_text(size = 11),
          legend.title  = element_text(size = 14, face = "bold"),
          legend.text   = element_text(size = 12))
}

pdf(here("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA.pdf"),
    width = 15, height = 11)
print(ra_map(ra_rank, cc_rank,
             "ASVs grouped by the rank order of their fraction means (top three orders + Other)"))
print(ra_map(ra_frac, cc_frac,
             "ASVs grouped by filter fraction (>50 % of their abundance in one fraction)"))
print(ra_map(ra_clus, cc_clus,
             "ASVs grouped by the processing cluster where they are most abundant"))
dev.off()

# ==============================================================================
# Section 4: Boxplot — group x depth bin
# ==============================================================================

ra_box <- function(res, palette, subtitle) {
  ggplot(res$df, aes(x = depth_type, y = RA_pct)) +
    geom_boxplot(aes(fill = group_lab), outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.18, size = 1.1, alpha = 0.5, colour = "black") +
    scale_fill_manual(values = group_palette(res, palette), guide = "none") +
    facet_wrap(~ group_lab, nrow = 1) +
    labs(x = "Depth bin", y = "Summed relative abundance (%)",
         subtitle = subtitle) +
    theme_bw(base_size = 15) +
    theme(strip.text    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 13),
          axis.title    = element_text(size = 17, face = "bold"),
          axis.text     = element_text(size = 13))
}

kw_by_group <- function(res, grouping) {
  do.call(rbind, lapply(levels(res$df$group_lab), function(g) {
    d  <- res$df[res$df$group_lab == g, ]
    kt <- kruskal.test(RA_pct ~ depth_type, data = d)
    data.frame(grouping = grouping, group = g,
               chi_sq = round(unname(kt$statistic), 3),
               df = unname(kt$parameter),
               p_value = signif(kt$p.value, 4), row.names = NULL)
  }))
}

kw <- rbind(kw_by_group(ra_rank, "fraction_rank_order"),
            kw_by_group(ra_frac, "filter_fraction"),
            kw_by_group(ra_clus, "processing_cluster"))

message("     Kruskal-Wallis, summed RA ~ depth within each group:")
for (i in seq_len(nrow(kw))) {
  message(sprintf("       %-18s %-22s chi-sq = %6.3f, df = %d, p = %s",
                  kw$grouping[i], kw$group[i], kw$chi_sq[i], kw$df[i],
                  format(kw$p_value[i])))
}

ra_summary <- dplyr::bind_rows(
  dplyr::mutate(ra_rank$df, grouping = "fraction_rank_order"),
  dplyr::mutate(ra_frac$df, grouping = "filter_fraction"),
  dplyr::mutate(ra_clus$df, grouping = "processing_cluster")) %>%
  dplyr::group_by(grouping, group_lab, depth_type) %>%
  dplyr::summarise(n = dplyr::n(),
                   median = round(median(RA_pct), 2),
                   mean   = round(mean(RA_pct), 2),
                   sd     = round(sd(RA_pct), 2),
                   min    = round(min(RA_pct), 2),
                   max    = round(max(RA_pct), 2),
                   .groups = "drop") %>%
  dplyr::left_join(kw[, c("grouping", "group", "chi_sq", "df", "p_value")],
                   by = c("grouping", "group_lab" = "group"))

write.csv(ra_summary,
          here("output", "survey", "maps", "survey_ASVgroup_RA_summary.csv"),
          row.names = FALSE)

pdf(here("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA_boxplot.pdf"),
    width = 11, height = 7)
print(ra_box(ra_rank, cc_rank,
             "ASVs grouped by the rank order of their fraction means (top three orders + Other)"))
print(ra_box(ra_frac, cc_frac,
             "ASVs grouped by filter fraction (>50 % of their abundance in one fraction)"))
print(ra_box(ra_clus, cc_clus,
             "ASVs grouped by the processing cluster where they are most abundant"))
dev.off()

message("S26_processing_cluster_asv_map.R: done.")
message("  PDF: output/survey/maps/ASGARD_survey_ASVgroup_RA.pdf")
message("  PDF: output/survey/maps/ASGARD_survey_ASVgroup_RA_boxplot.pdf")
message("  CSV: output/survey/maps/survey_ASVgroup_RA.csv")
message("  CSV: output/survey/maps/survey_ASVgroup_RA_summary.csv")
