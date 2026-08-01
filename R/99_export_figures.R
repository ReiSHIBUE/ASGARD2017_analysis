### 99_export_figures.R
### ASGARD 2017 — Export the manuscript figures as publication-resolution PNGs
### 論文用の図を高解像度PNGとして書き出す
###
### Rasterises the specific page of each figure PDF into output/figures_png/.
### Run after the pipelines have produced the PDFs:
###   Rscript R/99_export_figures.R          # 300 dpi (default)
###   FIG_DPI=600 Rscript R/99_export_figures.R
###
### 各図のPDFページをPNGに変換する。パイプライン実行後に走らせること。
###
### Uses pdftools if installed, otherwise the pdftoppm command line tool
### (poppler-utils). pdftools が無い場合は pdftoppm を使う。

library(here)

DPI     <- as.integer(Sys.getenv("FIG_DPI", "300"))
OUT_DIR <- here("output", "figures_png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# name = output file stem; pdf = source; page = page to rasterise
FIGURES <- list(
  list(name = "fig03_survey_heatmap",
       pdf  = c("output", "survey", "heatmaps", "ASGARD_hm_survey_16S_11clusters.pdf"),
       page = 4),
  list(name = "fig04_survey_map_clusters",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_map_11clusters.pdf"),
       page = 2),
  list(name = "fig04b_survey_map_divisions",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_map_11clusters_detail.pdf"),
       page = 3),
  list(name = "fig05_TS_survey_clusters",
       pdf  = c("output", "TS plot", "TS_diagram_11clusters.pdf"),
       page = 2),
  list(name = "fig15_processing_heatmap",
       pdf  = c("output_p", "heatmaps", "ASGARD_hm_processing_5000over.pdf"),
       page = 5),
  list(name = "fig16_ternary_cluster_and_class",
       pdf  = c("output_p", "ternary", "ASGARD_ternary_3pages.pdf"),
       page = 9),
  list(name = "fig16A_ternary_clusters",
       pdf  = c("output_p", "ternary", "ASGARD_ternary_3pages.pdf"),
       page = 7),
  list(name = "fig16B_ternary_classes",
       pdf  = c("output_p", "ternary", "ASGARD_ternary_3pages.pdf"),
       page = 8),
  list(name = "fig17_processing_map_by_depth",
       pdf  = c("output_p", "maps", "ASGARD_processing_map_4clusters_by_depth.pdf"),
       page = 1),
  list(name = "fig18_TS_processing_clusters",
       pdf  = c("output_p", "TS", "TS_diagram_4clusters.pdf"),
       page = 2),
  list(name = "new01_shannon_by_cluster",
       pdf  = c("output_p", "alpha_diversity", "ASGARD_shannon_processing.pdf"),
       page = 1),
  list(name = "new02_shannon_by_filter",
       pdf  = c("output_p", "alpha_diversity", "ASGARD_shannon_processing.pdf"),
       page = 2),
  list(name = "new03_survey_map_ASVgroup_RA_by_rank_order",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA.pdf"),
       page = 1),
  list(name = "new03b_survey_map_ASVgroup_RA_by_fraction",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA.pdf"),
       page = 2),
  list(name = "new03c_survey_map_ASVgroup_RA_by_cluster",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA.pdf"),
       page = 3),
  list(name = "new04_ASVgroup_RA_boxplot_by_rank_order",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA_boxplot.pdf"),
       page = 1),
  list(name = "new04b_ASVgroup_RA_boxplot_by_fraction",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA_boxplot.pdf"),
       page = 2),
  list(name = "new04c_ASVgroup_RA_boxplot_by_cluster",
       pdf  = c("output", "survey", "maps", "ASGARD_survey_ASVgroup_RA_boxplot.pdf"),
       page = 3),
  list(name = "new05_fraction_rank_order",
       pdf  = c("output_p", "fraction_order", "ASGARD_fraction_rank_order.pdf"),
       page = 1),
  list(name = "new06_fraction_rank_abundance_hist_log",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_rank_abundance_hist.pdf"),
       page = 1),
  list(name = "new06b_fraction_rank_abundance_hist_by_fraction_log",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_rank_abundance_hist.pdf"),
       page = 2),
  list(name = "new07_fraction_rank_abundance_hist_4throot",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_rank_abundance_hist.pdf"),
       page = 3),
  list(name = "new07b_fraction_rank_abundance_hist_by_fraction_4throot",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_rank_abundance_hist.pdf"),
       page = 4),
  list(name = "new08_fraction_logratio_volcano_by_rank",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_logratio_volcano.pdf"),
       page = 1),
  list(name = "new08b_fraction_logratio_volcano_by_class",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_logratio_volcano.pdf"),
       page = 2),
  list(name = "new08c_fraction_logratio_volcano_survey_weighted",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_logratio_volcano.pdf"),
       page = 3),
  list(name = "new09_whole_water_recovery_by_quadrant",
       pdf  = c("output_p", "fraction_order",
                "ASGARD_fraction_logratio_volcano.pdf"),
       page = 4)
)

# ==============================================================================
# Rasteriser: pdftools if available, otherwise pdftoppm
# ==============================================================================

has_pdftools <- requireNamespace("pdftools", quietly = TRUE)
pdftoppm     <- Sys.which("pdftoppm")

if (!has_pdftools && !nzchar(pdftoppm)) {
  stop("Need either the pdftools package or the pdftoppm command line tool ",
       "(poppler-utils) to export PNGs.")
}

render_page <- function(pdf_path, page, out_png, dpi) {
  if (has_pdftools) {
    pdftools::pdf_convert(pdf = pdf_path, format = "png", pages = page,
                          dpi = dpi, filenames = out_png, verbose = FALSE)
  } else {
    stem <- sub("\\.png$", "", out_png)
    system2(pdftoppm, c("-f", page, "-l", page, "-r", dpi, "-png",
                        shQuote(pdf_path), shQuote(stem)))
    # pdftoppm appends the page number: stem-<page>.png
    produced <- list.files(dirname(stem),
                           pattern = paste0("^", basename(stem), "-0*", page, "\\.png$"),
                           full.names = TRUE)
    if (length(produced) == 1) file.rename(produced, out_png)
  }
  file.exists(out_png)
}

# ==============================================================================
# Export
# ==============================================================================

message("Exporting figures at ", DPI, " dpi to output/figures_png/")

results <- lapply(FIGURES, function(fig) {
  src <- do.call(here, as.list(fig$pdf))
  out <- file.path(OUT_DIR, paste0(fig$name, ".png"))

  if (!file.exists(src)) {
    message("  MISSING  ", fig$name, "  (", paste(fig$pdf, collapse = "/"), ")")
    return(data.frame(figure = fig$name, ok = FALSE, note = "source PDF not found"))
  }

  n_pages <- if (has_pdftools) pdftools::pdf_info(src)$pages else NA_integer_
  if (!is.na(n_pages) && fig$page > n_pages) {
    message("  MISSING  ", fig$name, "  (page ", fig$page, " > ", n_pages, ")")
    return(data.frame(figure = fig$name, ok = FALSE, note = "page out of range"))
  }

  ok <- render_page(src, fig$page, out, DPI)
  size_mb <- if (ok) round(file.size(out) / 1048576, 2) else NA_real_
  message(sprintf("  %-8s %-42s %s MB", if (ok) "OK" else "FAILED",
                  fig$name, format(size_mb)))
  data.frame(figure = fig$name, ok = ok,
             note = paste0(basename(src), " p", fig$page, ", ", size_mb, " MB"))
})

results <- do.call(rbind, results)
write.csv(results, file.path(OUT_DIR, "figure_manifest.csv"), row.names = FALSE)

message("\n99_export_figures.R: ", sum(results$ok), " / ", nrow(results),
        " figures exported at ", DPI, " dpi.")
message("  Directory: output/figures_png/")
message("  Manifest:  output/figures_png/figure_manifest.csv")
