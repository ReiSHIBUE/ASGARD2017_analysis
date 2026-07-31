### 00_run_all.R
### ASGARD 2017 — Run all pipelines with corrected metadata (xlsx)
### 修正済みメタデータ (xlsx) を使用して全パイプラインを実行
###
### Usage:
###   source("R/00_run_all.R")     # from the repository root
###
### Output goes to output/ and output_p/, the same place the other runners
### write to. (Earlier versions redirected to output2/ by wrapping pdf() and
### sink(), but write.csv()/png()/ggsave() were not redirected, so results were
### split across two trees. If you want to keep a previous run, copy output/ and
### output_p/ aside before sourcing this file.)
###
### NOTE: map scripts (P04, S06, S16, S17, S23) need a Stadia Maps API key in
### .Renviron:  STADIA_MAPS_KEY=your_key_here

library(here)
library(readxl)

# ==============================================================================
# Step 1: Source original setup (loads ASV tables, taxonomy, output dirs, etc.)
# ==============================================================================

source(here("R", "00_setup.R"), echo = FALSE)

# ==============================================================================
# Step 2: Override meta_denovo_2 with corrected xlsx metadata
# ==============================================================================

message("\nOverriding meta_denovo_2 with corrected xlsx file...")

meta_denovo_2_new <- read_excel(
  here("data", "raw", "ASGARD_metaRDS_denovo_2.xlsx")
)

# The xlsx has an extra first column (...1) containing the original rownames.
# Convert to data.frame with proper rownames to match the original RDS format.
meta_denovo_2 <- as.data.frame(meta_denovo_2_new)
rownames(meta_denovo_2) <- meta_denovo_2[[1]]
meta_denovo_2 <- meta_denovo_2[, -1]  # drop the rowname column

message("  New meta_denovo_2: ", nrow(meta_denovo_2), " x ", ncol(meta_denovo_2))

# ==============================================================================
# Step 3: Script list
# 変数名は .run_ 接頭辞付き（スクリプトが i / script などを上書きするため）
# Loop variables are prefixed with .run_ because several scripts assign to
# short names like `i` and `script` in the global environment.
# ==============================================================================

.run_scripts <- c(
  # processing pipeline
  "R/P01_data_prep.R",
  "R/P02_ternary_plots.R",
  "R/P03_heatmaps_16S.R",
  "R/P04_maps.R",
  "R/P05_beta_diversity_pcoa.R",
  "R/P06_dbrda.R",
  "R/P07_18S_heatmaps.R",
  "R/P08_esv_heatmap.R",
  "R/P09_cluster_summary.R",
  "R/P11_env_boxplots.R",
  "R/P12_genus_bubble.R",
  "R/P13_ternary_v2.R",
  "R/P14_station_region_summary.R",
  "R/TS_diagram_p.R",
  # survey pipeline
  "R/S01_data_prep.R",
  "R/S02_heatmaps_16S.R",
  "R/S03_seq_depth.R",
  "R/S04_beta_diversity_pcoa.R",
  "R/S05_dbrda.R",
  "R/S06_maps.R",
  "R/S07_alpha_diversity.R",
  "R/S08_taxonomy.R",
  "R/S09_network.R",
  "R/S10_permanova.R",
  "R/S11_crosstable.R",
  "R/S12_indval.R",
  "R/S13_indval_18S.R",
  "R/S14_sampling_period.R",
  "R/S15_wSW_misclassification.R",
  "R/S16_NO3_summary.R",
  "R/S17_pca_env.R",
  "R/S18_cluster_specific_ASVs.R",
  "R/S19_env_boxplots.R",
  "R/S20_stratification.R",
  "R/S21_genus_bubble.R",
  "R/S22_dbo3_waffle.R",
  "R/S23_dbo3_envmap.R",
  "R/S25_transect_compare.R",
  "R/TS_diagram.R",
  "R/heatmap_watermass.R",
  # cross-pipeline: P10 also compares processing stations against the survey
  # stations, so it needs objects from both pipelines and runs last.
  "R/P10_water_mass.R",
  # rebuild the comparison caches from the session, then the report that reads
  # them (previously the caches were committed artifacts with no source script)
  "R/98_build_caches.R",
  "R/comparison_report.R"
)

# ==============================================================================
# Step 4: Run
# ==============================================================================

.run_ok  <- logical(length(.run_scripts))
.run_msg <- character(length(.run_scripts))

for (.run_i in seq_along(.run_scripts)) {
  .run_this <- .run_scripts[.run_i]
  message("\n", strrep("=", 60))
  message("Running: ", .run_this)
  message(strrep("=", 60))

  .run_res <- tryCatch({
    source(here(.run_this), echo = FALSE)
    list(ok = TRUE, msg = "")
  }, error = function(e) {
    message("ERROR in ", .run_this, ": ", conditionMessage(e))
    list(ok = FALSE, msg = conditionMessage(e))
  })

  # a failing script can leave a graphics device or sink open
  while (dev.cur() > 1) dev.off()
  while (sink.number() > 0) sink()

  .run_ok[.run_i]  <- .run_res$ok
  .run_msg[.run_i] <- .run_res$msg
}

# ==============================================================================
# Step 5: Report — do not claim success when scripts failed
# 失敗したスクリプトがある場合は成功と報告しない
# ==============================================================================

message("\n", strrep("=", 60))
message(sprintf("%d / %d scripts completed.", sum(.run_ok), length(.run_ok)))

if (any(!.run_ok)) {
  message("\nFAILED:")
  for (.run_i in which(!.run_ok)) {
    message("  ", .run_scripts[.run_i], ": ", .run_msg[.run_i])
  }
  message("\nOutputs from the failed scripts (and anything downstream of them) ",
          "are missing or stale.")
} else {
  message("All pipelines complete. Check output/ and output_p/ for PDFs.")
}
message(strrep("=", 60))

invisible(data.frame(script = .run_scripts, ok = .run_ok, message = .run_msg))
