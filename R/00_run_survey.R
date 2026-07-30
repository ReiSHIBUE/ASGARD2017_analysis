### 00_run_survey.R
### ASGARD 2017 Survey Site Analysis — Master Run Script
### マスター実行スクリプト（サーベイサイト）
###
### Run this script to execute the full survey pipeline in sequence.
### このスクリプトを実行すると、全サーベイパイプラインが順番に実行されます。
###
### Usage / 使い方:
###   setwd("~/Desktop/ASGARD2017_analysis")
###   source("R/00_run_survey.R")
###
### 00_setup.R loads all raw RDS files and builds upstream objects automatically.
### 00_setup.R がRDSファイルを読み込み、上流オブジェクトを自動的に構築します。
###
### NOTE: Maps require a Stadia Maps API key stored in .Renviron:
###   STADIA_MAPS_KEY=your_key_here

library(here)

# Create output subdirectories / 出力サブディレクトリを作成
source(here("R", "00_dirs.R"))
ensure_output_dirs()

scripts <- c(
  "R/00_setup.R",
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
  "R/TS_diagram.R"
  # heatmap_watermass.R needs the processing objects too — see 00_run_all.R
)

# 各スクリプトのエラーを記録して継続する / Record per-script errors and continue,
# so one failure (e.g. a missing Stadia Maps key) does not abort the whole run.
.run_ok  <- logical(length(scripts))
.run_msg <- character(length(scripts))

for (.run_i in seq_along(scripts)) {
  .run_this <- scripts[.run_i]
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

  while (dev.cur() > 1) dev.off()
  while (sink.number() > 0) sink()

  .run_ok[.run_i]  <- .run_res$ok
  .run_msg[.run_i] <- .run_res$msg
}

message("\n", strrep("=", 60))
message(sprintf("%d / %d scripts completed.", sum(.run_ok), length(.run_ok)))
if (any(!.run_ok)) {
  message("\nFAILED:")
  for (.run_i in which(!.run_ok)) message("  ", scripts[.run_i], ": ", .run_msg[.run_i])
  message("\nOutputs from the failed scripts (and anything downstream) are ",
          "missing or stale.")
} else {
  message("Survey pipeline complete. Check output/survey/ for PDFs.")
}
message(strrep("=", 60))
