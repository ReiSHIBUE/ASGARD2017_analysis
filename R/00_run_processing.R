### 00_run_processing.R
### ASGARD 2017 Processing Site Analysis — Master Run Script
### マスター実行スクリプト（プロセシングサイト）
###
### Run this script to execute the full processing pipeline in sequence.
### このスクリプトを実行すると、全プロセシングパイプラインが順番に実行されます。
###
### Usage / 使い方:
###   setwd("~/Desktop/ASGARD2017_analysis")
###   source("R/00_run_processing.R")
###
### 00_setup.R loads all raw RDS files and builds upstream objects automatically.
### 00_setup.R がRDSファイルを読み込み、上流オブジェクトを自動的に構築します。

library(here)

# Create output subdirectories / 出力サブディレクトリを作成
# The processing scripts write to output_p/, not output/.
source(here("R", "00_dirs.R"))
ensure_output_dirs()

scripts <- c(
  "R/00_setup.R",
  "R/P01_data_prep.R",
  "R/P02_ternary_plots.R",
  "R/P03_heatmaps_16S.R",
  "R/P04_maps.R",
  "R/P05_beta_diversity_pcoa.R",
  "R/P06_dbrda.R",
  "R/P07_18S_heatmaps.R",
  "R/P08_esv_heatmap.R",
  "R/P09_cluster_summary.R",
  "R/P10_water_mass.R",
  "R/P11_env_boxplots.R",
  "R/P12_genus_bubble.R",
  "R/P13_ternary_v2.R",
  "R/P14_station_region_summary.R",
  "R/TS_diagram_p.R"
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
  message("Processing pipeline complete. Check output_p/ for PDFs.")
}
message(strrep("=", 60))
