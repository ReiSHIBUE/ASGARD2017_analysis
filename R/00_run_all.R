### 00_run_all.R
### ASGARD 2017 Processing Site Analysis — Master Run Script
### マスター実行スクリプト
###
### Run this script to execute the full pipeline in sequence.
### このスクリプトを実行すると、全パイプラインが順番に実行されます。
###
### BEFORE RUNNING: the following objects must already exist in your R session
### (loaded by your upstream de novo / ODV processing script):
###   - seqtab_16Sprop
###   - seqtab_16Smat
###   - seqtab_16Spropcol_mra
###   - seqtab_filt
###   - shorternames
###   - meta_denovo_2
###
### Usage / 使い方:
###   setwd("~/Desktop/ASGARD2017_analysis")
###   source("R/00_run_all.R")

library(here)

scripts <- c(
  "R/01_data_prep.R",
  "R/02_ternary_plots.R",
  "R/03_heatmaps_16S.R",
  "R/04_maps.R",
  "R/05_beta_diversity_pcoa.R",
  "R/06_dbrda.R",
  "R/07_18S_heatmaps.R",
  "R/08_esv_heatmap.R"
)

for (script in scripts) {
  message("\n", strrep("=", 60))
  message("Running: ", script)
  message(strrep("=", 60))
  source(here(script), echo = FALSE)
}

message("\n", strrep("=", 60))
message("Pipeline complete. Check output/ for PDFs.")
message(strrep("=", 60))
