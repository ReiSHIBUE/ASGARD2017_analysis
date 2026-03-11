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
for (d in c("output/survey/heatmaps", "output/survey/maps",
            "output/survey/beta_diversity", "output/survey/dbrda",
            "output/survey/alpha_diversity", "output/survey/taxonomy",
            "output/survey/network")) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

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
  "R/S10_permanova.R"
)

for (script in scripts) {
  message("\n", strrep("=", 60))
  message("Running: ", script)
  message(strrep("=", 60))
  source(here(script), echo = FALSE)
}

message("\n", strrep("=", 60))
message("Survey pipeline complete. Check output/survey/ for PDFs.")
message(strrep("=", 60))
