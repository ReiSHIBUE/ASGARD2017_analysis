### 00_dirs.R
### ASGARD 2017 — Output directory setup / 出力ディレクトリの作成
###
### The analysis scripts write PDFs and CSVs with pdf()/write.csv() without
### creating their parent directories, and git does not track empty directories.
### In a fresh checkout that makes the very first plotting call fail (e.g.
### P02_ternary_plots.R -> output_p/ternary/), which then cascades into
### "object 'clusnum_p' not found" for every later script.
###
### 解析スクリプトは出力ディレクトリを作らずにPDF/CSVを書き出すため、
### 新しいクローンでは最初の描画で失敗する。ここで全ディレクトリを作成する。
###
### Sourced by 00_setup.R, so every entry point (00_run_all.R,
### 00_run_processing.R, 00_run_survey.R, or a manual session) gets the
### directories before any script runs.

OUTPUT_DIRS <- c(
  # survey pipeline (S01–S25)
  "output",
  "output/heatmaps",
  "output/comparison",
  "output/comparison/cache",
  "output/survey",
  "output/survey/alpha_diversity",
  "output/survey/beta_diversity",
  "output/survey/bubble",
  "output/survey/crosstable",
  "output/survey/dbrda",
  "output/survey/env_table",
  "output/survey/heatmaps",
  "output/survey/IndVal",
  "output/survey/maps",
  "output/survey/network",
  "output/survey/NO3 conc",
  "output/survey/stratification",
  "output/survey/taxonomy",
  "output/survey/transect_compare",
  "output/survey/waffle",
  "output/TS plot",
  # processing pipeline (P01–P14)
  "output_p",
  "output_p/beta_diversity",
  "output_p/boxplots",
  "output_p/bubble",
  "output_p/cluster_summary",
  "output_p/dbrda",
  "output_p/env_boxplots",
  "output_p/heatmaps",
  "output_p/maps",
  "output_p/ternary",
  "output_p/treewalk",
  "output_p/TS"
)

#' Create every output directory used by the pipelines.
#' 全出力ディレクトリを作成する（既存の場合は何もしない）。
ensure_output_dirs <- function(base = here::here(), dirs = OUTPUT_DIRS) {
  for (d in dirs) {
    dir.create(file.path(base, d), showWarnings = FALSE, recursive = TRUE)
  }
  invisible(dirs)
}
