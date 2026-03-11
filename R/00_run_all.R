### 00_run_all.R
### ASGARD 2017 — Combined Master Run Script (Processing + Survey)
### マスター実行スクリプト（プロセシング + サーベイ）
###
### Runs both pipelines sequentially. 00_setup.R is sourced once and shared.
### 両パイプラインを順番に実行します。00_setup.R は一度だけ読み込まれます。
###
### Usage / 使い方:
###   setwd("~/Desktop/ASGARD2017_analysis")
###   source("R/00_run_all.R")
###
### To run each pipeline independently:
###   source("R/00_run_processing.R")   # processing stations only
###   source("R/00_run_survey.R")       # survey stations only

library(here)
source(here("R/00_run_processing.R"), echo = FALSE)
source(here("R/00_run_survey.R"),     echo = FALSE)
