### S14_sampling_period.R
### ASGARD 2017 Survey — Sampling Period by Cluster
### クラスター別サンプリング時期
###
### REQUIRES (from S01, S02):
###   meta_asgard   - metadata 181 samples (must have date column)
###   clusnum10     - 10-cluster assignments
###
### OUTPUT:
###   output/survey/sampling_period_by_cluster.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Parse dates / 日付のパース
# ==============================================================================

# Set locale to C for English month abbreviations
Sys.setlocale("LC_TIME", "C")

df <- meta_asgard
df$cluster <- factor(clusnum10[rownames(df)], levels = as.character(1:10))
df$date_parsed <- as.POSIXct(df$date, format = "%b %d %Y %H:%M")
df$date_only <- as.Date(df$date_parsed)

# ==============================================================================
# Section 2: Summarise sampling period per cluster / クラスター別サンプリング期間の集計
# ==============================================================================

sampling_summary <- df %>%
  filter(!is.na(date_only)) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    start = min(date_only),
    end = max(date_only),
    span_days = as.integer(max(date_only) - min(date_only)),
    n_unique_dates = length(unique(date_only)),
    dates = paste(format(sort(unique(date_only)), "%b %d"), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(cluster)

sampling_summary$start <- format(sampling_summary$start, "%b %d")
sampling_summary$end <- format(sampling_summary$end, "%b %d")

# ==============================================================================
# Section 3: Print results / 結果表示
# ==============================================================================

for (i in seq_len(nrow(sampling_summary))) {
  s <- sampling_summary[i, ]
  message(sprintf("Cluster %2s (n=%2d): %s to %s (%d days span, %d unique dates)",
    s$cluster, s$n, s$start, s$end, s$span_days, s$n_unique_dates))
  message(sprintf("  Dates: %s\n", s$dates))
}

# ==============================================================================
# Section 4: Save CSV / CSV保存
# ==============================================================================

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)
write.csv(sampling_summary, here("output", "survey", "sampling_period_by_cluster.csv"), row.names = FALSE)

message("S14_sampling_period.R: done.")
message("  CSV: output/survey/sampling_period_by_cluster.csv")
