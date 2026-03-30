### S16_NO3_summary.R
### ASGARD 2017 Survey — NO3 Concentration Summary
### NO3濃度のサマリー（全体・クラスター別・サンプル別）
###
### REQUIRES (from S01, S02):
###   meta_asgard   - metadata 181 samples
###   clusnum10     - 10-cluster assignments
###
### OUTPUT:
###   output/survey/NO3_summary_overall.csv
###   output/survey/NO3_summary_by_cluster.csv
###   output/survey/NO3_summary_by_sample.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Prepare data / データ準備
# ==============================================================================

df <- meta_asgard
df$cluster <- factor(clusnum10[rownames(df)], levels = as.character(1:10))
df$sample_id <- rownames(df)

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 2: Overall summary / 全体のサマリー
# ==============================================================================

overall <- data.frame(
  n = sum(!is.na(df[["NO3(uM)"]])),
  n_NA = sum(is.na(df[["NO3(uM)"]])),
  mean = round(mean(df[["NO3(uM)"]], na.rm = TRUE), 2),
  sd = round(sd(df[["NO3(uM)"]], na.rm = TRUE), 2),
  median = round(median(df[["NO3(uM)"]], na.rm = TRUE), 2),
  Q1 = round(quantile(df[["NO3(uM)"]], 0.25, na.rm = TRUE), 2),
  Q3 = round(quantile(df[["NO3(uM)"]], 0.75, na.rm = TRUE), 2),
  min = round(min(df[["NO3(uM)"]], na.rm = TRUE), 2),
  max = round(max(df[["NO3(uM)"]], na.rm = TRUE), 2)
)

message("=== Overall NO3 ===")
message("  mean: ", overall$mean, " uM (SD: ", overall$sd, ")")
message("  median: ", overall$median, " (Q1=", overall$Q1, ", Q3=", overall$Q3, ")")
message("  range: ", overall$min, " - ", overall$max, " uM")
message("  NA: ", overall$n_NA)

write.csv(overall, here("output", "survey", "NO3_summary_overall.csv"), row.names = FALSE)

# ==============================================================================
# Section 3: Summary by cluster / クラスター別サマリー
# ==============================================================================

by_cluster <- df %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean = round(mean(`NO3(uM)`, na.rm = TRUE), 2),
    sd = round(sd(`NO3(uM)`, na.rm = TRUE), 2),
    median = round(median(`NO3(uM)`, na.rm = TRUE), 2),
    Q1 = round(quantile(`NO3(uM)`, 0.25, na.rm = TRUE), 2),
    Q3 = round(quantile(`NO3(uM)`, 0.75, na.rm = TRUE), 2),
    min = round(min(`NO3(uM)`, na.rm = TRUE), 2),
    max = round(max(`NO3(uM)`, na.rm = TRUE), 2),
    n_NA = sum(is.na(`NO3(uM)`)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean))

message("\n=== NO3 by cluster ===")
print(as.data.frame(by_cluster))

write.csv(by_cluster, here("output", "survey", "NO3_summary_by_cluster.csv"), row.names = FALSE)

# ==============================================================================
# Section 4: Per-sample values / サンプル別
# ==============================================================================

by_sample <- df %>%
  select(sample_id, station, depth_type, depth_m, cluster, `NO3(uM)`) %>%
  arrange(cluster, desc(`NO3(uM)`))

write.csv(by_sample, here("output", "survey", "NO3_summary_by_sample.csv"), row.names = FALSE)

message("\nS16_NO3_summary.R: done.")
message("  CSVs: NO3_summary_overall.csv, NO3_summary_by_cluster.csv, NO3_summary_by_sample.csv")
