### S13_indval_18S.R
### ASGARD 2017 Survey — Indicator Species Analysis (IndVal) for 18S
### 18S 指標種解析（IndVal）
###
### REQUIRES (from S01, S02, 00_setup):
###   clusnum10     - 10-cluster assignments
###   fullnamemat   - taxonomy matrix (from 00_setup)
###   data/esv_relabund_by_station_depth.tsv
###
### OUTPUT:
###   output/survey/indval_18S_results.csv

library(tidyverse)
library(indicspecies)
library(here)

# ==============================================================================
# Section 1: Load 18S data / 18Sデータ読み込み
# ==============================================================================

euk_esv <- read.table(here("data", "esv_relabund_by_station_depth.tsv"),
                      sep = "\t", header = TRUE, row.names = 1)
meta_cols <- c("station", "depth_m", "station_type", "filter")
taxon_cols <- setdiff(colnames(euk_esv), meta_cols)

# Keep only survey samples with clusnum10
survey_in <- intersect(names(clusnum10), rownames(euk_esv))
message("Survey samples with 18S data: ", length(survey_in))

euk_mat <- euk_esv[survey_in, taxon_cols]
for (col in taxon_cols) euk_mat[[col]] <- as.numeric(euk_mat[[col]])
euk_mat <- euk_mat[, colSums(euk_mat, na.rm = TRUE) > 0]
message("18S ESVs: ", ncol(euk_mat))

cluster_vec <- clusnum10[survey_in]

# ==============================================================================
# Section 2: IndVal analysis / IndVal解析
# ==============================================================================

set.seed(42)
indval_18S <- multipatt(euk_mat, cluster_vec, func = "IndVal.g",
                        control = how(nperm = 999))

# ==============================================================================
# Section 3: Extract results / 結果抽出
# ==============================================================================

sig <- indval_18S$sign
sig$ESV <- rownames(sig)
sig_only <- sig %>% filter(p.value < 0.05) %>% arrange(p.value)

s_cols <- grep("^s[.]", colnames(sig_only), value = TRUE)
sig_only$indicator_cluster <- apply(sig_only[, s_cols], 1, function(x) {
  paste(which(x == 1), collapse = "+")
})

# Map taxonomy
esv_lookup <- data.frame(
  ESV = fullnamemat[, "ESV"],
  Phylum = fullnamemat[, "Phylum"],
  Class = fullnamemat[, "Class"],
  Order = fullnamemat[, "Order"],
  Family = fullnamemat[, "Family"],
  Genus = fullnamemat[, "Genus"],
  stringsAsFactors = FALSE
)
sig_only <- sig_only %>% left_join(esv_lookup, by = "ESV")

single <- sig_only %>% filter(!grepl("[+]", indicator_cluster))

message("Total significant (p < 0.05): ", nrow(sig_only))
message("Single-cluster indicators: ", nrow(single))

# ==============================================================================
# Section 4: Print results / 結果表示
# ==============================================================================

for (cl in as.character(1:10)) {
  cl_ind <- single %>% filter(indicator_cluster == cl) %>% arrange(p.value)
  message(sprintf("\nCluster %s (%d indicators):", cl, nrow(cl_ind)))
  if (nrow(cl_ind) > 0) {
    for (j in seq_len(nrow(cl_ind))) {
      message(sprintf("  %-12s %-20s %-20s stat=%.3f p=%.3f",
        cl_ind$ESV[j], cl_ind$Class[j], cl_ind$Genus[j],
        cl_ind$stat[j], cl_ind$p.value[j]))
    }
  } else {
    message("  (none)")
  }
}

# ==============================================================================
# Section 5: Save / 保存
# ==============================================================================

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)
write.csv(sig_only, here("output", "survey", "indval_18S_results.csv"), row.names = FALSE)

message("\nS13_indval_18S.R: done.")
message("  CSV: output/survey/indval_18S_results.csv")
