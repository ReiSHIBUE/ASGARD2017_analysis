### P10_water_mass.R
### ASGARD 2017 Processing — Water-mass classification + survey comparison
### 水塊分類（Danielson 2020）と survey データとの station 比較
###
### REQUIRES (from 00_setup.R, P01–P03 for processing,
###            and S01–S02 for survey comparison):
###   meta_asgard_p2  - processing metadata 78 samples
###   clusnum_p       - 4-cluster assignment
###   meta_asgard     - survey metadata 181 samples
###   clusnum11       - 11-cluster assignment
###   hier_levels_11  - cluster name order
###
### OUTPUT:
###   output_p/cluster_summary/processing_cluster_watermass_counts.csv
###   output_p/cluster_summary/processing_cluster_watermass_percent.csv
###   output_p/cluster_summary/processing_cluster_watermass_test_results.csv
###   output_p/cluster_summary/processing_cluster_watermass_ranking.csv
###   output_p/cluster_summary/processing_cluster_watermass_primary_secondary.csv
###   output_p/cluster_summary/processing_mww_station_depth.csv
###   output_p/cluster_summary/station_processing_vs_survey.csv

library(dplyr)
library(tidyr)
library(here)

dir.create(here("output_p", "cluster_summary"),
           showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Section 1: Water mass classifier (Danielson et al. 2020, Table 2)
# ==============================================================================

classify_wm <- function(temp, sal) {
  wm <- rep(NA_character_, length(temp))
  for (i in seq_along(temp)) {
    s <- sal[i]; t <- temp[i]
    if (is.na(s) | is.na(t)) { wm[i] <- "NA"; next }
    if (s < 30.8 & t >= 3)                            { wm[i] <- "wCW"
    } else if (s < 30.8 & t < 3)                      { wm[i] <- "IMW & cCW"
    } else if (s >= 30.8 & s < 32.5 & t >= 0 & t < 3) { wm[i] <- "cSW"
    } else if (s >= 30.8 & s < 33.4 & t >= 3)         { wm[i] <- "wSW"
    } else if (s >= 32.5 & s < 33.8 & t >= 0 & t < 3) { wm[i] <- "AnW"
    } else if (s >= 30.8 & s < 33.8 & t >= -1 & t < 0){ wm[i] <- "MWW"
    } else if (s >= 30.8 & s < 33.8 & t < -1)         { wm[i] <- "WW"
    } else if ((s >= 33.4 & t >= 3) | (s >= 33.8 & t < 3)) {
      wm[i] <- "AtlW & BBW"
    } else { wm[i] <- "Unclassified" }
  }
  wm
}

# ==============================================================================
# Section 2: Cluster × water mass (processing, k=4)
# ==============================================================================

meta_p <- meta_asgard_p2[names(clusnum_p), ]
meta_p$cluster   <- factor(clusnum_p)
meta_p$watermass <- classify_wm(meta_p$temp, meta_p$salinity)

ct_p  <- table(meta_p$cluster, meta_p$watermass)
pct_p <- round(prop.table(ct_p, 1) * 100, 1)

write.csv(as.data.frame.matrix(ct_p),
          here("output_p", "cluster_summary",
               "processing_cluster_watermass_counts.csv"))
write.csv(as.data.frame.matrix(pct_p),
          here("output_p", "cluster_summary",
               "processing_cluster_watermass_percent.csv"))

# --- Association test: cluster x water mass (matches survey S11) ---
# NA / Unclassified columns are dropped, as in the survey analysis.
ct_p_test <- ct_p[, !colnames(ct_p) %in% c("NA", "Unclassified"), drop = FALSE]
set.seed(1)
chi_wm_p    <- suppressWarnings(chisq.test(ct_p_test, simulate.p.value = TRUE, B = 9999))
fisher_wm_p <- fisher.test(ct_p_test, simulate.p.value = TRUE, B = 9999)

sig_star <- function(p) ifelse(p <= 0.001, "***",
                        ifelse(p <= 0.01,  "**",
                        ifelse(p <= 0.05,  "*", "ns")))

wm_test_p <- data.frame(
  test      = c("Chi-squared (Monte Carlo, B=9999)", "Fisher exact (Monte Carlo, B=9999)"),
  statistic = c(round(unname(chi_wm_p$statistic), 2), NA),
  p_value   = c(chi_wm_p$p.value, fisher_wm_p$p.value),
  sig       = c(sig_star(chi_wm_p$p.value), sig_star(fisher_wm_p$p.value)),
  stringsAsFactors = FALSE)

write.csv(wm_test_p,
          here("output_p", "cluster_summary",
               "processing_cluster_watermass_test_results.csv"),
          row.names = FALSE)

# Primary / secondary water masses
ranking_p <- meta_p %>%
  filter(!watermass %in% c("NA", "Unclassified")) %>%
  count(cluster, watermass) %>%
  group_by(cluster) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(cluster, -n) %>%
  mutate(rank = row_number()) %>%
  ungroup()

primary_secondary_p <- ranking_p %>%
  filter(rank <= 2) %>%
  mutate(role = ifelse(rank == 1, "Primary", "Secondary")) %>%
  select(cluster, role, watermass, n, pct)

write.csv(ranking_p,
          here("output_p", "cluster_summary",
               "processing_cluster_watermass_ranking.csv"),
          row.names = FALSE)
write.csv(primary_secondary_p,
          here("output_p", "cluster_summary",
               "processing_cluster_watermass_primary_secondary.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 3: MWW samples — station × depth breakdown (processing only)
# ==============================================================================

mww_p <- meta_p %>% filter(watermass == "MWW")
mww_p_tab <- mww_p %>%
  count(station, depth_type) %>%
  arrange(station, depth_type)

write.csv(mww_p_tab,
          here("output_p", "cluster_summary",
               "processing_mww_station_depth.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 4: Station comparison — processing vs survey
# ==============================================================================

# This section is the only part of P10 that needs the *survey* objects. When the
# processing pipeline is run on its own (00_run_processing.R) they do not exist,
# so skip the comparison instead of aborting the script.
# このセクションのみ survey オブジェクトを必要とするため、無い場合はスキップする。
if (!exists("meta_asgard")) {
  message("  meta_asgard not found — skipping processing vs survey station ",
          "comparison (run S01_data_prep.R first to include it).")
} else {

stn_p <- meta_p %>%
  group_by(station) %>%
  summarise(n_processing = n(),
            lat_p = round(mean(lat, na.rm = TRUE), 2),
            lon_p = round(mean(lon, na.rm = TRUE), 2),
            .groups = "drop")

stn_s <- meta_asgard %>%
  group_by(station) %>%
  summarise(n_survey = n(),
            lat_s = round(mean(lat, na.rm = TRUE), 2),
            lon_s = round(mean(lon, na.rm = TRUE), 2),
            .groups = "drop")

stn_compare <- full_join(stn_p, stn_s, by = "station") %>%
  mutate(in_processing = !is.na(n_processing),
         in_survey     = !is.na(n_survey),
         shared        = in_processing & in_survey) %>%
  arrange(desc(shared), desc(in_processing),
          -coalesce(n_processing, 0L))

write.csv(stn_compare,
          here("output_p", "cluster_summary",
               "station_processing_vs_survey.csv"),
          row.names = FALSE)

}  # end of survey-comparison section

# ==============================================================================
# Section 5: Log summary
# ==============================================================================

message("\nP10_water_mass.R: done.")
message(sprintf("  MWW samples in processing: %d / 78", nrow(mww_p)))
if (exists("stn_compare")) {
  message(sprintf("  Stations shared: %d  | proc-only: %d  | survey-only: %d",
                  sum(stn_compare$shared),
                  sum(stn_compare$in_processing & !stn_compare$in_survey),
                  sum(!stn_compare$in_processing & stn_compare$in_survey)))
}
message("  CSV: output_p/cluster_summary/processing_cluster_watermass_counts.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_watermass_percent.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_watermass_test_results.csv")
message(sprintf("  Fisher (cluster x water mass): p = %.4g %s",
                fisher_wm_p$p.value, sig_star(fisher_wm_p$p.value)))
message("  CSV: output_p/cluster_summary/processing_cluster_watermass_ranking.csv")
message("  CSV: output_p/cluster_summary/processing_cluster_watermass_primary_secondary.csv")
message("  CSV: output_p/cluster_summary/processing_mww_station_depth.csv")
message("  CSV: output_p/cluster_summary/station_processing_vs_survey.csv")
