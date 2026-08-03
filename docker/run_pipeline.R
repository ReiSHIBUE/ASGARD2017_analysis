### run_pipeline.R
### Container / CI entry point: run both pipelines and report a hard pass/fail.
### コンテナ・CI用エントリポイント：両パイプラインを実行し合否を判定する
###
### Usage:  Rscript docker/run_pipeline.R
###
### Exits 0 only when every script completed, except the map scripts when no
### STADIA_MAPS_KEY is configured (they cannot work without an API key).
### Writes output/repro/script_status.csv and output/repro/key_numbers.txt.

library(here)

# Scripts that cannot run without a Stadia Maps API key.
MAP_SCRIPTS <- c("R/P04_maps.R", "R/S06_maps.R", "R/S16_NO3_summary.R",
                 "R/S23_dbo3_envmap.R")

has_key <- nzchar(Sys.getenv("STADIA_MAPS_KEY"))
if (!has_key) {
  message("STADIA_MAPS_KEY is not set — map scripts are expected to fail:")
  message("  ", paste(MAP_SCRIPTS, collapse = ", "))
}

t0 <- Sys.time()

status <- source(here("R", "00_run_all.R"))$value

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# ==============================================================================
# Record the run
# ==============================================================================

dir.create(here("output", "repro"), showWarnings = FALSE, recursive = TRUE)
write.csv(status, here("output", "repro", "script_status.csv"), row.names = FALSE)

sink(here("output", "repro", "key_numbers.txt"))
cat("=== RUN ===\n")
cat("elapsed (min):                            ", elapsed, "\n")
cat("scripts completed:                        ",
    sum(status$ok), "/", nrow(status), "\n")
cat("STADIA_MAPS_KEY set:                      ", has_key, "\n")

cat("\n=== SETUP ===\n")
cat("16S prokaryote ASVs:                      ", ncol(seqtab_16Smat), "\n")
cat("samples with >=1000 16S reads:            ", nrow(seqtab_16Smat), "\n")

cat("\n=== SURVEY ===\n")
if (exists("asgard_filtered")) {
  cat("survey samples:                           ", nrow(asgard_filtered), "\n")
  cat("survey ASVs after abundance filter:       ", ncol(asgard_filtered), "\n")
}
if (exists("clusnum11")) {
  cat("\nsurvey cluster sizes (11 clusters):\n")
  print(table(clusnum11))
}

cat("\n=== PROCESSING ===\n")
if (exists("asgard_filtered_p_hm2")) {
  cat("processing samples (>=5000 reads):        ", nrow(asgard_filtered_p_hm2), "\n")
  cat("processing ASVs after abundance filter:   ", ncol(asgard_filtered_p_hm2), "\n")
}
if (exists("clusnum_p")) {
  cat("\nprocessing cluster sizes (4 clusters):\n")
  print(table(clusnum_p))
}

cat("\n=== SESSION ===\n")
print(sessionInfo())
sink()

# ==============================================================================
# Pass / fail
# ==============================================================================

expected_failures <- if (has_key) character(0) else MAP_SCRIPTS
unexpected <- status$script[!status$ok & !(status$script %in% expected_failures)]

message("\n", strrep("=", 60))
message(sprintf("%d / %d scripts completed in %s min.",
                sum(status$ok), nrow(status), elapsed))

if (length(unexpected)) {
  message("\nUNEXPECTED FAILURES:")
  for (s in unexpected) {
    message("  ", s, ": ", status$message[status$script == s])
  }
  message(strrep("=", 60))
  quit(status = 1L)
}

if (any(!status$ok)) {
  message("Only the expected map-script failures occurred (no API key).")
}
message("Pipeline OK.")
message(strrep("=", 60))
