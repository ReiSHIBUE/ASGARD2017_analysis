### 00_run_all.R
### ASGARD 2017 — Run all pipelines with corrected metadata (xlsx)
### 修正済みメタデータ (xlsx) を使用して全パイプラインを実行し、output/ に出力
###
### Usage:
###   setwd("~/Desktop/ASGARD2017_analysis")
###   source("R/00_run_all.R")

library(here)
library(readxl)

# ==============================================================================
# Step 1: Source original setup (loads ASV tables, taxonomy, etc.)
# ==============================================================================

source(here("R", "00_setup.R"), echo = FALSE)

# ==============================================================================
# Step 2: Override meta_denovo_2 with corrected xlsx metadata
# ==============================================================================

message("\nOverriding meta_denovo_2 with corrected xlsx file...")

meta_denovo_2_new <- read_excel(
  "~/Desktop/ASGARD_metaRDS_denovo_2.xlsx"
)

# The xlsx has an extra first column (...1) containing the original rownames.
# Convert to data.frame with proper rownames to match the original RDS format.
meta_denovo_2 <- as.data.frame(meta_denovo_2_new)
rownames(meta_denovo_2) <- meta_denovo_2[[1]]
meta_denovo_2 <- meta_denovo_2[, -1]  # drop the rowname column

message("  New meta_denovo_2: ", nrow(meta_denovo_2), " x ", ncol(meta_denovo_2))

# ==============================================================================
# Step 3: Override pdf() and sink() to redirect output/ -> output2/
# ==============================================================================

original_pdf  <- pdf
original_sink <- sink

pdf <- function(file = if (onefile) "Rplots.pdf" else "Rplot%03d.pdf", ...,
                onefile = TRUE) {
  if (is.character(file)) {
    file <- sub("/output/", "/output2/", file, fixed = TRUE)
  }
  original_pdf(file = file, ..., onefile = onefile)
}

sink <- function(file = NULL, ...) {
  if (!is.null(file) && is.character(file)) {
    file <- sub("/output/", "/output2/", file, fixed = TRUE)
  }
  original_sink(file = file, ...)
}

# ==============================================================================
# Step 4: Create output2/ subdirectories (in case they don't exist)
# ==============================================================================

output2_dirs <- c(
  "output2/heatmaps", "output2/maps", "output2/ternary",
  "output2/beta_diversity", "output2/dbrda", "output2/boxplots",
  "output2/survey/heatmaps", "output2/survey/maps",
  "output2/survey/beta_diversity", "output2/survey/dbrda",
  "output2/survey/alpha_diversity", "output2/survey/taxonomy",
  "output2/survey/network"
)

for (d in output2_dirs) {
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# Step 5: Run processing pipeline (P01-P08)
# ==============================================================================

processing_scripts <- c(
  "R/P01_data_prep.R",
  "R/P02_ternary_plots.R",
  "R/P03_heatmaps_16S.R",
  "R/P04_maps.R",
  "R/P05_beta_diversity_pcoa.R",
  "R/P06_dbrda.R",
  "R/P07_18S_heatmaps.R",
  "R/P08_esv_heatmap.R"
)

for (script in processing_scripts) {
  message("\n", strrep("=", 60))
  message("Running: ", script)
  message(strrep("=", 60))
  tryCatch(
    source(here(script), echo = FALSE),
    error = function(e) message("ERROR in ", script, ": ", conditionMessage(e))
  )
}

# ==============================================================================
# Step 6: Run survey pipeline (S01-S10)
# ==============================================================================

survey_scripts <- c(
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

for (script in survey_scripts) {
  message("\n", strrep("=", 60))
  message("Running: ", script)
  message(strrep("=", 60))

  # S04 adds a "Sample" column to meta_asgard, which causes
  # rownames_to_column(meta_asgard, var = "Sample") in S07 to fail
  # with a duplicate column name error. Remove it before S07.
  if (script == "R/S07_alpha_diversity.R" && "Sample" %in% colnames(meta_asgard)) {
    meta_asgard <- meta_asgard[, colnames(meta_asgard) != "Sample", drop = FALSE]
  }

  tryCatch(
    source(here(script), echo = FALSE),
    error = function(e) message("ERROR in ", script, ": ", conditionMessage(e))
  )
}

# ==============================================================================
# Step 7: Restore original pdf() and sink()
# ==============================================================================

pdf  <- original_pdf
sink <- original_sink

message("\n", strrep("=", 60))
message("All pipelines complete with corrected metadata.")
message("Output saved to: output2/")
message(strrep("=", 60))

