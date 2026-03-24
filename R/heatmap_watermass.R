### heatmap_watermass.R
### ASGARD 2017 — 16S Heatmaps with Water Mass RowSideColors
### 水塊別色付き16Sヒートマップ
###
### REQUIRES (from P01/P02/S01):
###   asgard_filtered_p_hm2 - processing ASV matrix 78×221
###   meta_asgard_p2        - processing metadata 78 samples
###   asv_rgb2              - ternary col-side colours (length 221)
###   asgard_frtprop        - survey fourth-root matrix 181×258
###   meta_asgard           - survey metadata 181 samples
###
### OUTPUT:
###   output/heatmaps/ASGARD_hm_processing_watermass.pdf
###   output/survey/heatmaps/ASGARD_hm_survey_watermass.pdf

library(gplots)
library(viridis)
library(vegan)
library(here)

# ==============================================================================
# Section 1: 水塊分類関数 / Water mass classifier (Danielson et al. 2020, Table 2)
# ==============================================================================

classify_watermass <- function(temp, sal) {
  wm <- rep(NA_character_, length(temp))
  for (i in seq_along(temp)) {
    s <- sal[i]; t <- temp[i]
    if (is.na(s) | is.na(t)) { wm[i] <- "NA"; next }
    if (s < 30.8  & t >= 3)                        { wm[i] <- "wCW"
    } else if (s < 30.8  & t < 3)                  { wm[i] <- "IMW & cCW"
    } else if (s >= 30.8 & s < 32.5 & t >= 0 & t < 3) { wm[i] <- "cSW"
    } else if (s >= 30.8 & s < 33.4 & t >= 3)      { wm[i] <- "wSW"
    } else if (s >= 32.5 & s < 33.8 & t >= 0 & t < 3) { wm[i] <- "AnW"
    } else if (s >= 30.8 & s < 33.8 & t >= -1 & t < 0) { wm[i] <- "MWW"
    } else if (s >= 30.8 & s < 33.8 & t < -1)      { wm[i] <- "WW"
    } else if ((s >= 33.4 & t >= 3) | (s >= 33.8 & t < 3)) { wm[i] <- "AtlW & BBW"
    } else { wm[i] <- "Unclassified" }
  }
  wm
}

# 水塊カラーパレット / Water mass colour palette
wm_colors <- c(
  "wCW"         = "#E53935",
  "IMW & cCW"   = "#1E88E5",
  "cSW"         = "#43A047",
  "wSW"         = "#FB8C00",
  "AnW"         = "#8E24AA",
  "MWW"         = "#00ACC1",
  "WW"          = "#1A237E",
  "AtlW & BBW"  = "#6D4C41",
  "NA"          = "gray50",
  "Unclassified" = "gray80"
)

# ==============================================================================
# Section 2: Processing heatmap (78 samples)
# ==============================================================================

wm_p <- classify_watermass(meta_asgard_p2$temp, meta_asgard_p2$salinity)
names(wm_p) <- rownames(meta_asgard_p2)
wm_rgb_p <- wm_colors[wm_p]

message("Processing water mass counts:")
print(table(wm_p))

dir.create(here("output", "heatmaps"), showWarnings = FALSE, recursive = TRUE)
pdf(file = here::here("output", "heatmaps", "ASGARD_hm_processing_watermass.pdf"),
    width = 20, height = 20)

heatmap.2(
  (asgard_filtered_p_hm2)^.25,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = wm_rgb_p,
  ColSideColors = asv_rgb2,
  margins       = c(15, 15),
  scale         = "none",
  main          = "Processing 16S — Water mass (Danielson 2020)",
  trace         = "none",
  cexCol        = 0.2,
  cexRow        = 0.2,
  labRow        = meta_asgard_p2$side
)

dev.off()
message("Processing heatmap saved: output/heatmaps/ASGARD_hm_processing_watermass.pdf")

# ==============================================================================
# Section 3: Survey heatmap (181 samples)
# ==============================================================================

wm_s <- classify_watermass(meta_asgard$temp, meta_asgard$salinity)
names(wm_s) <- rownames(meta_asgard)
wm_rgb_s <- wm_colors[wm_s]

message("Survey water mass counts:")
print(table(wm_s))

dir.create(here("output", "survey", "heatmaps"), showWarnings = FALSE, recursive = TRUE)
pdf(file = here::here("output", "survey", "heatmaps", "ASGARD_hm_survey_watermass.pdf"),
    width = 20, height = 20)

asgard_frtmat <- asgard_frtprop[, colSums(asgard_frtprop) > 0]

heatmap.2(
  asgard_frtmat,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = wm_rgb_s,
  margins       = c(15, 15),
  scale         = "none",
  main          = "Survey 16S — Water mass (Danielson 2020)",
  trace         = "none",
  cexCol        = 0.2,
  cexRow        = 0.3
)

dev.off()
message("Survey heatmap saved: output/survey/heatmaps/ASGARD_hm_survey_watermass.pdf")
