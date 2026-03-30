### S15_wSW_misclassification.R
### ASGARD 2017 Survey — Identifying wSW-classified samples with high NO3
### wSWに分類されたが高NO3を持つサンプルの特定
###
### This script identifies samples that are classified as wSW by Danielson (2020)
### T-S definitions but have NO3 >= 10 uM, suggesting they are actually
### Anadyr Water (AnW) that has been warmed above 3°C.
### This is consistent with Hirawake et al. (2021) who showed that T-S
### classification can misidentify water masses during warm anomalies.
###
### REQUIRES (from S01, S02, S11):
###   meta_asgard    - metadata 181 samples
###   clusnum10      - 10-cluster assignments
###   classify_watermass() - from S11
###
### OUTPUT:
###   output/survey/crosstable/wSW_high_NO3_samples.csv

library(tidyverse)
library(here)

# ==============================================================================
# Section 1: Water mass classification (same as S11)
# ==============================================================================

classify_watermass <- function(temp, sal) {
  wm <- rep(NA_character_, length(temp))
  for (i in seq_along(temp)) {
    s <- sal[i]; t <- temp[i]
    if (is.na(s) | is.na(t))                              { wm[i] <- "NA"; next }
    if (s < 30.8  & t >= 3)                               { wm[i] <- "wCW"
    } else if (s < 30.8  & t < 3)                         { wm[i] <- "IMW & cCW"
    } else if (s >= 30.8 & s < 32.5 & t >= 0 & t < 3)    { wm[i] <- "cSW"
    } else if (s >= 30.8 & s < 33.4 & t >= 3)             { wm[i] <- "wSW"
    } else if (s >= 32.5 & s < 33.8 & t >= 0 & t < 3)    { wm[i] <- "AnW"
    } else if (s >= 30.8 & s < 33.8 & t >= -1 & t < 0)   { wm[i] <- "MWW"
    } else if (s >= 30.8 & s < 33.8 & t < -1)             { wm[i] <- "WW"
    } else if ((s >= 33.4 & t >= 3) | (s >= 33.8 & t < 3)) { wm[i] <- "AtlW & BBW"
    } else { wm[i] <- "Unclassified" }
  }
  wm
}

# ==============================================================================
# Section 2: Identify wSW samples with NO3 >= 10 uM
# ==============================================================================

df <- meta_asgard
df$cluster <- clusnum10[rownames(df)]
df$water_mass <- classify_watermass(df$temp, df$salinity)

wsw_high_no3 <- df %>%
  filter(water_mass == "wSW", `NO3(uM)` >= 10) %>%
  select(station, depth_type, depth_m, temp, salinity, `NO3(uM)`, `Sil(uM)`,
         `PO4(uM)`, cluster, water_mass) %>%
  arrange(desc(`NO3(uM)`))

message("wSW-classified samples with NO3 >= 10 uM: ", nrow(wsw_high_no3))
message("  out of total wSW samples: ", sum(df$water_mass == "wSW", na.rm = TRUE))

# ==============================================================================
# Section 3: Summary by cluster
# ==============================================================================

message("\nCluster distribution of misclassified samples:")
print(table(wsw_high_no3$cluster))

message("\nTemperature range: ", round(min(wsw_high_no3$temp), 2), " - ",
        round(max(wsw_high_no3$temp), 2), " C")
message("  These samples are just above the 3C AnW/wSW boundary")
message("  If T < 3C, they would be classified as AnW")

# ==============================================================================
# Section 4: Save results
# ==============================================================================

dir.create(here("output", "survey", "crosstable"), showWarnings = FALSE, recursive = TRUE)
write.csv(wsw_high_no3, here("output", "survey", "crosstable", "wSW_high_NO3_samples.csv"),
          row.names = TRUE)

message("\nS15_wSW_misclassification.R: done.")
message("  CSV: output/survey/crosstable/wSW_high_NO3_samples.csv")
message("  ", nrow(wsw_high_no3), " samples identified as potentially misclassified")
message("  Consistent with Hirawake et al. (2021) PiO 197:102641")
