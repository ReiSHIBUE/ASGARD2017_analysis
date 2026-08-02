### P14_station_region_summary.R
### ASGARD 2017 Processing — Station / region / DBO3 summaries
### 海域(Bering/Chukchi)・ステーション・DBO3 の集計と検定
###
### REQUIRES (from 00_setup.R, P01, P03):
###   asgard_filtered_p_hm2 - ASV count matrix 78 x 221
###   meta_asgard_p2        - metadata 78 samples (lat, station, depth_type)
###   clusnum_p             - 4-cluster assignment
###   shorternames          - ESV short labels
###   fullnameboot          - taxonomy + bootstrap df
###
### OUTPUT (output_p/cluster_summary/):
###   processing_region_depth_counts.csv         - region x depth counts
###   processing_cluster_region_depth_chisq.csv  - cluster x depth / region chi-sq
###   processing_filter_depth_counts.csv          - filter (size fraction) x depth counts
###   processing_filter_depth_chisq.csv           - filter x depth chi-sq / Fisher
###   processing_station_classRA.csv             - per-station class RA
###   processing_stationgroup_classRA.csv        - per station-group class RA
###   processing_DBO3_ASV_abundance.csv          - DBO3 per-ASV abundance

library(dplyr)
library(tidyr)
library(here)

dir.create(here("output_p", "cluster_summary"),
           showWarnings = FALSE, recursive = TRUE)

REGION_CUTOFF <- 65.7   # deg N: south = Bering, north = Chukchi

strip_boot <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)

# --- Class mapping per ASV ----------------------------------------------------
m   <- asgard_filtered_p_hm2
idx <- match(colnames(m), shorternames)
asv_class  <- strip_boot(fullnameboot$Class[idx]);  names(asv_class)  <- colnames(m)
asv_family <- strip_boot(fullnameboot$Family[idx]); names(asv_family) <- colnames(m)
asv_genus  <- strip_boot(fullnameboot$Genus[idx]);  names(asv_genus)  <- colnames(m)

RA <- m / rowSums(m) * 100   # per-sample relative abundance (%)

meta_p  <- meta_asgard_p2[rownames(m), ]
station <- meta_p$station
region  <- ifelse(is.na(meta_p$lat), "Unknown",
                  ifelse(meta_p$lat < REGION_CUTOFF, "Bering", "Chukchi"))
depth   <- as.character(meta_p$depth_type)
depth[is.na(depth) | depth == ""] <- "Unknown"
cl      <- as.character(clusnum_p[rownames(m)])

# ==============================================================================
# Section 1: region x depth counts (all 78 samples)
# ==============================================================================
region_f <- factor(region, levels = c("Bering", "Chukchi", "Unknown"))
depth_f  <- factor(depth,  levels = c("surf", "mid", "bottom", "Unknown"))
rd_tab   <- table(region = region_f, depth = depth_f)
write.csv(as.data.frame.matrix(addmargins(rd_tab)),
          here("output_p", "cluster_summary", "processing_region_depth_counts.csv"),
          row.names = TRUE)

# ==============================================================================
# Section 2: cluster x depth and cluster x region chi-square / Fisher
# (valid samples only; excludes the coordinate-less blank BOX_6_26)
# ==============================================================================
keep <- region != "Unknown" & depth != "Unknown"
cl_k <- factor(cl[keep], levels = c("1", "2", "3", "4"))

t_depth  <- table(cluster = cl_k, depth  = factor(depth[keep],  levels = c("surf","mid","bottom")))
t_region <- table(cluster = cl_k, region = factor(region[keep], levels = c("Bering","Chukchi")))

chi_d <- suppressWarnings(chisq.test(t_depth));  fis_d <- fisher.test(t_depth)
chi_r <- suppressWarnings(chisq.test(t_region)); fis_r <- fisher.test(t_region)

chisq_df <- data.frame(
  test             = c("cluster x depth", "cluster x region"),
  n                = c(sum(t_depth), sum(t_region)),
  chisq            = round(c(chi_d$statistic, chi_r$statistic), 3),
  df               = c(chi_d$parameter, chi_r$parameter),
  chisq_p          = round(c(chi_d$p.value, chi_r$p.value), 4),
  fisher_p         = round(c(fis_d$p.value, fis_r$p.value), 4),
  any_expected_lt5 = c(any(chi_d$expected < 5), any(chi_r$expected < 5))
)
write.csv(chisq_df,
          here("output_p", "cluster_summary", "processing_cluster_region_depth_chisq.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 2b: filter (size fraction) x depth chi-square / Fisher
# Do the 0.2 / 3 / 20 um fractions differ in their depth distribution?
# ==============================================================================
sig_star <- function(p) ifelse(p <= 0.001, "***",
                        ifelse(p <= 0.01,  "**",
                        ifelse(p <= 0.05,  "*", "ns")))

filt_raw <- as.character(meta_p$filter)   # "0.2 um", "3 um", "20 um" (with unit suffix)
filt_lab <- ifelse(grepl("^0\\.2", filt_raw), "0.2um",
             ifelse(grepl("^3",    filt_raw), "3um",
             ifelse(grepl("^20",   filt_raw), "20um", "Unknown")))

keep_fd  <- depth != "Unknown" & filt_lab != "Unknown"
filt_k   <- factor(filt_lab[keep_fd], levels = c("0.2um", "3um", "20um"))
depth_fd <- factor(depth[keep_fd],    levels = c("surf", "mid", "bottom"))
t_fd     <- table(filter = filt_k, depth = depth_fd)

# contingency table (counts)
write.csv(as.data.frame.matrix(addmargins(t_fd)),
          here("output_p", "cluster_summary", "processing_filter_depth_counts.csv"),
          row.names = TRUE)

chi_fd <- suppressWarnings(chisq.test(t_fd))
fis_fd <- fisher.test(t_fd)

filter_depth_df <- data.frame(
  test             = "filter x depth",
  n                = sum(t_fd),
  chisq            = round(unname(chi_fd$statistic), 3),
  df               = unname(chi_fd$parameter),
  chisq_p          = round(chi_fd$p.value, 4),
  fisher_p         = round(fis_fd$p.value, 4),
  sig              = sig_star(fis_fd$p.value),
  any_expected_lt5 = any(chi_fd$expected < 5),
  row.names        = NULL
)
write.csv(filter_depth_df,
          here("output_p", "cluster_summary", "processing_filter_depth_chisq.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 3: per-class RA per station and per station group
# ==============================================================================
targets <- c("Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia")
classRA <- sapply(targets, function(cls_i)
  rowSums(RA[, names(asv_class)[asv_class == cls_i], drop = FALSE]))

base_df <- data.frame(station = station, classRA, check.names = FALSE) %>%
  filter(!is.na(station), station != "")

# per station
station_df <- base_df %>%
  group_by(station) %>%
  summarise(n = n(),
            Alpha       = round(mean(Alphaproteobacteria), 1),
            Gamma       = round(mean(Gammaproteobacteria), 1),
            Bacteroidia = round(mean(Bacteroidia), 1),
            .groups = "drop") %>%
  mutate(Sum3 = Alpha + Gamma + Bacteroidia) %>%
  arrange(desc(n))
write.csv(station_df,
          here("output_p", "cluster_summary", "processing_station_classRA.csv"),
          row.names = FALSE)

# station groups (combine sub-stations, e.g. DBO3.3 + DBO3.8A -> DBO3)
grp_order <- c("DBO3", "CBW", "CL", "IL", "DBO2", "CBE", "CNL", "Incubation")
grp <- dplyr::case_when(
  grepl("^DBO3", base_df$station)                    ~ "DBO3",
  grepl("^DBO2", base_df$station)                    ~ "DBO2",
  grepl("^CBW",  base_df$station)                    ~ "CBW",
  grepl("^CBE",  base_df$station)                    ~ "CBE",
  grepl("^CNL",  base_df$station)                    ~ "CNL",
  grepl("^CL",   base_df$station)                    ~ "CL",
  grepl("^IL",   base_df$station)                    ~ "IL",
  grepl("incubation", base_df$station, ignore.case = TRUE) ~ "Incubation",
  TRUE ~ NA_character_)

group_df <- base_df %>%
  mutate(group = grp) %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(n_samples  = n(),
            n_stations  = dplyr::n_distinct(station),
            Alpha       = round(mean(Alphaproteobacteria), 1),
            Gamma       = round(mean(Gammaproteobacteria), 1),
            Bacteroidia = round(mean(Bacteroidia), 1),
            .groups = "drop") %>%
  mutate(Sum3 = Alpha + Gamma + Bacteroidia,
         group = factor(group, levels = grp_order)) %>%
  arrange(group)
write.csv(group_df,
          here("output_p", "cluster_summary", "processing_stationgroup_classRA.csv"),
          row.names = FALSE)

# ==============================================================================
# Section 4: DBO3 per-ASV abundance (family-genus level)
# ==============================================================================
is_dbo3 <- grepl("^DBO3", station)
dbo3_df <- data.frame(
  ASV         = sub(".*;\\s*", "", colnames(m)),
  Class       = asv_class,
  Family      = asv_family,
  Genus       = asv_genus,
  mean_RA_pct = round(colMeans(RA[is_dbo3, , drop = FALSE]), 3),
  row.names   = NULL
) %>% arrange(desc(mean_RA_pct))
write.csv(dbo3_df,
          here("output_p", "cluster_summary", "processing_DBO3_ASV_abundance.csv"),
          row.names = FALSE)

message("\nP14_station_region_summary.R: done.")
message("  region x depth, cluster chi-sq, station/group class RA, DBO3 ASV abundance written")
message(sprintf("  Bering/Chukchi (cutoff %.1f N): %d / %d  (Unknown %d)",
                REGION_CUTOFF, sum(region == "Bering"),
                sum(region == "Chukchi"), sum(region == "Unknown")))
message(sprintf("  cluster x depth:  Fisher p = %.4f | cluster x region: Fisher p = %.4f",
                fis_d$p.value, fis_r$p.value))
message(sprintf("  filter x depth:   Fisher p = %.4f %s (chisq p = %.4f)",
                fis_fd$p.value, sig_star(fis_fd$p.value), chi_fd$p.value))
