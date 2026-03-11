### S01_data_prep.R
### ASGARD 2017 Survey Site Analysis — Data Preparation
### データ準備スクリプト（サーベイサイト）
###
### REQUIRES (from 00_setup.R):
###   seqtab_filt           - raw count matrix (2193 × 18535)
###   seqtab_16Smat         - 16S prokaryote counts, ≥1000-read samples (1633 × 3076)
###   seqtab_16Sprop        - proportion table (1633 × 3076)
###   seqtab_16Spropcol_mra - mean relative abundance per ASV (length 3076)
###   meta_denovo_2         - metadata (1163 × 45)
###   shorternames          - "Order; Family; Genus; ESV_N" label vector (length 18535)
###
### PRODUCES (consumed by S02–S10):
###   asgard_filtered       - 181-sample ASV proportion matrix (181 × 258)
###   meta_asgard           - metadata for 181 survey samples (181 × 45)
###   asgard_seqcount       - integer count matrix for 181 samples (181 × 3076)
###   asgard_vec            - sample row names, length 181
###   asgard_frtprop        - fourth-root transformed proportion matrix (181 × 258)
###
### Survey samples are 0.2 µm (free-living) only because station_type=="S"
### filters to survey stations which only have 0.2 µm samples.
### サーベイサイトは0.2 µmフィルターのみ（自由生活型微生物）。

library(tidyverse)
library(janitor)
library(here)
library(vegan)

# ==============================================================================
# Section 1: 全サンプル含むdfを作る / Build full-sample data frame
# ==============================================================================

seqtab_16Spropcol_mradec <- sort(seqtab_16Spropcol_mra, decreasing = TRUE)
seqtab_16Spropcol_3076   <- seqtab_16Sprop[, names(seqtab_16Spropcol_mradec)] # 1633×3076

ODV.otu.pick_full        <- match(colnames(seqtab_16Spropcol_3076), colnames(seqtab_filt))
colnames(seqtab_16Spropcol_3076) <- shorternames[ODV.otu.pick_full]

seqtab_ODV_full   <- merge(meta_denovo_2, seqtab_16Spropcol_3076, by = 0, all.x = FALSE)
# 960×3122

# Survey stationsのみ抽出 / Keep survey stations (0.2 µm only)
seqtab_ODV_Survey_full <- seqtab_ODV_full[seqtab_ODV_full$station_type == "S", ]
# 689×3122

# ASGARD2017のみ抽出 / Filter for ASGARD2017 project
ASGARD <- seqtab_ODV_Survey_full %>%
  filter(project == "ASGARD2017") # 181×3122

# ==============================================================================
# Section 2: 181-sample subset
# サンプルIDベクトルとmetadataを取り出す / Extract sample IDs and metadata
# ==============================================================================

asgard_vec      <- ASGARD$Row.names  # length 181
meta_asgard     <- meta_denovo_2[asgard_vec, ] # 181×45

seqtab_prop_asgard <- seqtab_16Sprop[asgard_vec, ] # 181×3076

# ==============================================================================
# Section 3: 最小存在量フィルター / Apply minimum abundance cutoff
# max relative abundance > 0.001 AND present in >2 samples
# ==============================================================================

mincutoff      <- apply(seqtab_prop_asgard, 2, max) > 0.001
asgard_filtered <- seqtab_prop_asgard[, mincutoff]
asgard_filtered <- asgard_filtered[, colSums(asgard_filtered > 0) > 2]
# dim: 181×258

# ESV短縮名を列名に適用 / Assign shorter ESV names to columns
ODV.otu.pick_s      <- match(colnames(asgard_filtered), colnames(seqtab_filt))
colnames(asgard_filtered) <- shorternames[ODV.otu.pick_s]

# ==============================================================================
# Section 4: カウント行列と4乗根変換 / Count matrix and fourth-root transform
# ==============================================================================

common_asv      <- intersect(asgard_vec, rownames(seqtab_16Smat))
asgard_seqcount <- seqtab_16Smat[common_asv, ] # 181×3076

# 4乗根変換 / Fourth-root transformation (reduces dominance of abundant ASVs)
asgard_frtprop  <- asgard_filtered^0.25  # 181×258

message("S01_data_prep.R: done.")
message("  asgard_filtered:  ", nrow(asgard_filtered), " x ", ncol(asgard_filtered))
message("  meta_asgard:      ", nrow(meta_asgard), " x ", ncol(meta_asgard))
message("  asgard_seqcount:  ", nrow(asgard_seqcount), " x ", ncol(asgard_seqcount))
message("  asgard_frtprop:   ", nrow(asgard_frtprop), " x ", ncol(asgard_frtprop))
