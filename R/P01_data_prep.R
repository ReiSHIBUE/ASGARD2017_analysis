### 01_data_prep.R
### ASGARD 2017 Processing Site Analysis — Data Preparation
### データ準備スクリプト
###
### REQUIRES (must exist in R session before sourcing):
###   seqtab_16Sprop        - proportion table (samples × ASVs)
###   seqtab_16Smat         - count table (samples × ASVs)
###   seqtab_16Spropcol_mra - mean relative abundance vector per ASV
###   seqtab_filt           - raw filtered seqtab (for colnames / shorternames indexing)
###   shorternames          - short ESV label vector aligned to seqtab_filt columns
###   meta_denovo_2         - metadata table (samples × metadata cols)
###
### PRODUCES (used by later scripts):
###   asgard_processing        - 81-sample merged df (metadata + ASVs, 81 * 269)
###   asgard_processing2       - 78-sample merged df after ≥5000-read filter (78 * 267)
###   asgard_filtered_p        - filtered ASV proportion matrix, 81 samples
###   asgard_filtered_p_hm2    - matrix form of asgard_filtered_p2 (78 * 221)
###   meta_asgard_p            - metadata for 81 ASGARD processing samples
###   meta_asgard_p2           - metadata for 78 samples after read filter
###   asgard_vec_p             - sample row names, length 81
###   asgard_vec_p2            - sample row names, length 78
###   seqtab_asgard_16Srow     - row sums for 81 samples
###   seqtab_asgard_16Srow_2   - row sums for 78 samples

# ライブラリの読み込み / Load libraries
library(tidyverse)
library(janitor)
library(here)
library(gplots)   # heatmap.2()
library(viridis)
library(scales)   # rescale()
library(vegan)

# ==============================================================================
# Section 1: 全サンプル含むdfを作る / Build full-sample data frame
# ==============================================================================

# ASV列を平均相対存在量（降順）で並べ替える / Sort ASVs by mean relative abundance
seqtab_16Spropcol_mradec <- sort(seqtab_16Spropcol_mra, decreasing = TRUE)

seqtab_16Spropcol_3076 <- seqtab_16Sprop[, names(seqtab_16Spropcol_mradec)] # 1633*3076
# top100を抽出するmradec2ではなく、3076サンプルを指定するmradecを使う

ODV.otu.pick_full <- match(colnames(seqtab_16Spropcol_3076), colnames(seqtab_filt))
colnames(seqtab_16Spropcol_3076) <- shorternames[ODV.otu.pick_full] # colnames()= で列名変更

seqtab_ODV_full <- merge(meta_denovo_2, seqtab_16Spropcol_3076, by = 0, all.x = FALSE)
# 960*3122: 3122 = 3076 + 46 (meta_denovo_2は1163*45, Row.namesが含まれて+1)

# Processing station のみ抽出 / Keep only processing stations
seqtab_ODV_Processing_full <- seqtab_ODV_full[seqtab_ODV_full$station_type == "P", ]
# 241*3122  (S=689, P=241, E=30)

# ASGARD 2017 のみ抽出 / Filter for ASGARD2017 project
ASGARD_p <- seqtab_ODV_Processing_full %>%
  filter(project == "ASGARD2017") # 81*3122

# ==============================================================================
# Section 2: 81-sample subset
# サンプルIDベクトルとmetadataを取り出す / Extract sample IDs and metadata
# ==============================================================================

asgard_vec_p <- ASGARD_p$Row.names # length is 81
seqtab_prop_asgard_p <- seqtab_16Sprop[asgard_vec_p, ] # 81*3076

# ASGARDのmeta data抽出 / Extract ASGARD metadata
meta_asgard_p <- meta_denovo_2[asgard_vec_p, ] # 81*45

# 最小リードカットオフ適用 / Apply minimum abundance cutoff
mincutoff_p <- (apply(seqtab_prop_asgard_p, 2, max) > 0.001)
asgard_filtered_p <- seqtab_prop_asgard_p[, mincutoff_p]
asgard_filtered_p <- asgard_filtered_p[, colSums(asgard_filtered_p > 0) > 2] # 81*223

# ESV短縮名を列名に適用 / Assign shorter ESV names
ODV.otu.pick_p <- match(colnames(asgard_filtered_p), colnames(seqtab_filt))
colnames(asgard_filtered_p) <- shorternames[ODV.otu.pick_p]

# データフレーム化してmetadataと結合 / Convert to df and join with metadata
asgard_filtered_p <- as.data.frame(asgard_filtered_p)
asgard_filtered_p$Sample <- rownames(asgard_filtered_p)

meta_asgard_p$Sample <- rownames(meta_asgard_p)

# フィルターサイズに " µm" を付ける / Append " µm" to filter size labels
meta_asgard_p <- meta_asgard_p %>%
  mutate(filter = paste0(filter, " µm"))

asgard_processing <- left_join(meta_asgard_p, asgard_filtered_p, by = "Sample") # 81*269

# ==============================================================================
# Section 3: ≥5000リードフィルター後の78サンプルsubset
# 5000 reads minimum filter → 78-sample subset
# ==============================================================================

rownames(ASGARD_p) <- ASGARD_p$Row.names

common_asv <- intersect(rownames(ASGARD_p),
                        rownames(seqtab_16Smat)) # length is 81

asgard_seqcount <- seqtab_16Smat[common_asv, ] # 81*3076

seqtab_16Smatrow <- rowSums(seqtab_16Smat) # length is 1633

seqtab_16S_asgard_5000 <- seqtab_16Smatrow >= 5000 # logical vector, length 1633
seqtab_16Smat_asgard_5000 <- seqtab_16Smat[seqtab_16S_asgard_5000, ] # 1238*3076
seqtab_16Srow_asgard_5000 <- rowSums(seqtab_16Smat_asgard_5000)       # length 1238

common_asv2 <- intersect(rownames(ASGARD_p),
                         rownames(seqtab_16Smat_asgard_5000)) # length 78

asgard_seqcount2 <- seqtab_16Smat[common_asv2, ] # 78*3076
# プロポーション再計算 / Recompute proportions for 78-sample set
asgard_seqcount_16Sprop <- prop.table(asgard_seqcount2, 1) # 78*3076
asgard_vec_p2 <- rownames(asgard_seqcount2) # length 78
seqtab_prop_asgard_p2 <- asgard_seqcount_16Sprop[asgard_vec_p2, ] # 78*3076

mincutoff_p2 <- (apply(seqtab_prop_asgard_p2, 2, max) > 0.001)
asgard_filtered_p2 <- seqtab_prop_asgard_p2[, mincutoff_p2]
asgard_filtered_p2 <- asgard_filtered_p2[, colSums(asgard_filtered_p2 > 0) > 2] # 78*221

ODV.otu.pick_p2 <- match(colnames(asgard_filtered_p2), colnames(seqtab_filt))
colnames(asgard_filtered_p2) <- shorternames[ODV.otu.pick_p2]

asgard_filtered_p2 <- as.data.frame(asgard_filtered_p2)
asgard_filtered_p2$Sample <- rownames(asgard_filtered_p2) # 78*222
vec78 <- rownames(asgard_filtered_p2)

meta_asgard_p2 <- meta_asgard_p[vec78, ] # 78*46

asgard_processing2 <- left_join(meta_asgard_p2, asgard_filtered_p2, by = "Sample") # 78*267

# ==============================================================================
# Section 4: asgard_filtered_p_hm2 を作成する
# Create matrix for heatmap input
# ==============================================================================

asgard_filtered_p_hm2 <- asgard_filtered_p2 %>%
  select(-Sample) # 78*221

asgard_filtered_p_hm2 <- as.matrix(asgard_filtered_p_hm2) # 78*221 matrix

message("01_data_prep.R: done. Objects ready: asgard_processing (81), asgard_processing2 (78), asgard_filtered_p_hm2 (matrix 78x221).")
