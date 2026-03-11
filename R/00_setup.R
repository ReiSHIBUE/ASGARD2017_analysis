### 00_setup.R
### ASGARD 2017 Processing Site Analysis — Environment Setup
### 環境セットアップスクリプト
###
### Loads raw RDS files and builds all upstream objects required by 01_data_prep.R.
### このスクリプトはRDSファイルを読み込み、01_data_prep.R が必要とする
### 上流オブジェクトを全て構築します。
###
### DATA FILES (from data/raw/):
###   seqtab_filt.rds      - full ASV count table (2193 samples × 18535 ASVs)
###   table_list.rds       - primer/marker assignment per ASV column
###   meta_denovo_2.RDS    - sample metadata (1163 samples × 45 cols)
###   names_list.rds       - full taxonomic name strings (length 18535)
###   bootout_edit.rds     - bootstrap confidence values per taxonomic rank
###
### PRODUCES (consumed by 01_data_prep.R):
###   seqtab_filt          - raw count matrix (2193 × 18535)
###   seqtab_16Smat        - 16S prokaryote counts, ≥1000-read samples (1633 × 3076)
###   seqtab_16Sprop       - proportion table from seqtab_16Smat (1633 × 3076)
###   seqtab_16Spropcol_mra - mean relative abundance per ASV (length 3076)
###   meta_denovo_2        - metadata data frame (1163 × 45)
###   shorternames         - "Order; Family; Genus; ESV_N" label vector (length 18535)

library(here)

# ==============================================================================
# Section 1: RDSファイルの読み込み / Load raw RDS files
# ==============================================================================

message("Loading raw data files...")

seqtab_filt  <- readRDS(here("data", "raw", "seqtab_filt.rds"))
# 2193サンプル × 18535 ASVs (全マーカー含む)

table_list   <- readRDS(here("data", "raw", "table_list.rds"))
# 各ASV列がどのマーカー(16S_prokaryote, 18S等)に属するかを示すベクトル
table_list[is.na(table_list)] <- "unknown"

meta_denovo_2 <- readRDS(here("data", "raw", "meta_denovo_2.RDS"))
# サンプルメタデータ (station, depth, temp, sal, nutrients, etc.)

names_list   <- readRDS(here("data", "raw", "names_list.rds"))
# 分類名を含む名前リスト / Taxonomic name list (named list of character vectors)

bootout      <- readRDS(here("data", "raw", "bootout_edit.rds"))
# 各分類階層のブートストラップ信頼値 / Bootstrap confidence per rank

# ==============================================================================
# Section 2: 16S prokaryote ASVのみ抽出 / Subset to 16S prokaryote columns
# ==============================================================================

seqtab_prok <- table_list == "16S_prokaryote" # logical vector, length 18535
seqtab_mat  <- seqtab_filt[, seqtab_prok]     # 2193 × 3076

# ==============================================================================
# Section 3: リード数 ≥1000 のサンプルのみ残す / Keep samples with ≥1000 reads
# ==============================================================================

seqtab_matrow <- rowSums(seqtab_mat)           # total 16S reads per sample
seqtab_16S    <- seqtab_matrow >= 1000          # logical vector
seqtab_16Smat <- seqtab_mat[seqtab_16S, ]      # 1633 × 3076

# ==============================================================================
# Section 4: プロポーション表と平均相対存在量 / Proportion table & mean rel. abund.
# ==============================================================================

seqtab_16Sprop        <- prop.table(seqtab_16Smat, 1)            # 1633 × 3076
seqtab_16Spropcol     <- colSums(seqtab_16Sprop)                  # length 3076
seqtab_16Spropcol_mra <- seqtab_16Spropcol / nrow(seqtab_16Sprop) # mean rel. abund.

# ==============================================================================
# Section 5: shorternames の構築 / Build shorternames label vector
# Order; Family; Genus; ESV_N 形式のラベルを各ASVに付ける
# ==============================================================================

# 分類名の長い文字列を;で分割して行列にする / Split taxonomy strings into matrix
longernames  <- names_list$ref_dada2_silva_v132.fasta # length 18535
fullnamemat  <- do.call(rbind, strsplit(longernames, split = ";"))
colnames(fullnamemat) <- c(
  "root", "Domain", "Major Clade", "Superkingdom", "Kingdom",
  "Subkingdom", "Infrakingdom", "Superphylum", "Phylum", "Subphylum",
  "Infraphylum", "Superclass", "Class", "Subclass", "Infraclass",
  "Superorder", "Order", "Suborder", "Superfamily", "Family",
  "Subfamily", "Genus", "Accession", "ESV"
)

# ブートストラップ値と分類名を結合 / Combine taxonomy + bootstrap values
bosilva <- bootout$ref_dada2_silva_v132.fasta  # 18535 × 23 matrix
colnames(bosilva) <- c(
  "root", "Domain", "Major Clade", "Superkingdom", "Kingdom",
  "Subkingdom", "Infrakingdom", "Superphylum", "Phylum", "Subphylum",
  "Infraphylum", "Superclass", "Class", "Subclass", "Infraclass",
  "Superorder", "Order", "Suborder", "Superfamily", "Family",
  "Subfamily", "Genus", "Accession"
)

# "TaxonName (bootstrap)" 形式のdfを作る / Build "TaxonName (bootstrap)" df
fullnameboot <- bosilva
for (i in 1:nrow(bosilva)) {
  for (j in 1:ncol(bosilva)) {
    fullnameboot[i, j] <- paste0(fullnamemat[i, j], " (", bosilva[i, j], ")")
  }
}
fullnameboot       <- as.data.frame(fullnameboot)
fullnameboot$ESV   <- fullnamemat[, "ESV"]  # ESV列を追加 / Append ESV column

# "Order; Family; Genus; ESV_N" ラベルベクトル (length 18535)
shorternames <- paste(
  fullnameboot$Order,
  fullnameboot$Family,
  fullnameboot$Genus,
  fullnameboot$ESV,
  sep = "; "
)

message("00_setup.R: done.")
message("  seqtab_16Smat:        ", nrow(seqtab_16Smat), " x ", ncol(seqtab_16Smat))
message("  seqtab_16Sprop:       ", nrow(seqtab_16Sprop), " x ", ncol(seqtab_16Sprop))
message("  meta_denovo_2:        ", nrow(meta_denovo_2), " x ", ncol(meta_denovo_2))
message("  shorternames length:  ", length(shorternames))
