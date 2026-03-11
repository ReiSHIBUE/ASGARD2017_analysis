### S02_heatmaps_16S.R
### ASGARD 2017 Survey Site Analysis — 16S Heatmaps
### 16S ヒートマップスクリプト（サーベイサイト）
###
### REQUIRES (from S01_data_prep.R):
###   asgard_filtered  - 181×258 ASV proportion matrix
###   asgard_frtprop   - 181×258 fourth-root transformed matrix
###   meta_asgard      - metadata for 181 samples
###
### PRODUCES (consumed by S03–S10):
###   h1               - base heatmap object (provides row/col dendrograms)
###   h2               - heatmap.2 with row and col cluster side colours
###   clusnum          - sample cluster assignments, length 181, values 1–5
###   colclusnum       - ASV column cluster assignments, length 258, values 1–8
###   rsc              - 5-colour palette vector indexed by clusnum
###   colrsc           - 8-colour palette vector indexed by colclusnum
###
### OUTPUT:
###   output/survey/heatmaps/ASGARD_hm_survey_16S.pdf

library(gplots)
library(viridis)
library(scales)
library(vegan)

# ==============================================================================
# Section 1: ヒートマップ用行列の準備 / Prepare matrix for heatmap
# ==============================================================================

asgard_frtmat <- asgard_frtprop
asgard_frtmat <- asgard_frtmat[, colSums(asgard_frtmat) > 0] # 181×258

# ==============================================================================
# Section 2: h1 — デンドログラム抽出用ベースheatmap
# Base heatmap — used only to extract row and column dendrograms
# ==============================================================================

pdf(file = here::here("output", "survey", "heatmaps", "ASGARD_hm_survey_16S.pdf"),
    width = 20, height = 20)

h1 <- heatmap.2(
  asgard_frtmat,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  margins   = c(15, 15),
  scale     = "none",
  main      = "ASGARD Survey 16S — Bray/ward.D",
  trace     = "none",
  cexCol    = 0.2,
  cexRow    = 0.3
)

# ==============================================================================
# Section 3: クラスター割り当て / Extract cluster assignments
# 行: 5クラスター, 列: 8クラスター
# Rows: 5 clusters (samples), Cols: 8 clusters (ASVs)
# ==============================================================================

# 行クラスター (サンプル) / Row clusters (samples)
nclus      <- 5
oldclus    <- cutree(as.hclust(h1$rowDendrogram), k = nclus)
oldorder   <- unname(rle(oldclus[as.hclust(h1$rowDendrogram)$order])$values)
neworder   <- (1:nclus)
names(neworder) <- oldorder
clusnum    <- unname(neworder[as.character(oldclus)])
names(clusnum) <- names(oldclus) # length 181, values 1–5

# 列クラスター (ASVs) / Column clusters (ASVs)
colnclus      <- 8
cololdclus    <- cutree(as.hclust(h1$colDendrogram), k = colnclus)
cololdorder   <- unname(rle(cololdclus[as.hclust(h1$colDendrogram)$order])$values)
colneworder   <- (1:colnclus)
names(colneworder) <- cololdorder
colclusnum    <- unname(colneworder[as.character(cololdclus)])
names(colclusnum) <- names(cololdclus) # length 258, values 1–8

# カラーパレット / Colour palettes
rsc    <- hue_pal()(nclus)[clusnum]      # row side colours, length 181
colrsc <- plasma(colnclus)[colclusnum]   # col side colours, length 258

# ==============================================================================
# Section 4: h2 — クラスター色付きheatmap.2 / heatmap.2 with cluster side colours
# ==============================================================================

h2 <- heatmap.2(
  asgard_frtmat,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = rsc,
  ColSideColors = colrsc,
  Rowv          = h1$rowDendrogram,
  Colv          = h1$colDendrogram,
  margins       = c(15, 15),
  scale         = "none",
  main          = "ASGARD Survey 16S — Bray/ward.D (5 row clusters, 8 col clusters)",
  trace         = "none",
  cexCol        = 0.2,
  cexRow        = 0.3
)

# ==============================================================================
# Section 5: h5 — 列クラスター8のサブセット (particle-associated candidates)
# Subset heatmap for column cluster 8 — these ASVs have an interesting sub-structure
# ==============================================================================

cc8_names <- names(colclusnum)[colclusnum == 8]
asgard_frtmat_cc8 <- asgard_frtmat[, cc8_names, drop = FALSE]
asgard_frtmat_cc8 <- asgard_frtmat_cc8[
  rowSums(asgard_frtmat_cc8) > 0,
  colSums(asgard_frtmat_cc8) > 0,
  drop = FALSE
]

# 列クラスター8内のサブクラスター抽出 / Sub-cluster within col cluster 8
h5 <- heatmap.2(
  asgard_frtmat_cc8,
  distfun       = function(x) vegdist(x, method = "euclidean"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = rsc[rownames(asgard_frtmat_cc8)],
  Rowv          = h2$rowDendrogram,
  margins       = c(15, 15),
  scale         = "none",
  main          = "ASGARD Survey 16S — col cluster 8 subset",
  trace         = "none",
  cexCol        = 0.5,
  cexRow        = 0.3,
  labRow        = meta_asgard[rownames(asgard_frtmat_cc8), "depth_type"]
)

dev.off()

message("S02_heatmaps_16S.R: done.")
message("  clusnum length:    ", length(clusnum), " (", nclus, " row clusters)")
message("  colclusnum length: ", length(colclusnum), " (", colnclus, " col clusters)")
message("  PDF: output/survey/heatmaps/ASGARD_hm_survey_16S.pdf")
