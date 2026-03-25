### S02_heatmaps_16S.R
### ASGARD 2017 Survey Site Analysis — 16S Heatmaps (10 clusters)
### 16S ヒートマップスクリプト（サーベイサイト、10クラスター）
###
### REQUIRES (from S01_data_prep.R):
###   asgard_filtered  - 181×258 ASV proportion matrix
###   asgard_frtprop   - 181×258 fourth-root transformed matrix
###   meta_asgard      - metadata for 181 samples
###
### PRODUCES (consumed by S03–S10):
###   h1               - base heatmap object (provides row/col dendrograms)
###   h2               - heatmap.2 with 10-cluster row and 8-cluster col side colours
###   clusnum10        - sample cluster assignments, length 181, values 1–10
###   cc10             - 10-colour palette vector
###   rsc10            - row side colours indexed by clusnum10
###   colclusnum       - ASV column cluster assignments, length 258, values 1–8
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
# 行: 10クラスター (cutree直接), 列: 8クラスター
# ==============================================================================

# 行クラスター (サンプル) — 10クラスター / Row clusters — 10 clusters
nclus10       <- 10
oldclus10     <- cutree(as.hclust(h1$rowDendrogram), k = nclus10)
oldorder10    <- unname(rle(oldclus10[as.hclust(h1$rowDendrogram)$order])$values)
neworder10    <- 1:nclus10
names(neworder10) <- oldorder10
clusnum10     <- unname(neworder10[as.character(oldclus10)])
names(clusnum10) <- names(oldclus10) # length 181, values 1–10

# 列クラスター (ASVs) / Column clusters (ASVs)
colnclus      <- 8
cololdclus    <- cutree(as.hclust(h1$colDendrogram), k = colnclus)
cololdorder   <- unname(rle(cololdclus[as.hclust(h1$colDendrogram)$order])$values)
colneworder   <- (1:colnclus)
names(colneworder) <- cololdorder
colclusnum    <- unname(colneworder[as.character(cololdclus)])
names(colclusnum) <- names(cololdclus) # length 258, values 1–8

# カラーパレット / Colour palettes
cc10   <- hue_pal()(nclus10)
names(cc10) <- as.character(1:nclus10)
rsc10  <- cc10[as.character(clusnum10)]
names(rsc10) <- names(clusnum10)

colrsc <- plasma(colnclus)[colclusnum]   # col side colours, length 258

# ==============================================================================
# Section 4: h2 — 10クラスター色付きheatmap.2
# heatmap.2 with 10-cluster RowSideColors
# ==============================================================================

h2 <- heatmap.2(
  asgard_frtmat,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = rsc10[rownames(asgard_frtmat)],
  ColSideColors = colrsc,
  Rowv          = h1$rowDendrogram,
  Colv          = h1$colDendrogram,
  margins       = c(15, 15),
  scale         = "none",
  main          = "ASGARD Survey 16S — Bray/ward.D (10 row clusters, 8 col clusters)",
  trace         = "none",
  cexCol        = 0.2,
  cexRow        = 0.3
)

# station名ラベル付き / With station name labels
heatmap.2(
  asgard_frtmat,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = rsc10[rownames(asgard_frtmat)],
  ColSideColors = colrsc,
  Rowv          = h1$rowDendrogram,
  Colv          = h1$colDendrogram,
  margins       = c(15, 15),
  scale         = "none",
  main          = "ASGARD Survey 16S — 10 clusters (station names)",
  trace         = "none",
  cexCol        = 0.2,
  cexRow        = 0.3,
  labRow        = meta_asgard[rownames(asgard_frtmat), "station"]
)

dev.off()

message("S02_heatmaps_16S.R: done.")
message("  clusnum10 length:  ", length(clusnum10), " (", nclus10, " row clusters)")
message("  colclusnum length: ", length(colclusnum), " (", colnclus, " col clusters)")
message("  PDF: output/survey/heatmaps/ASGARD_hm_survey_16S.pdf")
