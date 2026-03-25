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

# カラーパレット (5クラスター) / Colour palettes (5 clusters)
rsc    <- hue_pal()(nclus)[clusnum]      # row side colours, length 181
colrsc <- plasma(colnclus)[colclusnum]   # col side colours, length 258

# ==============================================================================
# Section 4: 10クラスターへの分割 / Split into 10 clusters
# Cluster 1→2分割, 2→3分割, 3→2分割, 5→2分割, 4はそのまま
# ==============================================================================

c1_sub <- cutree(hclust(vegdist(asgard_frtprop[names(clusnum)[clusnum==1], colSums(asgard_frtprop[names(clusnum)[clusnum==1],])>0], "bray"), "ward.D"), k=2)
c2_sub <- cutree(hclust(vegdist(asgard_frtprop[names(clusnum)[clusnum==2], colSums(asgard_frtprop[names(clusnum)[clusnum==2],])>0], "bray"), "ward.D"), k=3)
c3_sub <- cutree(hclust(vegdist(asgard_frtprop[names(clusnum)[clusnum==3], colSums(asgard_frtprop[names(clusnum)[clusnum==3],])>0], "bray"), "ward.D"), k=2)
c5_sub <- cutree(hclust(vegdist(asgard_frtprop[names(clusnum)[clusnum==5], colSums(asgard_frtprop[names(clusnum)[clusnum==5],])>0], "bray"), "ward.D"), k=2)

# 仮ラベルで割り当て / Assign temporary labels
tmp10 <- rep(NA_character_, length(clusnum)); names(tmp10) <- names(clusnum)
tmp10[names(c1_sub)[c1_sub==1]] <- "o1a"; tmp10[names(c1_sub)[c1_sub==2]] <- "o1b"
tmp10[names(c2_sub)[c2_sub==1]] <- "o2a"; tmp10[names(c2_sub)[c2_sub==2]] <- "o2b"; tmp10[names(c2_sub)[c2_sub==3]] <- "o2c"
tmp10[names(c3_sub)[c3_sub==1]] <- "o3a"; tmp10[names(c3_sub)[c3_sub==2]] <- "o3b"
tmp10[names(clusnum)[clusnum==4]] <- "o4"
tmp10[names(c5_sub)[c5_sub==1]] <- "o5a"; tmp10[names(c5_sub)[c5_sub==2]] <- "o5b"

# デンドログラム順（下→上）で1-10を振り直す / Renumber by dendrogram order (bottom→top)
dend_order <- as.hclust(h1$rowDendrogram)$order
sample_names_ordered <- rownames(asgard_frtprop)[dend_order]
appearance_order <- unique(tmp10[sample_names_ordered])
remap10 <- setNames(as.character(1:10), appearance_order)

clusnum10 <- remap10[tmp10]; names(clusnum10) <- names(tmp10)

# 10クラスターカラーパレット / 10-cluster colour palette
cc10 <- hue_pal()(10); names(cc10) <- as.character(1:10)
rsc10 <- cc10[clusnum10]; names(rsc10) <- names(clusnum10)

# ==============================================================================
# Section 5: h2 — 10クラスター色付きheatmap.2
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

# h2b: station名ラベル付き / With station name labels
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
message("  clusnum length:    ", length(clusnum), " (", nclus, " row clusters)")
message("  clusnum10 length:  ", length(clusnum10), " (10 row clusters)")
message("  colclusnum length: ", length(colclusnum), " (", colnclus, " col clusters)")
message("  PDF: output/survey/heatmaps/ASGARD_hm_survey_16S.pdf")
