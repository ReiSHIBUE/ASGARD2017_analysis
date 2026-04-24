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
###   hier_names       - named vector mapping "1"–"10" to A1,A2,B1,B2,B3,C1–C5
###   hier_levels      - character vector of cluster names in dendrogram order
###   cc10             - 10-colour palette vector (named by hier_levels)
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

# 列クラスター (ASVs) — 6クラスター / Column clusters — 6 clusters
colnclus      <- 6
cololdclus    <- cutree(as.hclust(h1$colDendrogram), k = colnclus)
cololdorder   <- unname(rle(cololdclus[as.hclust(h1$colDendrogram)$order])$values)
colneworder   <- (1:colnclus)
names(colneworder) <- cololdorder
colclusnum    <- unname(colneworder[as.character(cololdclus)])
names(colclusnum) <- names(cololdclus) # length 258, values 1–6

# 階層命名 / Hierarchical naming based on binary dendrogram splits
# k=3: A, B, C.  Then each subtree split recursively:
#   A → A1(9), A2(15)
#   B → B1(17), B2 → B2a(17), B2b(26)
#   C → C1 → C1a(14), C1b(30); C2 → C2a(12), C2b → C2b1(12), C2b2(29)
hier_names <- c("1"="A1",  "2"="A2",   "3"="B1",   "4"="B2a",  "5"="B2b",
                "6"="C1a", "7"="C1b",  "8"="C2a",  "9"="C2b1", "10"="C2b2")
hier_levels <- unname(hier_names)  # dendrogram order

# カラーパレット / Colour palettes (A=赤系, B=緑系, C=青系)
cc10 <- c(
  "A1"   = "#E31A1C",  "A2"   = "#FF7F00",
  "B1"   = "#33A02C",  "B2a"  = "#B2DF8A",  "B2b"  = "#A6D854",
  "C1a"  = "#1F78B4",  "C1b"  = "#6A3D9A",  "C2a"  = "#A6CEE3",
  "C2b1" = "#CAB2D6",  "C2b2" = "#B3B3B3"
)
rsc10  <- cc10[hier_names[as.character(clusnum10)]]
names(rsc10) <- names(clusnum10)

colclus_colors <- c("1"="#FB9A99", "2"="#08519C", "3"="#20B2AA",
                     "4"="#DAA520", "5"="#DD3497", "6"="#252525")
colrsc <- colclus_colors[as.character(colclusnum)]
names(colrsc) <- names(colclusnum)  # ASV名で名前付け

# ==============================================================================
# Section 4: h2 — 10行クラスター + 6列クラスター色付きheatmap.2
# heatmap.2 with 10-cluster RowSideColors + 6-cluster ColSideColors
# ==============================================================================

h2 <- heatmap.2(
  asgard_frtmat,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = rsc10[rownames(asgard_frtmat)],
  ColSideColors = colrsc[colnames(asgard_frtmat)],
  Rowv          = h1$rowDendrogram,
  Colv          = h1$colDendrogram,
  margins       = c(15, 15),
  scale         = "none",
  main          = "ASGARD Survey 16S — 10 row clusters, 6 col clusters",
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
  ColSideColors = colrsc[colnames(asgard_frtmat)],
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

# ラベルなし版 / No row/col labels
heatmap.2(
  asgard_frtmat,
  distfun       = function(x) vegdist(x, method = "bray"),
  hclustfun     = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = rsc10[rownames(asgard_frtmat)],
  ColSideColors = colrsc[colnames(asgard_frtmat)],
  Rowv          = h1$rowDendrogram,
  Colv          = h1$colDendrogram,
  margins       = c(2, 2),
  scale         = "none",
  main          = "ASGARD Survey 16S — 10 row clusters, 6 col clusters",
  trace         = "none",
  labRow        = FALSE,
  labCol        = FALSE
)

dev.off()

# ==============================================================================
# Section 5: 11クラスター版 (C1b → C1b1 + C1b2)
# C1b(n=30) is split into C1b1(n=11) and C1b2(n=19) based on
# re-clustering of C samples at k=6
# ==============================================================================

# Cサンプルのサブツリーでk=6カット
hc_full <- as.hclust(h1$rowDendrogram)
c_samples <- names(clusnum10)[clusnum10 >= 6]
d_c <- as.dist(as.matrix(cophenetic(hc_full))[c_samples, c_samples])
hc_c <- hclust(d_c, method = "ward.D")
cuts6_c <- cutree(hc_c, k = 6)
c_order <- hc_c$order
ord6_c <- unname(rle(cuts6_c[hc_c$labels[c_order]])$values)
remap6_c <- setNames(1:6, ord6_c)
cuts6_c <- remap6_c[as.character(cuts6_c)]
names(cuts6_c) <- names(cutree(hc_c, k = 6))

# 11クラスター割り当て
hier_levels_11 <- c("A1", "A2", "B1", "B2a", "B2b", "C1a", "C1b1", "C1b2", "C2a", "C2b1", "C2b2")

clusnum11 <- character(length(clusnum10))
names(clusnum11) <- names(clusnum10)
for (s in names(clusnum10)) {
  cn <- hier_names[as.character(clusnum10[s])]
  if (cn == "C1b") {
    clusnum11[s] <- ifelse(cuts6_c[s] == 2, "C1b1", "C1b2")
  } else {
    clusnum11[s] <- cn
  }
}
clusnum11 <- factor(clusnum11, levels = hier_levels_11)

cc11 <- c(
  "A1"   = "#E31A1C",  "A2"   = "#FF7F00",
  "B1"   = "#33A02C",  "B2a"  = "#B2DF8A",  "B2b"  = "#A6D854",
  "C1a"  = "#1F78B4",  "C1b1" = "#6A3D9A",  "C1b2" = "#B15928",
  "C2a"  = "#A6CEE3",  "C2b1" = "#CAB2D6",  "C2b2" = "#B3B3B3"
)

rsc11 <- cc11[as.character(clusnum11)]
names(rsc11) <- names(clusnum11)

# 11クラスター版ヒートマップ
pdf(file = here::here("output", "survey", "heatmaps", "ASGARD_hm_survey_16S_11clusters.pdf"),
    width = 20, height = 20)

heatmap.2(asgard_frtmat,
  distfun = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col = viridis, Rowv = h1$rowDendrogram, Colv = h1$colDendrogram,
  margins = c(15, 15), scale = "none",
  main = "ASGARD Survey 16S — Bray/ward.D", trace = "none",
  cexCol = 0.2, cexRow = 0.3)

heatmap.2(asgard_frtmat,
  distfun = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col = viridis,
  RowSideColors = rsc11[rownames(asgard_frtmat)],
  ColSideColors = colrsc[colnames(asgard_frtmat)],
  Rowv = h1$rowDendrogram, Colv = h1$colDendrogram,
  margins = c(15, 15), scale = "none",
  main = "ASGARD Survey 16S — 11 row clusters, 6 col clusters", trace = "none",
  cexCol = 0.2, cexRow = 0.3)

heatmap.2(asgard_frtmat,
  distfun = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col = viridis,
  RowSideColors = rsc11[rownames(asgard_frtmat)],
  ColSideColors = colrsc[colnames(asgard_frtmat)],
  Rowv = h1$rowDendrogram, Colv = h1$colDendrogram,
  margins = c(15, 15), scale = "none",
  main = "ASGARD Survey 16S — 11 clusters (station names)", trace = "none",
  cexCol = 0.2, cexRow = 0.3,
  labRow = meta_asgard[rownames(asgard_frtmat), "station"])

heatmap.2(asgard_frtmat,
  distfun = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col = viridis,
  RowSideColors = rsc11[rownames(asgard_frtmat)],
  ColSideColors = colrsc[colnames(asgard_frtmat)],
  Rowv = h1$rowDendrogram, Colv = h1$colDendrogram,
  margins = c(2, 2), scale = "none",
  main = "ASGARD Survey 16S — 11 row clusters, 6 col clusters", trace = "none",
  labRow = FALSE, labCol = FALSE)

# ASV名（列名）表示版 / With ASV names as column labels
heatmap.2(asgard_frtmat,
  distfun = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col = viridis,
  RowSideColors = rsc11[rownames(asgard_frtmat)],
  ColSideColors = colrsc[colnames(asgard_frtmat)],
  Rowv = h1$rowDendrogram, Colv = h1$colDendrogram,
  margins = c(25, 5), scale = "none",
  main = "ASGARD Survey 16S — 11 clusters (ASV names)", trace = "none",
  cexCol = 0.15, cexRow = 0.15,
  labRow = FALSE)

dev.off()

message("  11-cluster heatmap: output/survey/heatmaps/ASGARD_hm_survey_16S_11clusters.pdf")
message("  clusnum11, cc11, rsc11, hier_levels_11 exported")

# ==============================================================================
# Section 6: Division C ヒートマップ (97 samples x 242 ASVs)
# Division C heatmap with 6 sub-clusters
# ==============================================================================

c_clusters <- c("C1a", "C1b1", "C1b2", "C2a", "C2b1", "C2b2")
c_samples <- names(clusnum11)[clusnum11 %in% c_clusters]

c_mat <- asgard_frtprop[c_samples, ]
c_mat <- c_mat[, colSums(c_mat) > 0]

c_rsc <- cc11[as.character(clusnum11[c_samples])]
names(c_rsc) <- c_samples

pdf(file = here::here("output", "survey", "heatmaps", "ASGARD_hm_divisionC.pdf"),
    width = 20, height = 15)

heatmap.2(
  as.matrix(c_mat),
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  RowSideColors = c_rsc[rownames(c_mat)],
  margins   = c(15, 15),
  scale     = "none",
  main      = "Division C - 16S heatmap (6 sub-clusters)",
  trace     = "none",
  cexCol    = 0.2,
  cexRow    = 0.3
)

heatmap.2(
  as.matrix(c_mat),
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  RowSideColors = c_rsc[rownames(c_mat)],
  margins   = c(15, 15),
  scale     = "none",
  main      = "Division C - 16S heatmap (station names)",
  trace     = "none",
  cexCol    = 0.2,
  cexRow    = 0.3,
  labRow    = meta_asgard[rownames(c_mat), "station"]
)

heatmap.2(
  as.matrix(c_mat),
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  RowSideColors = c_rsc[rownames(c_mat)],
  margins   = c(2, 2),
  scale     = "none",
  main      = "Division C - 16S heatmap",
  trace     = "none",
  labRow    = FALSE,
  labCol    = FALSE
)

dev.off()

message("  Division C heatmap: output/survey/heatmaps/ASGARD_hm_divisionC.pdf")

message("\nS02_heatmaps_16S.R: done.")
message("  clusnum10 length:  ", length(clusnum10), " (", nclus10, " row clusters)")
message("  hier_names: ", paste(hier_names, collapse = ", "))
message("  colclusnum length: ", length(colclusnum), " (", colnclus, " col clusters)")
message("  PDF: output/survey/heatmaps/ASGARD_hm_survey_16S.pdf")

