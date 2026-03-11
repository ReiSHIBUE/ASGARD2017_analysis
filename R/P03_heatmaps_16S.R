### 03_heatmaps_16S.R
### ASGARD 2017 Processing Site Analysis — 16S Heatmaps
### 16S ヒートマップスクリプト
###
### REQUIRES (from 01_data_prep.R and 02_ternary_plots.R):
###   asgard_filtered_p_hm2  - ASV proportion matrix, 78*221
###   meta_asgard_p2         - metadata for 78 samples
###   sample_rgb2            - row side colours by filter size (length 78)
###   asv_rgb2               - col side colours from ternary RGB (length 221)
###   zero_cols              - ASV column names absent in 0.2 µm fraction
###
### PRODUCES (used by later scripts):
###   h3           - heatmap.2 object (full 78-sample heatmap)
###   clusnum_p    - sample cluster assignments (length 78), values 1-4
###   rsc_p        - 4-colour palette vector for clusters
###   sample_rgb3  - row side colours aligned to clusnum_p (length 78)
###
### OUTPUT: output/ASGARD_hm_processing_5000over.pdf

library(gplots)
library(viridis)
library(vegan)

# ==============================================================================
# Section 1: 49 ASVのサブセット (3と20 µmに特異的なASVs)
# 49-ASV subset: particle-associated ASVs absent in 0.2 µm
# ==============================================================================

asgard_filtered_49asvs <- asgard_filtered_p_hm2[, zero_cols] # 78*49
asgard_filtered_49asvs2 <- asgard_filtered_49asvs[
  rowSums(asgard_filtered_49asvs) > 0,
  colSums(asgard_filtered_49asvs) > 0
]

meta_asgard_p2_49 <- meta_asgard_p2[rownames(asgard_filtered_49asvs2), ]

# ==============================================================================
# Section 2: ヒートマップの作成 / Generate heatmaps
# ==============================================================================

pdf(file = here::here("output", "heatmaps", "ASGARD_hm_processing_5000over.pdf"),
    width = 20, height = 20)

# h3: 全78サンプル、フィルターサイズ色付き / Full 78-sample heatmap coloured by filter size
h3 <- heatmap.2(
  (asgard_filtered_p_hm2)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col          = viridis,
  RowSideColors = sample_rgb2,
  ColSideColors = asv_rgb2,
  margins  = c(15, 15),
  scale    = "none",
  main     = "ASGARD_bray/ward.D2",
  trace    = "none",
  cexCol   = 0.2, # 列の文字サイズ / column label size
  cexRow   = 0.2, # 行の文字サイズ / row label size
  labRow   = meta_asgard_p2$side
)

# h4: 49 particle-associated ASVs のみ / 49-ASV particle-associated subset
h4 <- heatmap.2(
  (asgard_filtered_49asvs2)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col       = viridis,
  margins   = c(15, 15),
  scale     = "none",
  main      = "ASGARD_bray/ward.D2",
  trace     = "none",
  cexCol    = 0.8,
  cexRow    = 0.8,
  labRow    = meta_asgard_p2_49$side
)

# ==============================================================================
# Section 3: クラスター割り当て / Extract cluster assignments from h3
# ==============================================================================

nclus_p  <- 4
oldclus_p   <- cutree(as.hclust(h3$rowDendrogram), k = nclus_p)
oldorder_p  <- unname(rle(oldclus_p[as.hclust(h3$rowDendrogram)$order])$values)
neworder_p  <- (1:nclus_p)
names(neworder_p) <- oldorder_p
clusnum_p   <- unname(neworder_p[as.character(oldclus_p)])
names(clusnum_p) <- names(oldclus_p)
# clusnum_p: 整数ベクトル (1-4) / integer vector, length 78, values 1-4

# 4色パレット定義 / Define 4-colour palette
rsc_p <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
sample_rgb3 <- rsc_p[clusnum_p] # 行サイドカラー (クラスター色) / row-side colours by cluster

# h8: クラスター色付きで再描画、h3のデンドログラム固定 / Replot with cluster colours, fixed dendrogram
h8 <- heatmap.2(
  (asgard_filtered_p_hm2)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col          = viridis,
  RowSideColors = sample_rgb3,
  Rowv         = h3$rowDendrogram,
  ColSideColors = asv_rgb2,
  margins  = c(15, 15),
  scale    = "none",
  main     = "ASGARD_bray/ward.D2",
  trace    = "none",
  cexCol   = 0.2,
  cexRow   = 0.2,
  labRow   = meta_asgard_p2$side
)

# h9: サンプルIDラベル付き / Row labels = sample IDs
h9 <- heatmap.2(
  (asgard_filtered_p_hm2)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col          = viridis,
  RowSideColors = sample_rgb3,
  Rowv         = h3$rowDendrogram,
  ColSideColors = asv_rgb2,
  margins  = c(15, 15),
  scale    = "none",
  main     = "ASGARD_bray/ward.D2",
  trace    = "none",
  cexCol   = 0.2,
  cexRow   = 0.2,
  labRow   = meta_asgard_p2$Sample
)

dev.off()

message("03_heatmaps_16S.R: done. clusnum_p (length=", length(clusnum_p), "), PDF written.")
