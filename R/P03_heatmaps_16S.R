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

pdf(file = here::here("output_p", "heatmaps", "ASGARD_hm_processing_5000over.pdf"),
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

# ==============================================================================
# Section 4: メイン図 (h_main) — S02 風: クラスター名 overlay, タイトル/key 削除
# Main figure with S02-style cluster-name overlay (no title, no key)
# ==============================================================================

row_ord_p     <- as.hclust(h3$rowDendrogram)$order
samples_ord_p <- rownames(asgard_filtered_p_hm2)[row_ord_p]
clus_ord_p    <- as.character(clusnum_p[samples_ord_p])
cluster_centers_p <- sapply(c("1", "2", "3", "4"), function(cn) {
  idx <- which(clus_ord_p == cn)
  if (length(idx) == 0) return(NA)
  mean(idx)
})
nrow_h_p <- nrow(asgard_filtered_p_hm2)

heatmap.2(
  (asgard_filtered_p_hm2)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = sample_rgb3,
  Rowv          = h3$rowDendrogram,
  ColSideColors = asv_rgb2,
  margins   = c(2, 2),
  scale     = "none",
  main      = "",
  trace     = "none",
  labRow    = FALSE,
  labCol    = FALSE,
  key       = FALSE,
  add.expr  = {
    text(x = rep(-6, length(cluster_centers_p)),
         y = cluster_centers_p,
         labels = names(cluster_centers_p),
         col = "white",
         xpd = NA, adj = c(0.5, 0.5),
         srt = 90,
         cex = 3.5, font = 2)
  }
)

# ==============================================================================
# Section 5: 環境PCA + Top-down tree walk (Wilcoxon + BH)
# Top-down dendrogram walk: split sig if any retained PC has BH-adjusted p<0.05
# ==============================================================================

env_vars_p <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)",
                "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)", "depth_m")
env_df_p   <- meta_asgard_p2[, c("Sample", env_vars_p)]
env_complete_p <- env_df_p[complete.cases(env_df_p), ]
message("Tree walk — env complete cases: ", nrow(env_complete_p), " / 78")

env_mat_p <- as.matrix(env_complete_p[, env_vars_p])
log_vars_p <- c("NO3(uM)", "PO4(uM)", "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)")
for (v in log_vars_p) env_mat_p[, v] <- log1p(env_mat_p[, v])
env_scaled_p <- scale(env_mat_p)

pca_p_result <- prcomp(env_scaled_p, center = FALSE, scale. = FALSE)
rownames(pca_p_result$x) <- env_complete_p$Sample
var_expl_p   <- summary(pca_p_result)$importance["Proportion of Variance", ]
cum_var_p    <- cumsum(var_expl_p)
n_pcs_keep_p <- which(cum_var_p >= 0.90)[1]
pca_p_scores <- as.data.frame(pca_p_result$x[, 1:n_pcs_keep_p])
message(sprintf("  Retained PC1-PC%d (cum var = %.1f%%)",
                n_pcs_keep_p, 100 * cum_var_p[n_pcs_keep_p]))

# Dendrogram → leaves below each merge
hc_p <- as.hclust(h3$rowDendrogram)
leaves_below_p <- vector("list", nrow(hc_p$merge))
for (i in seq_len(nrow(hc_p$merge))) {
  l <- hc_p$merge[i, 1]; r <- hc_p$merge[i, 2]
  ll <- if (l < 0) -l else leaves_below_p[[l]]
  rr <- if (r < 0) -r else leaves_below_p[[r]]
  leaves_below_p[[i]] <- c(ll, rr)
}

# Test one node: returns sig flag + child sample sets
node_significant <- function(i) {
  l <- hc_p$merge[i, 1]; r <- hc_p$merge[i, 2]
  left_leaves  <- if (l < 0) -l else leaves_below_p[[l]]
  right_leaves <- if (r < 0) -r else leaves_below_p[[r]]
  left_samples  <- hc_p$labels[left_leaves]
  right_samples <- hc_p$labels[right_leaves]
  if (length(left_samples) < 5 || length(right_samples) < 5) return(FALSE)
  lin  <- left_samples  %in% rownames(pca_p_scores)
  rin  <- right_samples %in% rownames(pca_p_scores)
  if (sum(lin) < 5 || sum(rin) < 5) return(FALSE)
  pvals <- sapply(seq_len(n_pcs_keep_p), function(k) {
    tryCatch(wilcox.test(pca_p_scores[left_samples[lin], k],
                         pca_p_scores[right_samples[rin], k])$p.value,
             error = function(e) NA)
  })
  any(p.adjust(pvals, method = "BH") < 0.05, na.rm = TRUE)
}

# Per-node Wilcoxon + BH record for output
node_test_record <- function(i) {
  l <- hc_p$merge[i, 1]; r <- hc_p$merge[i, 2]
  left_leaves  <- if (l < 0) -l else leaves_below_p[[l]]
  right_leaves <- if (r < 0) -r else leaves_below_p[[r]]
  left_samples  <- hc_p$labels[left_leaves]
  right_samples <- hc_p$labels[right_leaves]
  if (length(left_samples) < 5 || length(right_samples) < 5) return(NULL)
  lin  <- left_samples  %in% rownames(pca_p_scores)
  rin  <- right_samples %in% rownames(pca_p_scores)
  if (sum(lin) < 5 || sum(rin) < 5) return(NULL)
  pvals <- sapply(seq_len(n_pcs_keep_p), function(k) {
    tryCatch(wilcox.test(pca_p_scores[left_samples[lin], k],
                         pca_p_scores[right_samples[rin], k])$p.value,
             error = function(e) NA_real_)
  })
  p_adj <- p.adjust(pvals, method = "BH")
  data.frame(node_id = i, height = hc_p$height[i],
             n_left = length(left_samples), n_right = length(right_samples),
             PC = paste0("PC", seq_len(n_pcs_keep_p)),
             p_value = round(pvals, 5),
             p_adj_BH = round(p_adj, 5),
             sig_BH = p_adj < 0.05)
}

treewalk_node_records <- do.call(rbind, lapply(seq_len(nrow(hc_p$merge)),
                                              node_test_record))

# Recursive top-down walk → list of accepted (= terminal) leaf-sets
treewalk_clusters <- list()
walk_node_topdown <- function(i) {
  l <- hc_p$merge[i, 1]; r <- hc_p$merge[i, 2]
  if (node_significant(i)) {
    if (l < 0) treewalk_clusters[[length(treewalk_clusters) + 1]] <<-
      hc_p$labels[-l]
    else walk_node_topdown(l)
    if (r < 0) treewalk_clusters[[length(treewalk_clusters) + 1]] <<-
      hc_p$labels[-r]
    else walk_node_topdown(r)
  } else {
    treewalk_clusters[[length(treewalk_clusters) + 1]] <<-
      hc_p$labels[leaves_below_p[[i]]]
  }
}
walk_node_topdown(nrow(hc_p$merge))

# Reorder cluster ids to follow dendrogram leaf order
leaf_order_p   <- hc_p$labels[hc_p$order]
cl_first_pos   <- sapply(treewalk_clusters,
                         function(s) min(match(s, leaf_order_p)))
treewalk_clusters <- treewalk_clusters[order(cl_first_pos)]
clusnum_tw <- integer(length(hc_p$labels))
names(clusnum_tw) <- hc_p$labels
for (k in seq_along(treewalk_clusters))
  clusnum_tw[treewalk_clusters[[k]]] <- k
clusnum_tw <- clusnum_tw[rownames(asgard_filtered_p_hm2)]
n_tw <- length(treewalk_clusters)
message(sprintf("  Tree walk yielded %d clusters", n_tw))
message("  Cluster sizes: ", paste(table(clusnum_tw), collapse = ", "))

# Save CSV outputs
dir.create(here::here("output_p", "treewalk"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(PC = paste0("PC", seq_along(var_expl_p)),
                     variance_pct  = round(var_expl_p * 100, 2),
                     cumulative_pct = round(cum_var_p * 100, 2),
                     retained = seq_along(var_expl_p) <= n_pcs_keep_p),
          here::here("output_p", "treewalk",
                     "processing_pca_summary.csv"),
          row.names = FALSE)
write.csv(round(as.data.frame(pca_p_result$rotation[, 1:n_pcs_keep_p]), 4),
          here::here("output_p", "treewalk",
                     "processing_pca_loadings.csv"))
write.csv(treewalk_node_records,
          here::here("output_p", "treewalk",
                     "processing_treewalk_node_results.csv"),
          row.names = FALSE)
write.csv(data.frame(Sample = names(clusnum_tw),
                     clusnum_p = unname(clusnum_p[names(clusnum_tw)]),
                     clusnum_treewalk = unname(clusnum_tw)),
          here::here("output_p", "treewalk",
                     "processing_treewalk_cluster_assignment.csv"),
          row.names = FALSE)
message("  CSV: output_p/treewalk/processing_pca_summary.csv")
message("  CSV: output_p/treewalk/processing_pca_loadings.csv")
message("  CSV: output_p/treewalk/processing_treewalk_node_results.csv")
message("  CSV: output_p/treewalk/processing_treewalk_cluster_assignment.csv")

# Palette (extend rsc_p if n_tw > 4)
rsc_tw_base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                 "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                 "#999999", "#66C2A5", "#FC8D62", "#8DA0CB")
rsc_tw <- rsc_tw_base[seq_len(n_tw)]
sample_rgb_tw <- rsc_tw[clusnum_tw]

# Centers for label overlay
row_ord_tw     <- as.hclust(h3$rowDendrogram)$order
samples_ord_tw <- rownames(asgard_filtered_p_hm2)[row_ord_tw]
clus_ord_tw    <- as.character(clusnum_tw[samples_ord_tw])
cluster_centers_tw <- sapply(as.character(seq_len(n_tw)), function(cn) {
  idx <- which(clus_ord_tw == cn)
  if (length(idx) == 0) return(NA)
  mean(idx)
})

# Page 6: Tree walk-derived clusters
heatmap.2(
  (asgard_filtered_p_hm2)^.25,
  distfun   = function(x) vegdist(x, method = "bray"),
  hclustfun = function(x) hclust(x, method = "ward.D"),
  col           = viridis,
  RowSideColors = sample_rgb_tw,
  Rowv          = h3$rowDendrogram,
  ColSideColors = asv_rgb2,
  margins   = c(2, 2),
  scale     = "none",
  main      = "",
  trace     = "none",
  labRow    = FALSE,
  labCol    = FALSE,
  key       = FALSE,
  add.expr  = {
    text(x = rep(-6, length(cluster_centers_tw)),
         y = cluster_centers_tw,
         labels = names(cluster_centers_tw),
         col = "white",
         xpd = NA, adj = c(0.5, 0.5),
         srt = 90,
         cex = 3.5, font = 2)
  }
)

dev.off()

message("03_heatmaps_16S.R: done. clusnum_p (length=", length(clusnum_p), "), PDF written.")
