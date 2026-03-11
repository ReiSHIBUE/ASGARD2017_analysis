### S09_network.R
### ASGARD 2017 Survey Site Analysis — ASV Co-occurrence Network
### 共起ネットワーク分析スクリプト（サーベイサイト）
###
### REQUIRES (from S01, S02):
###   asgard_frtprop   - fourth-root proportion matrix 181×258
###   meta_asgard      - metadata 181×45
###   clusnum          - sample cluster assignments (length 181)
###   colclusnum       - ASV column cluster assignments (length 258)
###
### PRODUCES:
###   asgard_cor_mat   - Spearman correlation matrix (258×258)
###   asgard_network   - igraph network object (edges: |r| > 0.6)
###
### OUTPUT:
###   output/survey/network/ASGARD_network_survey.pdf
###
### NOTE: igraph required. tSNE layout requires Rtsne:
###   install.packages(c("igraph", "Rtsne"))

library(vegan)
library(tidyverse)

if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph is required: install.packages('igraph')")
library(igraph)

# ==============================================================================
# Section 1: スピアマン相関行列 / Spearman correlation matrix between ASVs
# ==============================================================================

asgard_cor_mat <- cor(asgard_frtprop, method = "spearman")
hist(asgard_cor_mat, breaks = 100,
     main = "Distribution of Spearman Correlations — Survey ASVs",
     xlab = "Spearman r")

# ==============================================================================
# Section 2: 閾値適用・ネットワーク構築 / Apply threshold and build network
# |r| < 0.6 → 0; self-connections removed
# ==============================================================================

asgard_cor_thresh <- asgard_cor_mat
asgard_cor_thresh[abs(asgard_cor_thresh) < 0.6] <- 0
diag(asgard_cor_thresh) <- 0

asgard_network <- graph_from_adjacency_matrix(
  asgard_cor_thresh,
  mode     = "undirected",
  weighted = TRUE,
  diag     = FALSE
)

# 孤立ノードを削除 / Remove isolated vertices
asgard_network <- delete_vertices(asgard_network, degree(asgard_network) == 0)

# ネットワークレイアウトを計算 / Compute layout
asgard_network_layout <- layout_nicely(asgard_network)

# ==============================================================================
# Section 3: 頂点属性の設定 / Assign vertex colours by column cluster
# ==============================================================================

# colclusnum で列クラスター毎に色付け / Colour by ASV column cluster
asv_in_net    <- V(asgard_network)$name
net_colclusters <- colclusnum[asv_in_net]
net_colclusters[is.na(net_colclusters)] <- 0
vertex_colors <- plasma(8)[net_colclusters]
vertex_colors[net_colclusters == 0] <- "gray"

# depth_type が支配的な深度帯で色付けする場合の準備 / Node size by mean relative abundance
node_size <- sqrt(colMeans(asgard_frtprop)[asv_in_net]) * 30
node_size[is.na(node_size)] <- 1

# ==============================================================================
# Section 4: ネットワーク図を保存 / Write network plots to PDF
# ==============================================================================

pdf(file = here::here("output", "survey", "network", "ASGARD_network_survey.pdf"),
    width = 10, height = 10)

# layout_nicely を使った図 / Network with default layout
plot(asgard_network,
     layout           = asgard_network_layout,
     vertex.size      = node_size,
     vertex.color     = vertex_colors,
     vertex.frame.color = NA,
     vertex.label     = NA,
     edge.width       = abs(E(asgard_network)$weight) * 2,
     edge.color       = "gray70",
     main             = "ASGARD Survey ASV Co-occurrence Network (|r| > 0.6)")

# tSNE レイアウト (オプション) / Optional tSNE layout
if (requireNamespace("Rtsne", quietly = TRUE)) {
  tryCatch({
    library(Rtsne)
    asgard_tsne <- Rtsne::Rtsne(t(asgard_frtprop)^0.25,
                                 perplexity     = 10,
                                 dims           = 2,
                                 check_duplicates = FALSE,
                                 pca            = TRUE,
                                 max_iter       = 2500,
                                 verbose        = FALSE)
    tsne_coords <- asgard_tsne$Y
    rownames(tsne_coords) <- colnames(asgard_frtprop)

    tsne_layout <- tsne_coords[asv_in_net, ]

    plot(asgard_network,
         layout           = tsne_layout,
         vertex.size      = node_size,
         vertex.color     = vertex_colors,
         vertex.frame.color = NA,
         vertex.label     = NA,
         edge.width       = abs(E(asgard_network)$weight) * 2,
         edge.color       = "gray70",
         main             = "ASGARD Survey Network — tSNE Layout")
  }, error = function(e) {
    message("tSNE layout failed: ", conditionMessage(e))
  })
} else {
  message("Rtsne not installed — skipping tSNE layout.")
  message("  Install with: install.packages('Rtsne')")
}

dev.off()

message("S09_network.R: done.")
message("  Network: ", vcount(asgard_network), " vertices, ",
        ecount(asgard_network), " edges")
message("  PDF: output/survey/network/ASGARD_network_survey.pdf")
