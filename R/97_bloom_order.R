### 97_bloom_order.R
### ASGARD 2017 Survey — Bloom index and the cluster ordering derived from it
### ブルームインデックスと、それに基づくクラスター並び順
###
### The bloom index is the first component of the 9-variable environmental PCA,
### rescaled to [0, 1] and flipped so 0 = pre-bloom and 1 = post-bloom — the
### same definition S19_env_boxplots.R uses. Figures that order clusters by
### bloom stage (the faceted survey map and the T-S diagram) share the ordering
### computed here instead of each recomputing it.
###
### REQUIRES (from 00_setup.R, S01, S02):
###   meta_asgard, clusnum11, hier_levels_11
###
### PRODUCES:
###   bloom_index_df          - per-sample PC1, bloom index and cluster
###   cluster11_bloom         - named vector: mean bloom index per cluster
###   cluster11_bloom_order   - cluster names grouped by division (A, B, C) and,
###                             within each division, ordered by decreasing mean
###                             bloom index (post-bloom -> pre-bloom, i.e. the
###                             reverse of environmental PC1)
###   cluster11_bloom_order_flat - the same clusters ordered by bloom index only,
###                             ignoring division

library(here)

local({
  env_vars <- c("temp", "salinity", "DO", "NO3(uM)", "PO4(uM)", "Sil(uM)",
                "NH4(uM)", "FlECO-AFL(mg/m^3)", "depth_m")
  log_vars <- c("NO3(uM)", "PO4(uM)", "Sil(uM)", "NH4(uM)", "FlECO-AFL(mg/m^3)")

  df <- meta_asgard
  df$cluster11 <- factor(as.character(clusnum11[rownames(df)]),
                         levels = hier_levels_11)
  df$sample_id <- rownames(df)

  env_df       <- df[, c("sample_id", "cluster11", env_vars)]
  env_complete <- env_df[complete.cases(env_df), ]

  env_mat <- env_complete[, env_vars]
  for (v in log_vars) env_mat[[v]] <- log1p(env_mat[[v]])

  pca <- prcomp(scale(env_mat), center = FALSE, scale. = FALSE)

  rescale01 <- function(x) (x - min(x)) / (max(x) - min(x))

  out <- data.frame(
    sample_id   = env_complete$sample_id,
    cluster11   = env_complete$cluster11,
    PC1         = pca$x[, 1],
    bloom_index = 1 - rescale01(pca$x[, 1]),
    row.names   = NULL
  )

  means <- tapply(out$bloom_index, out$cluster11, mean)

  # division first (A, B, C), then descending bloom index within each division
  # まず division (A/B/C)、その中でブルームインデックス降順
  flat     <- names(sort(means, decreasing = TRUE))
  division <- sub("^(.).*", "\\1", flat)
  ordered  <- unlist(lapply(c("A", "B", "C"), function(d) flat[division == d]),
                     use.names = FALSE)

  bloom_index_df             <<- out
  cluster11_bloom            <<- means
  cluster11_bloom_order      <<- ordered
  cluster11_bloom_order_flat <<- flat
})

message("97_bloom_order.R: done.")
message("  facet order (division, then bloom index): ",
        paste(cluster11_bloom_order, collapse = " > "))
message("  bloom index alone:                       ",
        paste(cluster11_bloom_order_flat, collapse = " > "))
