### 98_build_caches.R
### ASGARD 2017 — Rebuild the comparison caches from the live session
### セッションのオブジェクトから比較用キャッシュを再作成する
###
### comparison_report.R, field_date_67stations.R and field_date_processing.R
### read output/comparison/cache/{survey,processing}_cache.rds, but no script in
### the repository wrote them — they were committed artifacts. That made those
### three scripts impossible to run from a bare checkout. This script rebuilds
### them from the objects the pipelines leave in the session.
###
### REQUIRES:
###   survey:     meta_asgard, clusnum11, asgard_seqcount, asgard_filtered,
###               cc11, hier_levels_11            (00_setup.R, S01, S02)
###   processing: meta_asgard_p2, clusnum_p, asgard_filtered_p_hm2,
###               asgard_vec_p, seqtab_16Smat     (00_setup.R, P01, P03)
###   optional:   comp_df_p7, rich_df_p, mapz     (P04 — needs a Stadia key)
###
### OUTPUT:
###   output/comparison/cache/survey_cache.rds
###   output/comparison/cache/processing_cache.rds

library(here)

CACHE_DIR <- here("output", "comparison", "cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Survey cache
# ==============================================================================

# S04 adds a "Sample" column to meta_asgard; drop it so the cache matches the
# metadata columns the consumers expect.
a_map <- meta_asgard[, colnames(meta_asgard) != "Sample", drop = FALSE]
a_map$cluster11 <- factor(as.character(clusnum11[rownames(a_map)]),
                          levels = hier_levels_11)
a_map$ID <- seq_len(nrow(a_map))

survey_cache <- list(
  a_map           = a_map,
  asgard_seqcount = asgard_seqcount,
  asgard_filtered = asgard_filtered,
  richness        = rowSums(asgard_seqcount > 0),
  richness258     = rowSums(asgard_filtered > 0),
  cc11            = cc11,
  hier_levels_11  = hier_levels_11
)

saveRDS(survey_cache, file.path(CACHE_DIR, "survey_cache.rds"))

# ==============================================================================
# Processing cache
# ==============================================================================

cc_p <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3")

processing_cache <- list(
  cc_p                  = cc_p,
  clusnum_p             = clusnum_p,
  meta_asgard_p2        = meta_asgard_p2,
  asgard_filtered_p_hm2 = asgard_filtered_p_hm2,
  # rebuilt the way P01 does it: `asgard_seqcount` itself cannot be used here
  # because S01 overwrites that name with the 181-sample survey matrix.
  # S01 が同名オブジェクトを上書きするため、P01 と同じ方法で作り直す。
  asgard_seqcount_p     = seqtab_16Smat[asgard_vec_p, ]
)

# P04_maps.R builds these, and needs a Stadia Maps key to run at all.
# P04 由来のオブジェクトはキーが無いと存在しないため、あるときだけ追加する。
for (nm in c("comp_df_p7", "rich_df_p", "mapz")) {
  if (exists(nm)) processing_cache[[nm]] <- get(nm)
}

saveRDS(processing_cache, file.path(CACHE_DIR, "processing_cache.rds"))

message("98_build_caches.R: done.")
message("  RDS: output/comparison/cache/survey_cache.rds (",
        length(survey_cache), " objects)")
message("  RDS: output/comparison/cache/processing_cache.rds (",
        length(processing_cache), " objects)")
