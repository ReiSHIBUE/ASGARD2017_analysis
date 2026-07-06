### comparison_report.R — living comparison report (Survey vs Processing)
### Dependency-free multi-page PDF (grid + gridExtra + ggplot2). Re-run to rebuild.
### TO ADD A SECTION: append a block in the "REPORT BODY" area and re-run.
### Reads cached objects from output/comparison/cache/ (regenerate with the
### rebuild scripts if missing). Run:  Rscript R/comparison_report.R
suppressMessages({library(tidyverse); library(vegan); library(grid); library(gridExtra); library(here)})
set.seed(1)

CACHE_DIR <- here::here("output", "comparison", "cache")
sc <- readRDS(file.path(CACHE_DIR, "survey_cache.rds"))
pc <- readRDS(file.path(CACHE_DIR, "processing_cache.rds"))
REPORT_DATE <- "2026-07-02"   # pass a fixed date (no Sys.time needed)
OUT <- here::here("output", "comparison", "ASGARD_survey_vs_processing_report.pdf")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

# ---------- page helpers ----------
text_page <- function(title, lines, sub = NULL) {
  grid.newpage()
  grid.text(title, x = 0.05, y = 0.94, just = c("left","top"),
            gp = gpar(fontsize = 22, fontface = "bold", col = "#123456"))
  y <- 0.86
  if (!is.null(sub)) { grid.text(sub, x = 0.05, y = 0.885, just = c("left","top"),
            gp = gpar(fontsize = 12, fontface = "italic", col = "grey30")); y <- 0.84 }
  for (ln in lines) {
    if (ln == "") { y <- y - 0.018; next }
    wrapped <- strwrap(ln, width = 105)
    for (w in wrapped) {
      grid.text(w, x = 0.06, y = y, just = c("left","top"), gp = gpar(fontsize = 12))
      y <- y - 0.026
    }
    y <- y - 0.010
  }
}
table_page <- function(title, df, fontsize = 10) {
  grid.newpage()
  grid.text(title, x = 0.05, y = 0.94, just = c("left","top"),
            gp = gpar(fontsize = 20, fontface = "bold", col = "#123456"))
  tg <- tableGrob(df, rows = NULL, theme = ttheme_default(base_size = fontsize))
  pushViewport(viewport(y = 0.46, height = 0.82)); grid.draw(tg); popViewport()
}

# ==========================================================================
# DATA / ANALYSES (recomputed from caches so the report is self-contained)
# ==========================================================================

## --- (A) 6 environmental variables: survey vs processing ---
vars6 <- c("temp","salinity","DO","FlECO-AFL(mg/m^3)","NO3(uM)","depth_m")
labs6 <- c("Temperature (C)","Salinity (PSU)","DO (umol/kg)","FlECO-AFL (mg/m3)","NO3 (uM)","Depth (m)")
surv <- sc$a_map; proc <- pc$meta_asgard_p2
proc_u <- proc[!duplicated(paste(proc$station, proc$depth_m)), ]
alpha6 <- 0.05/length(vars6)
env_long <- data.frame(); env_res <- data.frame()
for (i in seq_along(vars6)) {
  s <- as.numeric(surv[[vars6[i]]]); s <- s[!is.na(s)]
  p <- as.numeric(proc_u[[vars6[i]]]); p <- p[!is.na(p)]
  w <- wilcox.test(s, p)
  env_res <- rbind(env_res, data.frame(Variable = labs6[i],
    Med_Survey = round(median(s),2), Med_Proc = round(median(p),2),
    Survey_range = sprintf("%.1f-%.1f", min(s), max(s)),
    Proc_range = sprintf("%.1f-%.1f", min(p), max(p)),
    p = signif(w$p.value,3), sig = ifelse(w$p.value < alpha6, "*", "ns")))
  env_long <- rbind(env_long,
    data.frame(dataset="Survey", variable=labs6[i], value=s),
    data.frame(dataset="Processing", variable=labs6[i], value=p))
}
env_long$dataset  <- factor(env_long$dataset, levels=c("Survey","Processing"))
env_long$variable <- factor(env_long$variable, levels=labs6)
env_fig <- ggplot(env_long, aes(dataset, value, fill=dataset)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) + geom_jitter(width=0.2, size=0.4, alpha=0.3) +
  facet_wrap(~variable, scales="free_y", ncol=3) +
  scale_fill_manual(values=c(Survey="#1F78B4", Processing="#E31A1C"), guide="none") +
  labs(x=NULL, y=NULL, title="6 environmental variables: Survey vs Processing") +
  theme_bw(base_size=13) + theme(strip.text=element_text(face="bold"),
                                 plot.title=element_text(face="bold"))

## --- (B) 5-variable db-RDA: survey vs processing ---
vars5 <- c("temp","salinity","DO","NO3(uM)","FlECO-AFL(mg/m^3)")
run_dbrda <- function(comm, meta) {
  frt <- comm^0.25; frt <- frt[, colSums(frt) > 0]
  samp <- intersect(rownames(frt), rownames(meta)); env <- meta[samp, vars5]
  cc <- complete.cases(env); envs <- as.data.frame(scale(env[cc,]))
  colnames(envs) <- c("temp","salinity","DO","NO3","FlECO")
  D <- vegdist(frt[samp,][cc,], "bray")
  m <- capscale(D ~ temp+salinity+DO+NO3+FlECO, data=envs)
  ov <- anova.cca(m, permutations=999); mg <- as.data.frame(anova.cca(m, by="margin", permutations=999))
  list(n=sum(cc), asv=ncol(frt), con=m$CCA$tot.chi/m$tot.chi*100,
       r2=RsquareAdj(m)$adj.r.squared, F=ov$F[1], p=ov$`Pr(>F)`[1], mg=mg[1:5,])
}
S <- run_dbrda(sc$asgard_filtered, sc$a_map)
P <- run_dbrda(pc$asgard_filtered_p_hm2, pc$meta_asgard_p2)
dbrda_overall <- data.frame(
  Dataset=c("Survey","Processing"), n=c(S$n,P$n), ASVs=c(S$asv,P$asv),
  Constrained_pct=round(c(S$con,P$con),1), Adj_R2=round(c(S$r2,P$r2),3),
  Overall_F=round(c(S$F,P$F),2), Overall_p=c(S$p,P$p))
dbrda_var <- data.frame(
  Variable=c("temp","salinity","DO","NO3","FlECO"),
  Survey_F=round(S$mg$F,2), Survey_p=round(S$mg$`Pr(>F)`,3),
  Proc_F=round(P$mg$F,2), Proc_p=round(P$mg$`Pr(>F)`,3))
dbrda_var$Survey_sig <- ifelse(dbrda_var$Survey_p<0.05,"*","ns")
dbrda_var$Proc_sig   <- ifelse(dbrda_var$Proc_p<0.05,"*","ns")

# ==========================================================================
# REPORT BODY  (>>> ADD NEW SECTIONS BY APPENDING BLOCKS BEFORE dev.off() <<<)
# ==========================================================================
pdf(OUT, width = 11, height = 8.5)

## Cover
text_page("ASGARD 2017 — Survey vs Processing", sub = paste("Comparison report ·", REPORT_DATE),
  c("A living summary of differences between the survey (0.2 um, 11 clusters, 181 samples)",
    "and processing (0.2/3/20 um, 4 clusters, 78 samples) datasets.",
    "",
    "Contents:",
    "  1. Environmental variables (6-variable boxplot comparison)",
    "  2. db-RDA (5-variable constrained ordination)",
    "",
    "Note: processing environmental values are triplicated across filter fractions;",
    "they are de-duplicated to unique station x depth for fair environmental comparison."))

## Section 1 — environment
text_page("1. Environmental variables (6-var)",
  c("Q: Do survey and processing sample different environmental conditions?",
    "Method: Wilcoxon rank-sum per variable; processing de-duplicated to unique water",
    "samples (78 -> 32); Bonferroni alpha = 0.05/6 = 0.0083.",
    "",
    "Result: NO variable differs significantly after correction (all ns).",
    "  - Environments broadly OVERLAP between the two datasets.",
    "  - Survey spans a WIDER range of salinity and temperature (more diverse water masses).",
    "  - Temperature is closest to differing (survey slightly warmer)."))
print(env_fig)
table_page("1. Environmental variables — Wilcoxon (Bonferroni a=0.0083)", env_res, fontsize = 11)

## Section 2 — db-RDA
text_page("2. db-RDA (5-var constrained ordination)",
  c("Q: Which environmental variables structure the microbial community, and how much?",
    "Model: capscale(Bray-Curtis ~ temp + salinity + DO + NO3 + FlECO), 4th-root ASV props.",
    "",
    "Result:",
    "  - Both overall models are significant (p = 0.001); environment structures both.",
    "  - SURVEY: all 5 variables significant (esp. salinity & NO3); adj R2 ~ 0.15.",
    "  - PROCESSING: ONLY salinity significant; adj R2 ~ 0.10.",
    "  => Survey community tracks multiple gradients; processing is driven mainly by salinity.",
    "",
    "Caveats: marginal ANOVA is collinearity-sensitive; processing n is smaller (75 vs 174);",
    "CAP axis % are defined differently in S05 vs P06 and are not directly comparable."))
table_page("2. db-RDA — overall model", dbrda_overall, fontsize = 11)
table_page("2. db-RDA — per-variable marginal ANOVA (F, p)", dbrda_var, fontsize = 11)

## --- Section 3: Stations (survey vs processing) ---
stn_proc <- pc$meta_asgard_p2 %>% dplyr::group_by(Station = station) %>%
  dplyr::summarise(n_samples = dplyr::n(),
                   n_depths  = dplyr::n_distinct(depth_type),
                   .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_samples), Station) %>% as.data.frame()
stn_surv_tab <- as.data.frame(table(sc$a_map$station), stringsAsFactors = FALSE)
colnames(stn_surv_tab) <- c("Station", "n")
stn_surv_tab <- stn_surv_tab[order(-stn_surv_tab$n, stn_surv_tab$Station), ]
shared_stn <- sort(intersect(unique(sc$a_map$station), unique(pc$meta_asgard_p2$station)))
stn_proc$Shared <- ifelse(stn_proc$Station %in% shared_stn, "yes", "")

make_wide <- function(df, groups = 3) {
  n <- nrow(df); per <- ceiling(n / groups)
  parts <- lapply(seq_len(groups), function(g) {
    idx <- ((g - 1) * per + 1):min(g * per, n)
    d <- df[idx, , drop = FALSE]
    if (nrow(d) < per)
      d <- rbind(d, data.frame(Station = rep("", per - nrow(d)), n = rep("", per - nrow(d))))
    setNames(d, c(paste0("Station.", g), paste0("n.", g)))
  })
  do.call(cbind, parts)
}
stn_surv_wide <- make_wide(data.frame(Station = stn_surv_tab$Station,
                                      n = stn_surv_tab$n, stringsAsFactors = FALSE), 3)

text_page("3. Stations (survey vs processing)",
  c("Survey: 62 stations, 181 samples (1-7 samples/station; 0.2 um filter only).",
    "Processing: 12 stations, 78 samples (1-9 samples/station; station x depth x up to 3 filter fractions).",
    "",
    paste("Shared station names (6):", paste(shared_stn, collapse = ", ")),
    "  -> sampled in BOTH datasets; usable for direct same-site comparison.",
    "",
    "Tables: (3a) processing per-station counts; (3b) survey per-station counts (3-column layout)."))
table_page("3a. Processing - samples per station (12 stations)", stn_proc, fontsize = 11)
table_page("3b. Survey - samples per station (62 stations)", stn_surv_wide, fontsize = 8)

## --- Section 4: ASV counts (survey vs processing) ---
sf <- colnames(sc$asgard_filtered)
pf <- colnames(pc$asgard_filtered_p_hm2)
surv_obs <- sum(colSums(sc$asgard_seqcount   > 0) > 0)
proc_obs <- sum(colSums(pc$asgard_seqcount_p > 0) > 0)
asv_summary <- data.frame(
  Level = c("Full 16S pool (universe)", "Observed in dataset (count>0)",
            "Retained after filtering (clustering)"),
  Survey     = c(ncol(sc$asgard_seqcount),   surv_obs, length(sf)),
  Processing = c(ncol(pc$asgard_seqcount_p), proc_obs, length(pf)))
asv_overlap <- data.frame(
  Set = c("Shared", "Survey-only", "Processing-only", "Union"),
  ASVs = c(length(intersect(sf, pf)), length(setdiff(sf, pf)),
           length(setdiff(pf, sf)), length(union(sf, pf))))

text_page("4. ASV counts (survey vs processing)",
  c("Both datasets are drawn from the same 16S prokaryote ASV universe (3,076 ASVs).",
    "",
    "Observed (>=1 read): Survey 442, Processing 378.",
    "  - Survey observes more, but it also has more samples (181 vs 81) -> partly sampling effort.",
    "  - Processing additionally includes 3 and 20 um (particle-associated) fractions.",
    "",
    "Retained for clustering (max RA > 0.001 AND present in > 2 samples):",
    "  Survey 258, Processing 221; 178 shared (~59% of the 301-ASV union).",
    "  Survey has more unique retained ASVs (80) than processing (43)."))
table_page("4a. ASV counts by level", asv_summary, fontsize = 11)
table_page("4b. Retained-ASV overlap (258 vs 221)", asv_overlap, fontsize = 11)

## --- Section 5: Dataset summary (samples & ASVs) ---
overview <- data.frame(
  Metric = c("Samples (analysis set)", "Filter fractions", "Min 16S reads",
             "Stations", "16S ASV pool", "Observed ASVs (>=1 read)",
             "Retained ASVs (clustering)"),
  Survey = c(as.character(nrow(sc$a_map)), "0.2 um", "1000",
             as.character(length(unique(sc$a_map$station))),
             as.character(ncol(sc$asgard_seqcount)),
             as.character(sum(colSums(sc$asgard_seqcount > 0) > 0)),
             as.character(ncol(sc$asgard_filtered))),
  Processing = c(as.character(nrow(pc$asgard_filtered_p_hm2)), "0.2/3/20 um", "5000",
                 as.character(length(unique(pc$meta_asgard_p2$station))),
                 as.character(ncol(pc$asgard_seqcount_p)),
                 as.character(sum(colSums(pc$asgard_seqcount_p > 0) > 0)),
                 as.character(ncol(pc$asgard_filtered_p_hm2))),
  stringsAsFactors = FALSE)

text_page("5. Dataset summary — samples & ASVs",
  c("Consolidated overview of sample and ASV counts for the two datasets.",
    "",
    "ASV filtering (identical criteria for both): read counts converted to within-sample",
    "relative abundance; ASVs retained if max relative abundance > 0.001 (0.1%) AND",
    "present in > 2 samples.",
    "",
    "Datasets differ only in read threshold (survey >=1000, processing >=5000) and",
    "filter fractions sampled (survey 0.2 um only; processing 0.2/3/20 um)."))
table_page("5. Samples & ASVs — survey vs processing", overview, fontsize = 12)

## --- Section 6: Clustering procedure ---
text_page("6. Clustering procedure (survey & processing)",
  c("Common: community structure was summarised by hierarchical clustering. Bray-Curtis",
    "dissimilarities among samples were computed from the fourth-root-transformed ASV",
    "relative-abundance matrices (survey 181 x 258; processing 78 x 221), and samples were",
    "clustered using Ward's method (ward.D).",
    "",
    "Survey: the dendrogram was cut at k = 10 and organised into a three-level hierarchy with",
    "three top-level divisions (A, B, C). Cluster C1b was further subdivided by re-clustering",
    "the division-C samples (Ward, k = 6), yielding 11 final clusters (A1, A2, B1, B2a, B2b,",
    "C1a, C1b1, C1b2, C2a, C2b1, C2b2).",
    "",
    "Processing: the dendrogram was cut at k = 4; the four clusters were grouped into two",
    "divisions reflecting size-fraction behaviour (free-living = cluster 1; particle-associated",
    "= clusters 2-4).",
    "",
    "Environmental validation (processing only): a top-down tree walk compared, at each node,",
    "the two child groups (each >= 5 samples) on the principal components of the scaled",
    "environmental variables (PCs explaining >= 90% of variance) using Wilcoxon rank-sum tests",
    "with Benjamini-Hochberg correction; a node was split only if at least one PC differed",
    "significantly (adjusted p < 0.05). This yielded a single cluster (no environmentally",
    "significant split at the root)."))

## --- Section 7: Statistical tests inventory (comprehensive) ---
stat_tests <- data.frame(
  Test = c(
    "Kruskal-Wallis", "Wilcoxon rank-sum", "t-test (Welch/Student)",
    "Fisher exact", "Chi-squared", "Tukey HSD",
    "PERMANOVA (adonis2)", "PERMDISP (betadisper)", "Mantel test",
    "db-RDA perm. ANOVA", "Indicator species (IndVal)",
    "BH / FDR correction", "VIF (diagnostic)"),
  Scripts = c(
    "S03, S06, S07, P04", "P03, S07, S18, S20, S25", "S06, S11, S20, S25",
    "S06, S11, S20, P14", "S06, S11, P14", "S10, S17, P09",
    "S10, S17, P09", "S10, S17, P09", "S10",
    "S05, P06", "S12, S13",
    "P03, S06, S07, S10, S17, S18, P09", "S05"),
  Applied.to = c(
    "reads / alpha-div / richness / component-RA vs cluster or depth",
    "tree-walk env-PC; pairwise alpha-div; cluster-specific ASV; strat.; transect",
    "group means (e.g. delta-sigma ~ sea; transect shares)",
    "contingency tables (depth/region x cluster/station)",
    "contingency tables (simulated p-values)",
    "post-hoc for PERMDISP (within-group dispersion)",
    "community composition ~ cluster / environment",
    "homogeneity of multivariate dispersion",
    "community distance ~ environmental distance",
    "constrained ordination significance (overall/margin/axis)",
    "cluster indicator ASVs (16S = S12, 18S = S13)",
    "multiple-comparison correction (pairwise / tree-walk tests)",
    "collinearity check for db-RDA variable selection"),
  stringsAsFactors = FALSE)

text_page("7. Statistical tests — where each was used",
  c("Comprehensive inventory of statistical tests across the survey (S*) and",
    "processing (P*) scripts.",
    "",
    "Notes:",
    "  - Tukey's HSD is used ONLY as the post-hoc for PERMDISP (betadisper).",
    "  - BH (Benjamini-Hochberg) is the multiple-comparison correction; the depth",
    "    Kruskal-Wallis tests (S06, P04) instead use Bonferroni.",
    "  - PERMANOVA, db-RDA and Mantel use 999 permutations; PERMDISP is NOT permutation-",
    "    based (overall = ANOVA on betadisper; pairwise = Tukey's HSD).",
    "  - VIF is a collinearity diagnostic, not a significance test."))
table_page("7. Statistical tests inventory", stat_tests, fontsize = 9)

## --- Section 8: Cluster x depth association (unified Fisher test) ---
fisher_cxd <- function(clust, depth, seedval) {
  d <- data.frame(clust = as.character(clust),
                  depth = factor(as.character(depth), levels = c("surf", "mid", "bottom")),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$clust) & !is.na(d$depth), ]
  ct <- table(d$clust, d$depth)
  set.seed(seedval)
  fe <- fisher.test(ct, simulate.p.value = TRUE, B = 9999)
  cs <- suppressWarnings(chisq.test(ct, simulate.p.value = TRUE, B = 9999))
  list(n = nrow(d), k = nrow(ct), fp = fe$p.value, cp = cs$p.value)
}
S_cxd <- fisher_cxd(sc$a_map$cluster11, sc$a_map$depth_type, 1)
P_cxd <- fisher_cxd(clusnum_p_vec <- pc$clusnum_p[rownames(pc$meta_asgard_p2)],
                    pc$meta_asgard_p2$depth_type, 1)
cxd_tab <- data.frame(
  Dataset = c("Survey", "Processing"),
  n = c(S_cxd$n, P_cxd$n),
  Clusters = c(S_cxd$k, P_cxd$k),
  Fisher_p = c(signif(S_cxd$fp, 3), signif(P_cxd$fp, 3)),
  ChiSq_p  = c(signif(S_cxd$cp, 3), signif(P_cxd$cp, 3)),
  Significant = c(ifelse(S_cxd$fp < 0.05, "yes", "no"),
                  ifelse(P_cxd$fp < 0.05, "yes", "no")),
  stringsAsFactors = FALSE)

text_page("8. Cluster x depth association (unified test)",
  c("Q: Is cluster membership associated with sampling depth (surf/mid/bottom)?",
    "",
    "Both variables are categorical, so a contingency-table test is used (not Kruskal-Wallis).",
    "Unified method applied identically to both datasets: Fisher's exact test with Monte-Carlo",
    "simulation (B = 9,999); chi-squared (Monte-Carlo) reported alongside for reference.",
    "",
    "This harmonises the survey (S06) and processing (P14) tests, which previously differed",
    "slightly (survey used Monte-Carlo Fisher + per-cluster BH; processing used standard Fisher)."))
table_page("8. Cluster x depth — Fisher's exact (Monte-Carlo, B=9999)", cxd_tab, fontsize = 11)

## --- Section 9: Environmental testing of cluster splits (survey vs processing) ---
split_methods <- data.frame(
  Aspect = c("Primary clusters", "Env-node analysis", "Per-PC test",
             "Multiple-test correction", "Purpose", "Outcome"),
  Survey = c("Bray/Ward, cutree k=10 -> 11",
             "per-PC test at Division-C dendrogram nodes (S17)",
             "pairwise PERMANOVA (adonis2) on 1-D PC distance",
             "BH across PC1-PC9 within each node",
             "characterize which env axis separates a fixed split",
             "identifies the driving PC for each C sub-split"),
  Processing = c("Bray/Ward, cutree k=4",
                 "top-down tree walk (P03)",
                 "Wilcoxon rank-sum on PC scores",
                 "BH across retained PCs at each node",
                 "decide whether to split (stopping rule)",
                 "no env-significant split at root -> 1 cluster"),
  stringsAsFactors = FALSE)

text_page("9. Environmental testing of cluster splits — why survey and processing differ",
  c("Both pipelines DEFINE primary clusters identically (Bray-Curtis, Ward's method,",
    "cutree). They differ in the SECONDARY analysis that links cluster splits to",
    "environmental PCs (env PCA of 9 variables, PC1-PC9).",
    "",
    "Survey (S17): clusters are already fixed, so the per-PC test is descriptive -- for",
    "each existing split within Division C it asks which environmental PC best separates",
    "the two sub-clusters. PERMANOVA was used for consistency with the survey's",
    "PERMANOVA/PERMDISP framework (S10, S17) and because Division C is large (97 samples).",
    "",
    "Processing (P03): the tree walk is a decision rule -- it walks the dendrogram top-down",
    "and splits only where the two child groups differ environmentally, using a simple",
    "two-sample Wilcoxon rank-sum on PC scores, well suited to a binary split decision at",
    "the smaller processing sample size (78).",
    "",
    "In short: survey = 'which PC characterizes this (fixed) split?' (PERMANOVA);",
    "processing = 'is this split environmentally justified?' (Wilcoxon stopping rule).",
    "Both correct across PCs with Benjamini-Hochberg. (The divergence also partly reflects",
    "that the two pipelines were developed separately; they can be harmonised on request.)"))
table_page("9. Cluster-split environmental testing — method comparison", split_methods, fontsize = 9)

## --- Section 9: Environmental testing of cluster splits (survey vs processing) ---
split_methods <- data.frame(
  Aspect = c("Primary clusters", "Env test", "PCs used", "Per-PC test",
             "Correction", "Split rule", "Purpose", "Outcome"),
  Survey = c("Bray/Ward, cutree k=10 -> 11",
             "per-PC test at fixed Division-C splits",
             "all 9 (PC1-PC9)",
             "pairwise PERMANOVA (adonis2) + PERMDISP",
             "BH across PC1-PC9 per node",
             "none (all splits tested; descriptive)",
             "which env axis drives a fixed split",
             "driving PC identified per C sub-split"),
  Processing = c("Bray/Ward, cutree k=4 (FL/PA)",
                 "top-down tree walk (recursive)",
                 ">=90% cumulative-variance PCs",
                 "Wilcoxon rank-sum",
                 "BH across retained PCs per node",
                 "split if >=1 PC significant (p<0.05)",
                 "is a split env-justified?",
                 "single group (no env split at root)"),
  stringsAsFactors = FALSE)

text_page("9. Cluster-split environmental testing — survey vs processing",
  c("Both pipelines define primary clusters identically (Bray-Curtis, Ward, cutree). They",
    "differ in how each cluster split was tested against the environmental ordination.",
    "",
    "Survey (S17): clusters are fixed; each split within Division C is tested per-PC by",
    "pairwise PERMANOVA + PERMDISP (all 9 PCs, BH per node). Descriptive - identifies which",
    "environmental axis drives an already-defined split.",
    "",
    "Processing (P03): a top-down tree walk with a stopping rule; per-PC Wilcoxon rank-sum +",
    "BH (PCs up to 90% variance); a node is split only if >=1 PC is significant. A decision",
    "rule - asks whether a split is environmentally justified.",
    "",
    "Why different: survey Division C is large (97 samples) and embedded in a PERMANOVA/",
    "PERMDISP framework, so a multivariate permutation test was natural; processing is",
    "smaller (78) and needed a binary split rule, for which a simple two-sample Wilcoxon is",
    "well suited. Notably, the processing tree walk retained a SINGLE group (no env-significant",
    "split at the root), whereas survey Division C showed environmentally structured sub-splits.",
    "(The divergence also partly reflects independent development; the two can be harmonised.)"))

text_page("9a. Methods text — survey (Division C)",
  c(paste(
    "For the survey dataset, the finer structure within Division C was interpreted against an",
    "environmental ordination. The nine environmental variables (temperature, salinity, DO,",
    "NO3, PO4, silicate, NH4, fluorescence and depth) were standardized and ordinated by PCA,",
    "and all nine components (PC1-PC9) were retained and interpreted from their loadings (e.g.,",
    "PC1: nutrients, salinity and depth; PC2: DO versus temperature). The Division-C samples",
    "were re-clustered (Ward's method), and at each split of the resulting dendrogram (C1 vs C2,",
    "C1a vs C1b, C1b1 vs C1b2, C2a vs C2b, and C2b1 vs C2b2) the two descendant sub-clusters",
    "were compared along each PC separately: for every PC a one-dimensional Euclidean distance",
    "was computed from the sample PC scores and tested by pairwise PERMANOVA (adonis2, 999",
    "permutations), with PERMDISP (betadisper) used to check homogeneity of dispersion. Within",
    "each node, the p-values across the nine PCs were adjusted with the Benjamini-Hochberg",
    "procedure, and the principal component(s) showing the most significant separation were",
    "regarded as the environmental axis driving that split.")))

text_page("9b. Methods text — processing",
  c(paste(
    "For the processing dataset, the community dendrogram was cut at k = 4 and the four",
    "clusters were grouped into two size-fraction divisions: free-living (Cluster 1) and",
    "particle-associated (Clusters 2-4). The environmental separability of the processing",
    "dendrogram was assessed separately using a top-down tree walk: at each node the two child",
    "groups (each >= 5 samples) were compared on the environmental principal components (those",
    "explaining >= 90% of the cumulative variance) using Wilcoxon rank-sum tests with",
    "Benjamini-Hochberg correction, and a node was split only where at least one PC differed",
    "significantly (adjusted p < 0.05). This procedure retained a single group, indicating that",
    "-- unlike the environmentally structured subdivisions within survey Division C -- the",
    "processing clusters were not separated by a significant environmental gradient at the top",
    "level.")))

table_page("9c. Cluster-split env testing — method comparison", split_methods, fontsize = 9)

## >>> ADD FUTURE SECTIONS HERE (text_page / print(ggplot) / table_page) <<<

dev.off()
message("Report written: ", OUT)
message("Pages built. env n(survey)=", nrow(surv), " db-RDA S/P n=", S$n, "/", P$n)
