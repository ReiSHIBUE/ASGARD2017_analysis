### summary_html.R
### ASGARD 2017 — Cluster characteristics summary table (HTML)
###
### Builds an HTML summary of the survey (11) and processing (4) clusters:
###   primary/secondary water mass, bloom stage, Alpha/Gamma/Bacteroidia %,
###   characteristic ASV, and db-RDA significant variables.
###
### Values are the collated results from:
###   - water mass:  output_p/cluster_summary/processing_cluster_watermass_primary_secondary.csv
###                  output/survey/crosstable/cluster_watermass_proportions_11clusters.csv
###   - bloom index: PCA of 9 env variables (S19 / P11 recipe), median per cluster
###   - class %:     per-sample relative abundance by Class, averaged per cluster
###   - char. ASV:   ASV with the largest (in-cluster mean RA - out-of-cluster mean RA)
###   - db-RDA:      overall model significant variables (S05 / P06)
###
### OUTPUT:
###   output/cluster_characteristics/cluster_summary.html
###   output/cluster_characteristics/survey_11cluster_summary.csv
###   output/cluster_characteristics/processing_4cluster_summary.csv
###
### Usage:  Rscript R/summary_html.R

setwd("/Users/shibuerei/Desktop/ASGARD2017_analysis")
outdir <- "output/cluster_characteristics"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- data ----
surv <- data.frame(
  cluster        = c("A1","A2","B1","B2a","B2b","C1a","C1b1","C1b2","C2a","C2b1","C2b2"),
  primary_WM     = c("wSW (78)","cSW (47)","wSW (53)","cSW (35)","cSW (54)","wSW (79)",
                     "wSW (100)","wSW (63)","wCW (67)","wSW/cSW (42/42)","cSW (59)"),
  secondary_WM   = c("AnW (22)","wCW (33)","AnW (47)","wSW (29)","wSW (35)","AnW (21)",
                     "-","AnW (26)","cSW (17)","-","wSW (28)"),
  bloom_index    = c(0.38,0.59,0.23,0.59,0.49,0.16,0.38,0.16,0.67,0.58,0.49),
  bloom_stage    = c("pre","post","pre","post","mid","pre","pre","pre","post","post","mid"),
  Alpha_pct      = c(34,25,19,20,15,16,25,21,33,15,17),
  Gamma_pct      = c(55,68,38,47,43,36,29,35,36,43,42),
  Bacteroidia_pct= c(6,3,36,31,37,40,40,36,25,38,37),
  characteristic_ASV = c("Nitrincolaceae; uncultured","SAR86 clade","Thioglobaceae; SUP05_cluster",
                         "SAR86 clade","Nitrincolaceae; uncultured","Flavobacteriaceae; uncultured",
                         "Flavobacteriaceae; Ulvibacter","Thioglobaceae; SUP05_cluster",
                         "Clade_Ia (SAR11)","Flavobacteriaceae; NS5_marine_group",
                         "Nitrincolaceae; uncultured"),
  dbRDA_significant = "temp, salinity, DO, NO3, FlECO (all)",
  stringsAsFactors = FALSE)

proc <- data.frame(
  cluster        = c("1 (FL)","2 (PA)","3 (PA)","4 (PA)"),
  primary_WM     = c("AnW (48)","wSW (48)","AnW (40)","AnW (47)"),
  secondary_WM   = c("wSW (20)","MWW (13)","MWW (20)","cSW (21)"),
  bloom_index    = c(0.39,0.58,0.44,0.46),
  bloom_stage    = c("pre","post","mid","mid"),
  Alpha_pct      = c(20,21,28,21),
  Gamma_pct      = c(37,28,36,27),
  Bacteroidia_pct= c(39,42,25,48),
  characteristic_ASV = c("Nitrincolaceae; uncultured","Flavobacteriaceae; Ulvibacter",
                         "Moraxellaceae; Psychrobacter","Flavobacteriaceae; uncultured"),
  dbRDA_significant = "salinity only",
  stringsAsFactors = FALSE)

write.csv(surv, file.path(outdir, "survey_11cluster_summary.csv"),     row.names = FALSE)
write.csv(proc, file.path(outdir, "processing_4cluster_summary.csv"),  row.names = FALSE)

# ---- HTML helpers ----
esc <- function(x) { x <- gsub("&","&amp;",x); x <- gsub("<","&lt;",x); gsub(">","&gt;",x) }
stage_cell <- function(s) sprintf('<td class="%s">%s</td>', s, s)

table_html <- function(df) {
  head <- paste0(
    "<tr><th>Cluster</th><th>Primary WM (%)</th><th>Secondary WM (%)</th>",
    "<th>Bloom index</th><th>Bloom stage</th>",
    "<th>Alpha (%)</th><th>Gamma (%)</th><th>Bacteroidia (%)</th>",
    "<th>Characteristic ASV (Family; Genus)</th><th>db-RDA significant</th></tr>")
  rows <- apply(df, 1, function(r) {
    paste0("<tr><td>", esc(r["cluster"]), "</td><td>", esc(r["primary_WM"]),
           "</td><td>", ifelse(r["secondary_WM"]=="-","&ndash;",esc(r["secondary_WM"])),
           "</td><td class=\"num\">", r["bloom_index"], "</td>", stage_cell(r["bloom_stage"]),
           "<td class=\"num\">", r["Alpha_pct"], "</td><td class=\"num\">", r["Gamma_pct"],
           "</td><td class=\"num\">", r["Bacteroidia_pct"], "</td><td>", esc(r["characteristic_ASV"]),
           "</td><td>", esc(r["dbRDA_significant"]), "</td></tr>")
  })
  paste0("<table>\n", head, "\n", paste(rows, collapse = "\n"), "\n</table>")
}

css <- '
  body { font-family: -apple-system, "Helvetica Neue", Arial, sans-serif; margin: 24px; color: #222; }
  h1 { font-size: 20px; } h2 { font-size: 17px; margin-top: 28px; }
  table { border-collapse: collapse; width: 100%; font-size: 13px; margin-top: 8px; }
  th, td { border: 1px solid #ccc; padding: 6px 8px; text-align: left; vertical-align: top; }
  th { background: #123456; color: #fff; font-weight: bold; }
  tr:nth-child(even) { background: #f5f7fa; }
  td.num { text-align: right; }
  .pre { color: #1565C0; font-weight: bold; } .post { color: #C62828; font-weight: bold; }
  .mid { color: #6A1B9A; font-weight: bold; }
  .note { font-size: 12px; color: #555; margin-top: 10px; line-height: 1.5; }'

note <- '
<div class="note">
<b>Notes on methods:</b><br>
&bull; <b>Bloom stage</b>: from the median bloom index (= 1 &minus; rescaled PC1 of nine standardized environmental variables). pre &lt; 0.45, mid 0.45&ndash;0.55, post &gt; 0.55.<br>
&bull; <b>Alpha / Gamma / Bacteroidia (%)</b>: per-sample relative abundance of each class within the filtered ASV set, averaged across the samples of each cluster.<br>
&bull; <b>Characteristic ASV</b>: the ASV with the largest difference between its mean relative abundance inside vs outside the cluster.<br>
&bull; <b>db-RDA significant</b>: variables significant in the overall distance-based RDA (not cluster-specific).<br>
&bull; WM = water mass (AnW = Anadyr Water; wSW/cSW = warm/cold Shelf Water; wCW = warm Coastal Water; MWW = Modified Winter Water).
</div>'

html <- paste0(
  "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n",
  "<title>ASGARD 2017 — Cluster characteristics summary</title>\n<style>", css, "</style>\n</head>\n<body>\n",
  "<h1>ASGARD 2017 — Cluster characteristics summary</h1>\n",
  "<h2>Survey (11 clusters)</h2>\n", table_html(surv), "\n",
  "<h2>Processing (4 clusters)</h2>\n", table_html(proc), "\n",
  note, "\n</body>\n</html>\n")

writeLines(html, file.path(outdir, "cluster_summary.html"))
cat("DONE:\n  ", file.path(outdir, "cluster_summary.html"), "\n",
    "  ", file.path(outdir, "survey_11cluster_summary.csv"), "\n",
    "  ", file.path(outdir, "processing_4cluster_summary.csv"), "\n", sep = "")
