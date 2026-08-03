# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Pipeline

**Run both pipelines** (processing + survey, corrected xlsx metadata):
```r
source("R/00_run_all.R")   # from the repository root
```

**Run only the processing pipeline:**
```r
source("R/00_run_processing.R")
```

**Run only the survey pipeline:**
```r
source("R/00_run_survey.R")
```

**Run a single script** (objects from earlier scripts must already be in the session):
```r
library(here)
source(here("R/P03_heatmaps_16S.R"))
```

**Run only setup + one script** (e.g., to re-run data prep):
```r
source(here("R/00_setup.R"))
source(here("R/P01_data_prep.R"))
```

Output PDFs are written to `output/` (survey, shared) and `output_p/` (processing).
`R/00_dirs.R` holds the canonical directory list and is sourced by `00_setup.R`,
so every entry point creates the directories before any script writes to them.
The runners record per-script errors and print a failure summary rather than
always reporting success.

## Required R Packages

**Both pipelines:** `tidyverse`, `vegan`, `gplots`, `viridis`, `ggtern`, `ggmap`, `ggrepel`, `pheatmap`, `RColorBrewer`, `janitor`, `scales`, `here`, `readxl`, `fastcluster`

**Survey pipeline only:** `igraph`, `Rtsne` (optional tSNE layout), `indicspecies`, `sf`, `ggnewscale`, `gridExtra`, `waffle` (GitHub: `hrbrmstr/waffle`)

For maps (P04, S06, S16, S17, S23): `ggmap` requires a Stadia Maps API key stored in `.Renviron`:
```
STADIA_MAPS_KEY=your_key_here
```

## Architecture

The pipeline is **session-state based** — each script sources into a shared R session and leaves named objects in the global environment for the next script to consume. There is no package or module system; all inter-script communication is via named R objects.

### Processing pipeline (P01–P14)

| Script | Inputs (objects) | Key outputs |
|--------|-----------------|-------------|
| `00_setup.R` | RDS files in `data/raw/` | `seqtab_16Smat`, `seqtab_16Sprop`, `meta_denovo_2`, `shorternames` |
| `P01_data_prep.R` | setup outputs | `asgard_processing` (81×269), `asgard_processing2` (78×267), `asgard_filtered_p_hm2` (78×221 matrix), `meta_asgard_p2` |
| `P02_ternary_plots.R` | data_prep outputs | `ternary_prop_color2`, `asv_rgb2`, `sample_rgb2`, `zero_cols` (ASVs absent in 0.2 µm) |
| `P03_heatmaps_16S.R` | ternary outputs | `h3` (heatmap object), `clusnum_p` (cluster vector 1–4, length 78), `sample_rgb3` |
| `P04_maps.R` | cluster outputs + `zero_cols` | `mapz` (Stadia basemap) |
| `P05_beta_diversity_pcoa.R` | cluster outputs | `asgard_pcoa_df_p` (78×52) |
| `P06_dbrda.R` | PCoA outputs | `asgard_dbrda_model_p`, ANOVA results, `asgard_dbrda_merged_p` (75 samples) |
| `P07_18S_heatmaps.R` | cluster + TSV files | `asgard_euk_class_hm_filtered` (74×66), `sample_rgb4`, `h3_74` |
| `P08_esv_heatmap.R` | 18S outputs + TSV file | `h11` |
| `P09`–`P14`, `TS_diagram_p.R` | cluster outputs | cluster/water-mass/station summaries, env boxplots, bubble and ternary figures, T–S diagram |

### Survey pipeline (S01–S25)

| Script | Inputs (objects) | Key outputs |
|--------|-----------------|-------------|
| `00_setup.R` | (shared with processing) | (same as above) |
| `S01_data_prep.R` | setup outputs | `asgard_filtered` (181×258), `meta_asgard` (181×45), `asgard_frtprop` (181×258) |
| `S02_heatmaps_16S.R` | data_prep outputs | `h1`/`h2` (heatmaps), `clusnum10` (10 row clusters), `clusnum11` (11 named clusters A1–C2b2), `colclusnum` (6 col clusters), `cc11`, `hier_levels_11` |
| `S03_seq_depth.R` | data_prep outputs | seq depth plots, rarefaction curves |
| `S04_beta_diversity_pcoa.R` | cluster outputs | `asgard_pcoa_df`, distance matrices |
| `S05_dbrda.R` | PCoA outputs | `asgard_dbrda_model`, ANOVA results |
| `S06_maps.R` | cluster outputs | `mapz_survey` (Stadia basemap) |
| `S07_alpha_diversity.R` | data_prep outputs | Shannon, Simpson, Chao1 + Kruskal-Wallis |
| `S08_taxonomy.R` | cluster outputs | stacked bar charts, waffle charts |
| `S09_network.R` | cluster outputs | `asgard_network` (igraph, Spearman \|r\| > 0.6) |
| `S10_permanova.R` | distance matrices | adonis2, betadisper, Mantel test results |
| `S11`–`S25`, `TS_diagram.R`, `heatmap_watermass.R` | cluster outputs | cross-tables, IndVal, env PCA, boxplots, stratification, bubble/waffle figures, T–S diagram |

### Key design decisions

- **Two pipelines share `00_setup.R`**: processing stations (P, 3 filter sizes) use `P01–P14`; survey stations (S, 0.2 µm only) use `S01–S25`. No object name collisions (processing uses `_p` suffixes).
- **Two sample subsets in processing pipeline**: 81-sample set (≥1000 reads) and 78-sample set (≥5000 reads). Most downstream analysis uses the 78-sample set.
- **ASV filtering**: max relative abundance > 0.001 AND present in >2 samples.
- **Fourth-root transformation** (`^0.25`) is applied to proportion matrices before distance computation and heatmap display. Note `asgard_frtprop` (survey) is *already* transformed in `S01` — do not raise it to `0.25` again.
- **Clustering method** is `ward.D` (not `ward.D2`) everywhere.
- **`set.seed(42)`** is called in every script that runs permutation tests or a stochastic layout, so p-values are stable across runs.
- **`zero_cols`**: set of ASV names absent in the 0.2 µm (free-living) fraction — identifies particle-associated taxa. Computed in `P02`, used in `P03` and `P04`.
- **`clusnum11`**: the survey clustering used throughout. `S02` cuts the row dendrogram at k = 10 and then splits `C1b` into `C1b1`/`C1b2`, giving the 11 named clusters (A1, A2, B1, B2a, B2b, C1a, C1b1, C1b2, C2a, C2b1, C2b2). Column clusters: 6 ASV assemblages.
- **`clusnum_p`**: integer vector (values 1–4) assigning each of the 78 samples to a Bray-Curtis / Ward hierarchical cluster. Produced in `P03` and propagated to all later scripts.
- **`P10_water_mass.R`** needs both pipelines: its station-comparison section uses the survey object `meta_asgard` and is skipped (with a message) when only the processing pipeline has been run.
- **`h3$rowDendrogram`**: the dendrogram from the primary 16S heatmap is re-used as `Rowv` in later heatmaps to keep sample ordering consistent.
- **74-sample subset** in P07–P08: 4 samples dropped due to NA values in the 18S data.

### Raw data (`data/raw/`)

| File | Contents |
|------|----------|
| `seqtab_filt.rds` | Full ASV count matrix (2193 samples × 18535 ASVs, all markers) |
| `table_list.rds` | Marker assignment per ASV column (`16S_prokaryote`, `18S`, etc.) |
| `meta_denovo_2.RDS` | Sample metadata (1163 samples × 45 cols: station, depth, temp, sal, nutrients…) |
| `names_list.rds` | Full taxonomic name strings (SILVA v132, length 18535) |
| `bootout_edit.rds` | Bootstrap confidence values per taxonomic rank (18535 × 23) |

Tab-separated files in `data/` (`class_relabund_by_station_depth.tsv`, `phylum_relabund_by_station_depth.tsv`, `esv_relabund_by_station_depth.tsv`) are pre-computed 18S relative abundance tables read by P07 and P08.
