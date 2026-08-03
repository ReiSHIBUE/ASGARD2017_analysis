# 🌊 ASGARD2017 Analysis

R pipeline for 16S/18S amplicon analysis of microbial communities sampled during the ASGARD 2017 cruise. Two parallel pipelines analyze **processing stations** (3 filter size fractions) and **survey stations** (0.2 µm only).

## 🔭 Overview

The pipeline takes DADA2-processed ASV tables and sample metadata as input and produces:

- 🔺 Ternary diagrams showing ASV distributions across filter size fractions
- 🔥 Hierarchical clustering heatmaps (16S prokaryote and 18S eukaryote)
- 🗺️ Geographic maps of sample clusters and particle-associated taxa
- 📐 Beta diversity ordination (PCoA) and constrained ordination (db-RDA)
- 📦 Boxplots of environmental variables by microbial community cluster
- 🌿 Alpha diversity (Shannon, Simpson, Chao1) and taxonomy visualizations
- 🕸️ ASV co-occurrence networks (Spearman correlation)

## 📦 Requirements

**R packages:**
```r
install.packages(c("tidyverse", "vegan", "gplots", "viridis", "ggtern",
                   "ggmap", "ggrepel", "pheatmap", "RColorBrewer",
                   "janitor", "scales", "here", "igraph", "Rtsne",
                   "readxl", "indicspecies", "sf", "ggnewscale",
                   "gridExtra", "fastcluster"))
remotes::install_github("hrbrmstr/waffle")   # not on CRAN
```

🗺️ Maps require a [Stadia Maps](https://stadiamaps.com/) API key stored in `.Renviron`:
```
STADIA_MAPS_KEY=your_key_here
```
Without a key, `P04`, `S06`, `S16`, `S17` and `S23` fail; the rest of the
pipeline still runs.

🎲 Permutation-based scripts (`S06`, `S09`–`S11`, `S20`, `P09`, `S12`, `S13`,
`S17`) call `set.seed(42)`, so reported p-values are reproducible across runs.

## 🚀 Usage

Run both pipelines (uses the corrected xlsx metadata):
```r
source("R/00_run_all.R")   # from the repository root
```

Run only processing or survey (uses the RDS metadata):
```r
source("R/00_run_processing.R")
source("R/00_run_survey.R")
```

📄 Output PDFs are written to `output/` (survey) and `output_p/` (processing).
The runners create every output subdirectory before running, and report which
scripts failed instead of always printing "complete".

To re-run individual scripts, earlier scripts must have already been sourced in the same R session:
```r
library(here)
source(here("R/00_setup.R"))
source(here("R/P01_data_prep.R"))
source(here("R/P03_heatmaps_16S.R"))  # skipping P02 would fail — objects missing
```

## 🧬 Processing Pipeline (P01–P14)

| Script | Description | Output |
|--------|-------------|--------|
| `00_setup.R` | 📥 Loads raw RDS files, subsets to 16S prokaryote ASVs, builds taxonomy labels | — |
| `P01_data_prep.R` | 🔧 Filters to processing stations, applies abundance cutoffs, creates 81- and 78-sample subsets | — |
| `P02_ternary_plots.R` | 🔺 ASV distributions across filter fractions; identifies particle-associated ASVs | `output_p/ternary/` |
| `P03_heatmaps_16S.R` | 🔥 Bray-Curtis / Ward.D clustering heatmaps; assigns samples to 4 clusters | `output_p/heatmaps/` |
| `P04_maps.R` | 🗺️ Geographic maps of clusters and per-ESV presence/absence | `output_p/maps/` |
| `P05_beta_diversity_pcoa.R` | 📐 PCoA ordination + boxplots by cluster | `output_p/beta_diversity/` |
| `P06_dbrda.R` | 🌡️ db-RDA constrained by temp, salinity, DO, NO₃, fluorescence | `output_p/dbrda/` |
| `P07_18S_heatmaps.R` | 🦠 18S eukaryote class/phylum heatmaps and boxplots | `output_p/heatmaps/` |
| `P08_esv_heatmap.R` | 🧫 ESV-level relative abundance heatmap (74 samples) | `output_p/heatmaps/` |
| `P09_cluster_summary.R` | 📋 Per-cluster geographic + environmental summary | `output_p/cluster_summary/` |
| `P10_water_mass.R` | 🌊 Water-mass classification (Danielson 2020); station comparison vs survey | `output_p/cluster_summary/` |
| `P11_env_boxplots.R` | 📦 Environmental boxplots + bloom index (4 clusters) | `output_p/env_boxplots/` |
| `P12_genus_bubble.R` | 🫧 Genus bubble plots (4 clusters) | `output_p/bubble/` |
| `P13_ternary_v2.R` | 🔺 Improved ternary plots (class / cluster / RGB) | `output_p/ternary/` |
| `P14_station_region_summary.R` | 🗺️ Station / region / DBO3 summaries and tests | `output_p/cluster_summary/` |
| `TS_diagram_p.R` | 🌡️ T–S diagram with water-mass boxes (4 clusters) | `output_p/TS/` |

## 🔬 Survey Pipeline (S01–S25)

| Script | Description | Output |
|--------|-------------|--------|
| `S01_data_prep.R` | 🔧 Filters to survey stations (181 samples, 0.2 µm only) | — |
| `S02_heatmaps_16S.R` | 🔥 Clustering heatmap: 10 row clusters cut from the dendrogram, split into the 11 named clusters (A1–C2b2), and 6 ASV column clusters | `output/survey/heatmaps/` |
| `S03_seq_depth.R` | 📊 Sequencing depth dot plot + rarefaction curves | `output/survey/alpha_diversity/` |
| `S04_beta_diversity_pcoa.R` | 📐 PCoA ordination (Bray/Jaccard/Euclidean) + boxplots | `output/survey/beta_diversity/` |
| `S05_dbrda.R` | 🌡️ db-RDA constrained by temp, salinity, DO, NO₃, fluorescence | `output/survey/dbrda/` |
| `S06_maps.R` | 🗺️ Geographic maps of survey station clusters | `output/survey/maps/` |
| `S07_alpha_diversity.R` | 🌿 Shannon, Simpson, Chao1 + Kruskal-Wallis tests | `output/survey/alpha_diversity/` |
| `S08_taxonomy.R` | 🧬 Stacked bar charts by Order + waffle charts by cluster | `output/survey/taxonomy/` |
| `S09_network.R` | 🕸️ Co-occurrence network (Spearman \|r\| > 0.6) + tSNE layout | `output/survey/network/` |
| `S10_permanova.R` | 📏 PERMANOVA, PERMDISP, Mantel tests | `output/survey/beta_diversity/` |
| `S11_crosstable.R` | 🔀 Cluster × water-mass cross-tables and Fisher tests | `output/survey/crosstable/` |
| `S12_indval.R` / `S13_indval_18S.R` | 🎯 Indicator species analysis (16S / 18S) | `output/survey/` |
| `S14_sampling_period.R` | 📅 Sampling period per cluster | `output/survey/` |
| `S15_wSW_misclassification.R` | 🧪 wSW-classified samples with high NO₃ | `output/survey/env_table/` |
| `S16_NO3_summary.R` | 🧂 NO₃ summary + map | `output/survey/NO3 conc/` |
| `S17_pca_env.R` | 📉 Environmental PCA biplots + PC-score maps | `output/survey/beta_diversity/`, `output/survey/maps/` |
| `S18_cluster_specific_ASVs.R` | 🧬 Cluster-specific ASVs (Wilcoxon) | `output/survey/IndVal/` |
| `S19_env_boxplots.R` | 📦 Environmental boxplots + bloom index (11 clusters) | `output/survey/beta_diversity/` |
| `S20_stratification.R` | 🌊 Water-column stratification profiles | `output/survey/stratification/` |
| `S21_genus_bubble.R` | 🫧 Top-100 genus bubble plot | `output/survey/bubble/` |
| `S22_dbo3_waffle.R` | 🧇 Class composition at DBO3 stations | `output/survey/waffle/` |
| `S23_dbo3_envmap.R` | 🗺️ Chl-a and NO₃ maps at DBO3 stations | `output/survey/maps/` |
| `S25_transect_compare.R` | 🔬 Class composition compared across transects | `output/survey/transect_compare/` |
| `TS_diagram.R` | 🌡️ T–S diagram with water-mass boxes (11 clusters) | `output/TS plot/` |
| `heatmap_watermass.R` | 🔥 Heatmaps coloured by water mass | `output/heatmaps/` |

## 💾 Data

Raw data files are in `data/raw/` (RDS format, from DADA2 + SILVA v132 taxonomy):

| File | Description |
|------|-------------|
| `seqtab_filt.rds` | 🧮 ASV count matrix (2193 samples × 18535 ASVs) |
| `table_list.rds` | 🏷️ Marker type per ASV (`16S_prokaryote`, `18S`, etc.) |
| `meta_denovo_2.RDS` | 📋 Sample metadata (station, depth, temperature, salinity, nutrients, etc.) |
| `names_list.rds` | 🌿 Full SILVA taxonomic strings |
| `bootout_edit.rds` | 📊 Bootstrap confidence values per taxonomic rank |

Pre-computed 18S relative abundance tables (TSV) used by P07–P08 are in `data/`.

## 🗂️ Repository Structure

```
R/
  00_setup.R               # 📥 shared environment setup (also creates output dirs)
  00_dirs.R                # 📁 canonical list of output directories
  00_run_all.R             # 🚀 master run script (both pipelines)
  00_run_processing.R      # 🚀 processing pipeline runner
  00_run_survey.R          # 🚀 survey pipeline runner
  P01–P14_*.R              # 🔬 processing pipeline scripts
  S01–S25_*.R              # 🔬 survey pipeline scripts
  TS_diagram*.R, heatmap_watermass.R, comparison_report.R, summary_html.R
data/
  raw/                     # 💾 input RDS files
  *.tsv                    # 🦠 18S relative abundance tables
output/                    # 📄 survey PDFs + shared figures
  survey/                  # 📄 survey PDFs
output_p/                  # 📄 processing PDFs
```
