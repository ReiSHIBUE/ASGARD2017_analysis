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
                   "janitor", "scales", "here", "igraph", "Rtsne"))
```

🗺️ Maps require a [Stadia Maps](https://stadiamaps.com/) API key stored in `.Renviron`:
```
STADIA_MAPS_KEY=your_key_here
```

## 🚀 Usage

Run both pipelines:
```r
setwd("~/Desktop/ASGARD2017_analysis")
source("R/00_run_all.R")
```

Run only processing or survey:
```r
source("R/00_run_processing.R")
source("R/00_run_survey.R")
```

📄 Output PDFs are written to `output/` and `output/survey/`.

To re-run individual scripts, earlier scripts must have already been sourced in the same R session:
```r
library(here)
source(here("R/00_setup.R"))
source(here("R/P01_data_prep.R"))
source(here("R/P03_heatmaps_16S.R"))  # skipping P02 would fail — objects missing
```

## 🧬 Processing Pipeline (P01–P08)

| Script | Description | Output |
|--------|-------------|--------|
| `00_setup.R` | 📥 Loads raw RDS files, subsets to 16S prokaryote ASVs, builds taxonomy labels | — |
| `P01_data_prep.R` | 🔧 Filters to processing stations, applies abundance cutoffs, creates 81- and 78-sample subsets | — |
| `P02_ternary_plots.R` | 🔺 ASV distributions across filter fractions; identifies particle-associated ASVs | `output/ternary/` |
| `P03_heatmaps_16S.R` | 🔥 Bray-Curtis / Ward clustering heatmaps; assigns samples to 4 clusters | `output/heatmaps/` |
| `P04_maps.R` | 🗺️ Geographic maps of clusters and per-ESV presence/absence | `output/maps/` |
| `P05_beta_diversity_pcoa.R` | 📐 PCoA ordination + boxplots by cluster | `output/beta_diversity/` |
| `P06_dbrda.R` | 🌡️ db-RDA constrained by temp, salinity, DO, NO₃, fluorescence | `output/dbrda/` |
| `P07_18S_heatmaps.R` | 🦠 18S eukaryote class/phylum heatmaps and boxplots | `output/heatmaps/` |
| `P08_esv_heatmap.R` | 🧫 ESV-level relative abundance heatmap (74 samples) | `output/heatmaps/` |

## 🔬 Survey Pipeline (S01–S10)

| Script | Description | Output |
|--------|-------------|--------|
| `S01_data_prep.R` | 🔧 Filters to survey stations (181 samples, 0.2 µm only) | — |
| `S02_heatmaps_16S.R` | 🔥 Clustering heatmap (5 row clusters, 8 col clusters) | `output/survey/heatmaps/` |
| `S03_seq_depth.R` | 📊 Sequencing depth dot plot + rarefaction curves | `output/survey/alpha_diversity/` |
| `S04_beta_diversity_pcoa.R` | 📐 PCoA ordination (Bray/Jaccard/Euclidean) + boxplots | `output/survey/beta_diversity/` |
| `S05_dbrda.R` | 🌡️ db-RDA constrained by temp, salinity, DO, NO₃, fluorescence | `output/survey/dbrda/` |
| `S06_maps.R` | 🗺️ Geographic maps of survey station clusters | `output/survey/maps/` |
| `S07_alpha_diversity.R` | 🌿 Shannon, Simpson, Chao1 + Kruskal-Wallis tests | `output/survey/alpha_diversity/` |
| `S08_taxonomy.R` | 🧬 Stacked bar charts by Order + waffle charts by cluster | `output/survey/taxonomy/` |
| `S09_network.R` | 🕸️ Co-occurrence network (Spearman \|r\| > 0.6) + tSNE layout | `output/survey/network/` |
| `S10_permanova.R` | 📏 PERMANOVA, PERMDISP, Mantel tests | `output/survey/beta_diversity/` |

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
  00_setup.R               # 📥 shared environment setup
  00_run_all.R             # 🚀 master run script (both pipelines)
  00_run_processing.R      # 🚀 processing pipeline runner
  00_run_survey.R          # 🚀 survey pipeline runner
  P01–P08_*.R              # 🔬 processing pipeline scripts
  S01–S10_*.R              # 🔬 survey pipeline scripts
data/
  raw/                     # 💾 input RDS files
  *.tsv                    # 🦠 18S relative abundance tables
output/                    # 📄 processing PDFs (git-ignored)
  survey/                  # 📄 survey PDFs (git-ignored)
```
