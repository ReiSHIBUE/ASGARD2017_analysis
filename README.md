# ASGARD2017 Processing Site Analysis

R pipeline for 16S/18S amplicon analysis of microbial communities sampled across processing stations during the ASGARD 2017 cruise. Analyzes size-fractionated samples (0.2, 3, and 20 µm filters) to characterize free-living vs. particle-associated prokaryotic and eukaryotic communities.

## Overview

The pipeline takes DADA2-processed ASV tables and sample metadata as input and produces:

- Ternary diagrams showing ASV distributions across filter size fractions
- Hierarchical clustering heatmaps (16S prokaryote and 18S eukaryote)
- Geographic maps of sample clusters and particle-associated taxa
- Beta diversity ordination (PCoA) and constrained ordination (db-RDA)
- Boxplots of environmental variables by microbial community cluster

## Requirements

**R packages:**
```r
install.packages(c("tidyverse", "vegan", "gplots", "viridis", "ggtern",
                   "ggmap", "ggrepel", "pheatmap", "RColorBrewer",
                   "janitor", "scales", "here"))
```

Maps require a [Stadia Maps](https://stadiamaps.com/) API key:
```r
ggmap::register_stadiamaps("YOUR_API_KEY")
```

## Usage

Set your working directory to the repo root and source the master script:

```r
setwd("~/Desktop/ASGARD2017_analysis")
source("R/00_run_all.R")
```

Output PDFs are written to `output/`.

To re-run individual scripts, earlier scripts must have already been sourced in the same R session (objects are passed via the global environment):

```r
library(here)
source(here("R/00_setup.R"))
source(here("R/01_data_prep.R"))
source(here("R/03_heatmaps_16S.R"))  # skipping 02 would fail — objects missing
```

## Pipeline Scripts

| Script | Description | Output |
|--------|-------------|--------|
| `00_setup.R` | Loads raw RDS files, subsets to 16S prokaryote ASVs, filters to ≥1000-read samples, builds taxonomy labels | — |
| `01_data_prep.R` | Filters to ASGARD processing stations, applies abundance cutoffs, creates 81- and 78-sample subsets | — |
| `02_ternary_plots.R` | Plots ASV distributions across filter fractions; identifies particle-associated ASVs (`zero_cols`) | ternary plots (displayed) |
| `03_heatmaps_16S.R` | Bray-Curtis / Ward hierarchical clustering heatmaps; assigns samples to 4 clusters (`clusnum_p`) | `ASGARD_hm_processing_5000over.pdf` |
| `04_maps.R` | Geographic maps of sample clusters and per-ESV presence/absence | `processing_map.pdf`, `maps_pa.pdf` |
| `05_beta_diversity_pcoa.R` | PCoA (Bray-Curtis, Jaccard, Euclidean); boxplots of all variables by cluster | `asgard_boxplots_processing.pdf` |
| `06_dbrda.R` | Distance-based RDA constrained by temp, salinity, DO, NO₃, fluorescence | — |
| `07_18S_heatmaps.R` | 18S eukaryote class- and phylum-level heatmaps and boxplots by cluster | `ASGARD_hm_processing_18S.pdf`, boxplot PDFs |
| `08_esv_heatmap.R` | ESV-level relative abundance heatmap (74 samples shared with 18S) | `ASGARD_hm_processing_esv_relabund.pdf` |

## Data

Raw data files are in `data/raw/` (RDS format, from DADA2 + SILVA v132 taxonomy):

| File | Description |
|------|-------------|
| `seqtab_filt.rds` | ASV count matrix (2193 samples × 18535 ASVs) |
| `table_list.rds` | Marker type per ASV (`16S_prokaryote`, `18S`, etc.) |
| `meta_denovo_2.RDS` | Sample metadata (station, depth, temperature, salinity, nutrients, etc.) |
| `names_list.rds` | Full SILVA taxonomic strings |
| `bootout_edit.rds` | Bootstrap confidence values per taxonomic rank |

Pre-computed 18S relative abundance tables (TSV) used by scripts 07–08 are in `data/`.

## Repository Structure

```
R/
  00_setup.R           # environment setup
  00_run_all.R         # master run script
  01_data_prep.R       # filtering and subsetting
  02_ternary_plots.R   # size-fraction ternary diagrams
  03_heatmaps_16S.R    # 16S hierarchical clustering
  04_maps.R            # geographic maps
  05_beta_diversity_pcoa.R
  06_dbrda.R
  07_18S_heatmaps.R
  08_esv_heatmap.R
data/
  raw/                 # input RDS files
  *.tsv                # 18S relative abundance tables
output/                # generated PDFs (git-ignored)
archive/               # original monolithic script (ASGARD2017_p.R)
```
