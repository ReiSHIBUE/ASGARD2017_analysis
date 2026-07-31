### install_packages.R
### Package installation for the ASGARD2017 analysis container.
### コンテナ用パッケージインストールスクリプト
###
### Pinned to the Posit Package Manager snapshot set in the Dockerfile, so the
### same image build produces the same package versions.

options(
  repos = c(CRAN = Sys.getenv("CRAN_REPO", "https://cloud.r-project.org")),
  Ncpus = max(1L, parallel::detectCores())
)

cran_pkgs <- c(
  # core
  "here", "tidyverse", "readxl", "janitor", "scales",
  # ecology / stats
  "vegan", "fastcluster", "indicspecies", "Rtsne", "igraph",
  # plotting
  "gplots", "viridis", "viridisLite", "RColorBrewer", "ggrepel", "ggtern",
  "ggnewscale", "gridExtra", "pheatmap", "png",
  # maps
  "ggmap", "sf", "rnaturalearth", "rnaturalearthdata",
  # build helper for the GitHub-only package below
  "remotes"
)

install.packages(cran_pkgs)

# waffle is not on CRAN (S08_taxonomy.R, S22_dbo3_waffle.R). Installed from a
# source tarball rather than install_github() so the build does not depend on
# an authenticated GitHub API call.
# waffle は CRAN に無いためソースtarballから導入
remotes::install_url(
  "https://github.com/hrbrmstr/waffle/archive/refs/heads/master.tar.gz",
  upgrade = "never"
)

# Fail the build loudly if anything is missing, rather than at pipeline runtime.
missing <- setdiff(
  c(cran_pkgs, "waffle"),
  rownames(installed.packages())
)
if (length(missing)) {
  stop("packages failed to install: ", paste(missing, collapse = ", "))
}

message("All ", length(cran_pkgs) + 1L, " packages installed.")
