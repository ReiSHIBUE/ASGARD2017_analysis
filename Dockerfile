# ASGARD2017 analysis — reproducible runtime
#
# Build:  docker build -t asgard2017 .
# Run:    docker run --rm -v "$PWD/output:/work/output" \
#                    -v "$PWD/output_p:/work/output_p" \
#                    -e STADIA_MAPS_KEY=... asgard2017
#
# The image contains R and every package the pipeline needs, but no repository
# data — mount or copy the repo in at run time (the GitHub Actions workflow
# mounts the checkout at /work).

FROM rocker/r-ver:4.5.1

# Pin CRAN to a dated Posit Package Manager snapshot so rebuilds resolve the
# same package versions. Linux binaries make this install minutes, not hours.
ENV CRAN_REPO=https://p3m.dev/cran/__linux__/noble/2026-07-01
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# System libraries: sf/units/proj (maps), curl/ssl (ggmap, remotes),
# fontconfig/freetype (ggplot text), libglpk (igraph), libpng (png).
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libgdal-dev \
      libgeos-dev \
      libproj-dev \
      libudunits2-dev \
      libglpk-dev \
      libpng-dev \
      libfontconfig1-dev \
      libfreetype6-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libtiff-dev \
      libjpeg-dev \
      git \
    && rm -rf /var/lib/apt/lists/*

COPY docker/install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm /tmp/install_packages.R

WORKDIR /work

# Default: run both pipelines and fail if a script that should have worked did not.
CMD ["Rscript", "docker/run_pipeline.R"]
