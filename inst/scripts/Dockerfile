# =============================================================================
# Reproducible environment for cyRAVEN
# =============================================================================
#
# WHY A PINNED BASE, NOT r-base:latest
#   The analysis calls uwot (stochastic embedding) and derives every gate from a
#   kernel density estimate. Both are sensitive to their own version and to the
#   BLAS underneath them: the same --seed on a different uwot gives a visually
#   similar but numerically different embedding, and a different density
#   implementation can move a gate threshold, which changes the population
#   frequencies you would publish. `latest` would silently re-roll all of that on
#   a rebuild months from now.
#
#   rocker/r-ver:4.4.3 pins R itself AND sets a dated CRAN snapshot (Posit Package
#   Manager), so install.packages() resolves to the versions current at that R
#   release rather than today's. Bioconductor is pinned separately below, because
#   flowCore comes from Bioconductor and follows its own release cycle.
#
# VERSIONS THIS IMAGE REPRODUCES (the environment the pipeline was validated in):
#   R 4.4.3 . flowCore 2.18.0 . uwot 0.2.4 . data.table 1.18.4
#   ggplot2 4.0.3 . patchwork 1.3.2 . viridisLite 0.4.2 . scales 1.4.0
#   yaml 2.3.12 . optparse 1.8.2 . hexbin 1.28.5
#
# BUILD (from the package root -- the directory holding DESCRIPTION):
#   docker build -f inst/scripts/Dockerfile -t cyraven:1.0.0 .
#
# RUN -- mount the data directory and write results back out. Nothing is baked
# into the image: no FCS file, no patient table, no cohort name.
#   docker run --rm \
#     -v "$PWD:/data:ro" \
#     -v "$PWD/results:/results" \
#     cyraven:1.0.0 \
#     --dir /data --recursive \
#     --config /data/config_cohorts.yaml \
#     --sample-map /data/sample_map.csv \
#     --patient-table /data/patient_table.csv \
#     --outdir /results
#
#   On Windows PowerShell replace "$PWD" with "${PWD}".
#
# ITERATING WITHOUT REBUILDING THE IMAGE
#   The command above runs the copy of cyRAVEN installed at build time. Set
#   CYRAVEN_SOURCE to a mounted checkout and the entrypoint loads that source
#   tree instead, so routine edits need no rebuild:
#     -v "$PWD:/src:ro" -e CYRAVEN_SOURCE=/src
#   See entrypoint.sh. A rebuild is still needed if the source starts requiring a
#   package this image does not have.
# -----------------------------------------------------------------------------

FROM rocker/r-ver:4.4.3

LABEL org.opencontainers.image.title="cyRAVEN" \
      org.opencontainers.image.description="Gated flow-cytometry UMAP pipeline (FCS -> gates -> shared UMAP -> figures/tables)" \
      org.opencontainers.image.source="https://github.com/bhagesh-h/cyRAVEN" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

# Bioconductor release matching R 4.4.x. Pinning this is what keeps flowCore from
# drifting on a rebuild.
ENV BIOC_VERSION=3.20 \
    # uwot/data.table are OpenMP-threaded. Left unset they grab every core on the
    # host, which on a shared machine is antisocial and on a laptop causes thermal
    # throttling that makes the run slower, not faster. Override at run time with
    # `-e OMP_NUM_THREADS=8` or the pipeline's own --threads.
    OMP_NUM_THREADS=4 \
    # Fail on the missing-locale warnings that otherwise appear mid-run and look
    # like data problems.
    LANG=C.UTF-8

# System libraries. Each one is here because a specific R package needs it to
# compile; this is not a generic "build-essential everything" block.
#   libxml2      -> XML, pulled in by the Bioconductor stack
#   libcurl/ssl  -> curl, httr: how BiocManager fetches packages
#   libpng/jpeg/tiff, libfreetype/harfbuzz/fribidi -> ragg/systemfonts: PNG output
#   libhdf5      -> rhdf5-family dependencies of the flow stack
#   zlib/bz2/lzma-> compressed FCS payloads
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff5-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libhdf5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libglpk-dev \
        libfontconfig1-dev \
        fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cyraven

# Copy ONLY the dependency declaration before installing anything. The dependency
# layer is the expensive one (~15-25 min of compilation); keeping the analysis
# code out of it means editing R/ does not rebuild every dependency.
COPY DESCRIPTION inst/scripts/install_deps.R /opt/cyraven/

# install_deps.R reads Imports/Suggests straight out of DESCRIPTION rather than
# repeating the list, so the image can never install a different set of packages
# than the package declares.
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
                BiocManager::install(version=Sys.getenv('BIOC_VERSION'), ask=FALSE, update=FALSE)" \
 && Rscript /opt/cyraven/install_deps.R

# Now the package itself, installed from source so library(cyRAVEN) works and the
# CLI resolves through system.file().
COPY . /opt/cyraven/src
RUN R CMD INSTALL --no-multiarch --with-keep.source /opt/cyraven/src \
 && cp /opt/cyraven/src/inst/scripts/entrypoint.sh /opt/cyraven/entrypoint.sh \
 && chmod +x /opt/cyraven/entrypoint.sh

# Verify the image actually works at BUILD time, so a broken image fails here
# instead of an hour into someone's analysis. --help exercises the full option
# parser and every dependency the CLI touches.
RUN Rscript "$(Rscript -e "cat(system.file('scripts','cyraven.R',package='cyRAVEN'))")" --help > /dev/null \
 && Rscript -e "for (p in c('cyRAVEN','flowCore','uwot','data.table','ggplot2', \
                            'patchwork','viridisLite','scales','yaml','optparse','hexbin')) \
                  cat(sprintf('%-12s %s\n', p, as.character(packageVersion(p))))"

# Results are written here by default; mount over it to get them on the host.
VOLUME ["/results"]
WORKDIR /work

# ENTRYPOINT resolves which copy of the package to run at container-start time
# (see entrypoint.sh): the installed one by default, or a mounted source tree
# when CYRAVEN_SOURCE is set. Either way `docker run <image> --help` and
# `docker run <image> --dir /data ...` behave like the CLI itself.
ENTRYPOINT ["/opt/cyraven/entrypoint.sh"]
CMD ["--help"]
