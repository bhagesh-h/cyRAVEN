#!/usr/bin/env Rscript
# =============================================================================
# Install everything cyRAVEN declares, from the declaration itself
# =============================================================================
#
# WHY IT PARSES DESCRIPTION INSTEAD OF CARRYING ITS OWN LIST: a second list is a
# second thing to forget. If this file named the packages directly, adding an
# Import would build an image that installs one set and loads another, and the
# failure would surface at run time, inside the container, in the middle of an
# analysis. Reading Imports/Suggests out of DESCRIPTION makes that impossible by
# construction: there is one declaration and this script obeys it.
#
# WHY BiocManager AND NOT install.packages(): flowCore (and FlowSOM, if the
# clustering cross-check is wanted) live on Bioconductor, which install.packages()
# cannot see. BiocManager::install() resolves CRAN and Bioconductor together and
# honours the pinned BIOC_VERSION the image sets, so a rebuild months from now
# gets the same release rather than whatever is current.
#
# Usage:
#   Rscript install_deps.R              # DESCRIPTION in the working directory
#   Rscript install_deps.R path/to/DESCRIPTION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
desc_path <- if (length(args)) args[[1]] else "DESCRIPTION"
if (!file.exists(desc_path) && file.exists(file.path(desc_path, "DESCRIPTION")))
  desc_path <- file.path(desc_path, "DESCRIPTION")
if (!file.exists(desc_path))
  stop("DESCRIPTION not found at: ", desc_path, call. = FALSE)

dcf <- read.dcf(desc_path)

# Suggests is installed too, deliberately. They are optional to the PACKAGE but
# not to this image: without optparse there is no command line, without knitr the
# vignettes do not build, and without testthat the build-time self-check cannot
# run. An image that can only exercise half the package is not reproducible.
fields <- c("Depends", "Imports", "Suggests")
raw <- unlist(lapply(intersect(fields, colnames(dcf)), function(f) dcf[1, f]))

pkgs <- unlist(strsplit(paste(raw, collapse = ","), ","))
pkgs <- trimws(sub("\\(.*", "", pkgs))          # drop version constraints
pkgs <- pkgs[nzchar(pkgs)]
pkgs <- setdiff(pkgs, c("R", rownames(installed.packages(priority = "base"))))
pkgs <- unique(pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

message("installing ", length(pkgs), " package(s) declared in ", desc_path, ":")
message("  ", paste(sort(pkgs), collapse = ", "))

BiocManager::install(pkgs, ask = FALSE, update = FALSE)

# Verify rather than trust. BiocManager::install() warns on a package it could
# not install; it does not stop. Exiting non-zero here is what turns that warning
# into a failed image build instead of a container that dies later.
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("\nFAILED to install: ", paste(missing, collapse = ", "))
  quit(status = 1L)
}
message("\nall declared dependencies present")
