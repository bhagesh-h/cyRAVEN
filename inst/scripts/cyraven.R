#!/usr/bin/env Rscript
# =============================================================================
# cyRAVEN -- command-line front end
# =============================================================================
#
# This file is deliberately thin. Every line of analysis lives in the package,
# where it is documented, namespaced and tested; this script only turns command
# line arguments into the list that run_cyraven() expects.
#
# WHY THE SPLIT MATTERS. As one 5,000-line script the analysis could only be
# invoked one way -- as a process, with flags, writing files. Nothing could call a
# single step, substitute a different threshold rule, or test a figure in
# isolation, because none of it had a name outside the file. Moving the code into
# a package and leaving this behind keeps the command line exactly as it was
# while making the same work available to any R session.
#
#   Rscript cyraven.R --dir data/ --outdir results/ ...
#
# is equivalent to
#
#   library(cyRAVEN)
#   run_cyraven(list(dir = "data/", outdir = "results/", ...))
#
# Locate it after installing with:
#   system.file("scripts", "cyraven.R", package = "cyRAVEN")

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE))
    stop("the command-line front end needs the 'optparse' package:\n",
         "  install.packages(\"optparse\")", call. = FALSE)
  library(cyRAVEN)
  library(optparse)
})

opt <- parse_args(
  OptionParser(option_list = cyRAVEN:::build_option_list(),
               description = paste(
                 "Gated multi-sample flow cytometry UMAP.",
                 "Reads FCS files, derives gates from each sample's own marker",
                 "densities, scores populations from a declarative spec, embeds,",
                 "and tests differences between cohorts on per-sample values.")),
  args = commandArgs(trailingOnly = TRUE))

# A non-zero exit status matters here: this is run from shell scripts, Makefiles,
# CI and container entrypoints, all of which decide what to do next from it. A
# failure that exits 0 is worse than a failure.
ok <- withCallingHandlers(
  tryCatch({ run_cyraven(opt); TRUE },
           error = function(e) {
             message("\nERROR: ", conditionMessage(e))
             FALSE
           }),
  warning = function(w) {
    message("WARNING: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })

quit(status = if (isTRUE(ok)) 0L else 1L)
