#!/usr/bin/env Rscript
# =============================================================================
# Install flowCore with the fix for a $SPILLOVER keyword held as a matrix
# =============================================================================
#
# THE BUG. flowCore's getSpilloverFromFlowFrame() tests
#   if (class(mySpill) != "matrix")
# A matrix has had two classes since R 4.0, c("matrix", "array"), so the test
# has length two, and since R 4.2 that is an error: "the condition has length >
# 1". read.FCS() stores $SPILLOVER as a matrix, so the function fails on every
# compensated file it is given. It is reached when flowCore evaluates a
# Gating-ML gate declared compensation-ref="FCS", which is what applying
# cyRAVEN's --export-gates file with flowUtils does. cyRAVEN itself never calls
# it: it reads the keyword and compensates in maybe_compensate().
#
# THE FIX is upstream's own, one line: `if (!is(mySpill, "matrix"))` (RGLab
# flowCore commit 903b5b3, issue #272, 2024-08-22). It reached Bioconductor
# devel (flowCore 2.25) and no release: 2.18.0 (Bioconductor 3.20, the one this
# image pins), 2.20.0 and 2.22.1 still carry the bug.
#
# WHAT THIS DOES. Downloads the flowCore source of the pinned Bioconductor
# release, replaces that line, records the change in DESCRIPTION
# (`cyRAVENPatch`) and installs it over the stock build. The version number is
# left alone, so the patch shows in packageDescription("flowCore")$cyRAVENPatch
# rather than as a version no Bioconductor release has. A source that already
# carries the fix is left as it is; a source with neither line stops the build,
# because then this is no longer the code it was written for.
#
# Usage:
#   Rscript patch_flowcore.R            # into the default library
#   Rscript patch_flowcore.R path/lib   # into another library, for testing
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[[1]] else .libPaths()[1]
bad  <- 'if (class(mySpill) != "matrix")'
good <- 'if (!is(mySpill, "matrix"))'

tmp <- tempfile("flowCore_"); dir.create(tmp)
dl <- utils::download.packages("flowCore", destdir = tmp, type = "source",
                               repos = BiocManager::repositories())
if (!nrow(dl)) stop("could not download the flowCore source", call. = FALSE)
utils::untar(dl[1, 2], exdir = tmp)
src <- file.path(tmp, "flowCore")
version <- read.dcf(file.path(src, "DESCRIPTION"), fields = "Version")[1, 1]

f <- file.path(src, "R", "eval-methods.R")
x <- readLines(f)
i <- grep(bad, x, fixed = TRUE)
if (!length(i)) {
  if (any(grepl(good, x, fixed = TRUE))) {
    message("flowCore ", version, " already carries the fix; left as installed")
    quit(save = "no")
  }
  stop("flowCore ", version, ": neither the faulty nor the fixed line found in ",
       "R/eval-methods.R; review this patch against the new source", call. = FALSE)
}
x[i] <- sub(bad, good, x[i], fixed = TRUE)
writeLines(x, f)

d <- read.dcf(file.path(src, "DESCRIPTION"))
d <- cbind(d, cyRAVENPatch = paste(
  "getSpilloverFromFlowFrame(): !is(mySpill, \"matrix\") for",
  "class(mySpill) != \"matrix\" (upstream 903b5b3, issue #272)"))
write.dcf(d, file.path(src, "DESCRIPTION"))

message("installing flowCore ", version, " with the $SPILLOVER fix into ", lib)
utils::install.packages(src, lib = lib, repos = NULL, type = "source")

# Verify in a fresh R, as install_deps.R does: a failed install only warns.
chk <- sprintf(paste0(
  ".libPaths(c('%s', .libPaths())); ",
  "f <- deparse(get('getSpilloverFromFlowFrame', asNamespace('flowCore'))); ",
  "if (!any(grepl('is(mySpill, \"matrix\")', f, fixed = TRUE))) quit(status = 1)"), lib)
if (system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(chk))) != 0L) {
  message("FAILED: the installed flowCore does not carry the fix")
  quit(status = 1L)
}
message("flowCore ", version, " installed with the $SPILLOVER fix")
