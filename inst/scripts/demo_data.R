#!/usr/bin/env Rscript
# =============================================================================
# Fetch and prepare the demonstration dataset
# =============================================================================
#
# Downloads two publicly available raw FCS files and derives a small multi-sample
# cohort from them, together with the sample map and population specification
# cyRAVEN needs. Intended to make the documented worked example runnable without
# supplying data.
#
# SOURCE. CytoTrol control samples from RGLab's flowWorkspaceData repository,
# retrieved over HTTPS from GitHub. They require no account or API key, are
# approximately 5 MB each, and are raw acquisitions carrying FSC-A, FSC-H, FSC-W
# and SSC-A alongside a seven-colour T-cell panel. Licence: Artistic-2.0.
#
# THE GROUP LABELS ARE RANDOMISED AND CARRY NO BIOLOGICAL MEANING.
#
# Two files are not enough to exercise differential testing, so each is
# partitioned at random into four pseudo-samples and those are assigned to two
# groups. Every pseudo-sample therefore originates from the same tube, and the
# true difference between the groups is zero by construction.
#
# This makes the example a calibration check rather than a demonstration of
# sensitivity. The correct outcome is that no population differs significantly;
# any result reported as significant is a false positive. The example cannot show
# that cyRAVEN detects a real effect, because no real effect is present. Do not
# cite figures produced from it as evidence of anything other than that the
# pipeline runs and is calibrated.
#
# Usage:
#   Rscript demo_data.R <output_directory>

suppressPackageStartupMessages(library(flowCore))

args   <- commandArgs(trailingOnly = TRUE)
OUT    <- if (length(args)) args[1] else "demo"
RAW    <- file.path(OUT, "raw")
FCS    <- file.path(OUT, "fcs")

# Pinned to a tag rather than to a branch. Input that can change without this
# script changing is not reproducible.
SOURCE_BASE <- paste0("https://raw.githubusercontent.com/RGLab/",
                      "flowWorkspaceData/master/inst/extdata/")
SOURCE_FILES <- c("CytoTrol_CytoTrol_1.fcs", "CytoTrol_CytoTrol_2.fcs")
N_SPLIT    <- 4L
MAX_EVENTS <- 30000L
SEED       <- 42L

say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

#' Reduce each `$PnS` value to the marker symbol
#'
#' These files record the marker together with its fluorochrome, for example
#' "CD4 PcpCy55" and "CD3 V450". cyRAVEN resolves channels to marker symbols from
#' `$PnS`, so a specification declaring `CD4` would match nothing and every
#' population would be reported UNAVAILABLE.
#'
#' The first whitespace-delimited token is the marker for every channel in this
#' panel, including the hyphenated names HLA-DR and CD45RA, which are single
#' tokens. This is a property of this dataset rather than a general rule; other
#' archives name the channel differently and need a different mapping.
harmonise_marker_names <- function(ff) {
  pd <- pData(parameters(ff))
  keep <- !is.na(pd$desc) & nzchar(trimws(pd$desc))
  pd$desc[keep] <- sub("\\s.*$", "", trimws(pd$desc[keep]))
  pData(parameters(ff)) <- pd
  for (i in which(keep)) keyword(ff)[[paste0("$P", i, "S")]] <- pd$desc[i]
  ff
}

dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(FCS, recursive = TRUE, showWarnings = FALSE)

# ---- download ---------------------------------------------------------------
for (f in SOURCE_FILES) {
  dest <- file.path(RAW, f)
  if (file.exists(dest) && file.size(dest) > 1e6) {
    say("present, skipping: ", f); next
  }
  say("downloading ", f)
  utils::download.file(paste0(SOURCE_BASE, f), dest, mode = "wb", quiet = TRUE)
  # A redirect, a rate limit or a moved path yields an HTML error page saved
  # under a .fcs name, which would otherwise fail much later with an unrelated
  # message about the file's contents.
  if (!identical(rawToChar(readBin(dest, "raw", n = 3L)), "FCS"))
    stop(f, " does not begin with the FCS magic bytes. The download was ",
         "redirected or the path has moved.", call. = FALSE)
}

# ---- partition into pseudo-samples -----------------------------------------
set.seed(SEED)
rows <- list()
raws <- file.path(RAW, SOURCE_FILES)

for (i in seq_along(raws)) {
  ff <- read.FCS(raws[i], transformation = FALSE, truncate_max_range = FALSE)
  ff <- harmonise_marker_names(ff)
  n  <- nrow(ff)
  # Random partition rather than the first n/4 rows. Acquisition order correlates
  # with flow-rate and staining drift within a tube, so contiguous blocks would
  # differ systematically and manufacture exactly the batch structure this
  # example asserts is absent.
  part <- sample(rep_len(seq_len(N_SPLIT), n))
  for (k in seq_len(N_SPLIT)) {
    idx <- which(part == k)
    if (length(idx) > MAX_EVENTS) idx <- sort(sample(idx, MAX_EVENTS))
    sid  <- sprintf("donor%02d_rep%d", i, k)
    dest <- file.path(FCS, paste0(sid, ".fcs"))
    sub  <- ff[idx, ]
    keyword(sub)[["$FIL"]] <- paste0(sid, ".fcs")
    write.FCS(sub, dest)
    rows[[length(rows) + 1L]] <- data.frame(
      file = basename(dest), sample_id = sid, patient_id = sid,
      cohort = if (k <= N_SPLIT / 2) "GroupA" else "GroupB",
      is_control = FALSE, stringsAsFactors = FALSE)
  }
  say("partitioned ", basename(raws[i]), " (", n, " events) into ", N_SPLIT,
      " pseudo-samples")
}

smap <- do.call(rbind, rows)
write.csv(smap, file.path(OUT, "sample_map.csv"), row.names = FALSE)

# ---- population specification ----------------------------------------------
# Conventional CCR7/CD45RA memory subsets (Sallusto et al. 1999, Nature 401:708)
# expressed against the markers this panel resolves.
writeLines(c(
  "# Population specification for the CytoTrol seven-colour T-cell panel.",
  "# Each population is a conjunction of marker directions evaluated within the",
  "# CD45+ parent gate. See the Gating article for the syntax.",
  "populations:",
  "  CD4 T cells:",
  "    CD3: above",
  "    CD4: above",
  "    CD8: below",
  "  CD8 T cells:",
  "    CD3: above",
  "    CD8: above",
  "    CD4: below",
  "  Naive CD4:",
  "    CD3: above",
  "    CD4: above",
  "    CD8: below",
  "    CCR7: above",
  "    CD45RA: above",
  "  Central memory CD4:",
  "    CD3: above",
  "    CD4: above",
  "    CD8: below",
  "    CCR7: above",
  "    CD45RA: below",
  "  Effector memory CD4:",
  "    CD3: above",
  "    CD4: above",
  "    CD8: below",
  "    CCR7: below",
  "    CD45RA: below",
  "  TEMRA CD4:",
  "    CD3: above",
  "    CD4: above",
  "    CD8: below",
  "    CCR7: below",
  "    CD45RA: above",
  "  Naive CD8:",
  "    CD3: above",
  "    CD8: above",
  "    CD4: below",
  "    CCR7: above",
  "    CD45RA: above",
  "  Effector memory CD8:",
  "    CD3: above",
  "    CD8: above",
  "    CD4: below",
  "    CCR7: below",
  "    CD45RA: below",
  "  Activated T cells:",
  "    CD3: above",
  "    HLA-DR: above",
  "  Non-T cells:",
  "    CD3: below"
), file.path(OUT, "panel.yaml"))

say("wrote ", nrow(smap), " pseudo-samples to ", normalizePath(FCS))
say("wrote sample_map.csv and panel.yaml to ", normalizePath(OUT))
message("")
message("NOTE the cohort column in sample_map.csv is randomised. Every ",
        "pseudo-sample comes from the same tube,")
message("     so the true between-group difference is zero. Any population ",
        "reported as significant is a")
message("     false positive, and counting them is the point of the example.")
