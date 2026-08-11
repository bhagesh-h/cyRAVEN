#!/usr/bin/env Rscript
# =============================================================================
# Prepare the demonstration cohort
# =============================================================================
#
# Writes a runnable multi-sample cohort, together with the sample map and
# population specification cyRAVEN needs, so the worked example in the README can
# be reproduced without supplying data.
#
# SOURCE. The graft-versus-host disease dataset distributed with flowCore
# (Brinkman et al. 2007, Biol Blood Marrow Transplant 13:691). Peripheral blood
# was drawn from allogeneic blood and marrow transplant recipients at successive
# visits and stained with a four-colour myeloid panel. No download is required:
# the data ship inside a package cyRAVEN already depends on, which makes the
# example reproducible offline and immune to a repository moving or expiring.
# Licence: Artistic-2.0, as flowCore.
#
# WHY THIS COHORT. It exercises the parts of the pipeline that a single-tube
# example cannot:
#
#   five patients          donor-level variation is real, so the sample-level
#                          aggregation has something to aggregate over
#   seven visits           a genuine acquisition batch variable, not a
#                          constructed one
#   two GvHD grades        a real clinical contrast to test between groups
#   CD45 in the panel      the parent gate is live leukocytes as designed,
#                          rather than falling back to all scatter-gated events
#   a Time channel         acquisition stability can be assessed
#   2,205 to 66,105 events per file, so thresholds resolve from density rather
#                          than falling back to a quantile
#
# WHAT IT IS NOT. There is no healthy control arm. Both groups are transplant
# recipients and the contrast is GvHD grade 1 against grade 3, so a difference
# between them is a difference in disease severity and not between disease and
# health. The grades are also unbalanced across patients, which section 4 of the
# worked example takes up.
#
# TWO PREPARATION STEPS, AND WHY EACH IS NEEDED.
#
# $PnS carries the marker together with its fluorochrome, "CD15 FITC" and "CD45
# PE". cyRAVEN resolves channels to marker symbols from $PnS, so a specification
# declaring CD15 would match nothing and every population would be reported
# UNAVAILABLE. The first whitespace-delimited token is the marker for every
# channel in this panel. This is a property of this acquisition, not a general
# rule; other archives name channels differently and need a different mapping.
#
# The file also carries FL2-A, the pulse area of the same detector that FL2-H
# records the height of. It is the only area channel present, and cyRAVEN keeps
# area channels in preference to height precisely so that a marker appearing as
# both is not counted twice. Left in place, that rule would select FL2-A alone
# and discard the other three markers. It is dropped here, which leaves a
# consistently height-only file: read_fcs_resolved() then falls back to height
# for all four markers and says so, because height and area are different
# measurements and thresholds from one are not interchangeable with the other.
#
# Nothing else is altered. Event values are written as flowCore ships them.
#
# Usage:
#   Rscript demo_data.R <output_directory>

suppressPackageStartupMessages(library(flowCore))

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args)) args[1] else "demo"
FCS  <- file.path(OUT, "fcs")

say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

dir.create(FCS, recursive = TRUE, showWarnings = FALSE)

utils::data("GvHD", package = "flowCore", envir = environment())
gv <- get("GvHD", envir = environment())
pd <- pData(gv)

#' Reduce each `$PnS` to the marker symbol and drop the redundant area channel.
#' See the header for why each step is required.
prepare_frame <- function(ff, sid) {
  # Subset by name rather than by a logical vector: flowFrame indexing is
  # documented for character and numeric columns, and setdiff() is also a no-op
  # on a frame that does not carry the channel.
  ff <- ff[, setdiff(colnames(ff), "FL2-A")]
  pdat <- pData(parameters(ff))
  d <- as.character(pdat$desc)
  has <- !is.na(d) & nzchar(trimws(d))
  d[has] <- sub("\\s.*$", "", trimws(d[has]))
  pdat$desc <- d
  pData(parameters(ff)) <- pdat
  # Indices are those of the SUBSET frame, so the keywords stay consistent with
  # the parameter block after FL2-A is removed.
  for (i in which(has)) keyword(ff)[[paste0("$P", i, "S")]] <- d[i]
  keyword(ff)[["$FIL"]] <- paste0(sid, ".fcs")
  ff
}

# ---- write one FCS per sample ----------------------------------------------
rows <- list()
for (i in seq_along(gv)) {
  sid <- rownames(pd)[i]
  ff  <- prepare_frame(gv[[i]], sid)
  suppressWarnings(write.FCS(ff, file.path(FCS, paste0(sid, ".fcs"))))
  rows[[length(rows) + 1L]] <- data.frame(
    file       = paste0(sid, ".fcs"),
    sample_id  = sid,
    patient_id = paste0("patient", pd$Patient[i]),
    # The clinical contrast. Grades are as recorded in the source data.
    cohort     = paste0("GvHD grade ", as.character(pd$Grade[i])),
    # The acquisition batch. Visits are successive draws from the same patient,
    # so a visit is a distinct acquisition occasion rather than a label invented
    # for this example.
    visit      = paste0("visit", pd$Visit[i]),
    days       = pd$Days[i],
    is_control = FALSE,
    stringsAsFactors = FALSE)
}
smap <- do.call(rbind, rows)
write.csv(smap, file.path(OUT, "sample_map.csv"), row.names = FALSE)

# ---- population specification ----------------------------------------------
# A four-colour myeloid panel supports lineage-level populations and no more.
# Each definition below is a conventional one for these markers; none is tuned
# to this cohort, which is what lets the diagnostics contradict it.
writeLines(c(
  "# Population specification for the GvHD myeloid panel (CD15, CD45, CD14,",
  "# CD33). Every population is a conjunction of marker directions evaluated",
  "# within the CD45+ parent gate. See the Gating article for the syntax.",
  "#",
  "# These are conventional lineage definitions for a four-colour myeloid tube.",
  "# They are declared before the data are examined so that the phenotype,",
  "# threshold-drift and cluster-concordance outputs can contradict them.",
  "populations:",
  "  Granulocytes:",
  "    CD45: above",
  "    SSC-A: above",
  "    CD15: above",
  "  Monocytes:",
  "    CD45: above",
  "    CD14: above",
  "    CD15: below",
  "  CD33 positive myeloid:",
  "    CD45: above",
  "    CD33: above",
  "  Lymphocytes:",
  "    CD45: above",
  "    SSC-A: below",
  "    CD33: below",
  "    CD14: below",
  "  Myeloid marker negative:",
  "    CD45: above",
  "    CD14: below",
  "    CD15: below",
  "    CD33: below",
  "",
  "# Marker intensity read WITHIN an already-defined population, rather than",
  "# used to define one. Scope is declared per block, and the blocks are split so",
  "# that no marker is reported inside a gate its own threshold helped draw:",
  "# testing CD15 within a CD15-positive population returns 100 percent in every",
  "# sample, zero variance and an undefined p-value, which measures the",
  "# definition rather than the biology.",
  "functional_blocks:",
  "  CD33 on gated myeloid subsets:",
  "    markers: [CD33]",
  "    populations: [Granulocytes, Monocytes]",
  "  maturation markers on CD33 positive cells:",
  "    markers: [CD14, CD15]",
  "    populations: [CD33 positive myeloid]",
  "",
  "# A derived quantity: the ratio of two declared populations, tested with the",
  "# same statistics as any abundance. The granulocyte-to-lymphocyte ratio is a",
  "# conventional inflammatory index and is not derivable from either frequency",
  "# alone, because both are percentages of the same parent.",
  "ratios:",
  "  gran_lymph:",
  "    label: \"Granulocyte:lymphocyte ratio\"",
  "    numerator: Granulocytes",
  "    denominator: Lymphocytes"
), file.path(OUT, "panel.yaml"))

say("wrote ", nrow(smap), " sample(s) to ", normalizePath(FCS))
say("wrote sample_map.csv and panel.yaml to ", normalizePath(OUT))
message("")
message("cohort: ", paste(sprintf("%s n=%d", names(table(smap$cohort)),
                                  table(smap$cohort)), collapse = "; "))
message("patients: ", length(unique(smap$patient_id)),
        "   acquisition batches (visits): ", length(unique(smap$visit)))
message("")
message("NOTE both groups are transplant recipients. The contrast is GvHD grade,")
message("     not disease against health, and the grades are unbalanced across")
message("     patients. Read confounding_diagnostics.csv before interpreting any")
message("     between-group difference.")
