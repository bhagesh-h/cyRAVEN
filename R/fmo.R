# SECTION 2c -- FLUORESCENCE-MINUS-ONE REFERENCE CONTROLS
# =============================================================================
#
# WHY THIS FILE EXISTS. Control handling before this was a single unstained tube
# per panel: one distribution, used as the reference for every marker that failed
# to resolve a density minimum. That is the right control for asking "where does
# autofluorescence end", and the wrong one for asking "where does this marker's
# background end in this panel".
#
# The difference is spillover. In a stained sample, a channel's negative
# population sits higher and wider than in an unstained tube, because every other
# fluorochrome in the panel contributes to it. An unstained control cannot show
# that, so a cut anchored to it sits too low and calls spillover positive. A
# fluorescence-minus-one control is the same full panel with one reagent left
# out: its distribution in that channel IS the negative population under the
# spreading the real samples experience.
#
# WHAT THIS PACKAGE ADDS TO THE PRACTICE. Placing a gate on an FMO is routine and
# is also where practitioners disagree most, because the FMO shows a continuum
# and the analyst still picks a point on it. cyRAVEN picks by the same rule it
# uses everywhere (a fixed high quantile of the control), and then reports the
# quantity nobody usually has: how far the FMO-anchored cut sits from the cut the
# sample's own density implies, in units of that threshold's own uncertainty.
#
# That number is the useful one. A valley within its own uncertainty of the FMO
# is corroborated by an independent experiment. A valley three uncertainties
# above the FMO is calling real signal background. Neither fact is obtainable
# from either control alone, and the disagreement is what a reviewer should see.
#
# HOW A CONTROL IS DECLARED. Two new sample-map columns, both optional:
#
#   fmo_for        comma-separated markers this file is the FMO for. A file may
#                  serve several channels.
#   control_group  which samples this control applies to. Absent, a control
#                  applies to every sample in its panel; present, only to samples
#                  carrying the same value. Reagent lots and instrument settings
#                  change between batches, and an FMO acquired in one batch is
#                  not the negative population of another.

#' Parse the `fmo_for` column into a marker-to-file map
#'
#' @param smap the sample map, or NULL
#' @return a data.frame with one row per (file, marker) the file controls for,
#'   carrying `sample_id`, `marker` and `control_group`; NULL when no FMO is
#'   declared
#' @export
parse_fmo_map <- function(smap) {
  if (is.null(smap) || !"fmo_for" %in% names(smap)) return(NULL)
  v <- trimws(as.character(smap$fmo_for))
  keep <- which(!is.na(v) & nzchar(v))
  if (!length(keep)) return(NULL)
  sid <- smap$sample_id %||% smap$well %||% smap$file
  grp <- if ("control_group" %in% names(smap))
    trimws(as.character(smap$control_group)) else rep(NA_character_, nrow(smap))
  rows <- lapply(keep, function(i) {
    mk <- trimws(strsplit(v[i], ",")[[1]])
    mk <- mk[nzchar(mk)]
    if (!length(mk)) return(NULL)
    data.frame(sample_id = as.character(sid[i]), marker = mk,
               control_group = grp[i], stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Which FMO file, if any, controls a given marker for a given sample
#'
#' A control with no `control_group` applies everywhere. A control with one
#' applies only to samples sharing it, and a sample whose group has no control
#' for that marker falls back to a group-less control if one exists.
#'
#' @param fmo_map from [parse_fmo_map()]
#' @param sample_id the sample being gated
#' @param marker the marker being thresholded
#' @param group_of named character vector mapping sample_id to control group
#' @return the controlling sample_id, or NA
#' @export
fmo_for_sample <- function(fmo_map, sample_id, marker, group_of = NULL) {
  if (is.null(fmo_map) || !nrow(fmo_map)) return(NA_character_)
  cand <- fmo_map[fmo_map$marker == marker, , drop = FALSE]
  if (!nrow(cand)) return(NA_character_)
  # A control never controls itself: its own channel is the one left out.
  cand <- cand[cand$sample_id != sample_id, , drop = FALSE]
  if (!nrow(cand)) return(NA_character_)
  g <- if (!is.null(group_of)) unname(group_of[as.character(sample_id)]) else NA_character_
  if (!is.na(g)) {
    hit <- cand[!is.na(cand$control_group) & cand$control_group == g, , drop = FALSE]
    if (nrow(hit)) return(hit$sample_id[1])
  }
  hit <- cand[is.na(cand$control_group) | !nzchar(cand$control_group), , drop = FALSE]
  if (nrow(hit)) return(hit$sample_id[1])
  NA_character_
}

#' Distance between a derived cut and its FMO-anchored equivalent
#'
#' The diagnostic the feature exists for. Positive means the sample's own density
#' put the cut above the control; negative means below.
#'
#' HOW TO READ `distance_in_u`. It is the gap expressed in units of that
#' threshold's own standard uncertainty, from `threshold_uncertainty.csv`. Within
#' about one, the two methods agree to the precision either can claim, and the
#' derived cut is corroborated by an independent experiment. Beyond about three,
#' they disagree by more than either can explain and one of them is wrong: a
#' derived cut far above the FMO is discarding real signal, and one far below is
#' calling spillover positive.
#'
#' @param thr_all the thresholds table, carrying `sample_id`, `marker`,
#'   `threshold` and `source`
#' @param fmo_thresholds data.frame of `sample_id`, `marker`, `fmo_threshold`,
#'   `fmo_sample`
#' @param unc optional `thresholds` element of [run_gate_uncertainty()], used to
#'   scale the distance
#' @param agree_at distance in uncertainties within which the two agree
#' @param disagree_at distance beyond which they are reported as disagreeing
#' @return a data.frame, or NULL
#' @export
fmo_agreement <- function(thr_all, fmo_thresholds, unc = NULL,
                          agree_at = 1, disagree_at = 3) {
  if (is.null(thr_all) || !nrow(thr_all)) return(NULL)
  if (is.null(fmo_thresholds) || !nrow(fmo_thresholds)) return(NULL)
  m <- merge(thr_all[, intersect(c("sample_id", "panel", "marker", "threshold",
                                   "source"), names(thr_all))],
             fmo_thresholds, by = c("sample_id", "marker"))
  if (!nrow(m)) return(NULL)

  u <- rep(NA_real_, nrow(m))
  if (!is.null(unc) && nrow(unc) && "u_combined" %in% names(unc)) {
    k <- match(paste(m$sample_id, m$marker, sep = "\r"),
               paste(unc$sample_id, unc$marker, sep = "\r"))
    u <- unc$u_combined[k]
  }
  d <- m$threshold - m$fmo_threshold
  din <- ifelse(is.finite(u) & u > 0, d / u, NA_real_)

  out <- data.frame(
    sample_id = m$sample_id,
    panel = m$panel %||% NA_character_,
    marker = m$marker,
    derived_threshold = round(m$threshold, 4),
    derived_source = m$source,
    fmo_threshold = round(m$fmo_threshold, 4),
    fmo_sample = m$fmo_sample,
    distance = round(d, 4),
    threshold_u = round(u, 4),
    distance_in_u = round(din, 2),
    stringsAsFactors = FALSE)
  out$verdict <- ifelse(
    !is.finite(out$distance_in_u), "no uncertainty available",
    ifelse(abs(out$distance_in_u) <= agree_at, "corroborated by the FMO",
    ifelse(abs(out$distance_in_u) >= disagree_at,
           ifelse(out$distance_in_u > 0,
                  "derived cut well above the FMO: signal may be discarded",
                  "derived cut well below the FMO: spillover may be called positive"),
           "differs from the FMO, within explanation")))
  out[order(-abs(replace(out$distance_in_u, !is.finite(out$distance_in_u), -1))), ,
      drop = FALSE]
}
