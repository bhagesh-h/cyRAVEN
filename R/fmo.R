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
#'   carrying `sample_id`, `marker`, `control_group`, and `control_kind`
#'   (`minus_one` where the file leaves out a single reagent, `minus_multiple`
#'   where it leaves out several) with `n_markers_in_file`; NULL when no FMO
#'   is declared
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

  # FLUORESCENCE MINUS ONE, OR MINUS MANY? The distinction is not bookkeeping.
  # An FMO is the full panel with ONE reagent left out, so its distribution in
  # that channel is the negative population under exactly the spreading the real
  # samples experience (Roederer 2001, Cytometry 45:194). A tube with twelve
  # reagents left out is a different control: the eleven others are not spilling
  # into the channel either, so its negative is NARROWER than the true FMO, and a
  # cut anchored to it sits TOO LOW and calls spillover positive.
  #
  # The direction of that bias is knowable even when its size is not, so the kind
  # is recorded per file and travels into the threshold `source` string. It is
  # still a far better reference than an unstained tube or a quantile of the
  # parent; it is just not the control the FMO literature describes, and a reader
  # comparing two runs should be able to see which they have.
  per_file <- table(out$sample_id)
  out$control_kind <- ifelse(per_file[out$sample_id] > 1L,
                             "minus_multiple", "minus_one")
  out$n_markers_in_file <- as.integer(per_file[out$sample_id])
  out
}

#' Place a cut on a control channel, and refuse a control that cannot be one
#'
#' Two things go wrong with a control-anchored cut, and both were observed on
#' real minus-multiple data before this existed.
#'
#' A FAR QUANTILE RIDES AN ARTEFACT TAIL. The negative population's upper edge
#' is what the control is for, and `q = 0.995` finds it only when the control
#' contains nothing but that population. A control tube carrying about one
#' percent of bright events -- aggregates, dead cells, a residual stained
#' population -- puts the 99.5th percentile inside that tail rather than at the
#' edge of the negative. On the cohort this was written for, one channel's
#' control sat at 1.04 at the 95th percentile and 5.60 at the 99.5th: the cut
#' moved by four and a half units of a transformed scale, and landed above the
#' full stain's own 99th percentile, so nothing in any sample would have been
#' called positive. [tail_threshold()] is immune to this, because it locates the
#' negative mode and measures the spread of the half that cannot contain
#' positives. Taking whichever of the two is lower keeps the quantile where the
#' control is clean -- for a Gaussian negative the quantile is the lower of the
#' pair, so nothing changes -- and falls back to the robust estimate only where
#' the quantile has left the negative population altogether.
#'
#' THE CONTROL IS NOT ACTUALLY MISSING THE REAGENT. A minus control must be
#' dimmer in the channel it leaves out. When it is brighter, something is wrong
#' with the tube, the panel sheet, or the unmixing, and anchoring to it raises
#' the cut and deletes a real population. This is not a hypothetical: of twelve
#' channels declared on one minus-multiple tube, five were no dimmer than the
#' full stain and one was brighter at every quantile. Comparing the two negative
#' modes, in units of the sample's own spread, catches that before it reaches a
#' frequency table. A refused control is reported, not silently dropped, and the
#' threshold falls back to the sample's own data.
#'
#' @param control_x transformed control values in this channel.
#' @param sample_x transformed sample values in the same channel, used only to
#'   check that the control is the dimmer of the two.
#' @param q the quantile taken where the control is clean. Default `0.995`.
#' @param k passed to [tail_threshold()].
#' @param max_mode_shift how far above the sample's negative mode the control's
#'   may sit, in units of the sample's own robust spread, before the control is
#'   refused. Default `1`.
#' @return a list carrying `threshold` (`NA_real_` where the control is
#'   refused), `verdict`, and the two candidate cuts.
#' @export
control_cut <- function(control_x, sample_x = NULL, q = 0.995, k = 3.5,
                        max_mode_shift = 1) {
  out <- function(th, v, qc = NA_real_, rc = NA_real_)
    list(threshold = th, verdict = v, quantile_cut = qc, robust_cut = rc)
  cx <- control_x[is.finite(control_x)]
  if (length(cx) < 100L) return(out(NA_real_, "too few control events"))

  qc <- as.numeric(stats::quantile(cx, q, na.rm = TRUE))
  rc <- tail_threshold(cx, k = k)

  # IS THE CONTROL REALLY THE DIMMER TUBE? Both modes come from the same kernel,
  # and the margin is the sample's own spread, so the test does not care what
  # scale or transform is in force.
  if (!is.null(sample_x)) {
    sx <- sample_x[is.finite(sample_x)]
    if (length(sx) >= 200L) {
      mode_of <- function(z) {
        d <- stats::density(z, n = 1024, adjust = 2)
        d$x[which.max(d$y)]
      }
      mc <- mode_of(cx); ms <- mode_of(sx)
      lower <- sx[sx <= ms]
      ss <- stats::mad(c(lower, 2 * ms - lower), constant = 1.4826)
      if (!is.finite(ss) || ss <= 0) ss <- stats::sd(sx)
      if (is.finite(ss) && ss > 0 && mc > ms + max_mode_shift * ss)
        return(out(NA_real_, "control brighter than the sample", qc, rc))
    }
  }

  th <- if (is.finite(rc)) min(qc, rc) else qc
  if (!is.finite(th)) return(out(NA_real_, "no cut could be placed", qc, rc))
  out(th, if (is.finite(rc) && rc < qc) "tail trimmed" else "ok", qc, rc)
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
#'   `fmo_sample`, and optionally `fmo_verdict` from [control_cut()]. A row whose
#'   `fmo_threshold` is `NA` is a control that was refused, and is reported as
#'   such rather than compared.
#' @param unc optional `thresholds` element of [run_gate_uncertainty()], used to
#'   scale the distance. A cut the control supplied has no uncertainty of its
#'   own, because the events it came from are in a separate tube, so those rows
#'   carry `distance` but not `distance_in_u`. The scaled verdict is available
#'   for the rows where the sample's own cut was the one applied.
#' @param agree_at distance in uncertainties within which the two agree
#' @param disagree_at distance beyond which they are reported as disagreeing
#' @return a data.frame, or NULL
#' @export
fmo_agreement <- function(thr_all, fmo_thresholds, unc = NULL,
                          agree_at = 1, disagree_at = 3) {
  if (is.null(thr_all) || !nrow(thr_all)) return(NULL)
  if (is.null(fmo_thresholds) || !nrow(fmo_thresholds)) return(NULL)
  m <- merge(thr_all[, intersect(c("sample_id", "panel", "marker", "threshold",
                                   "source", "derived_threshold",
                                   "derived_source"), names(thr_all))],
             fmo_thresholds, by = c("sample_id", "marker"))
  if (!nrow(m)) return(NULL)

  # COMPARE THE RIGHT PAIR. Where the control replaced the data-derived cut,
  # `threshold` IS the control's cut, and subtracting one from the other gives
  # zero for every such row -- the control measured against itself. The cut the
  # sample's own density implied is carried alongside for exactly this reason;
  # it is used where present, and `threshold` only where the control was not the
  # one that won.
  d_thr <- if ("derived_threshold" %in% names(m))
    ifelse(is.finite(m$derived_threshold), m$derived_threshold, m$threshold)
    else m$threshold
  d_src <- if ("derived_source" %in% names(m))
    ifelse(!is.na(m$derived_source), m$derived_source, m$source) else m$source

  u <- rep(NA_real_, nrow(m))
  if (!is.null(unc) && nrow(unc) && "u_combined" %in% names(unc)) {
    k <- match(paste(m$sample_id, m$marker, sep = "\r"),
               paste(unc$sample_id, unc$marker, sep = "\r"))
    u <- unc$u_combined[k]
  }
  d <- d_thr - m$fmo_threshold
  din <- ifelse(is.finite(u) & u > 0, d / u, NA_real_)

  out <- data.frame(
    sample_id = m$sample_id,
    panel = m$panel %||% NA_character_,
    marker = m$marker,
    derived_threshold = round(d_thr, 4),
    derived_source = d_src,
    applied_threshold = round(m$threshold, 4),
    applied_source = m$source,
    fmo_threshold = round(m$fmo_threshold, 4),
    fmo_sample = m$fmo_sample,
    fmo_verdict = m$fmo_verdict %||% NA_character_,
    distance = round(d, 4),
    threshold_u = round(u, 4),
    distance_in_u = round(din, 2),
    stringsAsFactors = FALSE)
  out$verdict <- ifelse(
    !is.finite(out$fmo_threshold), "control refused, not used",
    ifelse(
    !is.finite(out$distance_in_u), "no uncertainty available",
    ifelse(abs(out$distance_in_u) <= agree_at, "corroborated by the FMO",
    ifelse(abs(out$distance_in_u) >= disagree_at,
           ifelse(out$distance_in_u > 0,
                  "derived cut well above the FMO: signal may be discarded",
                  "derived cut well below the FMO: spillover may be called positive"),
           "differs from the FMO, within explanation"))))
  out[order(-abs(replace(out$distance_in_u, !is.finite(out$distance_in_u), -1))), ,
      drop = FALSE]
}
