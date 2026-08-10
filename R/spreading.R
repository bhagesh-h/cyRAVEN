# SECTION 2d -- SPILLOVER SPREADING
# =============================================================================
#
# WHY THIS FILE EXISTS. cyRAVEN places every cut at a density minimum. Spreading
# is the one optical effect that removes density minima without moving anything
# else, so it is the most common reason a marker that should resolve does not.
#
# The mechanism. Compensation and unmixing subtract a fluorochrome's expected
# contribution to another detector, and the subtraction is correct on average.
# What it cannot undo is the photon-counting variance that came with the
# contribution: a channel receiving spillover from a bright neighbour has a
# negative population that is correctly CENTRED and abnormally WIDE. Two modes
# that were separable become one, the valley fills in, and `resolve_threshold()`
# reports `quantile_fallback` with no indication of why.
#
# Without this table the run says "no minimum was found" and stops there. With
# it, a `needs_review` marker that also sits at the top of the received-spreading
# list has an explanation, and the fix is a panel design decision rather than a
# gating one. On spectral platforms this is the effect that limits panel size in
# practice rather than detector count (Mage et al. 2026, Cytometry A 108:e70044).
#
# HOW IT IS MEASURED HERE, AND HOW THAT DIFFERS FROM THE SPILLOVER SPREADING
# MATRIX. The published SSM is computed from single-stain controls, where the
# contribution of one fluorochrome can be isolated. cyRAVEN does not have those:
# it is given stained samples. What it can measure is the observable consequence
# in the data it does have. For each ordered pair of markers, the spread of the
# receiver's NEGATIVE population is compared between cells that are negative for
# the source and cells that are positive for it. Restricting to the receiver's
# own negatives is what makes it a measure of spreading rather than of biology:
# genuine co-expression moves the receiver's positive cells, not the width of its
# negatives.
#
# This is a ranking, not a calibration. It says which channel pair is costing the
# most resolution in this panel on these samples. It does not give an SSM value
# comparable to a published one, and it is not a substitute for acquiring proper
# single-stain controls.

#' Spreading received by each marker from each other marker
#'
#' @param tmat transformed marker matrix for one sample
#' @param thr named threshold vector for the same markers
#' @param parent logical mask of the cells to use, normally the parent gate
#' @param min_cells smallest group of cells that can support a spread estimate
#' @return a data.frame with one row per ordered (source, receiver) pair, or NULL
#' @export
spreading_pairs <- function(tmat, thr, parent = NULL, min_cells = 200L) {
  if (is.null(tmat) || !ncol(tmat)) return(NULL)
  mk <- intersect(colnames(tmat), names(thr))
  mk <- mk[is.finite(thr[mk])]
  if (length(mk) < 2L) return(NULL)
  if (is.null(parent)) parent <- rep(TRUE, nrow(tmat))
  if (sum(parent) < min_cells * 2L) return(NULL)

  rows <- list()
  for (rec in mk) {
    yv <- tmat[, rec]
    # The receiver's own negative population, which is the distribution whose
    # width spreading inflates. Positive cells are excluded because a real
    # signal there is biology, not optics.
    neg <- parent & is.finite(yv) & yv < thr[[rec]]
    if (sum(neg) < min_cells * 2L) next
    for (src in setdiff(mk, rec)) {
      xv <- tmat[, src]
      lo <- neg & is.finite(xv) & xv <  thr[[src]]
      hi <- neg & is.finite(xv) & xv >= thr[[src]]
      if (sum(lo) < min_cells || sum(hi) < min_cells) next
      s_lo <- stats::mad(yv[lo], na.rm = TRUE)
      s_hi <- stats::mad(yv[hi], na.rm = TRUE)
      if (!is.finite(s_lo) || !is.finite(s_hi) || s_lo <= 0) next
      rows[[length(rows) + 1L]] <- data.frame(
        source = src, receiver = rec,
        n_source_neg = sum(lo), n_source_pos = sum(hi),
        spread_source_neg = round(s_lo, 4),
        spread_source_pos = round(s_hi, 4),
        spreading = round(s_hi - s_lo, 4),
        spreading_ratio = round(s_hi / s_lo, 3),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Spreading across a cohort, and which markers it is costing
#'
#' Aggregates [spreading_pairs()] over samples and joins the result to how each
#' receiver's threshold was obtained. A marker that receives a great deal of
#' spreading AND falls back to a quantile is the actionable finding: its cut is
#' unresolved for a reason the panel can fix.
#'
#' @param pops per-sample scoring results, carrying `tmat` and `thresholds`
#' @param gates per-sample gate objects
#' @param thr_all the thresholds table, used for the fallback rate per marker
#' @param panel_of named character vector mapping sample to panel
#' @param max_samples cap on samples contributing, since the ranking stabilises
#'   long before every sample is used
#' @param flag_ratio spreading ratio at or above which a pair is called
#'   substantial
#' @return list(pairs, receivers), each a data.frame, or NULL
#' @export
run_spreading_report <- function(pops, gates, thr_all = NULL, panel_of = NULL,
                                 max_samples = 8L, flag_ratio = 1.25) {
  if (is.null(pops) || !length(pops)) return(NULL)
  sids <- names(pops)
  if (length(sids) > max_samples) sids <- sids[seq_len(max_samples)]

  rows <- list()
  for (s in sids) {
    P <- pops[[s]]; g <- gates[[s]]
    if (is.null(P$tmat) || is.null(P$thresholds) || is.null(g)) next
    r <- tryCatch(spreading_pairs(P$tmat, P$thresholds, g$masks$cd45_pos),
                  error = function(e) NULL)
    if (is.null(r)) next
    r$sample_id <- s
    r$panel <- panel_of[[s]] %||% NA_character_
    rows[[length(rows) + 1L]] <- r
  }
  if (!length(rows)) return(NULL)
  pairs <- do.call(rbind, rows)

  key <- paste(pairs$panel, pairs$source, pairs$receiver, sep = "\r")
  agg <- do.call(rbind, lapply(split(pairs, key), function(d) data.frame(
    panel = d$panel[1], source = d$source[1], receiver = d$receiver[1],
    n_samples = nrow(d),
    median_spreading = round(stats::median(d$spreading, na.rm = TRUE), 4),
    median_ratio = round(stats::median(d$spreading_ratio, na.rm = TRUE), 3),
    stringsAsFactors = FALSE)))
  rownames(agg) <- NULL
  agg$substantial <- is.finite(agg$median_ratio) & agg$median_ratio >= flag_ratio
  agg <- agg[order(-replace(agg$median_ratio, !is.finite(agg$median_ratio), -1)), ]

  # Per receiver: how much it takes in total, from whom, and whether its own
  # threshold is one this package could not resolve.
  rkey <- paste(agg$panel, agg$receiver, sep = "\r")
  rec <- do.call(rbind, lapply(split(agg, rkey), function(d) {
    w <- d[which.max(replace(d$median_ratio, !is.finite(d$median_ratio), -1)), ]
    data.frame(panel = d$panel[1], receiver = d$receiver[1],
               n_sources = nrow(d),
               worst_source = w$source,
               worst_ratio = w$median_ratio,
               total_spreading = round(sum(pmax(d$median_spreading, 0), na.rm = TRUE), 4),
               n_substantial_sources = sum(d$substantial, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(rec) <- NULL

  rec$fallback_rate <- NA_real_
  if (!is.null(thr_all) && nrow(thr_all) && "source" %in% names(thr_all)) {
    fb <- do.call(rbind, lapply(
      split(thr_all, paste(thr_all$panel, thr_all$marker, sep = "\r")),
      function(d) data.frame(panel = d$panel[1], receiver = d$marker[1],
                             fallback_rate = mean(d$source == "quantile_fallback",
                                                  na.rm = TRUE),
                             stringsAsFactors = FALSE)))
    k <- match(paste(rec$panel, rec$receiver, sep = "\r"),
               paste(fb$panel, fb$receiver, sep = "\r"))
    rec$fallback_rate <- round(fb$fallback_rate[k], 3)
  }
  # The finding worth acting on: a marker whose cut this package could not
  # resolve AND which is receiving substantial spreading. The explanation for
  # the first is the second, and neither column says so on its own.
  rec$verdict <- ifelse(
    is.finite(rec$fallback_rate) & rec$fallback_rate > 0.5 &
      rec$n_substantial_sources > 0,
    "unresolved and heavily spread: a panel design problem, not a gating one",
    ifelse(rec$n_substantial_sources > 0, "receives substantial spreading",
           "no substantial spreading detected"))
  rec <- rec[order(-rec$total_spreading), ]

  list(pairs = agg, receivers = rec)
}
