# SECTION 2e -- BEAD CALIBRATION TO MESF OR ERF UNITS
# =============================================================================
#
# WHY THIS FILE EXISTS. Every intensity this package reports is in the
# instrument's own arbitrary units. Those are comparable within a run, and
# through `--baseline` within an instrument over time, and not at all between
# instruments or between studies. That is the ceiling on what conformance can
# do: it can tell you a cohort moved, and it cannot tell you whether the assay
# changed or the cytometer did, because both are expressed in the same
# uninterpretable scale.
#
# Calibration beads separate them. A bead set carries populations of known
# fluorophore equivalence, assigned in MESF (molecules of equivalent soluble
# fluorochrome) or ERF (equivalent reference fluorophores), traceable to primary
# standards. Acquiring one on the instrument that ran the samples gives a map
# from channel units to those units, and thresholds and medians expressed through
# it are portable (Wang et al., NIST; Cytometry A 91:540).
#
# WHAT IS FITTED. Detector response is linear in fluorophore count plus a
# background offset, so the map is a straight line in LINEAR units:
#
#     assigned_value = slope * channel_value + intercept
#
# fitted across the bead populations of that channel. The fit is reported with
# its R squared and its residuals, because a calibration nobody checked is worse
# than none: it converts an honest arbitrary number into a dishonest absolute
# one.
#
# WHAT THIS IS NOT. It does not correct spillover, it does not make two panels
# comparable, and it does not rescue a channel whose beads did not resolve. It
# converts units, for the channels where the fit holds, and says so per channel.
#
# WHY IT IS OPT-IN AND CANNOT BE OTHERWISE. It changes the scale every threshold
# and every median sits on. Every number in the run moves, by construction.

#' Locate bead populations in one channel
#'
#' Bead sets are discrete populations, so the peaks of the density are the
#' populations. Peaks are returned in ascending order, which is the order
#' assigned-value tables conventionally use.
#'
#' @param x channel values for the bead file, in linear instrument units
#' @param n_expected how many populations the bead set has
#' @param bins,smooth density resolution and smoothing
#' @param min_frac smallest share of events a population must hold
#' @param min_sep_log10 minimum separation between two populations, in decades.
#'   Bead populations are narrow, so histogram noise readily splits one of them
#'   into two adjacent maxima; without a separation rule the split pair occupies
#'   two of the expected slots and the brightest real population is dropped.
#' @return numeric vector of peak locations, ascending, or NULL
#' @export
bead_peaks <- function(x, n_expected, bins = 256L, smooth = 5, min_frac = 0.005,
                       min_sep_log10 = 0.15) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 500L) return(NULL)
  # Beads span decades, so peaks are found on a log scale where they are
  # comparably wide, then mapped back.
  lx <- log10(x)
  h <- graphics::hist(lx, breaks = bins, plot = FALSE)
  y <- h$counts
  if (smooth > 1) y <- as.numeric(stats::filter(y, rep(1 / smooth, smooth),
                                                sides = 2))
  y[!is.finite(y)] <- 0
  ctr <- h$mids
  tot <- sum(y)
  if (!is.finite(tot) || tot <= 0) return(NULL)

  # A local maximum that holds a meaningful share of the events.
  is_pk <- c(FALSE, y[-c(1, length(y))] > y[-c(length(y) - 1, length(y))] &
                    y[-c(1, 2)] <= y[-c(1, length(y))], FALSE)
  pk <- which(is_pk & y >= min_frac * tot)
  if (!length(pk)) return(NULL)

  # Merge maxima closer together than one population is wide, tallest first, so
  # a split population contributes one peak rather than two.
  pk <- pk[order(-y[pk])]
  kept <- integer(0)
  for (i in pk)
    if (!length(kept) || all(abs(ctr[i] - ctr[kept]) >= min_sep_log10))
      kept <- c(kept, i)

  # Then the n_expected tallest, restored to ascending order so they line up
  # with the assigned-value table.
  if (length(kept) > n_expected) kept <- kept[seq_len(n_expected)]
  10^sort(ctr[kept])
}

#' Fit a channel-units to assigned-units calibration from a bead acquisition
#'
#' @param bead_exprs the bead file's expression matrix, linear units
#' @param channel_cols named integer vector, marker or channel name to column
#' @param assigned data.frame with a `marker` column and one column per bead
#'   population, or a long form carrying `marker`, `peak` and `value`
#' @param min_r2 fit quality below which a channel is reported as not calibrated
#' @return a data.frame, one row per channel
#' @export
fit_bead_calibration <- function(bead_exprs, channel_cols, assigned,
                                 min_r2 = 0.98) {
  if (is.null(bead_exprs) || is.null(assigned) || !nrow(assigned)) return(NULL)
  names(assigned) <- tolower(names(assigned))
  if (!"marker" %in% names(assigned))
    stop("the assigned-value table needs a 'marker' column naming the channel",
         call. = FALSE)

  long <- if (all(c("peak", "value") %in% names(assigned))) {
    assigned[, c("marker", "peak", "value")]
  } else {
    vcols <- setdiff(names(assigned), "marker")
    do.call(rbind, lapply(seq_len(nrow(assigned)), function(i) data.frame(
      marker = assigned$marker[i], peak = seq_along(vcols),
      value = suppressWarnings(as.numeric(unlist(assigned[i, vcols]))),
      stringsAsFactors = FALSE)))
  }
  long <- long[is.finite(long$value), , drop = FALSE]
  if (!nrow(long)) return(NULL)

  rows <- list()
  for (mk in unique(long$marker)) {
    av <- sort(long$value[long$marker == mk])
    # `[[` on a named vector errors rather than returning NULL for an absent
    # name, and a marker the bead file does not carry is an ordinary state, not
    # a fault: the assigned-value table describes the bead lot, not this panel.
    j <- if (mk %in% names(channel_cols)) channel_cols[[mk]] else NA_integer_
    row <- data.frame(marker = mk, n_assigned = length(av), n_peaks = 0L,
                      slope = NA_real_, intercept = NA_real_, r_squared = NA_real_,
                      max_residual_pct = NA_real_,
                      verdict = "channel absent from the bead file",
                      stringsAsFactors = FALSE)
    if (!is.null(j) && !is.na(j) && j <= ncol(bead_exprs)) {
      pk <- bead_peaks(bead_exprs[, j], n_expected = length(av))
      row$n_peaks <- length(pk %||% numeric(0))
      if (!is.null(pk) && length(pk) >= 2L) {
        n <- min(length(pk), length(av))
        # Both ascending, matched by rank. A bead set whose dimmest population
        # is lost in background gives fewer peaks than assigned values, so the
        # BRIGHTEST n are matched and the count is reported.
        px <- utils::tail(pk, n); ay <- utils::tail(av, n)
        fit <- stats::lm(ay ~ px)
        pred <- stats::fitted(fit)
        r2 <- suppressWarnings(summary(fit)$r.squared)
        resid_pct <- 100 * max(abs(pred - ay) / pmax(ay, 1e-9))
        row$slope <- unname(stats::coef(fit)[2])
        row$intercept <- unname(stats::coef(fit)[1])
        row$r_squared <- round(r2, 5)
        row$max_residual_pct <- round(resid_pct, 2)
        row$verdict <- if (!is.finite(r2) || r2 < min_r2)
          "fit too poor to calibrate" else if (n < length(av))
          paste0("calibrated on the brightest ", n, " of ", length(av),
                 " populations") else "calibrated"
      } else {
        row$verdict <- "fewer than two bead populations resolved"
      }
    }
    rows[[length(rows) + 1L]] <- row
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Apply a fitted calibration to a sample's linear channel values
#'
#' Only channels whose fit was accepted are converted. A channel that was not
#' calibrated is returned untouched, and the caller must not present the two as
#' the same units.
#'
#' @param exprs a sample's expression matrix, linear units
#' @param channel_cols named integer vector, marker to column
#' @param calib the table from [fit_bead_calibration()]
#' @return list(exprs = converted matrix, applied = character vector of markers)
#' @export
apply_bead_calibration <- function(exprs, channel_cols, calib) {
  applied <- character(0)
  if (is.null(calib) || !nrow(calib)) return(list(exprs = exprs, applied = applied))
  ok <- calib[grepl("^calibrated", calib$verdict) & is.finite(calib$slope), ,
              drop = FALSE]
  for (i in seq_len(nrow(ok))) {
    j <- channel_cols[[ok$marker[i]]]
    if (is.null(j) || is.na(j) || j > ncol(exprs)) next
    exprs[, j] <- ok$intercept[i] + ok$slope[i] * exprs[, j]
    applied <- c(applied, ok$marker[i])
  }
  list(exprs = exprs, applied = applied)
}
