# SECTION 1b -- ACQUISITION-TIME QUALITY CONTROL
# =============================================================================
#
# WHY THIS FILE EXISTS. Every number this package reports for a sample is derived
# from the pooled events of that sample: one cofactor, one threshold per marker,
# one frequency per population. That is correct only if the instrument was doing
# the same thing throughout the acquisition.
#
# It frequently was not. A partial clog, an air bubble, a pressure change or a
# slow drift in laser power moves the event rate or the signal part-way through a
# file. A file that was two different instruments over its run produces one
# threshold that suits neither half, and nothing downstream can detect it: the
# threshold's own uncertainty widens a little, staining QC still passes, and the
# within-run peer check in threshold_scale_qc.csv compares the sample against
# other samples rather than against itself.
#
# WHAT IS MEASURED. The Time channel is binned into equal-width intervals, and
# two quantities are tracked across them:
#
#   event rate  -- events per bin. A clog is a drop, a bubble is a spike. Bins
#     are equal in TIME rather than equal in count, because a bin defined to hold
#     a fixed number of events cannot show that the rate changed.
#
#   channel location -- the median of each marker within the bin, on the analysis
#     scale. This catches drift and step changes that leave the rate untouched,
#     which is what a laser or detector problem looks like.
#
# Both are judged by robust z against the across-bin median, so a file whose
# signal legitimately varies a little is not flagged for it.
#
# WHAT IS DELIBERATELY NOT DONE. Nothing is removed by default. Removing events
# changes every count, every frequency and every threshold in the file, and doing
# that silently is the outcome this codebase treats as worst. The run reports the
# flagged windows and, separately, how far each population's frequency would move
# if they were excluded. A file with a visible clog that moves no population by
# more than its own gate uncertainty does not need re-acquiring, and that
# judgement is not one a QC tool should make on the analyst's behalf.
# `--drop-unstable-events` performs the removal, and is recorded in the manifest.
#
# RELATION TO PeacoQC AND flowCut. Both solve this problem well and both would be
# a Bioconductor dependency. The detector here is a binned robust location test,
# which is the shared core of those methods, written directly so the container
# stays as it is. Emmaneel et al. 2022, Cytometry A 101:325 is the reference for
# the approach and the comparison against flowAI, flowClean and flowCut.

#' Locate the Time channel in a read FCS file
#'
#' The reading stage keeps Time out of `marker_cols` and `scatter_cols` on
#' purpose, so it cannot drive a gate or an embedding. This finds it in the raw
#' matrix for the one use it does have.
#'
#' @param rd a list from [read_fcs_resolved()]
#' @return the column index, or NA when the file carries no Time channel
#' @export
find_time_column <- function(rd) {
  nm <- rd$channel_names %||% colnames(rd$exprs)
  if (is.null(nm)) return(NA_integer_)
  hit <- which(grepl("^time", nm, ignore.case = TRUE))
  if (!length(hit)) return(NA_integer_)
  as.integer(hit[1])
}

#' Flag acquisition intervals whose rate or signal departs from the file
#'
#' @param time numeric acquisition time per event, in the file's own units
#' @param mat matrix of channel values on the analysis scale, one column per
#'   channel, same number of rows as `time`
#' @param n_bins number of equal-width time intervals
#' @param mad_k robust z above which a bin is flagged
#' @param min_events_per_bin bins holding fewer events than this are not judged
#'   on their channel medians, because a median over a handful of events is
#'   noise; they are still judged on their rate, which is the point
#' @return list(bins = data.frame, flagged = logical per event, reason)
#' @export
detect_time_anomalies <- function(time, mat, n_bins = 40L, mad_k = 5,
                                  min_events_per_bin = 50L) {
  none <- function(reason) list(bins = NULL, flagged = rep(FALSE, length(time)),
                                reason = reason)
  time <- as.numeric(time)
  ok <- is.finite(time)
  if (sum(ok) < 200L) return(none("fewer than 200 events with a finite Time value"))
  rng <- range(time[ok])
  if (!is.finite(diff(rng)) || diff(rng) <= 0)
    return(none("Time does not advance; the channel is present but constant"))

  n_bins <- max(4L, as.integer(n_bins))
  brk <- seq(rng[1], rng[2], length.out = n_bins + 1L)
  bin <- .bincode(time, brk, include.lowest = TRUE)
  cnt <- tabulate(bin, nbins = n_bins)

  # Rate. A clog is a sustained drop and a bubble a spike; both are departures
  # from the file's own typical bin occupancy.
  rate_med <- stats::median(cnt)
  rate_mad <- max(stats::mad(cnt), 1)
  rate_z <- abs(cnt - rate_med) / rate_mad

  # Channel location. One robust z per channel per bin, reduced to the worst
  # channel, so a single misbehaving detector is not diluted by the others.
  chan_z <- matrix(0, nrow = n_bins, ncol = max(1L, ncol(mat)))
  worst <- rep(NA_character_, n_bins)
  if (!is.null(mat) && ncol(mat) > 0L) {
    for (j in seq_len(ncol(mat))) {
      v <- mat[, j]
      med_by <- vapply(seq_len(n_bins), function(b) {
        idx <- which(bin == b & is.finite(v))
        if (length(idx) < min_events_per_bin) NA_real_ else stats::median(v[idx])
      }, numeric(1))
      ref <- stats::median(med_by, na.rm = TRUE)
      # THE SCALE FLOOR, AND WHY IT IS NOT A SMALL CONSTANT. The natural
      # denominator is the spread of the bin medians, but for a channel whose
      # medians happen to agree closely that spread approaches zero and every
      # trivial difference becomes a large z. Flooring at an absolute epsilon
      # does not help, because "small" is not comparable between channels or
      # scales. The floor is therefore a fraction of the channel's OWN
      # event-level spread: a bin median must move by at least a twentieth of
      # the variation the channel shows anyway before it counts as a departure.
      ev_mad <- stats::mad(v, na.rm = TRUE)
      s <- max(stats::mad(med_by, na.rm = TRUE), 0.05 * ev_mad)
      # A channel with no event-level spread at all carries no information about
      # stability, so it contributes nothing rather than everything.
      z <- if (is.finite(s) && s > 0) abs(med_by - ref) / s else rep(0, n_bins)
      z[!is.finite(z)] <- 0
      chan_z[, j] <- z
    }
    cn <- colnames(mat) %||% paste0("ch", seq_len(ncol(mat)))
    worst <- cn[apply(chan_z, 1, which.max)]
  }
  chan_worst_z <- apply(chan_z, 1, max)

  flagged_bin <- (rate_z > mad_k) | (chan_worst_z > mad_k)
  bins <- data.frame(
    bin = seq_len(n_bins),
    time_from = brk[-length(brk)], time_to = brk[-1],
    n_events = cnt,
    rate_z = round(rate_z, 2),
    worst_channel = worst,
    worst_channel_z = round(chan_worst_z, 2),
    flagged = flagged_bin,
    stringsAsFactors = FALSE)

  list(bins = bins,
       flagged = !is.na(bin) & flagged_bin[pmax(bin, 1L)],
       reason = NA_character_)
}

#' Acquisition-time QC for one sample
#'
#' Runs before gating, so it uses log10 of the marker channels rather than the
#' fitted analysis transform. The test is a robust z on per-bin medians and both
#' scales are monotone, so the detection is equivalent; log10 needs no cofactor
#' and is therefore available at the point where a decision about the events is
#' still possible.
#'
#' @param rd a list from [read_fcs_resolved()]
#' @param tmat transformed marker matrix for the same events. When NULL, log10 of
#'   the marker columns is used.
#' @param n_bins,mad_k passed to [detect_time_anomalies()]
#' @return list(summary = one-row data.frame, bins, flagged)
#' @keywords internal
acquisition_qc_one <- function(rd, tmat = NULL, n_bins = 40L, mad_k = 5) {
  if (is.null(tmat) && length(rd$marker_cols %||% integer(0))) {
    tmat <- vapply(rd$marker_cols, function(j) log10(pmax(rd$exprs[, j], 1)),
                   numeric(nrow(rd$exprs)))
    colnames(tmat) <- names(rd$marker_cols)
  }
  tcol <- find_time_column(rd)
  base <- data.frame(
    sample_id = rd$sample_id, n_events = rd$n_events %||% NA_integer_,
    n_bins = NA_integer_, n_bins_flagged = NA_integer_,
    pct_events_flagged = NA_real_, worst_channel = NA_character_,
    max_rate_z = NA_real_, max_channel_z = NA_real_,
    verdict = "no Time channel", stringsAsFactors = FALSE)
  if (is.na(tcol)) return(list(summary = base, bins = NULL, flagged = NULL))

  res <- detect_time_anomalies(rd$exprs[, tcol], tmat, n_bins = n_bins,
                               mad_k = mad_k)
  if (is.null(res$bins)) {
    base$verdict <- res$reason
    return(list(summary = base, bins = NULL, flagged = NULL))
  }
  b <- res$bins
  nf <- sum(b$flagged)
  pct <- 100 * sum(res$flagged) / max(1L, length(res$flagged))
  base$n_bins <- nrow(b)
  base$n_bins_flagged <- nf
  base$pct_events_flagged <- round(pct, 3)
  base$max_rate_z <- max(b$rate_z, na.rm = TRUE)
  base$max_channel_z <- max(b$worst_channel_z, na.rm = TRUE)
  base$worst_channel <- if (nf) {
    w <- b$worst_channel[b$flagged]
    names(sort(table(w[!is.na(w)]), decreasing = TRUE))[1] %||% NA_character_
  } else NA_character_
  base$verdict <- if (!nf) "stable" else if (pct < 5) "unstable, minor" else "unstable"
  b$sample_id <- rd$sample_id
  list(summary = base, bins = b, flagged = res$flagged)
}

#' Run acquisition-time QC across a cohort
#'
#' @param reads named list of [read_fcs_resolved()] results
#' @param pops per-sample scoring results, used for the transformed matrix
#' @param n_bins,mad_k detector settings
#' @return list(summary, bins, flagged) where `flagged` is a named list of
#'   per-event logical vectors
#' @export
run_acquisition_qc <- function(reads, pops = NULL, n_bins = 40L, mad_k = 5) {
  if (is.null(reads) || !length(reads)) return(NULL)
  srows <- list(); brows <- list(); flags <- list()
  for (s in names(reads)) {
    rd <- reads[[s]]
    # The fitted transform is used when a scored run supplies it; otherwise
    # acquisition_qc_one() falls back to log10, which is what makes this runnable
    # before gating.
    tm <- pops[[s]]$tmat
    if (!is.null(tm) && nrow(tm) != nrow(rd$exprs)) tm <- NULL
    r <- tryCatch(acquisition_qc_one(rd, tm, n_bins = n_bins, mad_k = mad_k),
                  error = function(e) NULL)
    if (is.null(r)) next
    srows[[length(srows) + 1L]] <- r$summary
    if (!is.null(r$bins)) brows[[length(brows) + 1L]] <- r$bins
    if (!is.null(r$flagged)) flags[[s]] <- r$flagged
  }
  if (!length(srows)) return(NULL)
  list(summary = do.call(rbind, srows),
       bins = if (length(brows)) do.call(rbind, brows) else NULL,
       flagged = flags)
}

#' How far each population moves if the flagged windows are excluded
#'
#' The quantity that turns a QC flag into a decision. A file with a visible
#' anomaly that moves no population by more than its own gate uncertainty does
#' not need re-acquiring.
#'
#' @param tmat transformed marker matrix for one sample
#' @param thr named threshold vector
#' @param parent logical parent-gate mask
#' @param spec population specification
#' @param keep logical, TRUE for events to retain
#' @return data.frame with the reported and cleaned percentage per population
#' @export
frequency_delta_if_cleaned <- function(tmat, thr, parent, spec, keep) {
  if (is.null(tmat) || !length(keep) || all(keep)) return(NULL)
  pct_at <- function(p) {
    if (!sum(p)) return(NULL)
    hi <- derive_intermediate_bounds(tmat, thr, p, spec)
    sp <- score_populations(tmat, thr, p, spec, hi_thr = hi)
    vapply(sp$masks, function(m) 100 * sum(m) / sum(p), numeric(1))
  }
  a <- pct_at(parent)
  b <- pct_at(parent & keep)
  if (is.null(a) || is.null(b)) return(NULL)
  pops <- intersect(names(a), names(b))
  if (!length(pops)) return(NULL)
  data.frame(population = pops,
             pct_reported = round(unname(a[pops]), 4),
             pct_if_cleaned = round(unname(b[pops]), 4),
             pct_delta_if_cleaned = round(unname(b[pops] - a[pops]), 4),
             row.names = NULL, stringsAsFactors = FALSE)
}
