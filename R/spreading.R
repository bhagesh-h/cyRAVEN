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
# practice rather than detector count (Mage et al. 2026, Cytometry A 109:436).
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
# single-stain controls: where those are supplied, --single-stains computes the
# published SSM from them (spillover_spreading_matrix()). Its other limit is that
# source-positive cells are often another cell type, with its own
# autofluorescence, so part of what it measures can be cell type rather than
# optics.
#
# A SCREEN FOR COMPENSATION AND UNMIXING ERRORS, beside it. Compensation is
# correct when the receiver's median does not move as the source gets brighter
# (Roederer 2001); a residual trend, up or down, is what a wrong reference
# control or a wrong spillover value leaves. `trend_rho` is the rank
# correlation between source and receiver among cells positive for the source
# and negative for the receiver, the region an NxN plot is read in. It is
# judged only there, never among cells negative for both, where a diagonal is
# expected from autofluorescence differences between cell types (Roederer et
# al. 2016). A rank correlation is unchanged by the transform, so it reads the
# same on any scale. CompensAID (Olsman et al. 2026) screens the stained sample
# in the same spirit, with a validated rule: it segments each pair's
# source-positive population and flags the pair when the Secondary Stain Index
# of the last segment falls below -1. The errors either screen finds originate
# in the reference controls; their input is the sample.

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
      ih <- which(hi)
      if (length(ih) > 5000L) ih <- ih[unique(round(seq(1, length(ih), length.out = 5000L)))]
      rho <- suppressWarnings(stats::cor(xv[ih], yv[ih], method = "spearman"))
      # The receiver's width at a given source level: deviations from the
      # median of its source decile. A trend widens the marginal spread above
      # and leaves this one alone, so it separates spreading from a
      # compensation error that the marginal width cannot.
      dec <- findInterval(xv[ih], stats::quantile(xv[ih], seq(0.1, 0.9, 0.1)))
      s_cond <- stats::mad(yv[ih] - stats::ave(yv[ih], dec, FUN = stats::median))
      rows[[length(rows) + 1L]] <- data.frame(
        source = src, receiver = rec,
        n_source_neg = sum(lo), n_source_pos = sum(hi),
        spread_source_neg = round(s_lo, 4),
        spread_source_pos = round(s_hi, 4),
        spreading = round(s_hi - s_lo, 4),
        spreading_ratio = round(s_hi / s_lo, 3),
        spreading_ratio_given_source = round(s_cond / s_lo, 3),
        trend_rho = round(rho, 3),
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
#' @param trend_flag absolute rank correlation at or above which a pair is
#'   flagged as a possible compensation or unmixing error. Default 0.3, a
#'   conventional moderate correlation, not a validated cut. A negative trend
#'   on a pair whose receiver is `flag_ratio` times wider at a given source
#'   level is reported as `trend_unresolved` instead, since truncation at the
#'   receiver's cut produces one from spreading alone.
#' @return list(pairs, receivers), each a data.frame, or NULL
#' @export
run_spreading_report <- function(pops, gates, thr_all = NULL, panel_of = NULL,
                                 max_samples = 8L, flag_ratio = 1.25,
                                 trend_flag = 0.3) {
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
    median_ratio_given_source = round(stats::median(d$spreading_ratio_given_source,
                                                    na.rm = TRUE), 3),
    median_trend_rho = round(stats::median(d$trend_rho, na.rm = TRUE), 3),
    stringsAsFactors = FALSE)))
  rownames(agg) <- NULL
  agg$substantial <- is.finite(agg$median_ratio) & agg$median_ratio >= flag_ratio
  # A moderate rank correlation (Cohen's 0.3) in the median sample: the
  # receiver drifts with the source where it should not. A screen, to send the
  # reader to that pair's NxN plot, not a verdict.
  #
  # SPREADING ALONE PULLS THE TREND NEGATIVE. The trend is read among cells
  # below the receiver's cut, and where the source spreads into the receiver,
  # brighter source cells have the wider receiver distribution and lose more of
  # its upper half to that cut, so the median of what remains falls with the
  # source though compensation is correct. A negative trend on a pair whose
  # receiver is substantially wider at a given source level therefore cannot be
  # told from overcompensation here (`trend_unresolved`); a positive one shows
  # against that pull and is still flagged. The width at a given source level,
  # not `substantial`, decides: a trend widens the marginal spread by itself,
  # and judged on that the strongest overcompensation would be called
  # unresolved. Single-stain controls (`--single-stains`) decide such a pair.
  trend <- is.finite(agg$median_trend_rho) & abs(agg$median_trend_rho) >= trend_flag
  agg$trend_unresolved <- trend & agg$median_trend_rho < 0 &
    is.finite(agg$median_ratio_given_source) &
    agg$median_ratio_given_source >= flag_ratio
  agg$possible_compensation_error <- trend & !agg$trend_unresolved
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
               n_trend_flags = sum(d$possible_compensation_error, na.rm = TRUE),
               n_trend_unresolved = sum(d$trend_unresolved, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(rec) <- NULL

  rec$fallback_rate <- NA_real_
  if (!is.null(thr_all) && nrow(thr_all) && "source" %in% names(thr_all)) {
    fb <- do.call(rbind, lapply(
      split(thr_all, paste(thr_all$panel, thr_all$marker, sep = "\r")),
      function(d) data.frame(panel = d$panel[1], receiver = d$marker[1],
                             fallback_rate = mean(no_minimum(d$source),
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

#' Spillover spreading matrix from single-stain controls
#'
#' The published measure (Nguyen et al. 2013): for a single-stain control of
#' fluorochrome i, split into its positive population P and negative population
#' N on its primary detector, the spreading coefficient into detector j is
#'
#'     SS_ij = sqrt(max(0, sigma(P, j)^2 - sigma(N, j)^2)) / sqrt(F(P, i) - F(N, i))
#'
#' where sigma is the robust standard deviation (84.13th minus 50th percentile)
#' of the compensated values and F the median. The square root of the intensity
#' difference normalises out how bright the control was, so the coefficient is a
#' property of the fluorochrome-detector pair on this instrument, comparable
#' with published values and between panels.
#'
#' The primary detector of each control is the one where it is brightest (the
#' highest 99th percentile); the positive and negative populations are split at the density minimum on that
#' detector, on an arcsinh scale with cofactor 150, or at the median where none
#' is found.
#'
#' @param stains named list of reads from [read_fcs_resolved()], compensated.
#' @param min_events fewest events each of P and N must hold. Default 200.
#' @return data.frame, one row per control and receiving detector: `control`,
#'   `source_detector`, `receiver_detector`, `n_positive`, `n_negative`,
#'   `delta_median`, `spread_positive`, `spread_negative`, `ss`. NULL when no
#'   control can be used.
#' @export
spillover_spreading_matrix <- function(stains, min_events = 200L) {
  rsd <- function(x) as.numeric(stats::quantile(x, 0.8413, na.rm = TRUE) -
                                  stats::median(x, na.rm = TRUE))
  rows <- list()
  for (nm in names(stains)) {
    rd <- stains[[nm]]
    ex <- rd$exprs; mc <- rd$marker_cols
    if (is.null(ex) || length(mc) < 2L) next
    det <- rd$channel_names[unname(mc)] %||% names(mc)
    # The fluorochrome's own detector is where it is brightest: spillover into
    # any other detector is a fraction of the primary signal.
    score <- vapply(unname(mc), function(j)
      stats::quantile(ex[, j], 0.99, na.rm = TRUE, names = FALSE), numeric(1))
    pj <- unname(mc)[which.max(score)]
    x <- ex[, pj]
    v <- density_valley(asinh(x / 150))
    cut <- if (is.finite(v)) 150 * sinh(v) else stats::median(x, na.rm = TRUE)
    P <- is.finite(x) & x > cut; N <- is.finite(x) & x <= cut
    if (sum(P) < min_events || sum(N) < min_events) next
    dF <- stats::median(x[P]) - stats::median(x[N])
    if (!is.finite(dF) || dF <= 0) next
    for (k in seq_along(mc)) {
      j <- unname(mc)[k]
      if (j == pj) next
      sP <- rsd(ex[P, j]); sN <- rsd(ex[N, j])
      rows[[length(rows) + 1L]] <- data.frame(
        control = nm, source_detector = det[match(pj, unname(mc))],
        receiver_detector = det[k], n_positive = sum(P), n_negative = sum(N),
        delta_median = round(dF, 2), spread_positive = round(sP, 3),
        spread_negative = round(sN, 3),
        ss = round(sqrt(max(0, sP^2 - sN^2)) / sqrt(dF), 4),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
