# SECTION 3 -- DENSITY VALLEY AND THRESHOLD RESOLUTION
# =============================================================================

#' Find the deepest density valley between two modes
#'
#' WHAT: smooths a histogram, locates peaks above `peak_frac` of the maximum,
#' and returns the valley between the peak pair maximising depth * min(height).
#' Returns NA when the distribution is unimodal -- the honest answer, which then
#' routes threshold resolution to the unstained-control fallback.
#' WHY: for a bimodal marker the valley is the least-arbitrary cutoff available
#' and adapts to each sample's own staining intensity.
#' @param min_rel_depth minimum valley depth AS A FRACTION of the smaller of the
#'   two flanking peaks. This is the parameter that makes minority positive
#'   populations findable -- see the note on scoring below.
#' @param min_upper_frac the upper mode must hold at least this fraction of events,
#'   so a cut cannot be placed out in the extreme tail where a handful of events
#'   form a bump.
#' @param x A vector of values.
#' @param bins The bins. Default `220L`.
#' @param smooth The smooth. Default `4`.
#' @param peak_frac The peak frac. Default `0.02`.
#' @param min_gap_frac The min gap frac. Default `0.06`.
#' @param range_q The range q. Default `c(0.001, 0.999)`.
#' @param details when TRUE return a list carrying the cut plus the relative
#'   depth of the valley and why no cut was found, instead of the bare cut.
#'   The default returns exactly what it always has: one number or `NA_real_`.
#'   [threshold_uncertainty()] uses the list form, because how deep a valley is
#'   determines how far the cut can move under resampling, and recomputing the
#'   smoothed histogram elsewhere would let the two drift apart.
#' @export
density_valley <- function(x, bins = 220L, smooth = 4, peak_frac = 0.02,
                           min_gap_frac = 0.06, range_q = c(0.001, 0.999),
                           min_rel_depth = 0.30, min_upper_frac = 0.002,
                           details = FALSE) {
  none <- function(reason)
    if (details) list(cut = NA_real_, rel_depth = NA_real_, n_peaks = NA_integer_,
                      reason = reason) else NA_real_
  x <- x[is.finite(x)]
  if (length(x) < 500L) return(none("fewer than 500 finite values"))
  rg <- quantile(x, range_q, na.rm = TRUE)
  if (!all(is.finite(rg)) || diff(rg) <= 0) return(none("degenerate range"))
  h <- graphics::hist(x[x >= rg[1] & x <= rg[2]], breaks = seq(rg[1], rg[2], length.out = bins + 1L),
            plot = FALSE)
  y <- h$counts
  # moving-average smoothing
  k <- max(1L, as.integer(smooth))
  ys <- stats::filter(y, rep(1 / (2 * k + 1), 2 * k + 1), sides = 2)
  ys[is.na(ys)] <- y[is.na(ys)]
  ctr <- h$mids
  mx <- max(ys)
  pk <- which(c(FALSE, ys[-1] >= ys[-length(ys)]) &
              c(ys[-length(ys)] >= ys[-1], FALSE) & ys > peak_frac * mx)
  if (length(pk) < 2L) return(none("unimodal: fewer than two peaks"))
  min_gap <- min_gap_frac * diff(range(ctr))
  tot <- sum(ys[is.finite(ys)])
  best <- NULL
  for (i in seq_len(length(pk) - 1L)) for (j in (i + 1L):length(pk)) {
    a <- pk[i]; b <- pk[j]
    if (ctr[b] - ctr[a] < min_gap) next
    vi <- a + which.min(ys[a:b]) - 1L
    lo <- min(ys[a], ys[b])
    depth <- lo - ys[vi]
    if (depth <= 0) next
    # RELATIVE depth: how completely the two modes separate, independent of how
    # tall they are. The previous scoring used depth * min(height), which is an
    # ABSOLUTE quantity and therefore preferred whichever peak pair contained the
    # most events. That is wrong whenever the positive population is a minority:
    # in a sample where CD45+ cells are ~5% of events, a shallow ripple inside the
    # enormous negative peak scores higher than the clean, deep separation between
    # negative and positive, so the threshold lands INSIDE the negative
    # distribution and calls background positive. Because CD45 sits at the root of
    # the gate hierarchy, that one misplaced cut then corrupts every population
    # beneath it, and it does so silently -- the valley is real, just the wrong one.
    rel <- depth / lo
    if (rel < min_rel_depth) next
    # The upper mode must be a population, not a tail artefact.
    if (sum(ys[vi:length(ys)], na.rm = TRUE) < min_upper_frac * tot) next
    # Among candidates that genuinely separate, prefer the deepest separation and
    # break ties toward the wider peak spacing (a better-resolved pair).
    score <- rel + 0.05 * (ctr[b] - ctr[a]) / diff(range(ctr))
    if (is.null(best) || score > best$score)
      best <- list(score = score, cut = ctr[vi], rel_depth = rel)
  }
  if (is.null(best)) return(none("no peak pair separates enough"))
  if (details) list(cut = best$cut, rel_depth = best$rel_depth,
                    n_peaks = length(pk), reason = NA_character_) else best$cut
}

#' Resolve a threshold for one marker, recording where the value came from
#'
#' Priority (per the agreed contract):
#'   1. explicit config value          -- always wins
#'   2. density valley within parent   -- bimodal markers
#'   3. unstained-control quantile     -- unimodal markers with a control present
#'   4. quantile fallback              -- flagged needs_review
#' A valley falling BELOW the control-derived value is rejected: it lies inside
#' the background distribution and would call noise positive.
#' @param marker Marker name.
#' @param x_parent The x parent.
#' @param cfg_value The cfg value.
#' @param control_x The control x.
#' @param control_q The control q. Default `0.995`.
#' @param fallback_q The fallback q. Default `0.9`.
#' @export
resolve_threshold <- function(marker, x_parent, cfg_value = NULL,
                              control_x = NULL, control_q = 0.995,
                              fallback_q = 0.90) {
  if (!is.null(cfg_value) && is.finite(cfg_value))
    return(list(threshold = cfg_value, source = "config", needs_review = FALSE))
  ctrl <- if (!is.null(control_x) && length(control_x) > 100L)
    as.numeric(quantile(control_x, control_q, na.rm = TRUE)) else NA_real_
  v <- density_valley(x_parent)
  if (is.finite(v)) {
    if (is.finite(ctrl) && v < ctrl)
      return(list(threshold = ctrl, source = "control_q995_valley_rejected",
                  needs_review = FALSE))
    return(list(threshold = v, source = "valley", needs_review = FALSE))
  }
  if (is.finite(ctrl))
    return(list(threshold = ctrl, source = "control_q995", needs_review = FALSE))
  list(threshold = as.numeric(quantile(x_parent, fallback_q, na.rm = TRUE)),
       source = "quantile_fallback", needs_review = TRUE)
}

# =============================================================================
