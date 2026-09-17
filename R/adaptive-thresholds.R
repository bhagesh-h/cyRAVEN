# =============================================================================
# ADAPTIVE THRESHOLD SELECTION
# =============================================================================
#
# WHY THIS EXISTS. `density_valley()` answers one question well: where is the
# deepest separation between two modes. It answers with NA when the marker is
# unimodal, and `resolve_threshold()` then falls back to a fixed quantile of the
# parent. On a spectral panel that fallback is not rare. On the cohort this was
# written for, 65% of every sample-by-marker cut came back `quantile_fallback`,
# which means two thirds of the gates in the run were placed by the constant
# `fallback_q` rather than by the data, and the viability gate was skipped in all
# 20 samples because L/D never resolved.
#
# The diagnosis was not that the markers are unimodal. It is that the histogram
# `density_valley()` smooths -- 220 bins, a 9-bin moving average -- is badly
# under-smoothed for these data: L/D showed 32 "peaks" in one sample and CD45
# showed 10, nearly all of them noise. Widen the kernel and CD45 resolves to
# exactly two modes.
#
# WHAT THIS ADDS, and it is deliberately the openCyto vocabulary rather than a
# new one (Finak and others, PLoS Comput Biol 2014;10:e1003806):
#
#   * mindensity over a BANDWIDTH LADDER rather than one fixed smoothing
#   * tailgate, for a marker with one real mode and a positive tail
#   * Otsu, for a marker whose two classes are of comparable size
#
# and then one selector over all of them, so no marker needs configuring by hand.
#
# HOW A CANDIDATE IS JUDGED. Every method proposes a cut; each cut is scored the
# same way, by how deep a density gap it sits in. That is method-agnostic, needs
# no per-marker prior, and measures the property a gate is supposed to have. A
# cut through a real trough scores near 1; a cut through the middle of one
# population scores near 0.
#
# WHAT IT WILL NOT DO. It cannot find a population that was not stained. On this
# cohort one sample carries CD3 at 0.3% positive while its CD45 stains normally,
# and no threshold rule recovers T cells from a channel that holds none. The
# selector reports what is there.

#' How deep a density gap does a candidate cut sit in
#'
#' @param x transformed values for one marker within the parent gate.
#' @param cut candidate threshold.
#' @param adjust kernel bandwidth multiplier. Default `2`.
#' @return a number in `(0, 1)`: 1 is a cut at the floor of a clean trough, 0 a
#'   cut through the middle of a single mode. `-Inf` for a degenerate split.
#' @keywords internal
gate_gap_quality <- function(x, cut, adjust = 2) {
  if (!is.finite(cut)) return(-Inf)
  x <- x[is.finite(x)]
  if (length(x) < 200L) return(-Inf)
  d <- stats::density(x, n = 1024, adjust = adjust)
  i <- which.min(abs(d$x - cut))
  if (i < 10L || i > (length(d$x) - 10L)) return(-Inf)
  fl <- min(max(d$y[1:i]), max(d$y[i:length(d$y)]))
  if (!is.finite(fl) || fl <= 0) return(-Inf)
  p <- mean(x > cut)
  # A gate keeping essentially everything is not separating anything, whatever
  # the density happens to do at that point. Without this bound a ripple in the
  # left tail scores well and the cut lands below the whole distribution.
  if (p < 0.0005 || p > 0.98) return(-Inf)
  1 - d$y[i] / fl
}

#' Cut a robust spread above (or below) the dominant mode
#'
#' openCyto's `tailgate`. The spread is estimated by mirroring the half of the
#' distribution the positive tail cannot reach, so a long tail does not inflate
#' the width of the mode it is being measured against.
#'
#' @param x transformed values.
#' @param k how many robust standard deviations from the mode. Default `3.5`.
#' @param side `"upper"` for a minority-positive marker, `"lower"` for a marker
#'   whose dominant mode is itself the positive population.
#' @param adjust kernel bandwidth multiplier. Default `2`.
#' @return the cut.
#' @keywords internal
tail_threshold <- function(x, k = 3.5, side = c("upper", "lower"), adjust = 2) {
  side <- match.arg(side)
  x <- x[is.finite(x)]
  if (length(x) < 200L) return(NA_real_)
  d <- stats::density(x, n = 1024, adjust = adjust)
  m <- d$x[which.max(d$y)]
  half <- if (side == "upper") x[x <= m] else x[x >= m]
  s <- stats::mad(c(half, 2 * m - half), constant = 1.4826)
  if (!is.finite(s) || s <= 0) s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  if (side == "upper") m + k * s else m - k * s
}

#' Otsu's between-class variance threshold
#'
#' Parameter-free and direction-neutral, and right where the two classes are of
#' comparable size. It splits a single population when the positive fraction is
#' small, which is why it is one candidate among several rather than the rule.
#'
#' @param x transformed values.
#' @param nbins histogram resolution. Default `512L`.
#' @return the cut.
#' @keywords internal
otsu_threshold <- function(x, nbins = 512L) {
  x <- x[is.finite(x)]
  if (length(x) < 200L) return(NA_real_)
  rg <- stats::quantile(x, c(0.001, 0.999), na.rm = TRUE)
  if (!all(is.finite(rg)) || diff(rg) <= 0) return(NA_real_)
  xx <- x[x >= rg[1] & x <= rg[2]]
  h <- graphics::hist(xx, breaks = seq(rg[1], rg[2], length.out = nbins + 1L),
                      plot = FALSE)
  p <- h$counts / sum(h$counts)
  w0 <- cumsum(p); w1 <- 1 - w0
  mu <- cumsum(p * h$mids); muT <- sum(p * h$mids)
  h$mids[which.max((muT * w0 - mu)^2 / pmax(w0 * w1, 1e-12))]
}

#' Choose a threshold by trying every method and scoring the cuts alike
#'
#' @param x transformed values for one marker within the parent gate.
#' @param ladder smoothing windows for the mindensity sweep, widest last.
#' @param bins histogram resolution for the mindensity sweep. Default `400L`.
#' @param min_quality a candidate must sit in a trough at least this deep to be
#'   preferred over the tail rule. Default `0.05`.
#' @param fallback_q quantile used only when nothing else is usable.
#' @return `list(threshold, source, quality)`, where `source` names the method
#'   that won so a reader can see how each cut was placed.
#' @keywords internal
best_threshold <- function(x, ladder = c(4, 6, 8, 12, 16, 24, 32, 48, 64),
                           bins = 400L, min_quality = 0.05, fallback_q = 0.90) {
  x <- x[is.finite(x)]
  if (length(x) < 500L)
    return(list(threshold = NA_real_, source = "too_few_events", quality = NA_real_))

  cand <- list()
  add <- function(nm, cut) if (is.numeric(cut) && length(cut) == 1L && is.finite(cut))
    cand[[length(cand) + 1L]] <<- list(source = nm, cut = as.numeric(cut))

  for (s in ladder) {
    if ((2L * s + 1L) >= bins) next
    v <- tryCatch(density_valley(x, bins = bins, smooth = s, details = TRUE),
                  error = function(e) NULL)
    if (!is.null(v) && is.finite(v$cut)) add(paste0("valley_s", s), v$cut)
  }
  add("tailgate_upper", tryCatch(tail_threshold(x, side = "upper"), error = function(e) NA_real_))
  add("tailgate_lower", tryCatch(tail_threshold(x, side = "lower"), error = function(e) NA_real_))
  add("otsu",           tryCatch(otsu_threshold(x),                 error = function(e) NA_real_))

  if (length(cand)) {
    q <- vapply(cand, function(z) gate_gap_quality(x, z$cut), numeric(1))
    if (any(is.finite(q)) && max(q, na.rm = TRUE) > min_quality) {
      b <- cand[[which.max(q)]]
      return(list(threshold = b$cut, source = b$source, quality = max(q, na.rm = TRUE)))
    }
  }

  # NO REAL TROUGH ANYWHERE. That is a fact about the staining rather than a
  # failure to look hard enough, and the honest answer is a cut anchored to this
  # sample's own mode and spread. The quantile declares a fixed share positive
  # whatever the data look like, which is the fabrication this module exists to
  # remove, so it is the last resort and still reported as one.
  tu <- tryCatch(tail_threshold(x, side = "upper"), error = function(e) NA_real_)
  if (is.finite(tu)) {
    p <- mean(x > tu)
    if (p > 0.0002 && p < 0.60)
      return(list(threshold = tu, source = "tailgate", quality = NA_real_))
  }
  list(threshold = as.numeric(stats::quantile(x, fallback_q, na.rm = TRUE)),
       source = "quantile_fallback", quality = NA_real_)
}
