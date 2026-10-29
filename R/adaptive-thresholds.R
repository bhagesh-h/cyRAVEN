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
#   * a tail gate, for a marker with one real mode and a positive tail (the
#     openCyto name for the case; the rule here is the mode plus robust SDs,
#     not gate_tail's slope rule, see tail_threshold())
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
#' The robust form of the conventional positivity rule, the negative
#' population's centre plus k of its standard deviations. The centre is the
#' dominant mode. Its spread is estimated twice, each estimate immune to a
#' different contamination, and the smaller is used, since contamination only
#' ever widens a spread:
#'
#' * the MAD of the half the positive tail cannot reach, mirrored about the
#'   mode, which a positive tail or shoulder on the cut side cannot inflate;
#' * the half-width at half maximum of the mode's flank on the cut side,
#'   converted to a standard deviation (HWHM / sqrt(2 ln 2)) after removing the
#'   kernel's own variance, which a second negative population on the far side
#'   cannot inflate. Peak width at half maximum is cytometry's established
#'   measure of a population's spread (the half-peak CV).
#'
#' The second exists because the first fails where the negative cells are
#' themselves two populations with different background, as monocytes and
#' lymphocytes are on PBMC: the mirrored half takes in the second hump and the
#' cut lands past a small positive shoulder, or past every cell. On one sepsis
#' cohort the MAD alone put CD19 beyond every cell in 3 of 30 samples. For a
#' single normal negative the two estimates agree.
#'
#' NOT openCyto's `gate_tail`, which places the cut where the density's first
#' derivative, past the reference peak's steepest descent, falls below a
#' tolerance (cytoUtils `.cytokine_cutpoint`). On a normal negative mode the two
#' land near each other. Under `auto_tol` the tolerance is 1% of the steepest
#' slope, which a normal reaches 3.57 SD out. The default `tol = 0.01` is an
#' absolute slope, so where it cuts depends on the data's units and the
#' population's share: 3.1 SD out of a unit-variance normal holding every event,
#' 2.5 SD at SD 2, 3.6 SD at SD 0.5. The slope rule also
#' depends on the kernel bandwidth; tried on the sepsis cohort it called 55% of
#' one sample's cells Vd2+ at a bandwidth where it found the B cells, and was not
#' adopted.
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
  i <- which.max(d$y)
  m <- d$x[i]
  half <- if (side == "upper") x[x <= m] else x[x >= m]
  s <- stats::mad(c(half, 2 * m - half), constant = 1.4826)
  sf <- flank_sd(d, i, side)
  if (is.finite(sf) && sf > 0 && (!is.finite(s) || s <= 0 || sf < s)) s <- sf
  if (!is.finite(s) || s <= 0) s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  if (side == "upper") m + k * s else m - k * s
}

#' Standard deviation of a density mode from its half-maximum width on one side
#'
#' The distance from the mode to where the density falls to half its height,
#' interpolated between grid points, is the half-width at half maximum; for a
#' normal population it is sqrt(2 ln 2) standard deviations. A kernel density
#' estimate of a normal population with SD `s` is normal with variance
#' `s^2 + bw^2`, so the kernel's variance is subtracted.
#'
#' @param d a [stats::density()] result.
#' @param i index of the mode in `d$y`.
#' @param side `"upper"` measures the flank above the mode, `"lower"` below.
#' @return the SD, or NA where the density never falls to half on that side.
#' @keywords internal
flank_sd <- function(d, i, side = c("upper", "lower")) {
  side <- match.arg(side)
  idx <- if (side == "upper") i:length(d$y) else i:1L
  j <- idx[which(d$y[idx] < d$y[i] / 2)[1]]
  if (is.na(j)) return(NA_real_)
  j0 <- if (side == "upper") j - 1L else j + 1L
  xh <- d$x[j0] + (d$y[i] / 2 - d$y[j0]) * (d$x[j] - d$x[j0]) / (d$y[j] - d$y[j0])
  sk <- abs(xh - d$x[i]) / sqrt(2 * log(2))
  sqrt(max(sk^2 - d$bw^2, 0))
}

#' The tail gate of a marker with no density minimum, where it can be placed
#'
#' Above the dominant mode, openCyto's default side, never below zero signal (a
#' negative population ending below zero puts the cut at zero), and only where
#' it leaves under 60% of the parent above it (more means the dominant mode was
#' not a negative one). NOTHING above the cut is an
#' answer rather than a failure: a viability dye in a sample with no dead
#' cells, a marker no cell in this sample carries. Shared by the default path
#' of [resolve_threshold()] and the last resort of [best_threshold()], so the
#' two cannot disagree on where a tail goes.
#'
#' WHY THE SIDE IS NOT GUESSED. A single mode does not say whether it is the
#' marker's negative or positive population. Its distance from zero signal
#' does not either: granulocyte autofluorescence and the spillover spread of
#' bright cells put a NEGATIVE mode well above zero, and a per-marker cofactor
#' small enough to resolve the negatives stretches that offset into many
#' standard deviations. A version that cut below any mode lying 3.5 of its SDs
#' above zero took the lower side 17 times in two runs on one 24-sample whole-blood cohort
#' scored against manual gating, and was wrong every time it could be checked:
#' Vd1, Vd2 and CD11c, all rare, called positive on 90-99% of CD45+ cells,
#' emptying CD4 and CD8 T cells (defined Vd1- Vd2-), and a CD16 cut below the
#' bulk that took classical monocytes to 0% and non-classical to 12% (manual
#' 6.4% and 0.1%). openCyto takes the side from the analyst for this reason.
#' A marker carried by most of its parent with no density minimum gets a cut
#' above its mode, flagged `needs_review`, unless the analyst declares the
#' other side, as openCyto's `side` argument does: `thresholds: <marker>:
#' {tail: lower}` in the configuration. A control, a fixed cut or a `parent:`
#' in which the marker separates are the alternatives.
#'
#' @param x transformed values within the parent gate.
#' @param floor the transformed value of zero signal, `-Inf` for scatter.
#' @param side `"upper"` (default), or `"lower"` where the analyst has declared
#'   the dominant mode positive; the cut is then the mode minus 3.5 SD and must
#'   leave over 40% of the parent above it.
#' @return the cut, or NA where no tail gate can be placed.
#' @keywords internal
tail_gate_cut <- function(x, floor = -Inf, side = c("upper", "lower")) {
  side <- match.arg(side)
  x <- x[is.finite(x)]
  tg <- tryCatch(tail_threshold(x, side = side), error = function(e) NA_real_)
  if (!is.finite(tg)) return(NA_real_)
  # A negative population ending below zero signal puts the cut at zero
  # signal, not below it: no stained cell sits below zero.
  tg <- max(tg, floor)
  p <- mean(x > tg)
  # Above a negative mode under 60% may be positive; below a declared positive
  # mode over 40% must be, or the mode was not the positive one.
  if (if (side == "upper") p < 0.60 else p > 0.40) tg else NA_real_
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
#' @param floor no candidate at or below this value is considered; see
#'   [density_valley()]. Default `-Inf`.
#' @param prefer `"deepest"` (default) scores every candidate by the trough it
#'   sits in. `"lowest"`, for a gate that removes only the negative mode, takes
#'   the lowest density valley that sits in a real trough, and falls back to the
#'   scored choice when no valley does; see [density_valley()].
#' @param tail `"upper"` (default, also meaning nothing was declared) offers a
#'   tail gate on either side of the dominant mode and lets the trough decide;
#'   `"lower"`, the analyst's declaration that the dominant mode is positive
#'   (`thresholds: <marker>: {tail: lower}`), offers only the lower one and
#'   takes it as the last resort too, as the default path does. Without this
#'   the declaration was ignored under `--adaptive-gates`.
#' @return `list(threshold, source, quality)`, where `source` names the method
#'   that won so a reader can see how each cut was placed.
#' @keywords internal
best_threshold <- function(x, ladder = c(4, 6, 8, 12, 16, 24, 32, 48, 64),
                           bins = 400L, min_quality = 0.05, fallback_q = 0.90,
                           floor = -Inf, prefer = c("deepest", "lowest"),
                           tail = c("upper", "lower")) {
  prefer <- match.arg(prefer)
  tail <- match.arg(tail)
  x <- x[is.finite(x)]
  if (length(x) < 500L)
    return(list(threshold = NA_real_, source = "too_few_events", quality = NA_real_))

  cand <- list()
  # A candidate at or below `floor` (zero signal) is not a cut; see density_valley().
  add <- function(nm, cut) if (is.numeric(cut) && length(cut) == 1L && is.finite(cut) &&
                               cut > floor)
    cand[[length(cand) + 1L]] <<- list(source = nm, cut = as.numeric(cut))

  for (s in ladder) {
    if ((2L * s + 1L) >= bins) next
    v <- tryCatch(density_valley(x, bins = bins, smooth = s, floor = floor, prefer = prefer,
                                 details = TRUE),
                  error = function(e) NULL)
    if (!is.null(v) && is.finite(v$cut)) add(paste0("valley_s", s), v$cut)
  }
  if (tail == "upper")
    add("tailgate_upper", tryCatch(tail_threshold(x, side = "upper"), error = function(e) NA_real_))
  add("tailgate_lower", tryCatch(tail_threshold(x, side = "lower"), error = function(e) NA_real_))
  add("otsu",           tryCatch(otsu_threshold(x),                 error = function(e) NA_real_))

  if (length(cand)) {
    q <- vapply(cand, function(z) gate_gap_quality(x, z$cut), numeric(1))
    if (prefer == "lowest") {
      ok <- which(startsWith(vapply(cand, `[[`, "", "source"), "valley") &
                  is.finite(q) & q > min_quality)
      if (length(ok)) {
        b <- cand[[ok[which.min(vapply(cand[ok], `[[`, 0, "cut"))]]]
        return(list(threshold = b$cut, source = b$source, quality = q[[match(list(b), cand)]]))
      }
    }
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
  tu <- tail_gate_cut(x, floor, tail)
  if (is.finite(tu))
    return(list(threshold = tu, source = "tailgate", quality = NA_real_))
  # Floored and untied as on the default path: a quantile below zero signal
  # would call events positive that no stain put there.
  list(threshold = untie_cut(x, max(as.numeric(stats::quantile(x, fallback_q, na.rm = TRUE)),
                                    floor), 1 - fallback_q),
       source = "quantile_fallback", quality = NA_real_)
}
