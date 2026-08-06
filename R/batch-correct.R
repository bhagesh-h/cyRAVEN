# =============================================================================
# cyRAVEN -- batch correction, opt-in and guarded
# =============================================================================
#
# THE POSITION THIS PACKAGE TAKES, AND WHY IT NOW SHIPS A CORRECTOR ANYWAY.
# batch_mixing_report() measures batch structure and deliberately stops there,
# because where acquisition batch is confounded with the biological group --
# one cohort run on one set of days, which is how clinical cohorts are normally
# collected -- "removing the batch" and "removing the finding" are the same
# operation and no algorithm can tell them apart.
#
# That argument is against correcting BLINDLY, not against correcting. When
# batch and cohort are separable, an uncorrected batch effect is its own source
# of false positives. So correction is available here, and the guard that made
# the original position worth holding is enforced in code rather than left to
# the reader: correct_batch() computes Cramer's V between batch and group and
# REFUSES above `max_cramers_v` unless explicitly overridden. A user who wants
# to proceed anyway has to say so, and the manifest records that they did.
#
# WHY QUANTILE ALIGNMENT AND NOT HARMONY. Harmony moves cells in PCA space so
# that batch labels mix, which changes coordinates but not the intensities any
# threshold is applied to. Aligning the marker distributions themselves keeps
# the corrected values on the same scale as the uncorrected ones, so a gate, a
# median or a cluster centroid computed afterwards still means what it did.
#
# WHAT IT REACHES, PRECISELY. The correction is applied to the marker matrix
# assembled for the embedding, so it reaches the UMAP, the unsupervised
# clustering, the gate-vs-cluster agreement, and anything else read from
# cells_umap.csv. It does NOT reach the per-sample MFI table, the population
# frequencies, or the differential-abundance and differential-state tests: those
# are computed from each sample's own transformed matrix and its own thresholds,
# both derived before this point and both already batch-local by construction.
#
# That split is deliberate rather than an omission. A per-sample threshold is
# already immune to a between-batch shift, so correcting the intensities it was
# derived from would move the data out from under a gate that had adapted to it.
# What a batch effect genuinely distorts is the SHARED space -- one embedding
# built from every sample at once -- and that is what is corrected here.
#
# The method is the simple, stated case of what CytoNorm does: per batch and per
# marker, map that batch's quantiles onto the pooled reference distribution by
# monotone interpolation. Being monotone matters -- it cannot reorder cells
# within a batch, so it cannot invent a population or change which side of a
# gate a cell falls on relative to its own batch-mates. It corrects location and
# spread, and nothing else.


#' Cramer's V between two categorical variables
#'
#' Used here to quantify how far batch and biological group overlap. 0 is
#' independent, 1 is one perfectly determined by the other.
#'
#' @param a,b character or factor vectors of equal length
#' @return numeric in 0..1, or NA when either has a single level
#' @keywords internal
cramers_v <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- factor(a[ok]); b <- factor(b[ok])
  if (nlevels(a) < 2L || nlevels(b) < 2L) return(NA_real_)
  tb <- table(a, b)
  n <- sum(tb)
  chi <- suppressWarnings(stats::chisq.test(tb, correct = FALSE)$statistic)
  as.numeric(sqrt((chi / n) / max(1L, min(nlevels(a), nlevels(b)) - 1L)))
}


#' Align one marker's distribution across batches
#'
#' WHAT: for each batch, build the monotone map that carries that batch's
#' quantiles onto the pooled reference quantiles, and apply it.
#'
#' WHY MONOTONE INTERPOLATION AND NOT A SHIFT OR A Z-SCORE: a shift corrects the
#' location of a distribution and leaves its shape, so a batch whose negative
#' population is wider stays wider and its gate still lands somewhere else. A
#' z-score assumes the spread is meaningful on both sides of zero, which is
#' false for a bimodal marker where the two modes have different variances.
#' Matching quantiles corrects the shape too, and because the map is monotone it
#' preserves the ORDER of cells within a batch -- the property that stops a
#' correction from inventing structure.
#'
#' @param x numeric vector for one marker, across all cells
#' @param batch batch label per cell
#' @param probs quantile grid used to build the map
#' @return numeric vector, corrected
#' @keywords internal
align_quantiles <- function(x, batch, probs = seq(0, 1, length.out = 101)) {
  ok <- is.finite(x) & !is.na(batch)
  if (sum(ok) < 10L) return(x)
  ref <- stats::quantile(x[ok], probs, na.rm = TRUE, names = FALSE)
  out <- x
  for (b in unique(batch[ok])) {
    i <- ok & batch == b
    if (sum(i) < 10L) next
    q <- stats::quantile(x[i], probs, na.rm = TRUE, names = FALSE)
    # Ties in the batch quantiles (a marker that is constant over part of its
    # range) make approx() ambiguous; jitter-free deduplication keeps the map a
    # function without perturbing the data.
    keep <- !duplicated(q)
    if (sum(keep) < 2L) next
    out[i] <- stats::approx(q[keep], ref[keep], xout = x[i],
                            rule = 2, ties = "ordered")$y
  }
  out
}


#' Correct batch effects in the transformed marker matrix
#'
#' WHAT IT DOES: aligns each marker's distribution across batches, after
#' checking that batch and biological group are separable enough for that to be
#' meaningful.
#'
#' WHAT IT REFUSES TO DO: correct when Cramer's V between batch and group
#' exceeds `max_cramers_v`. At that point the two are close to the same
#' variable, and any correction removes the effect being looked for. Passing
#' `force = TRUE` proceeds and records the decision in the returned object so it
#' reaches the run manifest.
#'
#' This is descriptive of intensity only. It does not touch the gate hierarchy,
#' which is derived per sample and is already batch-local by construction.
#'
#' @param tmat numeric matrix, cells x markers, already transformed
#' @param batch batch label per cell
#' @param group biological group per cell, for the confounding check
#' @param markers markers to correct; defaults to all columns
#' @param max_cramers_v refusal threshold
#' @param force proceed despite confounding
#' @return list(tmat, cramers_v, corrected, markers, reason)
#' @export
correct_batch <- function(tmat, batch, group = NULL, markers = NULL,
                          max_cramers_v = 0.6, force = FALSE) {
  tmat <- as.matrix(tmat)
  if (is.null(markers)) markers <- colnames(tmat)
  markers <- intersect(markers, colnames(tmat))
  nb <- length(unique(batch[!is.na(batch)]))
  if (nb < 2L)
    return(list(tmat = tmat, cramers_v = NA_real_, corrected = FALSE,
                markers = character(0),
                reason = "only one batch level -- nothing to correct"))

  v <- if (!is.null(group)) cramers_v(batch, group) else NA_real_
  if (!is.na(v) && v > max_cramers_v && !isTRUE(force)) {
    log_msg("  batch correction REFUSED: Cramer's V(batch, group) = ",
            round(v, 3), " > ", max_cramers_v)
    log_msg("    At this overlap, removing the batch and removing the ",
            "biological effect are the same operation. Re-run with ",
            "--force-batch-correction to proceed anyway.")
    return(list(tmat = tmat, cramers_v = v, corrected = FALSE,
                markers = character(0),
                reason = sprintf("refused: batch confounded with group (V = %.3f)", v)))
  }

  for (m in markers) tmat[, m] <- align_quantiles(tmat[, m], batch)
  reason <- if (!is.na(v) && v > max_cramers_v)
    sprintf("FORCED despite confounding (V = %.3f)", v)
  else sprintf("corrected across %d batches (V = %s)", nb,
               if (is.na(v)) "not assessed" else sprintf("%.3f", v))
  log_msg("  batch correction: ", reason, ", ", length(markers), " marker(s)")
  list(tmat = tmat, cramers_v = v, corrected = TRUE, markers = markers,
       reason = reason)
}
