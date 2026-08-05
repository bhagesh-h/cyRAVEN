# SECTION 2 -- AUTO-DERIVED TRANSFORM
# =============================================================================

#' Derive the arcsinh cofactor from the data
#'
#' WHAT: for each marker, bisect the cofactor so the inter-quartile range of the
#' transformed background bulk equals `target_iqr`; the panel cofactor is the
#' median across markers.
#' WHY: the cofactor sets how much the near-zero region is expanded. Spectrally
#' unmixed data has a wide negative/background band; too small a cofactor (the
#' template's default of 5) expands that noise until it dominates all distances.
#' On the test batch this recovers ~150 for the baseline panel -- derived, not
#' assumed, so a different instrument or gain setting gets its own value.
#' @param ex Numeric expression matrix, events x channels.
#' @param marker_cols Named integer vector mapping marker symbol to column index.
#' @param target_iqr Interquartile range the transformed background is bisected towards. Default `0.5`.
#' @param lo Lower bound of the search interval. Default `1`.
#' @param hi Upper bound of the search interval. Default `5000`.
#' @param max_markers Ceiling on the number of markers used. Default `40L`.
#' @export
derive_cofactor <- function(ex, marker_cols, target_iqr = 0.5,
                            lo = 1, hi = 5000, max_markers = 40L) {
  mk <- names(marker_cols)
  if (length(mk) > max_markers) mk <- sample(mk, max_markers)
  cands <- vapply(mk, function(m) {
    x <- ex[, marker_cols[[m]]]
    q <- stats::quantile(x, c(0.05, 0.60), na.rm = TRUE)
    bg <- x[x > q[1] & x < q[2]]
    bg <- bg[is.finite(bg)]
    n <- length(bg)
    if (n < 100L) return(NA_real_)

    # ---- exact fast path for IQR(asinh(bg / mid)) ---------------------------
    #
    # THE COST THIS REMOVES. The bisection below runs 40 iterations per marker,
    # and the original body evaluated asinh() over the WHOLE background vector
    # on every one of them. `bg` is everything between the 5th and 60th
    # percentile, so at 300,000 events that is ~165,000 elements: 40 iterations
    # x 40 markers x 165k = 264 million transcendental calls PER SAMPLE, and
    # pooling the cofactor across 8 samples multiplied that by eight.
    #
    # WHY IT CAN BE SKIPPED WITHOUT APPROXIMATING ANYTHING. asinh(./mid) is
    # strictly increasing for mid > 0, so it maps order statistics to order
    # statistics: the k-th smallest of asinh(bg/mid) is asinh(k-th smallest of
    # bg / mid), whatever mid is. R's default (type 7) quantile is a linear
    # interpolation between two adjacent order statistics at positions that
    # depend only on n and p -- NOT on mid. So the four order statistics the IQR
    # needs can be selected ONCE, and each bisection step becomes four asinh
    # calls on scalars instead of 165,000 on a vector.
    #
    # This is an identity, not an approximation: verified bit-identical against
    # IQR(asinh(bg/mid)) across n = 101/5000/123457 and mid = 1/7.3/150/4999
    # (see tests/testthat/test-transform.R).
    #
    # partial = sorts only far enough to fix the four positions needed, which is
    # O(n) rather than the O(n log n) of a full sort.
    hpos <- function(p) { h <- (n - 1) * p + 1; j <- floor(h); c(j, h - j) }
    p25 <- hpos(0.25); p75 <- hpos(0.75)
    idx <- unique(pmin(n, c(p25[1], p25[1] + 1L, p75[1], p75[1] + 1L)))
    o <- sort(bg, partial = idx)[idx]
    names(o) <- as.character(idx)
    at <- function(k) o[[as.character(min(n, k))]]
    b25 <- c(at(p25[1]), at(p25[1] + 1L)); g25 <- p25[2]
    b75 <- c(at(p75[1]), at(p75[1] + 1L)); g75 <- p75[2]

    a <- lo; b <- hi
    for (i in seq_len(40L)) {
      mid <- sqrt(a * b)
      lo_q <- (1 - g25) * asinh(b25[1] / mid) + g25 * asinh(b25[2] / mid)
      hi_q <- (1 - g75) * asinh(b75[1] / mid) + g75 * asinh(b75[2] / mid)
      if (hi_q - lo_q > target_iqr) a <- mid else b <- mid
    }
    sqrt(a * b)
  }, numeric(1))
  cf <- median(cands, na.rm = TRUE)
  if (!is.finite(cf) || cf <= 0) { warning("cofactor derivation failed; using 150"); cf <- 150 }
  # The per-marker candidates are carried along rather than discarded: they are
  # what --per-marker-cofactor uses, and their spread is worth logging even when
  # the single pooled value is what gets applied.
  attr(cf, "per_marker") <- cands
  cf
}

#' Derive the panel cofactor from SEVERAL samples, not from whichever file was
#' read first
#'
#' WHY THIS EXISTS. The original call site derived the cofactor from
#' `reads[[sids[1]]]` alone -- one sample sets the arcsinh transform for the whole
#' panel, and therefore the scale on which every threshold, every marker median
#' and every UMAP distance in the run is computed. If that first file happens to
#' be weakly stained, an unstained control, or simply an outlier, the transform
#' is wrong for everyone else and nothing downstream can detect it. The choice of
#' "first" is alphabetical, so it is not even a considered sample -- it is a
#' filename.
#'
#' WHAT IT DOES INSTEAD: derives the cofactor independently on up to `max_samples`
#' samples and takes the median of those. The median, not the mean, because the
#' failure this guards against is precisely one aberrant sample, and a mean would
#' let it back in.
#'
#' WHY IT ALSO REPORTS THE SPREAD: if per-sample cofactors vary by more than
#' `warn_ratio`, the panel does not have one shared background and a single
#' cofactor is a compromise rather than a description. That is worth saying out
#' loud -- it usually means a gain change or an instrument setting drifted
#' mid-batch, which is the same thing the batch diagnostic looks for from the
#' other end.
#'
#' BEHAVIOUR CHANGE, stated plainly: this alters derived numbers relative to the
#' baseline whenever samples disagree. Pass --cofactor-from-first-sample to
#' restore the previous single-sample behaviour exactly.
#'
#' @param reads named list of read objects (each carrying exprs + marker_cols)
#' @param sids sample ids belonging to this panel
#' @param max_samples Ceiling on the number of samples used. Default `8L`.
#' @param warn_ratio Ratio of per-sample values above which a warning is emitted. Default `2`.
#' @param ... Passed to derive_cofactor().
#' @return single numeric cofactor, with attributes recording what it came from
#' @export
derive_cofactor_pooled <- function(reads, sids, max_samples = 8L,
                                   warn_ratio = 2.0, ...) {
  sids <- sids[!is.na(sids) & sids %in% names(reads)]
  if (!length(sids)) stop("derive_cofactor_pooled: no readable samples")
  # Evenly spaced across the sample order rather than the first N, so a batch
  # acquired in cohort order does not derive its transform from one cohort.
  if (length(sids) > max_samples)
    sids <- sids[unique(round(seq(1, length(sids), length.out = max_samples)))]
  vals <- vapply(sids, function(s) {
    v <- try(derive_cofactor(reads[[s]]$exprs, reads[[s]]$marker_cols, ...),
             silent = TRUE)
    if (inherits(v, "try-error")) NA_real_ else as.numeric(v)
  }, numeric(1))
  vals <- vals[is.finite(vals) & vals > 0]
  if (!length(vals)) {
    warning("pooled cofactor derivation failed on every sample; using 150")
    return(structure(150, n_samples = 0L, source = "fallback"))
  }
  cf <- stats::median(vals)
  rng <- range(vals)
  if (rng[1] > 0 && rng[2] / rng[1] > warn_ratio)
    log_msg("  NOTE per-sample cofactors span ", round(rng[1], 1), "-",
            round(rng[2], 1), " (ratio ", round(rng[2] / rng[1], 1),
            "x) across ", length(vals), " samples. One shared cofactor is a ",
            "compromise here \u2014 check for a gain change mid-batch, and see the ",
            "batch diagnostic.")
  structure(cf, n_samples = length(vals), per_sample = vals,
            source = "pooled_median")
}

#' Apply the spillover matrix if one is present
#' WHY: spectral instruments (e.g. Sony ID7000) write already-unmixed data and
#' carry no $SPILLOVER keyword. Compensating twice, or erroring because the
#' matrix is absent, are both wrong -- detect and report.
#' @param ff_exprs The ff exprs.
#' @param keywords FCS keyword list.
#' @export
maybe_compensate <- function(ff_exprs, keywords) {
  key <- intersect(c("$SPILLOVER", "SPILL", "$SPILL", "spillover"), names(keywords))
  if (!length(key)) {
    log_msg("  no spillover matrix in keywords \u2014 data is already compensated/unmixed")
    return(ff_exprs)
  }
  sp <- keywords[[key[1]]]
  if (!is.matrix(sp) || nrow(sp) == 0) {
    log_msg("  spillover keyword present but unusable \u2014 skipping compensation")
    return(ff_exprs)
  }
  # Index the matrix POSITIONALLY via its column names.
  #
  # WHY not sp[cn, cn]: an FCS spillover matrix is square and its row order is
  # its column order, so writers frequently emit colnames only and leave rownames
  # NULL. Character row indexing then fails with "subscript out of bounds" -- on a
  # file that is otherwise perfectly readable. Matching column names to positions
  # works whether or not rownames were written.
  cn <- intersect(colnames(sp), colnames(ff_exprs))
  if (length(cn) < 2) {
    log_msg("  spillover channels do not match data columns \u2014 skipping compensation")
    return(ff_exprs)
  }
  j <- match(cn, colnames(sp))
  m <- sp[j, j, drop = FALSE]
  if (nrow(m) != ncol(m)) {
    log_msg("  spillover matrix is not square \u2014 skipping compensation")
    return(ff_exprs)
  }
  inv <- try(solve(m), silent = TRUE)
  # A singular or near-singular matrix means the file's own compensation is
  # unusable; inverting it would amplify noise without bound. Report and proceed
  # uncompensated rather than emit numbers nobody can interpret.
  if (inherits(inv, "try-error")) {
    log_msg("  spillover matrix is singular \u2014 skipping compensation")
    return(ff_exprs)
  }
  log_msg("  applying spillover matrix (", length(cn), " channels)")
  ff_exprs[, cn] <- as.matrix(ff_exprs[, cn]) %*% inv
  ff_exprs
}

# =============================================================================
