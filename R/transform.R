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
#'
#' WHAT THE CRITERION AMOUNTS TO. For a marker whose events are mostly negative,
#' N(0, sigma), the background band is the 5th to 60th percentile and its IQR
#' is about 0.79 sigma, so an IQR of 0.5 after asinh(x / c) puts c near 1.6
#' sigma: the linear region of the transform just wider than the negative
#' population. That is the condition that keeps a symmetric negative population
#' one mode, which it stays exactly while c >= sigma, and the same idea as the
#' logicle width rule (Parks et al. 2006). It is a heuristic, not a published
#' criterion.
#'
#' PER MARKER, FROM THE NEGATIVE POPULATION. Detectors differ many-fold in how
#' wide their negatives are: on one spectral panel the 5th percentile of the
#' values below zero ran from 147 (CCR5) to 2,330 (CCR9). One cofactor then
#' splits the widest negatives into two peaks and compresses the narrowest, and
#' the published recommendation is a cofactor per channel (den Braanker et al.
#' 2021; flowVS, Azad et al. 2016). The bulk criterion above cannot be applied
#' per marker as it stands, because for a marker most cells carry (CD45, CD3)
#' the 5th-60th percentile band is the POSITIVE population, which gave CD45 a
#' cofactor of 3,256 on the same panel. So each marker's own negative population
#' is reconstructed from its values below zero, mirrored about zero, which
#' positive cells never reach; the same bisection is run on that. Returned as
#' attribute `per_marker_negative`, NA where fewer than 100 values sit below
#' zero (uncompensated or linearised log data, for instance).
#'
#' @param ex Numeric expression matrix, events x channels.
#' @param marker_cols Named integer vector mapping marker symbol to column index.
#' @param target_iqr Interquartile range the transformed background is bisected towards. Default `0.5`.
#' @param lo Lower bound of the search interval. Default `1`.
#' @param hi Upper bound of the search interval. Default `5000`.
#' @param max_markers Ceiling on the number of markers used. Default `40L`.
#' @param rows Logical or integer index of the events to estimate from, or
#'   NULL for every event. The pipeline passes the scatter- and singlet-gated
#'   cells, so debris does not set the background.
#' @export
derive_cofactor <- function(ex, marker_cols, target_iqr = 0.5,
                            lo = 1, hi = 5000, max_markers = 40L, rows = NULL) {
  mk <- names(marker_cols)
  if (length(mk) > max_markers) mk <- sample(mk, max_markers)
  col <- function(m) {
    x <- ex[, marker_cols[[m]]]
    if (is.null(rows)) x else x[rows]
  }
  cands <- vapply(mk, function(m) bisect_cofactor(col(m), target_iqr, lo, hi),
                  numeric(1))
  neg <- vapply(names(marker_cols), function(m) {
    x <- col(m)
    n <- x[is.finite(x) & x < 0]
    if (length(n) < 100L) return(NA_real_)
    bisect_cofactor(c(n, -n), target_iqr, lo, hi)
  }, numeric(1))
  cf <- median(cands, na.rm = TRUE)
  if (!is.finite(cf) || cf <= 0) { warning("cofactor derivation failed; using 150"); cf <- 150 }
  # The bulk per-marker candidates are carried along rather than discarded:
  # their spread is worth logging even when the single pooled value is what
  # gets applied. The negative-population ones are what --per-marker-cofactor
  # applies.
  attr(cf, "per_marker") <- cands
  attr(cf, "per_marker_negative") <- neg
  cf
}

#' Bisect one marker's cofactor to the background-IQR target
#'
#' The search [derive_cofactor()] runs per marker, on the 5th to 60th
#' percentile of `x`.
#' @param x numeric vector, linear units.
#' @param target_iqr,lo,hi see [derive_cofactor()].
#' @return the cofactor, or NA when fewer than 100 values fall in the band.
#' @keywords internal
bisect_cofactor <- function(x, target_iqr = 0.5, lo = 1, hi = 5000) {
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
}

#' The events the transform is estimated from: scatter-gated single cells
#'
#' Debris is most of the events in many files, two thirds on one cohort, and it
#' scatters little and fluoresces less, so a background estimated over every
#' event is the debris's and is narrower than the cells'. The scatter and
#' singlet gates need no fluorescence transform, so they can run first and the
#' cofactor can be estimated on the cells it will be applied to.
#'
#' @param rd A list from [read_fcs_resolved()].
#' @return logical vector over the events, or NULL when the gates cannot run.
#' @keywords internal
cofactor_rows <- function(rd) {
  r <- tryCatch(suppressWarnings(suppressMessages({
    sg <- derive_scatter_gate(rd$exprs, rd$scatter_cols)
    derive_singlet_band(rd$exprs, rd$scatter_cols, sg$mask)$mask
  })), error = function(e) NULL)
  if (is.null(r) || sum(r) < 500L) NULL else r
}

#' The samples a panel's transform is estimated from
#'
#' The biological samples, never the declared controls. A control is in the
#' panel to be compared against, not to set the scale everything is measured
#' on: an unstained tube drags the cofactor down, and a stained FMO or
#' full-stain tube is one donor's distribution counted again. Including them
#' also changes WHICH samples the pooled estimate draws, since it takes an
#' evenly spaced subset of the panel. On the 30-sample sepsis cohort adding two
#' control tubes moved the cofactor from 261 to 326 and every cut in every
#' sample with it, in markers the controls had nothing to do with. When a panel
#' holds nothing but controls they are all it has, and are used.
#' @param sids sample ids in the panel.
#' @param smap the sample map, or NULL.
#' @return the sample ids to estimate from.
#' @keywords internal
transform_samples <- function(sids, smap = NULL) {
  keep <- setdiff(sids, declared_controls(smap))
  if (length(keep)) keep else sids
}

#' Sample ids the sample sheet declares as controls
#'
#' Used wherever a cohort-wide quantity is pooled over samples -- the
#' transform, and the peer cuts and peer comparison of thresholds -- so that a
#' control tube is gated like a sample but never counted as one of the cohort.
#' @param smap the sample map, or NULL.
#' @return character vector of sample ids, empty when none are declared.
#' @keywords internal
declared_controls <- function(smap = NULL) {
  if (is.null(smap) || !all(c("sample_id", "is_control") %in% names(smap)))
    return(character(0))
  as.character(smap$sample_id[vapply(smap$is_control, isTRUE, logical(1))])
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
#' @param cells When `TRUE` (default) each sample's cofactor is estimated on its
#'   scatter-gated single cells (see [cofactor_rows()]); `FALSE` uses every
#'   event, as versions before 1.0.0's audit did.
#' @param ... Passed to derive_cofactor().
#' @return single numeric cofactor, with attributes recording what it came from:
#'   `per_marker_negative`, the median over samples of each marker's
#'   negative-population cofactor (see [derive_cofactor()]), among them.
#' @export
derive_cofactor_pooled <- function(reads, sids, max_samples = 8L,
                                   warn_ratio = 2.0, cells = TRUE, ...) {
  sids <- sids[!is.na(sids) & sids %in% names(reads)]
  if (!length(sids)) stop("derive_cofactor_pooled: no readable samples")
  # Evenly spaced across the sample order rather than the first N, so a batch
  # acquired in cohort order does not derive its transform from one cohort.
  if (length(sids) > max_samples)
    sids <- sids[unique(round(seq(1, length(sids), length.out = max_samples)))]
  per <- list()
  vals <- vapply(sids, function(s) {
    rows <- if (isTRUE(cells)) cofactor_rows(reads[[s]]) else NULL
    v <- try(derive_cofactor(reads[[s]]$exprs, reads[[s]]$marker_cols,
                             rows = rows, ...), silent = TRUE)
    if (inherits(v, "try-error")) return(NA_real_)
    per[[s]] <<- attr(v, "per_marker_negative")
    as.numeric(v)
  }, numeric(1))
  mk <- unique(unlist(lapply(per, names)))
  pm <- vapply(mk, function(m) {
    v <- unlist(lapply(per, function(p) if (m %in% names(p)) p[[m]] else NA_real_))
    v <- v[is.finite(v) & v > 0]
    if (length(v)) stats::median(v) else NA_real_
  }, numeric(1))
  vals <- vals[is.finite(vals) & vals > 0]
  if (!length(vals)) {
    warning("pooled cofactor derivation failed on every sample; using 150")
    return(structure(150, n_samples = 0L, source = "fallback", per_marker_negative = pm))
  }
  cf <- stats::median(vals)
  rng <- range(vals)
  if (rng[1] > 0 && rng[2] / rng[1] > warn_ratio)
    log_msg("  NOTE per-sample cofactors span ", round(rng[1], 1), "-",
            round(rng[2], 1), " (ratio ", round(rng[2] / rng[1], 1),
            "x) across ", length(vals), " samples. One shared cofactor is a ",
            "compromise here, check for a gain change mid-batch, and see the ",
            "batch diagnostic.")
  structure(cf, n_samples = length(vals), per_sample = vals,
            source = "pooled_median", per_marker_negative = pm,
            events = if (isTRUE(cells)) "scatter-gated single cells" else "all events")
}

#' Apply the spillover matrix if one is present
#'
#' WHY: spectral instruments (e.g. Sony ID7000) write already-unmixed data and
#' carry no $SPILLOVER keyword. Compensating twice, or erroring because the
#' matrix is absent, are both wrong -- detect and report.
#'
#' WHAT AN ABSENT KEYWORD MEANS. Nothing, on its own. `$SPILLOVER` is the only
#' standard carrier of a compensation matrix (FCS 3.1, Spidlen et al. 2010), so
#' a file without it is either unmixed or compensated already, or conventional
#' data nobody compensated. The keyword cannot tell those apart, and the result
#' says the status is unknown rather than guessing the reassuring answer. FCS 3.2
#' files can say it per parameter in `$PnTYPE` ("Unmixed_Fluorescence"), and
#' that is read where it is written.
#'
#' A MATRIX FROM ELSEWHERE. `external`, from `--spillover`, replaces the file's
#' own matrix: a matrix recomputed from single-stain controls, in AutoSpill or in
#' the acquisition software, is applied to every file instead of whatever the
#' files carry.
#'
#' @param ff_exprs The ff exprs.
#' @param keywords FCS keyword list.
#' @param external A spillover matrix from [read_spillover_file()], or NULL.
#' @return `ff_exprs`, compensated where a matrix applied, with attribute
#'   `compensation` stating what was done.
#' @export
maybe_compensate <- function(ff_exprs, keywords, external = NULL) {
  done <- function(x, status) structure(x, compensation = status)
  key <- intersect(c("$SPILLOVER", "SPILL", "$SPILL", "spillover"), names(keywords))
  if (!is.null(external)) {
    sp <- external
    src <- "supplied matrix (--spillover)"
  } else if (!length(key)) {
    types <- unlist(keywords[grepl("^\\$P[0-9]+TYPE$", names(keywords))])
    if (any(grepl("unmixed", types, ignore.case = TRUE))) {
      log_msg("  no spillover matrix; $PnTYPE declares the data unmixed")
      return(done(ff_exprs, "unmixed, declared by $PnTYPE"))
    }
    log_msg("  no spillover matrix in the keywords: compensation status UNKNOWN. ",
            "Treated as already compensated or unmixed, which is what spectral ",
            "instruments write; conventional data never compensated would need ",
            "--spillover")
    return(done(ff_exprs, "unknown: no $SPILLOVER, used as stored"))
  } else {
    sp <- keywords[[key[1]]]
    src <- "file $SPILLOVER"
  }
  if (!is.matrix(sp) || nrow(sp) == 0) {
    log_msg("  spillover keyword present but unusable, skipping compensation")
    return(done(ff_exprs, "not applied: $SPILLOVER unusable"))
  }
  # A supplied matrix names detectors as they are written ("FITC-A"); the data
  # carry them as alter.names makes them ("FITC.A"). Accept either spelling.
  if (!is.null(external) && !any(colnames(sp) %in% colnames(ff_exprs)))
    colnames(sp) <- make.names(colnames(sp))
  # Index the matrix POSITIONALLY via its column names.
  #
  # WHY not sp[cn, cn]: an FCS spillover matrix is square and its row order is
  # its column order, so writers frequently emit colnames only and leave rownames
  # NULL. Character row indexing then fails with "subscript out of bounds" -- on a
  # file that is otherwise perfectly readable. Matching column names to positions
  # works whether or not rownames were written.
  cn <- intersect(colnames(sp), colnames(ff_exprs))
  if (length(cn) < 2) {
    log_msg("  spillover channels do not match data columns, skipping compensation")
    return(done(ff_exprs, paste0("not applied: ", src, " names no data column")))
  }
  j <- match(cn, colnames(sp))
  m <- sp[j, j, drop = FALSE]
  if (nrow(m) != ncol(m)) {
    log_msg("  spillover matrix is not square, skipping compensation")
    return(done(ff_exprs, paste0("not applied: ", src, " not square")))
  }
  inv <- try(solve(m), silent = TRUE)
  # A singular or near-singular matrix means the file's own compensation is
  # unusable; inverting it would amplify noise without bound. Report and proceed
  # uncompensated rather than emit numbers nobody can interpret.
  if (inherits(inv, "try-error")) {
    log_msg("  spillover matrix is singular, skipping compensation")
    return(done(ff_exprs, paste0("not applied: ", src, " singular")))
  }
  log_msg("  applying ", src, " (", length(cn), " channels)")
  ff_exprs[, cn] <- as.matrix(ff_exprs[, cn]) %*% inv
  done(ff_exprs, paste0("compensated, ", src))
}

#' Read a spillover matrix from a CSV file
#'
#' The layout FlowJo and AutoSpill export and the `$SPILLOVER` keyword encodes:
#' a square matrix whose columns are the detectors, row i giving the fraction
#' of detector i's fluorochrome seen in each detector, with ones on the
#' diagonal. A first column holding the row names is accepted and dropped.
#'
#' @param path CSV path, or NULL.
#' @return numeric matrix with detector column names, or NULL when `path` is
#'   NULL or empty.
#' @export
read_spillover_file <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) == nrow(d) + 1L && !is.numeric(d[[1]])) d <- d[, -1, drop = FALSE]
  m <- suppressWarnings(as.matrix(sapply(d, as.numeric)))
  if (is.null(dim(m)) || nrow(m) != ncol(m) || any(!is.finite(m)))
    stop("--spillover: ", path, " is not a square numeric matrix with one ",
         "column per detector", call. = FALSE)
  if (any(abs(diag(m) - 1) > 1e-6))
    stop("--spillover: the diagonal of ", path, " is not all ones, so it is not ",
         "a spillover matrix. A compensation matrix (its inverse) is not accepted",
         call. = FALSE)
  colnames(m) <- names(d)
  m
}

# =============================================================================
