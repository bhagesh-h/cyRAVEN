# Explore mode: unsupervised discovery alongside the declared analysis.
#
# WHAT IT IS
#
# The rest of cyRAVEN starts from a specification: you declare populations, it
# scores them, and six diagnostics try to break the declaration. That design
# cannot find a population nobody declared, and it cannot look outside the parent
# gate, because both are consequences of the declaration.
#
# Explore mode is the complement. It takes every event in the file, every
# eligible channel, embeds and clusters without reference to any specification,
# and reports what it found. It is what cyCONDOR, FlowSOM and the rest of the
# unsupervised family do.
#
# WHAT IS DIFFERENT ABOUT DOING IT HERE
#
# An unsupervised tool hands you cluster 17 and a heatmap, and you squint at the
# colour scale to decide what it is. That inference is eyeballed, it is done on
# pooled data, and it has no per-sample calibration behind it.
#
# When the declared pipeline has already run, this module has things no
# standalone clusterer has:
#
#   - the transform, with a cofactor ESTIMATED from this panel rather than
#     assumed;
#   - a threshold for every marker IN EVERY SAMPLE, derived at that sample's own
#     density minimum;
#   - staining QC verdicts;
#   - the batch/group confounding verdict.
#
# So a cluster is named, not guessed at: for each marker, the fraction of the
# cluster's cells above THEIR OWN SAMPLE'S threshold. "CD19+ HLA-DR+ CD3- CD14-"
# with a number behind every call. And its abundance is tested per donor, with
# Benjamini-Hochberg and an effect size, carrying the confounding verdict -- so
# explore mode cannot quietly reintroduce the failure modes the declared pipeline
# exists to prevent.
#
# Standalone (`--explore` with no config) it degrades to the unsupervised default:
# a pooled cofactor, and 2-means on per-cluster medians for the QC gate, which is
# what an unsupervised tool would do anyway. Every such fallback is recorded in
# `explore_provenance.csv` so a reader can tell which mode produced a number.
#
# ISOLATION
#
# Everything lands in <outdir>/explore/. Not one existing output path is written
# to, and with --explore absent no function in this file is called. The claim is
# tested rather than asserted: tests/testthat/test-explore.R runs the pipeline
# with and without the flag and compares the file lists.
#
# FEATURES
#
# Every eligible channel, including scatter and viability. That is deliberate and
# it is the unsupervised default: it is what lets the QC gate find debris (CD45
# low), dead cells (viability high) and granulocytes (side scatter high), none of
# which any lineage antibody identifies. The cost is that some clusters split on
# scatter rather than on lineage. `--explore-markers` overrides the list and
# `--explore-exclude` drops channels from it.

#' Eligible explore features for a panel
#'
#' Every fluorescence marker the panel resolves, plus scatter, minus anything
#' excluded. Unlike [select_umap_features()] there is no lineage preference: the
#' point of explore mode is to assume nothing about which channels matter.
#'
#' @param panel_markers Character vector of marker names in the panel.
#' @param scatter_names Character vector of scatter channel names.
#' @param prefer Optional explicit list; when given it replaces the default
#'   entirely and only its intersection with what is available is used.
#' @param exclude Channels to drop.
#' @return Character vector of feature names.
#' @keywords internal
explore_features <- function(panel_markers, scatter_names = character(0),
                             prefer = NULL, exclude = NULL) {
  feats <- if (!is.null(prefer) && length(prefer)) {
    intersect(prefer, c(panel_markers, scatter_names))
  } else {
    c(panel_markers, scatter_names)
  }
  feats <- setdiff(feats, exclude %||% character(0))
  unique(feats[!is.na(feats) & nzchar(feats)])
}

#' Transformed matrix over every eligible channel for one sample
#'
#' Mirrors how the declared path builds `tmat`, with one difference: it is not
#' restricted to the channels a specification happens to reference. Fluorescence
#' goes through the panel transform; scatter is log10, because it is strictly
#' positive and spans decades and the arcsinh cofactor derived for fluorescence
#' is meaningless on it.
#'
#' @param rd One element of `reads`.
#' @param feats Feature names wanted.
#' @param tr Transform object from [make_transform()].
#' @return Numeric matrix, cells x features.
#' @keywords internal
explore_matrix <- function(rd, feats, tr) {
  fl <- intersect(feats, names(rd$marker_cols))
  sc <- intersect(feats, names(rd$scatter_cols))
  n <- nrow(rd$exprs)
  out <- NULL
  if (length(fl)) {
    out <- vapply(fl, function(m) tr$fn(rd$exprs[, rd$marker_cols[[m]]], m),
                  numeric(n))
    if (!is.matrix(out)) out <- matrix(out, nrow = n, dimnames = list(NULL, fl))
  }
  if (length(sc)) {
    sm <- vapply(sc, function(k) log10(pmax(rd$exprs[, rd$scatter_cols[[k]]], 1)),
                 numeric(n))
    if (!is.matrix(sm)) sm <- matrix(sm, nrow = n, dimnames = list(NULL, sc))
    out <- if (is.null(out)) sm else cbind(out, sm)
  }
  out
}

#' Split per-cluster medians into a low and a high mode
#'
#' 2-means on one column of per-cluster medians. Used only when no per-sample
#' threshold is available, i.e. standalone explore. Returns the indices of the
#' HIGH group and the boundary between them, or NULL when the split is not
#' meaningful (fewer than 3 clusters, or no separation).
#'
#' @param v Numeric vector, one median per cluster.
#' @param seed Integer seed.
#' @return list(high = integer, cut = numeric) or NULL.
#' @keywords internal
two_mode_split <- function(v, seed = 42L) {
  v <- as.numeric(v)
  ok <- is.finite(v)
  if (sum(ok) < 3L) return(NULL)
  if (diff(range(v[ok])) <= .Machine$double.eps) return(NULL)
  withr::with_seed(seed, {
    km <- try(stats::kmeans(matrix(v[ok], ncol = 1), centers = 2L, nstart = 10L),
              silent = TRUE)
  })
  if (inherits(km, "try-error")) return(NULL)
  hi_cluster <- which.max(km$centers[, 1])
  hi_idx <- which(ok)[km$cluster == hi_cluster]
  cut <- mean(c(max(v[setdiff(which(ok), hi_idx)], na.rm = TRUE),
                min(v[hi_idx], na.rm = TRUE)))
  list(high = hi_idx, cut = cut)
}

#' Fraction of each cluster's cells positive for each marker
#'
#' The heart of what this module adds over a standalone clusterer. Positivity is
#' judged against each cell's OWN SAMPLE'S threshold, so a brightly stained
#' sample and a dim one are not compared against a common cut. Falls back to the
#' pooled per-cluster median when a sample has no threshold for that marker.
#'
#' @param X Transformed matrix, cells x features.
#' @param cluster Integer cluster per cell.
#' @param sample_id Character sample id per cell.
#' @param thr_by_sample Named list: sample id -> named numeric vector of
#'   thresholds on the transformed scale. May be empty.
#' @return Numeric matrix, clusters x features, values in [0, 1], with attribute
#'   "source" naming how each feature was called.
#' @keywords internal
explore_positivity <- function(X, cluster, sample_id, thr_by_sample = list()) {
  feats <- colnames(X)
  ks <- sort(unique(cluster))
  # Called twice with differently typed cluster ids: raw integers for the coarse
  # QC pass, and already-prefixed "k1"/"k2" strings for the final clustering.
  # Prefixing unconditionally produced "kk1" rownames, and every later
  # pos[ks, ] lookup then failed with "subscript out of bounds".
  rn <- if (is.numeric(cluster)) paste0("k", ks) else as.character(ks)
  pos <- matrix(NA_real_, nrow = length(ks), ncol = length(feats),
                dimnames = list(rn, feats))
  src <- setNames(rep("pooled_median", length(feats)), feats)

  for (f in feats) {
    thr_vec <- rep(NA_real_, length(cluster))
    if (length(thr_by_sample)) {
      tv <- vapply(thr_by_sample, function(t) {
        # `[[` on a NAMED NUMERIC VECTOR raises "subscript out of bounds" for a
        # name that is absent, where the same call on a list returns NULL. Not
        # every feature has a threshold -- scatter channels never do -- so the
        # membership test has to come first.
        if (!f %in% names(t)) return(NA_real_)
        v <- t[[f]]
        if (is.null(v) || !is.finite(v)) NA_real_ else as.numeric(v)
      }, numeric(1))
      have <- names(tv)[is.finite(tv)]
      if (length(have)) {
        hit <- sample_id %in% have
        thr_vec[hit] <- tv[sample_id[hit]]
        src[[f]] <- if (all(hit)) "per_sample_threshold" else "per_sample_threshold_partial"
      }
    }
    # Anything with no per-sample threshold falls back to the pooled median of
    # that feature, which is what an unsupervised tool has available.
    if (anyNA(thr_vec)) {
      med <- stats::median(X[, f], na.rm = TRUE)
      thr_vec[is.na(thr_vec)] <- med
    }
    above <- X[, f] > thr_vec
    for (i in seq_along(ks)) {
      sel <- cluster == ks[i]
      pos[i, f] <- if (any(sel)) mean(above[sel], na.rm = TRUE) else NA_real_
    }
  }
  attr(pos, "source") <- src
  pos
}

#' Phenotype string for each cluster
#'
#' `CD19+ HLA-DR+ CD3- CD14-`, from the positivity fractions. A marker is called
#' `+` above `hi`, `-` below `lo`, and omitted in between -- an omitted marker is
#' an honest "this cluster is mixed for it", not a silent negative.
#'
#' @param pos Matrix from [explore_positivity()].
#' @param hi,lo Fractions bounding a positive and a negative call.
#' @param max_markers Most markers to name, highest |evidence| first.
#' @return Character vector, one phenotype per cluster.
#' @keywords internal
explore_phenotype <- function(pos, hi = 0.65, lo = 0.20, max_markers = 6L) {
  apply(pos, 1, function(r) {
    calls <- character(0)
    ev <- numeric(0)
    for (f in names(r)) {
      v <- r[[f]]
      if (!is.finite(v)) next
      if (v >= hi) { calls <- c(calls, paste0(f, "+")); ev <- c(ev, v) }
      else if (v <= lo) { calls <- c(calls, paste0(f, "-")); ev <- c(ev, 1 - v) }
    }
    if (!length(calls)) return("unresolved")
    o <- order(ev, decreasing = TRUE)
    paste(calls[o][seq_len(min(max_markers, length(calls)))], collapse = " ")
  })
}

#' Cluster-level QC gate
#'
#' cyCONDOR's idea: cluster coarsely first, then judge WHOLE CLUSTERS on their
#' marker profile, so no per-event cutoff is invented. Three calls, all derived
#' rather than hard-coded:
#'
#'   debris     leukocyte marker low
#'   dead       viability marker high
#'   saturated  top-percentile in most channels at once -- no real cell type is,
#'              so these are aggregates and doublets
#'
#' The improvement over doing it blind: when per-sample thresholds exist, the
#' leukocyte and viability calls use them, so "CD45 low" means "below this
#' sample's own CD45 cut" rather than "in the lower of two pooled modes".
#'
#' @param pos Positivity matrix from [explore_positivity()].
#' @param med Per-cluster medians, clusters x features.
#' @param sizes Cells per cluster.
#' @param leukocyte,viability Channel names, or NULL.
#' @param have_thresholds Whether `pos` was built from per-sample thresholds.
#' @param drop_dead Whether to drop the dead call.
#' @param sat_quantile,sat_min_fraction Saturation rule.
#' @param seed Integer seed.
#' @return data.frame with one row per cluster and a `call` column.
#' @keywords internal
explore_qc_gate <- function(pos, med, sizes, leukocyte = NULL, viability = NULL,
                            have_thresholds = FALSE, drop_dead = TRUE,
                            sat_quantile = 0.99, sat_min_fraction = 0.7,
                            seed = 42L) {
  ks <- rownames(pos)
  call <- rep("keep", length(ks))
  basis <- rep(NA_character_, length(ks))

  mark_low <- function(chan, label) {
    if (is.null(chan) || !chan %in% colnames(pos)) return(invisible(NULL))
    if (have_thresholds) {
      lowi <- which(pos[, chan] < 0.5)
      why <- "below the sample's own threshold in most of its cells"
    } else {
      sp <- two_mode_split(med[, chan], seed = seed)
      if (is.null(sp)) return(invisible(NULL))
      lowi <- setdiff(seq_along(ks), sp$high)
      why <- "in the low mode of a 2-means split on per-cluster medians"
    }
    sel <- lowi[call[lowi] == "keep"]
    call[sel] <<- label
    basis[sel] <<- why
    invisible(NULL)
  }
  mark_high <- function(chan, label) {
    if (is.null(chan) || !chan %in% colnames(pos)) return(invisible(NULL))
    if (have_thresholds) {
      hii <- which(pos[, chan] >= 0.5)
      why <- "above the sample's own threshold in most of its cells"
    } else {
      sp <- two_mode_split(med[, chan], seed = seed)
      if (is.null(sp)) return(invisible(NULL))
      hii <- sp$high
      why <- "in the high mode of a 2-means split on per-cluster medians"
    }
    sel <- hii[call[hii] == "keep"]
    call[sel] <<- label
    basis[sel] <<- why
    invisible(NULL)
  }

  # Saturation first: an aggregate is bright everywhere, including the leukocyte
  # marker, so calling it debris or dead would mislabel it.
  if (ncol(med) >= 3L) {
    tops <- vapply(colnames(med), function(f)
      med[, f] >= stats::quantile(med[, f], sat_quantile, na.rm = TRUE),
      logical(nrow(med)))
    if (!is.matrix(tops)) tops <- matrix(tops, nrow = nrow(med))
    frac_top <- rowMeans(tops, na.rm = TRUE)
    sat <- which(frac_top >= sat_min_fraction)
    call[sat] <- "saturated"
    basis[sat] <- sprintf("top %.0f%% in >= %.0f%% of channels at once",
                          100 * (1 - sat_quantile), 100 * sat_min_fraction)
  }

  mark_low(leukocyte, "debris")
  if (isTRUE(drop_dead)) mark_high(viability, "dead")

  data.frame(cluster = ks, cells = as.integer(sizes),
             pct = round(100 * sizes / sum(sizes), 3),
             call = call, basis = basis, stringsAsFactors = FALSE)
}
