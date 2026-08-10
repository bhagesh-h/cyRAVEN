# SECTION 10c -- PER-MARKER BATCH DRIFT
# =============================================================================
#
# WHY THIS EXISTS ALONGSIDE batch_mixing_report(). That function answers whether
# the batches mix in the shared embedding, and iLISI is the right statistic for
# it. What it cannot say is WHICH channel moved. A flagged embedding names
# nothing a wet-lab decision can act on; a flagged marker names a reagent lot, a
# detector voltage, or a laser that was serviced between runs.
#
# WHY NOT JUST REUSE stats_threshold_drift() WITH THE BATCH COLUMN. That is worth
# doing and the pipeline now does it, at no cost, because that function was
# always generic in its grouping variable. But a threshold is one number per
# sample, so it only sees drift that moves the CUT. A marker can shift its spread,
# grow a tail, or lose the separation between its modes while the density minimum
# between them stays exactly where it was. The distribution has to be compared
# directly to catch that.
#
# THE STATISTIC. One-dimensional Earth Mover's distance, the L1 distance between
# two empirical quantile functions. It is the natural distance for "how far would
# the cells have to move", it needs no binning choice, and it is in the units of
# the analysis scale so it can be reported as a distance rather than as a score.
# Because those units differ per marker, it is also reported divided by the
# marker's own pooled MAD, which is what makes two markers comparable.
#
# RNG. Subsampling consumes draws, so the stream is saved and restored on the way
# out. See uncertainty.R for why: run_cyraven() seeds once and the UMAP cell
# selection draws from that one stream.

#' One-dimensional Earth Mover's distance between two samples
#'
#' The mean absolute difference between the two empirical quantile functions,
#' evaluated on a common grid. Equivalent to the area between the two cumulative
#' distributions, and to the Wasserstein-1 distance.
#'
#' @param a,b numeric vectors
#' @param n quantile grid resolution
#' @return the distance in the units of `a` and `b`, or NA if either side is
#'   too small to describe a distribution
#' @export
emd_1d <- function(a, b, n = 512L) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 20L || length(b) < 20L) return(NA_real_)
  p <- (seq_len(n) - 0.5) / n
  qa <- stats::quantile(a, p, names = FALSE, type = 7L)
  qb <- stats::quantile(b, p, names = FALSE, type = 7L)
  mean(abs(qa - qb))
}

#' Per-marker distributional drift between acquisition batches
#'
#' For each marker, the largest Earth Mover's distance between any two batches,
#' reported in analysis units and scaled by the marker's pooled MAD so markers can
#' be compared with each other.
#'
#' WHAT A FLAG MEANS. A marker whose distribution differs between batches by an
#' appreciable fraction of its own spread was measured differently in them. That
#' is a statement about the assay, not about the donors, unless batch and study
#' group coincide -- which is the question `batch_mixing_report()` answers with
#' Cramer's V, and which should be read first. Where the two overlap, this table
#' cannot separate a reagent lot from the biology either.
#'
#' @param cells cell-level table carrying `sample_id`, the batch column, and one
#'   column per marker on the analysis scale
#' @param batch_col name of the batch column in `cells`
#' @param markers marker columns to test; defaults to every numeric column that
#'   is not structural
#' @param max_cells cells sampled per batch before the distances are computed;
#'   the quantile function stops moving long before the whole batch is used
#' @param min_cells smallest batch worth comparing
#' @param flag_at value of `emd_over_mad` at or above which a marker is flagged
#' @param seed seed for the local RNG stream, which is restored on exit
#' @return a data.frame ordered by `emd_over_mad`, or NULL
#' @export
marker_batch_drift <- function(cells, batch_col, markers = NULL,
                               max_cells = 20000L, min_cells = 200L,
                               flag_at = 0.5, seed = 42L) {
  if (is.null(cells) || !nrow(cells) || !batch_col %in% names(cells)) return(NULL)
  b <- as.character(cells[[batch_col]])
  ok <- !is.na(b) & nzchar(trimws(b))
  if (!any(ok)) return(NULL)
  cells <- cells[ok, , drop = FALSE]; b <- b[ok]

  if (is.null(markers)) {
    drop_pat <- "^FSC|^SSC|^Time|^Event|Width$|Height$|-H$|-W$|^umap|^cluster"
    num <- names(cells)[vapply(cells, is.numeric, logical(1))]
    markers <- setdiff(num[!grepl(drop_pat, num, ignore.case = TRUE)],
                       c("event_index", "sample_id"))
  }
  markers <- intersect(markers, names(cells))
  if (!length(markers)) return(NULL)

  keep_b <- names(which(table(b) >= min_cells))
  if (length(keep_b) < 2L) return(NULL)

  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  idx <- lapply(keep_b, function(bb) {
    i <- which(b == bb)
    if (length(i) > max_cells) sort(sample(i, max_cells)) else i
  })
  names(idx) <- keep_b

  pairs <- utils::combn(keep_b, 2L, simplify = FALSE)
  rows <- lapply(markers, function(mk) {
    x <- cells[[mk]]
    vals <- lapply(idx, function(i) x[i])
    # Pooled MAD over the batches actually compared, so the scaling matches the
    # comparison rather than the whole table.
    pooled <- unlist(vals, use.names = FALSE)
    mad <- stats::mad(pooled, na.rm = TRUE)
    d <- vapply(pairs, function(p) emd_1d(vals[[p[1]]], vals[[p[2]]]), numeric(1))
    if (!any(is.finite(d))) return(NULL)
    w <- which.max(replace(d, !is.finite(d), -Inf))
    med <- vapply(vals, function(v) stats::median(v, na.rm = TRUE), numeric(1))
    data.frame(
      marker = mk,
      n_batches = length(keep_b),
      n_cells_compared = length(pooled),
      emd_max = round(d[w], 4),
      worst_pair = paste(pairs[[w]], collapse = " vs "),
      pooled_mad = round(mad, 4),
      emd_over_mad = if (is.finite(mad) && mad > 0) round(d[w] / mad, 3) else NA_real_,
      medians_by_batch = paste(sprintf("%s=%.3f", names(med), med), collapse = "; "),
      stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$verdict <- ifelse(!is.finite(out$emd_over_mad), "not comparable",
                 ifelse(out$emd_over_mad >= flag_at, "differs between batches",
                        "consistent"))
  out[order(-replace(out$emd_over_mad, !is.finite(out$emd_over_mad), -Inf)), ,
      drop = FALSE]
}
