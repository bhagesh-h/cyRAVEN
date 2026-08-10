# =============================================================================
# Per-marker batch drift
#
# The claim under test: a marker that was measured differently between batches is
# named, and one that was not is left alone. Both halves matter -- a diagnostic
# that flags everything is as useless as one that flags nothing.
# =============================================================================

#' Cells in two batches. `shift` displaces CD4 in batch B only; `spread`
#' multiplies its scale in batch B without moving its centre.
drift_cells <- function(shift = 0, spread = 1, n = 3000L, seed = 5) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  data.frame(
    sample_id = rep(c("S1", "S2"), each = n),
    batch = rep(c("A", "B"), each = n),
    CD4 = c(stats::rnorm(n, 2, 1), stats::rnorm(n, 2 + shift, 1 * spread)),
    CD8 = stats::rnorm(2 * n, 1, 1),
    stringsAsFactors = FALSE)
}

test_that("the distance is zero for identical samples and grows with separation", {
  x <- stats::rnorm(4000)
  expect_equal(emd_1d(x, x), 0)
  # A pure location shift of d moves every quantile by d, so the distance is d.
  expect_equal(emd_1d(x, x + 2), 2, tolerance = 1e-8)
  expect_gt(emd_1d(x, x + 3), emd_1d(x, x + 1))
})

test_that("too few events returns NA rather than a distance", {
  expect_true(is.na(emd_1d(stats::rnorm(5), stats::rnorm(4000))))
  expect_true(is.na(emd_1d(numeric(0), numeric(0))))
})

test_that("a marker that shifted between batches is flagged and one that did not is not", {
  d <- marker_batch_drift(drift_cells(shift = 2), "batch", min_cells = 100L)
  expect_s3_class(d, "data.frame")
  rownames(d) <- d$marker

  expect_identical(d["CD4", "verdict"], "differs between batches")
  expect_identical(d["CD8", "verdict"], "consistent")
  # The shifted marker sorts first, because the table is ordered by the scaled
  # distance and a reader should meet the worst channel before the others.
  expect_identical(d$marker[1], "CD4")
  expect_gt(d["CD4", "emd_over_mad"], d["CD8", "emd_over_mad"])
  expect_identical(d["CD4", "worst_pair"], "A vs B")
})

test_that("drift in spread alone is caught, which a threshold test cannot see", {
  # THE REASON THIS EXISTS ALONGSIDE stats_threshold_drift(). Batch B has the same
  # centre and a much wider distribution. A density minimum between two modes, and
  # therefore a threshold, can sit in exactly the same place in both.
  d <- marker_batch_drift(drift_cells(shift = 0, spread = 3), "batch",
                          min_cells = 100L)
  rownames(d) <- d$marker
  expect_identical(d["CD4", "verdict"], "differs between batches")
  expect_identical(d["CD8", "verdict"], "consistent")
})

test_that("nothing is flagged when the batches agree", {
  d <- marker_batch_drift(drift_cells(shift = 0), "batch", min_cells = 100L)
  expect_true(all(d$verdict == "consistent"))
})

test_that("structural channels are never treated as markers", {
  cl <- drift_cells(shift = 2)
  cl$`FSC-A` <- stats::runif(nrow(cl))
  cl$umap_1 <- stats::runif(nrow(cl))
  cl$event_index <- seq_len(nrow(cl))
  d <- marker_batch_drift(cl, "batch", min_cells = 100L)
  expect_setequal(d$marker, c("CD4", "CD8"))
})

test_that("a single batch is not a comparison", {
  cl <- drift_cells()
  cl$batch <- "A"
  expect_null(marker_batch_drift(cl, "batch", min_cells = 100L))
  # Nor is a batch too small to describe a distribution.
  expect_null(marker_batch_drift(drift_cells(), "batch", min_cells = 1e6))
  expect_null(marker_batch_drift(drift_cells(), "not_a_column"))
  expect_null(marker_batch_drift(NULL, "batch"))
})

test_that("the drift scan leaves the RNG stream exactly where it was", {
  # Same hazard as the uncertainty code: this subsamples each batch, and
  # run_cyraven() seeds once for the whole run.
  cl <- drift_cells(shift = 2, n = 5000L)
  set.seed(123)
  before <- .Random.seed
  expected <- stats::runif(5)
  set.seed(123)
  invisible(marker_batch_drift(cl, "batch", max_cells = 1000L, min_cells = 100L))
  expect_identical(.Random.seed, before)
  expect_identical(stats::runif(5), expected)
})

test_that("the threshold test accepts a batch map as readily as a group map", {
  # This is the claim that made the threshold half of the feature free: the
  # function was always generic in its grouping vector.
  thr <- data.frame(
    sample_id = sprintf("S%02d", 1:8),
    marker = "CD4",
    threshold = c(1.0, 1.1, 0.9, 1.05, 3.0, 3.1, 2.9, 3.05),
    source = "valley", stringsAsFactors = FALSE)
  thr <- rbind(thr, transform(thr, marker = "CD8", threshold = rep(2, 8)))
  bmap <- stats::setNames(rep(c("run1", "run2"), each = 4), sprintf("S%02d", 1:8))

  out <- stats_threshold_drift(thr, bmap)
  rownames(out) <- out$marker
  expect_lt(out["CD4", "p_value"], 0.05)
  expect_gt(out["CD4", "gap_over_sd"], out["CD8", "gap_over_sd"])
})
