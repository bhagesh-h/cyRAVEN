# =============================================================================
# The arcsinh cofactor, and the exactness of its fast path
# =============================================================================
#
# derive_cofactor() bisects on IQR(asinh(bg / mid)). The implementation does NOT
# evaluate that expression: it exploits the fact that asinh is strictly
# increasing to select the four order statistics the interquartile range needs
# once, then interpolates them per iteration. These tests exist to hold that
# identity in place, because the optimisation is invisible at the call site and a
# well-meaning future edit could replace it with the "obvious" version and never
# notice the run got 100x slower — or, worse, replace it with an approximation
# and change every threshold in every result set.

test_that("the order-statistic identity is exact, not approximate", {
  # asinh(./m) maps order statistics to order statistics for any m > 0, and R's
  # type-7 quantile interpolates between positions that depend only on n and p.
  # So IQR(asinh(x/m)) is computable from x's order statistics alone.
  iqr_from_order_stats <- function(x, m) {
    n <- length(x); xs <- sort(x)
    at <- function(p) {
      h <- (n - 1) * p + 1; j <- floor(h); g <- h - j
      (1 - g) * asinh(xs[j] / m) + g * asinh(xs[min(n, j + 1)] / m)
    }
    at(0.75) - at(0.25)
  }
  withr::local_seed(1)
  for (n in c(101L, 5000L, 20011L)) {
    x <- abs(stats::rnorm(n, 300, 400))
    for (m in c(1, 7.3, 150, 4999)) {
      expect_identical(iqr_from_order_stats(x, m), stats::IQR(asinh(x / m)),
                       info = sprintf("n = %d, cofactor = %g", n, m))
    }
  }
})

test_that("derive_cofactor recovers a planted scale and is deterministic", {
  withr::local_seed(42)
  ex <- matrix(abs(stats::rnorm(20000 * 4, 0, 250)), ncol = 4,
               dimnames = list(NULL, c("A", "B", "C", "D")))
  mc <- stats::setNames(seq_len(4), colnames(ex))

  cf1 <- derive_cofactor(ex, mc)
  cf2 <- derive_cofactor(ex, mc)

  expect_true(is.finite(as.numeric(cf1)))
  expect_gt(as.numeric(cf1), 0)
  expect_identical(as.numeric(cf1), as.numeric(cf2))
  # The per-marker candidates are carried for --per-marker-cofactor.
  expect_length(attr(cf1, "per_marker"), 4L)
})

test_that("derive_cofactor scales with the data it is given", {
  # Doubling every intensity should not leave the cofactor unchanged: the
  # cofactor is what maps the background band onto a fixed IQR, so a wider band
  # needs a larger one. This is the property that makes deriving it worthwhile at
  # all, and a fast path that quietly returned a constant would pass every test
  # above but fail this one.
  withr::local_seed(7)
  base <- matrix(abs(stats::rnorm(20000 * 3, 0, 200)), ncol = 3,
                 dimnames = list(NULL, c("A", "B", "C")))
  mc <- stats::setNames(seq_len(3), colnames(base))
  small <- as.numeric(derive_cofactor(base, mc))
  large <- as.numeric(derive_cofactor(base * 4, mc))
  expect_gt(large, small * 2)
})

test_that("derive_cofactor_pooled resists one aberrant sample", {
  withr::local_seed(3)
  mk <- stats::setNames(seq_len(3), c("A", "B", "C"))
  mkfake <- function(scale) list(
    exprs = matrix(abs(stats::rnorm(8000 * 3, 0, 300 * scale)), ncol = 3,
                   dimnames = list(NULL, names(mk))),
    marker_cols = mk)

  # One badly-scaled file, three normal ones. The median must ignore the outlier;
  # deriving from the FIRST file alone — the behaviour this replaced — would take
  # the outlier as the truth for the whole panel.
  reads <- list(bad = mkfake(0.02), a = mkfake(1), b = mkfake(1), c = mkfake(1))
  from_first <- as.numeric(derive_cofactor(reads$bad$exprs, mk))
  pooled <- suppressMessages(derive_cofactor_pooled(reads, names(reads)))

  expect_identical(attr(pooled, "n_samples"), 4L)
  expect_gt(as.numeric(pooled), from_first * 2)
  expect_identical(attr(pooled, "source"), "pooled_median")
})

test_that("maybe_compensate leaves unmixed data alone and reports why", {
  ex <- matrix(stats::rnorm(300), ncol = 3,
               dimnames = list(NULL, c("A", "B", "C")))
  expect_message(out <- maybe_compensate(ex, list()), "no spillover matrix")
  expect_identical(out, ex)
})

test_that("maybe_compensate survives a spillover matrix with no rownames", {
  # FCS writers frequently emit column names only. Indexing the matrix by
  # character row name then fails with "subscript out of bounds" on a file that
  # is otherwise perfectly readable, so the implementation matches positionally.
  ex <- matrix(stats::rnorm(300), ncol = 3,
               dimnames = list(NULL, c("A", "B", "C")))
  sp <- diag(3); colnames(sp) <- c("A", "B", "C")   # rownames deliberately NULL
  expect_message(out <- maybe_compensate(ex, list(`$SPILLOVER` = sp)),
                 "applying spillover")
  expect_equal(unname(out), unname(ex), tolerance = 1e-12)
})
