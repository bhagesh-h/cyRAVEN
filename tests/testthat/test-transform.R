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
# notice the run got 100x slower, or, worse, replace it with an approximation
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
  # deriving from the FIRST file alone, the behaviour this replaced, would take
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
  # The status is reported as unknown, not as "already compensated": an absent
  # keyword cannot tell unmixed data from conventional data nobody compensated.
  expect_match(attr(out, "compensation"), "^unknown")
  attr(out, "compensation") <- NULL
  expect_identical(out, ex)
})

test_that("a supplied spillover matrix replaces the file's own", {
  ex <- cbind(`FITC.A` = c(1000, 0, 500), `PE.A` = c(250, 800, 125))
  own <- diag(2); colnames(own) <- c("FITC.A", "PE.A")
  # row i: the share of detector i's fluorochrome seen in each detector
  sup <- matrix(c(1, 0.25, 0, 1), 2, byrow = TRUE,
                dimnames = list(NULL, c("FITC-A", "PE-A")))
  out <- suppressMessages(maybe_compensate(ex, list(`$SPILLOVER` = own), external = sup))
  # PE-A receives 25% of FITC, which compensation removes
  expect_equal(unname(out[, "PE.A"]), c(0, 800, 0), tolerance = 1e-9)
  expect_match(attr(out, "compensation"), "supplied")
  f <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(as.data.frame(sup, check.names = FALSE), f, row.names = FALSE)
  expect_equal(unname(read_spillover_file(f)), unname(sup))
  bad <- sup; bad[1, 1] <- 0.5
  utils::write.csv(as.data.frame(bad, check.names = FALSE), f, row.names = FALSE)
  expect_error(read_spillover_file(f), "diagonal")
})

test_that("per-marker cofactors come from each marker's negative population", {
  # A marker almost every cell carries (CD45) and a mostly negative one with a
  # wide negative population. The bulk criterion takes CD45's background band
  # from its positive cells; the negative-population estimate does not.
  withr::local_seed(11)
  n <- 20000
  ex <- cbind(CD45 = c(rnorm(n * 0.98, 6000, 1500), rnorm(n * 0.02, 0, 100)),
              CCR9 = c(rnorm(n * 0.95, 0, 1200), rnorm(n * 0.05, 15000, 3000)))
  mc <- c(CD45 = 1L, CCR9 = 2L)
  cf <- derive_cofactor(ex, mc)
  neg <- attr(cf, "per_marker_negative")
  bulk <- attr(cf, "per_marker")
  # about 1.6 sigma of each negative population, so they track the twelve-fold
  # difference in negative spread
  expect_equal(neg[["CD45"]], 1.6 * 100, tolerance = 0.25)
  expect_equal(neg[["CCR9"]], 1.6 * 1200, tolerance = 0.25)
  # while the bulk band for CD45 is its positive cells, and says nothing about
  # its negatives
  expect_gt(abs(log(bulk[["CD45"]] / neg[["CD45"]])), log(3))
  tr <- make_transform("arcsinh", cofactor = 300, per_marker = neg)
  expect_equal(tr$fn(1000, "CCR9"), asinh(1000 / neg[["CCR9"]]))
  # a marker without its own cofactor, and scatter, use the panel value
  expect_equal(tr$fn(1000, "CD3"), asinh(1000 / 300))
  expect_equal(tr$inv(tr$fn(777, "CD45"), "CD45"), 777)
})

test_that("maybe_compensate survives a spillover matrix with no rownames", {
  # FCS writers frequently emit column names only. Indexing the matrix by
  # character row name then fails with "subscript out of bounds" on a file that
  # is otherwise perfectly readable, so the implementation matches positionally.
  ex <- matrix(stats::rnorm(300), ncol = 3,
               dimnames = list(NULL, c("A", "B", "C")))
  sp <- diag(3); colnames(sp) <- c("A", "B", "C")   # rownames deliberately NULL
  expect_message(out <- maybe_compensate(ex, list(`$SPILLOVER` = sp)),
                 "applying file \\$SPILLOVER")
  expect_identical(attr(out, "compensation"), "compensated, file $SPILLOVER")
  attr(out, "compensation") <- NULL
  expect_equal(unname(out), unname(ex), tolerance = 1e-12)
})

test_that("declared controls do not set the transform", {
  smap <- data.frame(sample_id = c("s1", "s2", "fmm", "fs"),
                     is_control = c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(transform_samples(c("s1", "fmm", "s2", "fs"), smap),
                   c("s1", "s2"))
  # A panel of nothing but controls keeps them: they are all it has.
  expect_identical(transform_samples(c("fmm", "fs"), smap), c("fmm", "fs"))
  # No sheet, or no is_control column, changes nothing.
  expect_identical(transform_samples(c("a", "b")), c("a", "b"))
  expect_identical(transform_samples(c("a", "b"), data.frame(sample_id = "a")),
                   c("a", "b"))
})

test_that("declared controls are read from the sheet, and only from it", {
  smap <- data.frame(sample_id = c("s1", "fmm"), is_control = c(FALSE, TRUE))
  expect_identical(declared_controls(smap), "fmm")
  expect_identical(declared_controls(NULL), character(0))
  expect_identical(declared_controls(data.frame(sample_id = "s1")), character(0))
})
