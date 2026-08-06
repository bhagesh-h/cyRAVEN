# =============================================================================
# Gate placement uncertainty, and the RNG discipline that keeps it additive
# =============================================================================

#' Two modes separated by `gap`, so the depth of the valley is known in advance.
planted_bimodal <- function(gap, n = 8000L, sd = 0.6, seed = 3) {
  withr::local_seed(seed)
  c(stats::rnorm(n, 0, sd), stats::rnorm(n, gap, sd))
}

test_that("density_valley(details = TRUE) reports the same cut it always did", {
  # The details argument was added so the uncertainty code could read the valley
  # depth without recomputing the histogram somewhere else. It must not have
  # changed a single threshold in the process.
  for (gap in c(3, 4, 6)) {
    x <- planted_bimodal(gap)
    expect_identical(density_valley(x), density_valley(x, details = TRUE)$cut)
  }
  # And on input with no valley, both forms still decline, with a reason. Two
  # cases, one for each way the function can refuse: too few events to judge, and
  # a distribution with no separating gap.
  expect_true(is.na(density_valley(stats::rnorm(100))))
  d1 <- density_valley(stats::rnorm(100), details = TRUE)
  expect_true(is.na(d1$cut))
  expect_match(d1$reason, "fewer than 500")

  withr::local_seed(1)
  flat <- stats::runif(20000)
  expect_true(is.na(density_valley(flat)))
  d2 <- density_valley(flat, details = TRUE)
  expect_true(is.na(d2$cut))
  expect_true(nzchar(d2$reason))
})

test_that("resampling exposes a threshold that only noise made bimodal", {
  # The estimator finds a spurious valley in a plain normal a fair fraction of
  # the time at a few thousand events. That is exactly what the bootstrap rate is
  # for: a cut that survives resampling is a population boundary, a cut that
  # comes and goes is histogram noise, and the reported threshold looks identical
  # in both cases.
  withr::local_seed(2)
  noise <- stats::rnorm(4000)
  real  <- planted_bimodal(6)
  u_noise <- threshold_uncertainty(noise, B = 60L, seed = 4)
  u_real  <- threshold_uncertainty(real,  B = 60L, seed = 4)
  expect_gt(u_real$valley_rate, u_noise$valley_rate)
})

test_that("a shallow valley carries more uncertainty than a deep one", {
  # THE CLAIM THE WHOLE FEATURE RESTS ON. Two populations, one cleanly separated
  # and one barely, must not be reported with the same precision.
  deep    <- threshold_uncertainty(planted_bimodal(7), B = 60L, seed = 5)
  shallow <- threshold_uncertainty(planted_bimodal(2.6), B = 60L, seed = 5)
  expect_true(is.finite(deep$u_combined))
  expect_true(is.finite(shallow$u_combined))
  expect_gt(shallow$u_combined, deep$u_combined)
  # The deep valley is also the one the bootstrap finds every time.
  expect_gt(deep$valley_rate, shallow$valley_rate - 1e-9)
  expect_gt(deep$rel_depth, shallow$rel_depth)
})

test_that("a declared threshold has zero uncertainty, not small uncertainty", {
  u <- threshold_uncertainty(planted_bimodal(5), source = "config")
  expect_identical(u$u_combined, 0)
  expect_identical(u$basis, "declared")
})

test_that("a quantile fallback reports how arbitrary it is", {
  # No valley exists here, so the cut is a quantile. The method component sweeps
  # WHICH quantile, and on a unimodal distribution that spread should not be
  # mistaken for a well-determined threshold.
  u <- threshold_uncertainty(stats::rnorm(6000), source = "quantile_fallback",
                             B = 60L, seed = 9)
  expect_true(is.finite(u$u_combined))
  expect_gt(u$u_method, 0)
  expect_match(u$basis, "fallback quantile")
})

test_that("a threshold from a control tube says so rather than guessing", {
  u <- threshold_uncertainty(planted_bimodal(5), source = "control_q995")
  expect_true(is.na(u$u_combined))
  expect_match(u$basis, "control tube")
})

test_that("the uncertainty machinery leaves the RNG stream exactly where it was", {
  # THE REGRESSION THIS FILE EXISTS FOR. run_cyraven() seeds once and STEP 6's
  # UMAP cell selection draws from that one stream. A step that consumed draws
  # without restoring the stream would change which cells are embedded and
  # silently redraw every UMAP figure in the run, with no error anywhere.
  set.seed(123)
  before_seed <- .Random.seed
  expected <- stats::runif(5)

  set.seed(123)
  invisible(threshold_uncertainty(planted_bimodal(5), B = 40L, seed = 77))
  expect_identical(.Random.seed, before_seed)
  expect_identical(stats::runif(5), expected)
})

test_that("uncertainty does not depend on marker order or on what ran before", {
  # Each threshold draws from a stream keyed on its own names, so reordering a
  # panel cannot change a published number.
  x <- planted_bimodal(4)
  a <- threshold_uncertainty(x, B = 40L, seed = 11)
  invisible(stats::runif(1000))
  b <- threshold_uncertainty(x, B = 40L, seed = 11)
  expect_identical(a$u_combined, b$u_combined)
})

test_that("frequency uncertainty grows with the uncertainty of the cut", {
  withr::local_seed(21)
  n <- 4000L
  pos <- rep(c(TRUE, FALSE), each = n / 2)
  tmat <- cbind(CD3 = ifelse(pos, stats::rnorm(n, 4, 0.6), stats::rnorm(n, 0, 0.6)),
                CD4 = stats::rnorm(n, 2, 1))
  thr <- c(CD3 = 2, CD4 = 2)
  spec <- list(`T cells` = list(CD3 = "above"))
  parent <- rep(TRUE, n)

  tight <- population_frequency_uncertainty(tmat, thr, parent, spec,
                                            c(CD3 = 0.05, CD4 = 0.05))
  loose <- population_frequency_uncertainty(tmat, thr, parent, spec,
                                            c(CD3 = 0.60, CD4 = 0.05))
  expect_equal(tight$per_population$pct_of_cd45_pos,
               loose$per_population$pct_of_cd45_pos)
  expect_gt(loose$per_population$u_pct_points, tight$per_population$u_pct_points)
  # A marker the population does not read contributes nothing to it.
  expect_false("CD4" %in% tight$budget$term)
})

test_that("annotating a comparison appends columns and changes none of them", {
  gs <- data.frame(population = c("A", "B"), median_reference = c(10, 5),
                   median_comparison = c(14, 5.2), p_value = c(0.01, 0.03),
                   stringsAsFactors = FALSE)
  uf <- data.frame(sample_id = rep(c("S1", "S2"), each = 2),
                   population = rep(c("A", "B"), 2),
                   u_pct_points = c(1, 4, 1, 4), qc_status = "pass",
                   stringsAsFactors = FALSE)
  out <- annotate_gate_uncertainty(gs, uf)

  expect_identical(out[, names(gs)], gs)
  expect_true(all(c("gate_u_pct_points", "difference_over_gate_u") %in% names(out)))
  # A is a four-point difference on a one-point uncertainty: resolvable.
  # B is a fifth of a point on a four-point uncertainty: not.
  expect_gt(out$difference_over_gate_u[out$population == "A"], 1)
  expect_lt(out$difference_over_gate_u[out$population == "B"], 1)
})

test_that("annotation with no uncertainty table still returns the table", {
  gs <- data.frame(population = "A", median_reference = 1, median_comparison = 2,
                   stringsAsFactors = FALSE)
  out <- annotate_gate_uncertainty(gs, NULL)
  expect_identical(out[, names(gs)], gs)
  expect_true(all(is.na(out$gate_u_pct_points)))
})
