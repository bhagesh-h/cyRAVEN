# =============================================================================
# Gate placement uncertainty, and the RNG discipline that keeps it additive
# =============================================================================

#' Two modes separated by `gap`, so the depth of the valley is known in advance.
#'
#' WHY THIS SAVES AND RESTORES BY HAND RATHER THAN CALLING withr::local_seed().
#' The RNG regression test below asserts that the package leaves .Random.seed
#' exactly as it found it. That assertion is only about the package if the
#' FIXTURE is stream-neutral too. withr 3.0.2's local_seed() restores the
#' generator's position but not its state vector, so a fixture built with it
#' leaks a changed stream into the test and the assertion fails while pointing at
#' the wrong code. Using the same guard the package uses keeps the test measuring
#' the package.
planted_bimodal <- function(gap, n = 8000L, sd = 0.6, seed = 3) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
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

# =============================================================================
# Counting uncertainty and detection limits
# =============================================================================

test_that("an unobserved population is not reported as certainly absent", {
  # THE REASON WILSON IS USED AND THE TEXTBOOK FORMULA IS NOT. sqrt(p(1-p)/n)
  # evaluates to exactly zero at k = 0, which claims perfect knowledge of a
  # population nobody saw. That is the regime the number is consulted in.
  z <- counting_uncertainty(0, 10000)
  expect_gt(z$u_counting_pct_points, 0)
  expect_identical(z$detection, "below LOD")

  # And it still behaves at the other boundary.
  full <- counting_uncertainty(10000, 10000)
  expect_gt(full$u_counting_pct_points, 0)
})

test_that("counting uncertainty agrees with the textbook value once counts are large", {
  # Wilson is used for its behaviour at the boundary, not because it disagrees
  # with the ordinary binomial standard error where that one is valid. If the two
  # parted company in the common case, the column would not be comparable with
  # anything else in the literature.
  n <- 200000; k <- 60000
  wald <- 100 * sqrt((k / n) * (1 - k / n) / n)
  got  <- counting_uncertainty(k, n)$u_counting_pct_points
  expect_equal(got, wald, tolerance = 1e-3)
})

test_that("relative precision improves with the number of events counted", {
  # The quantity a reader acts on is the uncertainty RELATIVE to the frequency.
  # In absolute percentage points a 50% population is the least certain of all,
  # which is correct and is not what "a rare population is unreliable" means.
  k <- c(10, 100, 1000, 10000)
  d <- counting_uncertainty(k, 1e6)
  rel <- d$u_counting_pct_points / (100 * k / 1e6)
  expect_true(all(diff(rel) < 0))

  # The absolute form peaks in the middle, which is the binomial behaving.
  a <- counting_uncertainty(c(1000, 500000, 999000), 1e6)$u_counting_pct_points
  expect_gt(a[2], a[1])
  expect_gt(a[2], a[3])
})

test_that("the limits scale with the denominator this run actually saw", {
  # --max-events-per-file lowers the parent count and therefore raises every
  # limit in proportion. Stating that is the point of expressing them as
  # percentages rather than as the fixed event counts they come from.
  a <- counting_uncertainty(100, 100000)
  b <- counting_uncertainty(100, 10000)
  expect_equal(b$lod_pct, a$lod_pct * 10)
  expect_equal(b$loq_pct, a$loq_pct * 10)
  expect_equal(a$lod_pct, 100 * 20 / 100000)
})

test_that("the detection verdict changes exactly at the declared event counts", {
  d <- counting_uncertainty(c(19, 20, 49, 50), 100000)
  expect_identical(d$detection,
                   c("below LOD", "detected, below LOQ",
                     "detected, below LOQ", "quantified"))
  # And the thresholds are arguments, not constants.
  d2 <- counting_uncertainty(c(19, 20), 100000, lod_events = 5L, loq_events = 10L)
  expect_identical(d2$detection, c("quantified", "quantified"))
})

test_that("impossible or missing counts return NA rather than a number", {
  d <- counting_uncertainty(c(-1, NA, 11, 5), 10)
  expect_true(all(is.na(d$u_counting_pct_points[1:3])))
  expect_true(all(is.na(d$detection[1:3])))
  expect_false(is.na(d$u_counting_pct_points[4]))
  expect_true(is.na(counting_uncertainty(5, 0)$u_counting_pct_points))
})

test_that("two populations behind the same sharp cut are told apart by their counts", {
  # THE CLAIM THE FEATURE RESTS ON. Before this, a population of a handful of
  # events and one of thousands, both separated by the same clean valley, were
  # reported with the same uncertainty, because displacing a well-separated cut
  # moves neither of them much. Their gate terms are still alike; what separates
  # them is how many cells there were to count.
  withr::local_seed(31)
  n <- 20000L
  common <- 4000L; rare <- 12L
  x <- c(stats::rnorm(common, 6, 0.4),      # abundant, well above the cut
         stats::rnorm(rare, 6, 0.4),        # scarce, same place
         stats::rnorm(n - common - rare, 0, 0.4))
  tmat <- cbind(A = x, B = c(rep(5, common), rep(-5, rare),
                             rep(-5, n - common - rare)))
  thr <- c(A = 3, B = 0)
  spec <- list(Common = list(A = "above", B = "above"),
               Rare   = list(A = "above", B = "below"))
  fu <- population_frequency_uncertainty(tmat, thr, rep(TRUE, n), spec,
                                         c(A = 0.02, B = 0.02))
  pp <- fu$per_population
  rownames(pp) <- pp$population

  expect_identical(pp["Common", "n_cells"], common)
  expect_identical(pp["Rare", "n_cells"], rare)
  expect_identical(pp["Common", "detection"], "quantified")
  expect_identical(pp["Rare", "detection"], "below LOD")

  # Relative to its own frequency the rare population is far less certain, which
  # is the statement the gate terms alone could not make.
  rel <- pp$u_total_pct_points / pp$pct_of_cd45_pos
  expect_gt(rel[pp$population == "Rare"], rel[pp$population == "Common"])
})

test_that("the counting term is appended and never folded into the gate term", {
  # ADDITIVITY. u_pct_points must keep meaning gate placement alone and must hold
  # the value it held before this feature existed, whatever the limits are set
  # to. Only the new columns may move.
  withr::local_seed(32)
  n <- 6000L
  pos <- rep(c(TRUE, FALSE), each = n / 2)
  tmat <- cbind(CD3 = ifelse(pos, stats::rnorm(n, 4, 0.6), stats::rnorm(n, 0, 0.6)))
  thr <- c(CD3 = 2)
  spec <- list(`T cells` = list(CD3 = "above"))
  parent <- rep(TRUE, n)

  # Roughly 3000 of the 6000 events clear the cut, so the second call brackets
  # that count between its two limits and the first leaves it above both.
  a <- population_frequency_uncertainty(tmat, thr, parent, spec, c(CD3 = 0.1))
  b <- population_frequency_uncertainty(tmat, thr, parent, spec, c(CD3 = 0.1),
                                        lod_events = 2000L, loq_events = 5000L)
  keep <- c("population", "pct_of_cd45_pos", "u_pct_points", "n_terms",
            "n_terms_missing")
  expect_identical(a$per_population[, keep], b$per_population[, keep])
  expect_identical(a$budget, b$budget)
  # Moving the limits moves only the verdict.
  expect_identical(a$per_population$detection, "quantified")
  expect_identical(b$per_population$detection, "detected, below LOQ")

  # The total is the quadrature sum, and never smaller than either part.
  pp <- a$per_population
  expect_equal(pp$u_total_pct_points,
               round(sqrt(pp$u_pct_points^2 + pp$u_counting_pct_points^2), 4))
  expect_gte(pp$u_total_pct_points, pp$u_pct_points)
})

test_that("the frequency itself is untouched by the counting analysis", {
  # The percentage is scored from the same masks whether the counts are read off
  # them or not, so it must be exactly count/parent, not a value recovered from a
  # rounded percentage.
  withr::local_seed(33)
  n <- 5000L
  tmat <- cbind(CD3 = stats::rnorm(n))
  fu <- population_frequency_uncertainty(tmat, c(CD3 = 0), rep(TRUE, n),
                                         list(P = list(CD3 = "above")),
                                         c(CD3 = 0.05))
  pp <- fu$per_population
  expect_identical(pp$n_parent_events, n)
  expect_equal(pp$pct_of_cd45_pos, 100 * pp$n_cells / n)
})

test_that("the comparison gains a stricter ratio without disturbing the first", {
  gs <- data.frame(population = c("A", "B"), median_reference = c(10, 5),
                   median_comparison = c(14, 5.2), stringsAsFactors = FALSE)
  uf <- data.frame(sample_id = rep(c("S1", "S2"), each = 2),
                   population = rep(c("A", "B"), 2),
                   u_pct_points = c(1, 4, 1, 4),
                   u_total_pct_points = c(2, 5, 2, 5),
                   qc_status = "pass", stringsAsFactors = FALSE)
  old <- annotate_gate_uncertainty(gs, uf[, setdiff(names(uf), "u_total_pct_points")])
  new <- annotate_gate_uncertainty(gs, uf)

  # The pre-existing columns are identical with and without the new one present.
  expect_identical(new$gate_u_pct_points, old$gate_u_pct_points)
  expect_identical(new$difference_over_gate_u, old$difference_over_gate_u)
  # The total-based ratio is the stricter of the two, always.
  expect_true(all(new$difference_over_total_u <= new$difference_over_gate_u))
  # A run predating the counting term leaves the new columns NA rather than
  # borrowing the gate value.
  expect_true(all(is.na(old$difference_over_total_u)))
})
