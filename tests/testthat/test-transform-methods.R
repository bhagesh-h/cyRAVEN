# =============================================================================
# Transform choice and batch correction
# =============================================================================

test_that("logicle parameters follow the automatic rule", {
  # w = (m - log10(t/|r|))/2 with r a low quantile of the negatives. Constructed
  # so the answer is checkable by hand rather than by re-running the function.
  withr::local_seed(3)
  x <- c(stats::rnorm(5000, -200, 60), stats::rnorm(5000, 20000, 4000))
  p <- logicle_params_one(x, m = 4.5)
  r <- stats::quantile(x[x < 0], 0.05, names = FALSE)
  expect_equal(p$t, max(x))
  expect_equal(p$w, (4.5 - log10(max(x) / abs(r))) / 2, tolerance = 1e-9)
  expect_gte(p$w, 0); expect_lte(p$w, 4.5 / 2)
})

test_that("a channel with no negative values degenerates to log, not an error", {
  # Nothing below zero means there is no negative population for a linear region
  # to rescue. w = 0 is the right answer, not a failure.
  p <- logicle_params_one(abs(stats::rnorm(1000, 5000, 500)), m = 4.5)
  expect_identical(p$w, 0)
})

test_that("the transform object applies the method it was built for", {
  x <- c(-500, 0, 100, 10000)
  a <- make_transform("arcsinh", cofactor = 150)
  expect_equal(a$fn(x), asinh(x / 150))
  expect_identical(a$method, "arcsinh")

  n <- make_transform("none")
  expect_identical(n$fn(x), x)

  lg <- make_transform("logicle", logicle = list(CD3 = logicle_params_one(x)))
  expect_length(lg$fn(x, "CD3"), 4L)
  # Monotone: a transform that reordered intensities would move cells across
  # gates for no reason connected to the data.
  expect_false(is.unsorted(lg$fn(sort(x), "CD3")))
})

test_that("a transform refuses the inputs it cannot honour", {
  expect_error(make_transform("arcsinh", cofactor = 0), "positive cofactor")
  expect_error(make_transform("arcsinh", cofactor = NULL), "positive cofactor")
  expect_error(make_transform("logicle"), "derive_logicle_pooled")
  # An unparameterised marker must stop rather than fall back to a default
  # scale: one channel silently on a different ruler is the failure the pooled
  # derivation exists to prevent.
  lg <- make_transform("logicle", logicle = list(CD3 = logicle_params_one(1:100)))
  expect_error(lg$fn(1:10, "CD4"), "no logicle parameters")
})

test_that("logicle keeps compensated negatives on scale", {
  # The property that motivates offering logicle at all: arcsinh with a large
  # derived cofactor compresses the negative population toward zero, while
  # logicle spreads it over its linear region.
  withr::local_seed(5)
  x <- c(stats::rnorm(4000, -300, 80), stats::rnorm(4000, 30000, 5000))
  lg <- make_transform("logicle", logicle = list(M = logicle_params_one(x)))
  neg <- x[x < 0]
  expect_gt(diff(range(lg$fn(neg, "M"))), 0.1)
})

test_that("Cramer's V is 1 when one variable determines the other", {
  b <- rep(c("run1", "run2"), each = 20)
  expect_equal(cramers_v(b, b), 1, tolerance = 1e-8)
  # Independent labels sit near zero; the exact value is sampling noise, so the
  # assertion is on the scale rather than a point value.
  withr::local_seed(2)
  expect_lt(cramers_v(b, sample(c("A", "B"), 40, TRUE)), 0.5)
  expect_true(is.na(cramers_v(rep("only", 10), rep(c("x", "y"), 5))))
})

test_that("quantile alignment removes a batch shift and preserves order", {
  withr::local_seed(11)
  n <- 3000
  batch <- rep(c("b1", "b2"), each = n / 2)
  x <- c(stats::rnorm(n / 2, 0, 1), stats::rnorm(n / 2, 2.5, 1))   # b2 shifted
  before <- abs(diff(tapply(x, batch, stats::median)))
  after  <- abs(diff(tapply(align_quantiles(x, batch), batch, stats::median)))
  expect_gt(before, 2)
  expect_lt(after, 0.15)

  # Monotone within batch: the map may not reorder cells, or it would invent
  # structure that was never measured.
  i <- which(batch == "b2")
  expect_identical(order(align_quantiles(x, batch)[i]), order(x[i]))
})

test_that("batch correction refuses when batch is confounded with group", {
  # Perfect confounding: each batch is exactly one group. Correcting here would
  # remove the biology along with the batch, so it must decline and say so.
  withr::local_seed(4)
  tm <- matrix(stats::rnorm(600), ncol = 3,
               dimnames = list(NULL, c("CD3", "CD4", "CD8")))
  batch <- rep(c("b1", "b2"), each = 100)
  group <- batch
  r <- suppressMessages(correct_batch(tm, batch, group))
  expect_false(r$corrected)
  expect_identical(r$tmat, tm)
  expect_match(r$reason, "^refused")
  expect_equal(r$cramers_v, 1, tolerance = 1e-8)

  # Overridable, and the override is recorded rather than silent.
  f <- suppressMessages(correct_batch(tm, batch, group, force = TRUE))
  expect_true(f$corrected)
  expect_match(f$reason, "FORCED")
})

test_that("batch correction proceeds when batch and group are separable", {
  withr::local_seed(6)
  tm <- matrix(stats::rnorm(1200), ncol = 3,
               dimnames = list(NULL, c("CD3", "CD4", "CD8")))
  batch <- rep(c("b1", "b2"), each = 200)
  group <- rep(rep(c("case", "ctrl"), each = 100), 2)   # crossed with batch
  tm[batch == "b2", ] <- tm[batch == "b2", ] + 1.5
  r <- suppressMessages(correct_batch(tm, batch, group))
  expect_true(r$corrected)
  expect_lt(abs(diff(tapply(r$tmat[, "CD3"], batch, stats::median))), 0.2)
})

test_that("a single batch level is a no-op, not an error", {
  tm <- matrix(stats::rnorm(60), ncol = 3,
               dimnames = list(NULL, c("a", "b", "c")))
  r <- correct_batch(tm, rep("one", 20))
  expect_false(r$corrected)
  expect_identical(r$tmat, tm)
  expect_match(r$reason, "only one batch")
})
