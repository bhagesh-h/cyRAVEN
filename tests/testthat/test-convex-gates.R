# =============================================================================
# Learned convex gates: is the geometry right, is the gradient right, and does
# it find a gate that is actually there?
# =============================================================================

# A target blob sitting on a DIAGONAL boundary. Chosen deliberately: an
# axis-aligned rectangle -- which is what the population spec can express --
# cannot separate this without either admitting non-targets or losing targets,
# so a test that passes here is testing the thing convex gates are for.
synth_gate_data <- function(n = 3000, seed = 11) {
  withr::local_seed(seed)
  a <- stats::rnorm(n); b <- stats::rnorm(n)
  X <- cbind(CD4 = a, CD8 = b, CD3 = stats::rnorm(n), CD19 = stats::rnorm(n))
  y <- as.integer(a + b > 1.6 & a - b > -1.2 & a - b < 1.2)
  list(X = X, y = y)
}

test_that("the analytic gradient matches a numeric one", {
  # The whole reason this package needs no automatic-differentiation dependency
  # is that gate_objective() carries its own derivative. If that derivative is
  # wrong the optimiser still returns something, and it is silently not the
  # minimum -- which is exactly the kind of failure that survives to publication.
  withr::local_seed(3)
  Z  <- matrix(stats::runif(400), ncol = 2L)
  y  <- as.integer(rowSums(Z) > 1)
  wt <- ifelse(y == 1, 3, 1)
  fw <- cbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  ct <- colMeans(Z[y == 1, , drop = FALSE])

  # 4 free planes (one angle each) + 8 offsets = 12 parameters.
  for (lam in c(0, 2.5)) {
    th <- c(stats::runif(4, 0, 2 * pi), stats::rnorm(8) * 0.3)
    f <- function(t) as.numeric(gate_objective(t, Z, y, wt, fw, ct, 40, lam))
    analytic <- attr(gate_objective(th, Z, y, wt, fw, ct, 40, lam), "gradient")
    h <- 1e-6
    numeric <- vapply(seq_along(th), function(i) {
      tp <- th; tm <- th; tp[i] <- tp[i] + h; tm[i] <- tm[i] - h
      (f(tp) - f(tm)) / (2 * h)
    }, numeric(1))
    expect_equal(analytic, numeric, tolerance = 1e-4)
  }
})

test_that("half-plane intersection returns the polygon it should", {
  # The unit square, given as four half-planes. Any correct implementation
  # returns its four corners; an implementation that mishandles the feasibility
  # test returns extra points that lie outside.
  W <- cbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  b <- c(0, 1, 0, 1)
  p <- halfplane_polygon(W, b)
  expect_identical(nrow(p), 4L)
  expect_equal(sort(unique(round(as.vector(p), 8))), c(0, 1))

  # An unbounded wedge has no polygon, and must say so rather than return a
  # partial one. This is why four normals are pinned to the PC axes.
  expect_null(halfplane_polygon(cbind(c(1, 0), c(0, 1)), c(0, 0)))
})

test_that("point-in-polygon agrees with the half-plane test it replaces", {
  withr::local_seed(5)
  W <- cbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1))
  b <- c(0, 1, 0, 1, -0.3)
  poly <- halfplane_polygon(W, b)
  pts <- matrix(stats::runif(2000, -0.2, 1.2), ncol = 2L)
  by_planes <- apply(pts, 1L, function(x) all(as.vector(crossprod(W, x)) + b > 0))
  by_poly   <- point_in_polygon(pts, poly)
  # Points within a whisker of an edge are allowed to differ: the two tests
  # break the boundary the same way only up to floating point.
  d <- apply(pts, 1L, function(x) min(abs(as.vector(crossprod(W, x)) + b)))
  expect_identical(by_poly[d > 1e-6], by_planes[d > 1e-6])
})

test_that("a diagonal population is recovered better than any rectangle can", {
  d <- synth_gate_data()
  g <- learn_convex_gate(d$X[, c("CD4", "CD8")], d$y, seed = 1)
  expect_false(is.null(g))
  expect_gt(g$metrics[["f1"]], 0.8)

  # The comparison that justifies the whole file: the best axis-aligned
  # rectangle, which is the most the population spec can express for this pair.
  # It is fitted generously -- the exact bounding box of the targets -- and the
  # convex gate must still beat it.
  tg <- d$X[d$y == 1, c("CD4", "CD8"), drop = FALSE]
  box <- apply(tg, 2L, range)
  in_box <- d$X[, "CD4"] >= box[1, 1] & d$X[, "CD4"] <= box[2, 1] &
            d$X[, "CD8"] >= box[1, 2] & d$X[, "CD8"] <= box[2, 2]
  expect_gt(g$metrics[["f1"]], gate_metrics(in_box, d$y)[["f1"]])
})

test_that("the polygon stays within the range of the data", {
  # REGRESSION. Parametrising the free normals as free VECTORS makes the
  # objective scale-degenerate: (w, b) and (cw, cb) describe the same half-plane
  # but sharpen the sigmoid, so on separable data the optimiser drives ||w|| to
  # infinity. The offsets grow with it until the PC-box planes are satisfied
  # everywhere, stop bounding anything, and the polygon escapes to coordinates
  # thousands of times the data range -- while still reporting a high F1, which
  # is what made it survive a first round of testing. Angles cannot do this.
  d <- synth_gate_data()
  for (sd in 1:3) {
    g <- learn_convex_gate(d$X[, c("CD4", "CD8")], d$y, seed = sd)
    expect_false(is.null(g))
    rg <- apply(d$X[, c("CD4", "CD8")], 2L, range)
    pad <- 0.2 * (rg[2, ] - rg[1, ])
    expect_true(all(g$polygon[, 1L] >= rg[1, 1] - pad[1] &
                    g$polygon[, 1L] <= rg[2, 1] + pad[1]))
    expect_true(all(g$polygon[, 2L] >= rg[1, 2] - pad[2] &
                    g$polygon[, 2L] <= rg[2, 2] + pad[2]))
  }
})

test_that("recall does not depend on the subsample ceiling", {
  # REGRESSION. n_target_total was captured before the max_cells subsample and
  # then used as the denominator for counts taken after it, so recall was
  # divided by the subsample ratio: exactly 0.5 at the default ceiling against
  # 40,000 cells. It reads as "the gate loses half the population" while
  # precision sits at 1.00, which looks like a conservative gate rather than an
  # arithmetic error -- the reason it survived the first round of tests.
  withr::local_seed(7)
  n <- 12000
  a <- stats::rnorm(n); b <- stats::rnorm(n)
  X <- cbind(M1 = a, M2 = b)
  y <- as.integer(a > 1.2 & b > 0.6)

  full <- learn_convex_gate(X, y, max_cells = 20000L, seed = 1)   # no subsample
  half <- learn_convex_gate(X, y, max_cells = 6000L,  seed = 1)   # 2x subsample
  expect_false(is.null(full)); expect_false(is.null(half))
  # Subsampling changes which cells are seen, so the gates differ a little. What
  # must NOT differ is the SCALE of recall: a factor-of-two gap is the bug.
  expect_lt(abs(full$metrics[["recall"]] - half$metrics[["recall"]]), 0.2)
  expect_gt(half$metrics[["recall"]], 0.7)
})

test_that("a population living in the far tail is still recovered", {
  # The realistic shape of an undescribed cluster: a corner of the marker space,
  # largely made of cells beyond the percentiles used to set the working scale.
  withr::local_seed(7)
  n <- 20000
  a <- stats::rnorm(n); b <- stats::rnorm(n)
  y <- as.integer(a > 1.8 & b > 0.9)
  g <- learn_convex_gate(cbind(M1 = a, M2 = b), y, seed = 1)
  expect_false(is.null(g))
  expect_gt(g$metrics[["recall"]], 0.7)
  expect_gt(g$metrics[["f1"]], 0.75)
})

test_that("held-out metrics are reported and are not the in-sample ones", {
  d <- synth_gate_data()
  g <- learn_convex_gate(d$X[, c("CD4", "CD8")], d$y, seed = 2)
  expect_true(all(c("precision", "recall", "f1") %in% names(g$metrics)))
  expect_false(is.null(g$metrics_insample))
  # Both are finite and in range; the point is that BOTH exist, so a reader can
  # see the optimism rather than being handed one number.
  expect_true(all(g$metrics[c("precision", "recall", "f1")] >= 0))
  expect_true(all(g$metrics[c("precision", "recall", "f1")] <= 1))
  expect_gte(g$n_test, 1L)
  expect_gt(g$n_train, g$n_test)
})

test_that("marker ranking puts the informative markers first", {
  d <- synth_gate_data()
  rk <- rank_gate_markers(d$X, d$y)
  expect_setequal(rk$marker[1:2], c("CD4", "CD8"))
  # And the noise markers score near zero rather than merely lower.
  expect_lt(max(rk$score[rk$marker %in% c("CD3", "CD19")]), 0.01)
})

test_that("recall is cumulative down a hierarchy, not per level", {
  d <- synth_gate_data()
  st <- explain_cluster(d$X, d$y, max_depth = 3L, seed = 4)
  expect_false(is.null(st))
  expect_true(all(diff(st$summary$cumulative_recall) <= 1e-9))   # monotone down
  expect_lte(max(st$summary$depth), 3L)
  # Each level's gate must actually shrink the surviving set.
  expect_true(all(diff(st$summary$cells_in_gate) < 0))
})

test_that("a population with no signal yields no strategy rather than a bad one", {
  withr::local_seed(9)
  X <- matrix(stats::rnorm(3000 * 3), ncol = 3L,
              dimnames = list(NULL, c("A", "B", "C")))
  y <- stats::rbinom(3000, 1L, 0.25)          # label independent of every marker
  st <- explain_cluster(X, y, max_depth = 3L, seed = 6)
  # It may fit one gate -- noise always admits some polygon -- but it must not
  # claim to separate anything.
  if (!is.null(st)) expect_lt(max(st$summary$holdout_f1, na.rm = TRUE), 0.6)
})

test_that("too few target cells returns NULL, not a fitted-to-nothing gate", {
  d <- synth_gate_data()
  y <- integer(length(d$y)); y[1:5] <- 1L
  expect_null(learn_convex_gate(d$X[, 1:2], y))
  expect_null(explain_cluster(d$X, y))
})

test_that("learning a gate restores the RNG stream", {
  d <- synth_gate_data()
  set.seed(77); before <- stats::runif(3)
  set.seed(77)
  invisible(learn_convex_gate(d$X[, c("CD4", "CD8")], d$y, seed = 3))
  expect_identical(before, stats::runif(3))
})

test_that("the strategy figure is written and is not blank", {
  d <- synth_gate_data()
  st <- explain_cluster(d$X, d$y, max_depth = 2L, seed = 8)
  f <- file.path(withr::local_tempdir(), "strategy.png")
  suppressMessages(fig_gate_strategy(d$X, d$y, st, f, label = "test population"))
  expect_true(file.exists(f))
  expect_gt(file.size(f), 5000)
})
