# =============================================================================
# Three fixes found by running cyRAVEN against cyCONDOR on a 24-sample cohort
#
# 1. Values above a channel's declared range ($PnR) are capped on reading, so a
#    stray float cannot set the top of the logicle scale.
# 2. No density-derived cut may sit at or below zero signal, so the bump that
#    compensated events form below zero cannot be mistaken for a negative mode.
# 3. Equalised embedding draws stop at the lower quartile of the per-sample
#    counts, so one small sample cannot shrink every other sample to its size.
# =============================================================================

test_that("values above $PnR are capped at it", {
  ex <- cbind(a = c(1, 50, 200, 9e7), b = c(1, 2, 3, 4))
  colnames(ex) <- c("FL1-A", "FL2-A")
  kw <- list(`$P1R` = "262144", `$P2R` = "1024")
  out <- clamp_to_declared_range(ex, kw)
  expect_equal(out[, 1], c(1, 50, 200, 262144))
  expect_equal(out[, 2], ex[, 2])
})

test_that("a channel with many values above $PnR is still capped", {
  # More than 1% over used to leave the channel alone, and one such file set the
  # logicle scale for a whole panel. It is capped now, and named in a warning.
  ex <- cbind(c(rep(100, 90), rep(5e7, 10)))
  colnames(ex) <- "FL1-A"
  out <- clamp_to_declared_range(ex, list(`$P1R` = "262144"))
  expect_equal(max(out), 262144)
})

test_that("a channel without a declared range is left alone", {
  ex <- cbind(c(1, 9e9))
  colnames(ex) <- "Time"
  expect_equal(clamp_to_declared_range(ex, list()), ex)
})

# A mostly negative marker on a logicle-like scale: a small bump of compensated
# events below zero, the large negative peak just above it, and a small positive
# population well above both.
floor_case <- function(seed = 1) {
  set.seed(seed)
  c(stats::rnorm(600, -0.9, 0.12),     # compensated events below zero
    stats::rnorm(20000, 0.9, 0.25),    # negative population
    stats::rnorm(900, 2.8, 0.20))      # positive population
}

test_that("without a floor the valley can land below the negative peak", {
  v <- density_valley(floor_case())
  expect_true(is.finite(v))
  expect_lt(v, 0.9)
})

test_that("a floor at zero signal moves the cut above the negative peak", {
  v <- density_valley(floor_case(), floor = 0)
  expect_true(is.finite(v))
  expect_gt(v, 1.5)
  expect_lt(v, 2.8)
})

test_that("resolve_threshold passes the floor through and reports it", {
  r <- resolve_threshold("Vd1", floor_case(), floor = 0)
  expect_identical(r$source, "valley")
  expect_gt(r$threshold, 1.5)
  expect_identical(r$floor, 0)
})

test_that("a manual or config cut is never floored", {
  expect_identical(resolve_threshold("Vd1", floor_case(), cfg_value = -1, floor = 0)$threshold, -1)
  ov <- list(threshold = -1, reason = "test")
  expect_identical(resolve_threshold("Vd1", floor_case(), override = ov, floor = 0)$threshold, -1)
})

test_that("best_threshold ignores candidates at or below the floor", {
  b <- best_threshold(floor_case(), floor = 0)
  expect_gt(b$threshold, 0)
})

test_that("equalised draws stop at the lower quartile, not the smallest sample", {
  n <- c(a = 393, b = 2000, c = 3000, d = 4000, e = 5000)
  take <- plan_subsample(n, cap = 25000)
  level <- floor(stats::quantile(n, 0.25, names = FALSE))
  expect_equal(unname(take), unname(pmin(level, n)))
  # the small sample gives everything it has, and does not shrink the others
  expect_equal(take[["a"]], 393L)
  expect_gt(take[["e"]], 393L)
})

test_that("--no-equalise takes every cell up to the cap", {
  n <- c(a = 393, b = 2000, c = 30000)
  expect_equal(unname(plan_subsample(n, cap = 25000, equalise = FALSE)),
               c(393L, 2000L, 25000L))
})

test_that("auto scaling is none under logicle and robust otherwise", {
  expect_identical(umap_scale_method(list(transform = "logicle")), "none")
  expect_identical(umap_scale_method(list(transform = "arcsinh")), "robust")
  expect_identical(umap_scale_method(list()), "robust")
  expect_identical(umap_scale_method(list(transform = "logicle", umap_scale = "zscore")), "zscore")
})

# Whole blood inside the leukocyte gate: a small, wide CD45-negative mode
# overlapping a large CD45-dim mode (granulocytes), and a large CD45-bright mode
# (lymphocytes) cleanly apart from it.
cd45_case <- function(seed = 2) {
  set.seed(seed)
  c(stats::rnorm(1500, 1.3, 0.35), stats::rnorm(12000, 2.4, 0.15),
    stats::rnorm(6000, 3.1, 0.12))
}

test_that("the deepest valley is the dim/bright one", {
  v <- density_valley(cd45_case(), floor = 0)
  expect_gt(v, 2.5)
})

test_that("prefer = 'lowest' keeps the dim cells above the cut", {
  v <- density_valley(cd45_case(), floor = 0, prefer = "lowest")
  expect_gt(v, 1.3)
  expect_lt(v, 2.2)
  r <- resolve_threshold("CD45", cd45_case(), floor = 0, prefer = "lowest")
  expect_lt(r$threshold, 2.2)
  b <- best_threshold(cd45_case(), floor = 0, prefer = "lowest")
  expect_lt(b$threshold, 2.2)
})

test_that("prefer = 'lowest' ignores a ripple in the negative mode's tail", {
  set.seed(2)
  x <- c(stats::rnorm(1500, 0.5, 0.15), stats::rnorm(12000, 2.4, 0.15),
         stats::rnorm(6000, 3.1, 0.12))
  expect_gt(density_valley(x, floor = 0, prefer = "lowest"), 0.6)
})

test_that("a CD45 cut through the CD45+ cells is detected against the debris", {
  set.seed(3)
  debris <- stats::rnorm(5000, 0.5, 0.2)
  parent <- c(stats::rnorm(6000, 2.4, 0.15), stats::rnorm(4000, 3.1, 0.12))
  x <- c(debris, parent)
  in_gate <- c(rep(FALSE, 5000), rep(TRUE, 10000))
  ref <- cd45_negative_reference(x, in_gate)
  expect_lt(ref, 1.2)
  # a cut between dim and bright excludes the 60% that are dim
  expect_true(cd45_cut_splits_positive(parent, 2.8, ref))
  # a cut between negative and dim does not
  expect_false(cd45_cut_splits_positive(parent, 1.5, ref))
  # too little debris to trust: no reference
  expect_true(is.na(cd45_negative_reference(x[1:100], in_gate[1:100])))
})

test_that("logicle parameters cover a marker some files lack", {
  # --panel-optional-markers: CD16 is missing from the first file and LAG3 from
  # the second. Both still get parameters, from the files that carry them.
  set.seed(4)
  mk_read <- function(markers) {
    ex <- matrix(stats::rnorm(600 * length(markers), 500, 800), ncol = length(markers))
    colnames(ex) <- markers
    list(exprs = ex, marker_cols = setNames(seq_along(markers), markers))
  }
  reads <- list(a = mk_read(c("CD3", "LAG3")), b = mk_read(c("CD3", "CD16")),
                c = mk_read(c("CD3", "CD16", "LAG3")))
  p <- derive_logicle_pooled(reads, names(reads))
  expect_setequal(names(p), c("CD3", "LAG3", "CD16"))
})

test_that("a population cut the sample could not measure is taken from its peers", {
  thr <- data.frame(
    sample_id = sprintf("S%02d", 1:8), panel = "P1", marker = "SSC-A",
    threshold = c(4.00, 4.02, 3.98, 4.05, 3.97, 4.01, 4.70, 5.30),
    source = c(rep("valley_s6", 5), "otsu", "quantile_fallback", "tailgate"),
    stringsAsFactors = FALSE)
  r <- cyRAVEN:::peer_cut_replacements(thr)
  expect_setequal(r$sample_id, c("S07", "S08"))
  expect_equal(unique(r$threshold), stats::median(thr$threshold[1:6]))
  expect_identical(r$old_source[r$sample_id == "S07"], "quantile_fallback")

  # A tail-rule cut close to the minima stays; so does every measured cut.
  thr$threshold[8] <- 4.10
  expect_identical(cyRAVEN:::peer_cut_replacements(thr)$sample_id, "S07")

  # Too few samples with a minimum: nothing to borrow from.
  thr$source[1:4] <- "tailgate"
  expect_null(cyRAVEN:::peer_cut_replacements(thr))
})

test_that("a borrowed cut reports no resampling uncertainty rather than a fake one", {
  u <- threshold_uncertainty(stats::rnorm(2000), source = "peer_median")
  expect_true(is.na(u$u_combined))
  expect_match(u$basis, "other samples")
})
