# =============================================================================
# Three fixes found by running cyRAVEN against cyCONDOR on a 24-sample cohort
#
# 1. Values above a channel's declared range ($PnR) are capped on reading under
#    --keep-margin-events (removed by default), so a stray float cannot set the
#    top of the logicle scale.
# 2. No density-derived cut may sit at or below zero signal, so the bump that
#    compensated events form below zero cannot be mistaken for a negative mode.
# 3. Equalised embedding draws stop at the lower quartile of the per-sample
#    counts, so one small sample cannot shrink every other sample to its size.
# =============================================================================

test_that("values above $PnR are capped at it under --keep-margin-events", {
  ex <- cbind(a = c(1, 50, 200, 9e7), b = c(1, 2, 3, 4))
  colnames(ex) <- c("FL1-A", "FL2-A")
  kw <- list(`$P1R` = "262144", `$P2R` = "1024", `$DATATYPE` = "F")
  sc <- linearise_channels(ex, kw)$scale
  out <- clamp_to_declared_range(ex, sc)
  expect_equal(out[, 1], c(1, 50, 200, 262144))
  expect_equal(out[, 2], ex[, 2])
})

test_that("a channel with many values above $PnR is still capped", {
  # More than 1% over used to leave the channel alone, and one such file set the
  # logicle scale for a whole panel. It is capped now, and named in a warning.
  ex <- cbind(c(rep(100, 90), rep(5e7, 10)))
  colnames(ex) <- "FL1-A"
  sc <- linearise_channels(ex, list(`$P1R` = "262144", `$DATATYPE` = "F"))$scale
  out <- clamp_to_declared_range(ex, sc)
  expect_equal(max(out), 262144)
})

test_that("a channel without a declared range is left alone", {
  ex <- cbind(c(1, 9e9))
  colnames(ex) <- "Time"
  sc <- linearise_channels(ex, list())$scale
  expect_equal(clamp_to_declared_range(ex, sc), ex)
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

test_that("auto scaling is none on a shared scale and the percentile range otherwise", {
  expect_identical(umap_scale_method(list(transform = "logicle")), "none")
  expect_identical(umap_scale_method(list(per_marker_cofactor = TRUE)), "none")
  expect_identical(umap_scale_method(list(transform = "arcsinh")), "range")
  expect_identical(umap_scale_method(list()), "range")
  expect_identical(umap_scale_method(list(transform = "logicle", umap_scale = "zscore")), "zscore")
  expect_identical(umap_scale_method(list(umap_scale = "robust")), "robust")
})

test_that("range scaling does not inflate a mostly negative marker's background", {
  # Median/MAD divides a 98%-negative marker by the width of its negatives, so
  # its background ends up as wide as a balanced marker's whole range.
  set.seed(3)
  rare <- c(rnorm(9800, 0, 0.1), rnorm(200, 4, 0.3))
  half <- c(rnorm(5000, 0, 0.1), rnorm(5000, 4, 0.3))
  m <- cbind(rare = rare, half = half)
  neg_sd <- function(z) c(sd(z[1:9800, 1]), sd(z[1:5000, 2]))
  r <- neg_sd(scale_features(m, "robust")$mat)
  g <- neg_sd(scale_features(m, "range")$mat)
  expect_gt(r[1] / r[2], 10)
  expect_lt(g[1] / g[2], 1.5)
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
  # Six of eight samples found a density minimum: the marker is bimodal here,
  # so both the fallback and the tail cut take the median of the minima.
  thr <- data.frame(
    sample_id = sprintf("S%02d", 1:8), panel = "P1", marker = "CD19",
    threshold = c(2.40, 2.42, 2.38, 2.45, 2.37, 2.41, 1.34, 3.24),
    source = c(rep("valley_s6", 5), "otsu", "quantile_fallback", "tailgate"),
    stringsAsFactors = FALSE)
  r <- cyRAVEN:::peer_cut_replacements(thr)
  expect_setequal(r$sample_id, c("S07", "S08"))
  expect_equal(unique(r$threshold), stats::median(thr$threshold[1:6]))
  expect_identical(r$old_source[r$sample_id == "S08"], "tailgate")

  # Without a majority of minima, as for a rare marker like Vd1, the tail rule
  # stands unless a cut is far from the others; a fallback is still replaced.
  thr$source[1:4] <- "tailgate"
  thr$threshold[8] <- 2.50
  expect_identical(cyRAVEN:::peer_cut_replacements(thr)$sample_id, "S07")
  thr$threshold[8] <- 0.80
  expect_setequal(cyRAVEN:::peer_cut_replacements(thr)$sample_id, c("S07", "S08"))

  # Measured cuts, and cuts from a control or an override, are never touched.
  thr$source <- c(rep("valley_s6", 6), "control_q995", "manual")
  expect_null(cyRAVEN:::peer_cut_replacements(thr))
})

test_that("a borrowed cut reports no resampling uncertainty rather than a fake one", {
  u <- threshold_uncertainty(stats::rnorm(2000), source = "peer_median")
  expect_true(is.na(u$u_combined))
  expect_match(u$basis, "other samples")
})

test_that("the total cap keeps every sample's name", {
  # Callers look each sample's draw up by name. The cap used to rebuild the
  # vector through pmax(1L, ...), which drops names, so a capped plan made every
  # lookup fail and explore mode stopped on the 35-sample demonstration cohort.
  n <- stats::setNames(c(17289L, 4530L, 6765L, 0L, 25515L), paste0("s", 1:5))
  take <- suppressMessages(plan_subsample(n, cap = 6000L, total_cap = 10000L))
  expect_identical(names(take), c("s1", "s2", "s3", "s5"))
  expect_lte(sum(take), 10000L)
  expect_true(is.na(take["s4"]))
})
