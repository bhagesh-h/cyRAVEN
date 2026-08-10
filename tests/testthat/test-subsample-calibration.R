# =============================================================================
# Rare-preserving subsampling, and bead calibration
#
# Both change numbers the previous version reported, which is why both are
# opt-in. The tests below check that they do what they claim AND that the
# default path is untouched.
# =============================================================================

#' A common population and a rare one, well separated in marker space.
rare_cohort <- function(n_common = 9500L, n_rare = 500L, seed = 9) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  X <- rbind(
    cbind(A = stats::rnorm(n_common, 0, 1), B = stats::rnorm(n_common, 0, 1)),
    cbind(A = stats::rnorm(n_rare, 8, 0.4), B = stats::rnorm(n_rare, 8, 0.4)))
  list(X = X, is_rare = rep(c(FALSE, TRUE), c(n_common, n_rare)))
}

test_that("inverse-density weights are larger in sparse regions", {
  d <- rare_cohort()
  w <- density_weights(d$X, k = 20L, n_ref = 1000L)
  expect_length(w, nrow(d$X))
  expect_true(all(w >= 0))
  expect_gt(stats::median(w[d$is_rare]), stats::median(w[!d$is_rare]))
})

test_that("the rare draw over-represents the sparse population", {
  # THE CLAIM THE FEATURE RESTS ON. Under a uniform draw a 5% population stays
  # at 5%, which at these sizes is too few cells to form a cluster.
  d <- rare_cohort()
  idx <- seq_len(nrow(d$X))
  n_take <- 1000L

  set.seed(1); unif <- sort(sample(idx, n_take))
  rare <- draw_subsample_rare(idx, d$X, n_take, seed = 1L)

  f_unif <- mean(d$is_rare[unif])
  f_rare <- mean(d$is_rare[rare$idx])
  expect_gt(f_rare, f_unif)
  expect_length(rare$idx, n_take)
  expect_length(rare$weight, n_take)
})

test_that("the sampling weight reweights the drawn cells back toward the truth", {
  # Without this the embedded set cannot be used for any quantity, because it is
  # a biased sample by design.
  d <- rare_cohort()
  idx <- seq_len(nrow(d$X))
  r <- draw_subsample_rare(idx, d$X, 2000L, seed = 2L)

  truth <- mean(d$is_rare)
  naive <- mean(d$is_rare[r$idx])
  weighted <- sum(r$weight * d$is_rare[r$idx]) / sum(r$weight)
  # The weighted estimate is closer to the cohort's real composition than the
  # unweighted one taken over the same drawn cells.
  expect_lt(abs(weighted - truth), abs(naive - truth))
})

test_that("asking for everything returns everything, unweighted", {
  d <- rare_cohort(n_common = 400L, n_rare = 100L)
  idx <- seq_len(nrow(d$X))
  r <- draw_subsample_rare(idx, d$X, 10000L)
  expect_identical(r$idx, idx)
  expect_true(all(r$weight == 1))
})

test_that("the rare draw borrows and returns the RNG stream", {
  d <- rare_cohort()
  set.seed(123)
  before <- .Random.seed
  expected <- stats::runif(5)
  set.seed(123)
  invisible(draw_subsample_rare(seq_len(nrow(d$X)), d$X, 500L, seed = 7L))
  expect_identical(.Random.seed, before)
  expect_identical(stats::runif(5), expected)
})

# ---- bead calibration -------------------------------------------------------

#' A bead file: `n_pop` discrete populations at known linear intensities.
synth_beads <- function(peaks = c(500, 5000, 50000, 200000), n_each = 3000L,
                        seed = 4) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  v <- unlist(lapply(peaks, function(p) p * 10^stats::rnorm(n_each, 0, 0.03)))
  cbind(FITC = v, OTHER = stats::runif(length(v), 1, 10))
}

test_that("bead populations are recovered in ascending order", {
  peaks <- c(500, 5000, 50000, 200000)
  b <- synth_beads(peaks)
  found <- bead_peaks(b[, "FITC"], n_expected = length(peaks))
  expect_length(found, length(peaks))
  expect_false(is.unsorted(found))
  # Within a few per cent of where they were planted.
  expect_true(all(abs(found - peaks) / peaks < 0.15))
})

test_that("a linear map to the assigned units is fitted and scored", {
  peaks <- c(500, 5000, 50000, 200000)
  b <- synth_beads(peaks)
  # Assigned values that are an exact linear function of the planted peaks, so
  # the fit must recover it.
  assigned <- data.frame(marker = "FITC", p1 = 2 * 500, p2 = 2 * 5000,
                         p3 = 2 * 50000, p4 = 2 * 200000, stringsAsFactors = FALSE)
  cal <- fit_bead_calibration(b, c(FITC = 1L, OTHER = 2L), assigned)
  row <- cal[cal$marker == "FITC", ]
  expect_identical(row$verdict, "calibrated")
  expect_gt(row$r_squared, 0.99)
  expect_equal(row$slope, 2, tolerance = 0.1)
})

test_that("a channel whose fit does not hold is left in instrument units", {
  # A calibration nobody checked is worse than none: it turns an honest
  # arbitrary number into a dishonest absolute one.
  b <- synth_beads(c(500, 5000, 50000, 200000))
  assigned <- data.frame(marker = "FITC", p1 = 1, p2 = 900, p3 = 20, p4 = 5,
                         stringsAsFactors = FALSE)
  cal <- fit_bead_calibration(b, c(FITC = 1L), assigned)
  expect_identical(cal$verdict, "fit too poor to calibrate")

  out <- apply_bead_calibration(b, c(FITC = 1L), cal)
  expect_identical(out$applied, character(0))
  expect_identical(out$exprs, b)
})

test_that("a channel absent from the bead file is reported, not guessed", {
  b <- synth_beads()
  assigned <- data.frame(marker = "NOT_PRESENT", p1 = 1, p2 = 2,
                         stringsAsFactors = FALSE)
  cal <- fit_bead_calibration(b, c(FITC = 1L), assigned)
  expect_identical(cal$verdict, "channel absent from the bead file")
  expect_identical(cal$n_peaks, 0L)
})

test_that("applying the calibration converts only the channels that passed", {
  peaks <- c(500, 5000, 50000, 200000)
  b <- synth_beads(peaks)
  assigned <- data.frame(marker = "FITC", p1 = 1000, p2 = 10000, p3 = 1e5,
                         p4 = 4e5, stringsAsFactors = FALSE)
  cal <- fit_bead_calibration(b, c(FITC = 1L, OTHER = 2L), assigned)

  sample_ex <- cbind(FITC = c(500, 5000), OTHER = c(1, 2))
  out <- apply_bead_calibration(sample_ex, c(FITC = 1L, OTHER = 2L), cal)
  expect_identical(out$applied, "FITC")
  # FITC doubled, OTHER untouched.
  expect_equal(out$exprs[, "FITC"], c(1000, 10000), tolerance = 0.05 * 10000)
  expect_identical(out$exprs[, "OTHER"], sample_ex[, "OTHER"])
})

test_that("a long-form assigned table is accepted as readily as a wide one", {
  peaks <- c(500, 5000, 50000, 200000)
  b <- synth_beads(peaks)
  long <- data.frame(marker = "FITC", peak = 1:4, value = 2 * peaks,
                     stringsAsFactors = FALSE)
  cal <- fit_bead_calibration(b, c(FITC = 1L), long)
  expect_identical(cal$verdict, "calibrated")
  expect_equal(cal$slope, 2, tolerance = 0.1)
})

test_that("an assigned table with no marker column is refused with a reason", {
  expect_error(fit_bead_calibration(synth_beads(), c(FITC = 1L),
                                    data.frame(p1 = 1, p2 = 2)),
               "marker")
})
