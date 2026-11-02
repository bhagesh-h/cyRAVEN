# =============================================================================
# Acquisition-time quality control
#
# The claim: a file whose event rate or signal moved part-way through the run is
# named, a file that ran cleanly is not, and nothing is removed unless asked.
# =============================================================================

#' A synthetic acquisition. `clog` empties a window of events; `drift` shifts one
#' channel over a window. Both leave the rest of the file untouched.
synth_run <- function(n = 20000L, clog = NULL, drift = NULL, seed = 5) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  time <- sort(stats::runif(n, 0, 1000))
  mat <- cbind(CD3 = stats::rnorm(n, 2, 0.5), CD4 = stats::rnorm(n, 1, 0.5))
  keep <- rep(TRUE, n)
  # A clog is a stretch of time during which almost nothing was recorded.
  if (!is.null(clog)) keep[time > clog[1] & time < clog[2]] <- stats::runif(
    sum(time > clog[1] & time < clog[2])) < 0.05
  # A drift moves one channel over a stretch without touching the rate.
  if (!is.null(drift)) {
    w <- time > drift[1] & time < drift[2]
    mat[w, "CD3"] <- mat[w, "CD3"] + drift[3]
  }
  list(time = time[keep], mat = mat[keep, , drop = FALSE])
}

#' A read-shaped fixture. Marker values are raised to instrument scale, as a
#' real file stores them, because acquisition_qc_one() judges the stored values.
fake_rd <- function(n = 20000L, with_time = TRUE, seed = 5) {
  s <- synth_run(n = n, seed = seed)
  ex <- cbind(10^(s$mat + 2), Time = s$time)
  nm <- c("CD3", "CD4", if (with_time) "Time" else "NotTime")
  colnames(ex) <- nm
  list(exprs = ex, marker_cols = c(CD3 = 1L, CD4 = 2L), channel_names = nm,
       n_events = nrow(ex), sample_id = "S1", file = "S1.fcs")
}

test_that("a clean acquisition is not flagged", {
  s <- synth_run()
  r <- detect_time_anomalies(s$time, s$mat)
  expect_s3_class(r$bins, "data.frame")
  expect_false(any(r$bins$flagged))
  expect_false(any(r$flagged))
})

test_that("a drop in event rate is flagged", {
  # THE CANONICAL CASE. A partial clog empties a stretch of the acquisition; the
  # signal is untouched, so only the rate test can see it.
  s <- synth_run(clog = c(400, 550))
  r <- detect_time_anomalies(s$time, s$mat)
  expect_true(any(r$bins$flagged))
  hit <- r$bins[r$bins$flagged, ]
  # The flagged interval overlaps the window the events were removed from.
  expect_true(any(hit$time_to > 400 & hit$time_from < 550))
  expect_gt(max(r$bins$rate_z), 5)
})

test_that("a shift in one channel is flagged even when the rate is steady", {
  # The complementary failure: the instrument kept delivering cells, but not the
  # same measurement. A rate test alone cannot see this.
  s <- synth_run(drift = c(300, 500, 3))
  r <- detect_time_anomalies(s$time, s$mat)
  expect_true(any(r$bins$flagged))
  expect_gt(max(r$bins$worst_channel_z), 5)
  # And it names the channel responsible rather than only the interval.
  expect_identical(unique(r$bins$worst_channel[r$bins$flagged]), "CD3")
})

test_that("a file with no usable Time channel says so rather than passing", {
  # Silence would be indistinguishable from a clean acquisition.
  expect_match(detect_time_anomalies(rep(1, 5000), NULL)$reason, "does not advance")
  expect_match(detect_time_anomalies(1:10, NULL)$reason, "fewer than 200")

  rd <- fake_rd(with_time = FALSE)
  out <- acquisition_qc_one(rd)
  expect_identical(out$summary$verdict, "no Time channel")
  expect_null(out$bins)
})

test_that("the Time column is found regardless of channel-name case", {
  rd <- fake_rd()
  expect_identical(find_time_column(rd), 3L)
  rd$channel_names <- c("CD3", "CD4", "TIME")
  expect_identical(find_time_column(rd), 3L)
  rd$channel_names <- c("CD3", "CD4", "Time.Stamp")
  expect_identical(find_time_column(rd), 3L)
  rd$channel_names <- c("CD3", "CD4", "SSC-A")
  expect_true(is.na(find_time_column(rd)))
})

test_that("the per-sample verdict scales with how much of the file is affected", {
  clean <- acquisition_qc_one(fake_rd())
  expect_identical(clean$summary$verdict, "stable")
  expect_identical(clean$summary$n_bins_flagged, 0L)

  s <- synth_run(clog = c(300, 700))
  ex <- cbind(10^(s$mat + 2), Time = s$time)
  colnames(ex) <- c("CD3", "CD4", "Time")
  rd <- list(exprs = ex, marker_cols = c(CD3 = 1L, CD4 = 2L),
             channel_names = colnames(ex), n_events = nrow(ex),
             sample_id = "S2", file = "S2.fcs")
  bad <- acquisition_qc_one(rd)
  expect_true(grepl("^unstable", bad$summary$verdict))
  expect_gt(bad$summary$n_bins_flagged, 0L)
  expect_gt(bad$summary$pct_events_flagged, 0)
})

test_that("the cohort driver returns one row per sample and per-event flags", {
  reads <- list(S1 = fake_rd(seed = 1), S2 = fake_rd(seed = 2))
  reads$S1$sample_id <- "S1"; reads$S2$sample_id <- "S2"
  out <- run_acquisition_qc(reads)
  expect_identical(nrow(out$summary), 2L)
  expect_setequal(out$summary$sample_id, c("S1", "S2"))
  for (s in names(out$flagged))
    expect_length(out$flagged[[s]], nrow(reads[[s]]$exprs))
  expect_null(run_acquisition_qc(NULL))
  expect_null(run_acquisition_qc(list()))
})

test_that("the impact table states what excluding the events would cost", {
  # This is what makes the flag actionable rather than alarming.
  withr::local_seed(11)
  n <- 4000L
  tmat <- cbind(CD3 = c(stats::rnorm(n / 2, 4, 0.5), stats::rnorm(n / 2, 0, 0.5)))
  thr <- c(CD3 = 2)
  spec <- list(`T cells` = list(CD3 = "above"))
  parent <- rep(TRUE, n)

  # Drop events that are all CD3-positive, so the frequency must fall.
  keep <- rep(TRUE, n); keep[1:500] <- FALSE
  d <- frequency_delta_if_cleaned(tmat, thr, parent, spec, keep)
  expect_s3_class(d, "data.frame")
  expect_identical(d$population, "T cells")
  expect_lt(d$pct_if_cleaned, d$pct_reported)
  expect_equal(d$pct_delta_if_cleaned, d$pct_if_cleaned - d$pct_reported,
               tolerance = 1e-6)

  # Nothing flagged means nothing to report, not a table of zeros.
  expect_null(frequency_delta_if_cleaned(tmat, thr, parent, spec, rep(TRUE, n)))
})

test_that("the detector consumes no randomness", {
  # Same discipline as the uncertainty code: run_cyraven() seeds once and the
  # UMAP cell selection draws from that stream.
  s <- synth_run(clog = c(400, 550))
  set.seed(123)
  before <- .Random.seed
  expected <- stats::runif(5)
  set.seed(123)
  invisible(detect_time_anomalies(s$time, s$mat))
  expect_identical(.Random.seed, before)
  expect_identical(stats::runif(5), expected)
})

test_that("a compensated channel centred on zero is judged without a kink", {
  # log10(max(x, 1)) pinned every value at or below 1 to 0, where the median
  # of a mostly negative compensated channel sits; the stored values have no
  # such kink.
  set.seed(5)
  n <- 40000L; time <- sort(stats::runif(n, 0, 1000))
  x <- c(stats::rnorm(n * 0.9, 0, 100), stats::rnorm(n * 0.1, 3000, 600))[sample(n)]
  rd <- function(v) list(exprs = cbind(M = v, Time = time), marker_cols = c(M = 1L),
                         channel_names = c("M", "Time"), n_events = n,
                         sample_id = "S", file = "S.fcs")
  expect_identical(acquisition_qc_one(rd(x))$summary$n_bins_flagged, 0L)
  w <- time > 600 & time < 700
  b <- acquisition_qc_one(rd(x + 80 * w))$bins
  hit <- b[b$flagged, ]
  expect_gt(nrow(hit), 0L)
  expect_true(all(hit$time_to > 600 & hit$time_from < 700))
})

test_that("PeacoQC removes a disturbed stretch of the acquisition", {
  skip_if_not_installed("PeacoQC")
  # A stable acquisition with a 10% stretch where every channel shifts up.
  withr::local_seed(21)
  n <- 30000L
  time <- seq(0, 300, length.out = n)
  mk <- function(m) c(stats::rnorm(n * 0.7, 0, 150), stats::rnorm(n * 0.3, m, m / 5))[sample(n)]
  ex <- cbind(`FSC-A` = stats::rnorm(n, 80000, 8000), `SSC-A` = stats::rnorm(n, 30000, 5000),
              CD3 = mk(5000), CD4 = mk(3000), CD8 = mk(4000), Time = time)
  bad <- time > 150 & time < 180
  ex[bad, c("CD3", "CD4", "CD8")] <- ex[bad, c("CD3", "CD4", "CD8")] + 3000
  rd <- list(exprs = ex, marker_cols = c(CD3 = 3L, CD4 = 4L, CD8 = 5L),
             scatter_cols = c(`FSC-A` = 1L, `SSC-A` = 2L),
             channel_names = colnames(ex), sample_id = "S1", file = "S1.fcs",
             n_events = n)
  r <- suppressMessages(peacoqc_one(rd))
  expect_length(r$good, n)
  expect_gt(mean(!r$good[bad]), 0.5)          # most of the disturbance removed
  expect_lt(mean(!r$good[!bad]), 0.2)         # most of the rest kept
  # PeacoQC skips its isolation tree below 150 time bins with a warning that
  # is suppressed; the skip has to be said in the summary instead.
  if (is.na(r$pct_it)) expect_match(r$reason, "isolation tree not run")
})
