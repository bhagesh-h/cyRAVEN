# =============================================================================
# Reading: log-amplified channels, margin events, and cuts on tied values
#
# Three faults found by auditing the reader against FCS 3.1 and the published
# QC tools. Log-amplified integer channels were analysed as channel numbers;
# every value above $PnR was capped whatever scale $PnR described, which on the
# demo cohort piled 6% of CD45 events onto one value in the middle of the data;
# and a quantile cut sitting on that pile selected nothing under a strict `>`.
# =============================================================================

test_that("integer log channels are linearised as flowCore does", {
  ex <- cbind(`FL1-H` = c(0, 256, 512, 1023), `FSC-H` = c(0, 100, 500, 1023))
  kw <- list(`$DATATYPE` = "I", `$P1E` = "4,0", `$P1R` = "1024",
             `$P2E` = "0,0", `$P2R` = "1024")
  l <- linearise_channels(ex, kw)
  # 10^(4 * x / 1024), f2 = 0 read as 1
  expect_equal(l$exprs[, 1], 10^(4 * c(0, 256, 512, 1023) / 1024))
  expect_equal(l$exprs[, 2], ex[, 2])
  expect_equal(l$scale$range, c(1e4, 1024))
  expect_identical(l$scale$stored[1], "log, linearised")
})

test_that("a log $PnE on floating-point data is not applied, and its $PnR not used", {
  # flowCore's GvHD frames written back out: linear floats to 10,000 under the
  # FACSCalibur keywords. FCS 3.1 forbids log storage of floats.
  ex <- cbind(`FL1-H` = c(1, 500, 9000))
  kw <- list(`$DATATYPE` = "F", `$P1E` = "4,0", `$P1R` = "1024")
  l <- linearise_channels(ex, kw)
  expect_equal(l$exprs, ex)
  expect_true(is.na(l$scale$range[1]))
  expect_true(is.na(l$scale$margin_at[1]))
  # so nothing is capped in the middle of the data
  expect_equal(clamp_to_declared_range(ex, l$scale), ex)
  expect_identical(margin_events(ex, l$scale, 1L)$n, 0L)
})

test_that("margin events are the ones at the top of a used channel", {
  ex <- cbind(`CD3-A` = c(10, 262143, 50, 7e8), `W-A` = c(262143, 1, 1, 1),
              `SSC-A` = c(5, 5, 262143, 5))
  kw <- list(`$DATATYPE` = "F", `$P1R` = "262144", `$P2R` = "262144",
             `$P3R` = "262144")
  sc <- linearise_channels(ex, kw)$scale
  m <- margin_events(ex, sc, used = c(1L, 3L))
  # events 2 and 4 in CD3, event 3 in SSC; the saturated width channel is not
  # in use, so event 1 is a cell
  expect_identical(which(m$event), c(2L, 3L, 4L))
  expect_identical(m$n, 3L)
  expect_identical(unname(m$by_channel[["CD3-A"]]), 2L)
})

test_that("lower-limit outliers are counted by flowAI's rule", {
  # flow_margin_check(neg_values = 1): among the negative values, those at or
  # below median - 3.5 * mad / 0.6745, with R's mad().
  withr::local_seed(21)
  x <- c(rnorm(5000, 0, 100), -5000, -4000)
  y <- rnorm(5002, 500, 50)
  ex <- cbind(`PE-A` = x, `FITC-A` = y)
  neg <- x[x < 0]
  lim <- -3.5 * stats::mad(neg) / 0.6745 + stats::median(neg)
  by <- lower_limit_outliers(ex, 1:2)
  expect_identical(names(by), "PE-A")
  expect_identical(unname(by[["PE-A"]]), sum(x <= lim))
  expect_gte(by[["PE-A"]], 2L)
  expect_identical(attr(by, "events"), sum(x <= lim))
  # a channel with no negative values has nothing below its lower limit
  expect_length(lower_limit_outliers(ex, 2L), 0L)
  # an event below the limit in two channels is one event
  ex2 <- cbind(ex, `APC-A` = x)
  by2 <- lower_limit_outliers(ex2, c(1L, 3L))
  expect_identical(sum(by2), 2L * sum(x <= lim))
  expect_identical(attr(by2, "events"), sum(x <= lim))
})

test_that("a quantile cut on a pile of tied values is moved off it", {
  x <- c(seq(0, 1, length.out = 85), rep(5, 15))
  q <- as.numeric(quantile(x, 0.9))
  expect_equal(q, 5)
  # a strict `>` at the pile selects nothing
  expect_equal(mean(x > q), 0)
  cut <- untie_cut(x, q, target_above = 0.10)
  # the pile is now above the cut: 15% positive, the nearer of 0% and 15%
  expect_equal(mean(x > cut), 0.15)
  # a cut touching no value is returned unchanged
  expect_identical(untie_cut(x, 0.123456), 0.123456)
})

test_that("a unimodal marker without a control gets the tail gate, not the 90th percentile", {
  set.seed(4)
  x <- c(rnorm(20000, 0, 0.3), rnorm(40, 2.2, 0.2))
  r <- resolve_threshold("CD25", x)
  expect_identical(r$source, "tailgate")
  expect_true(r$needs_review)
  # a constant 10% would be the fabricated answer
  expect_lt(mean(x > r$threshold), 0.02)
  # with a control present the control still wins, as before
  rc <- resolve_threshold("CD25", x, control_x = rnorm(5000, 0, 0.3))
  expect_identical(rc$source, "control_q995")
  expect_identical(rc$derived_source, "tailgate")
})

test_that("--fmo-sets-cut lets a true FMO lower a valley, and nothing else", {
  set.seed(7)
  x <- c(rnorm(15000, 0, 0.25), rnorm(6000, 1.6, 0.25))
  fmo <- rnorm(5000, 0, 0.25)
  ctrl <- as.numeric(quantile(fmo, 0.995))
  v <- resolve_threshold("CD25", x, control_x = fmo, control_kind = "fmo_q995")
  expect_identical(v$source, "valley")
  f <- resolve_threshold("CD25", x, control_x = fmo, control_kind = "fmo_q995",
                         fmo_sets_cut = TRUE)
  expect_identical(f$source, "fmo_q995_valley_lowered")
  expect_lt(f$threshold, v$threshold)
  expect_equal(f$derived_threshold, v$threshold)
  # a minus-multiple tube is permissive and never argues a cut down
  m <- resolve_threshold("CD25", x, control_x = fmo, control_kind = "fmm_q995",
                         fmo_sets_cut = TRUE)
  expect_identical(m$source, "valley")
})

test_that("a negative mode carrying background is not read as positive", {
  # Vd2 on CD45+ whole blood at its own per-marker cofactor: every cell carries
  # background (granulocyte autofluorescence, spillover spread), so the
  # dominant NEGATIVE mode sits many of its SDs above zero signal. Cutting below
  # it, as an earlier side test did, called 90-99% of CD45+ cells Vd2+.
  # A lower-background shoulder merging into the mode, with no trough.
  set.seed(21)
  x <- c(rnorm(18000, 2.8, 0.25), rnorm(2500, 2.2, 0.45))
  r <- resolve_threshold("Vd2", x, floor = 0)
  expect_identical(r$source, "tailgate")
  expect_lt(mean(x > r$threshold), 0.01)
  expect_gt(r$threshold, 2.8)
})

test_that("a dye with no dead cells gets a tail gate keeping everyone, not 90%", {
  set.seed(31)
  x <- rnorm(20000, 2.7, 0.3)
  r <- resolve_threshold("L/D", x, floor = 2.7)
  expect_identical(r$source, "tailgate")
  expect_gt(mean(x < r$threshold), 0.999)
})

test_that("--adaptive-gates falls back to the same tail gate as the default path", {
  # one mode far above zero, a shoulder down to zero and no trough anywhere:
  # the cut goes above the mode on both paths
  set.seed(23)
  x <- c(rnorm(16000, 4, 0.4), runif(4000, 0, 4))
  tg <- tail_gate_cut(x, floor = 0)
  expect_gt(tg, 4)
  # min_quality = Inf rejects every scored candidate, reaching the last resort
  b <- best_threshold(x, floor = 0, min_quality = Inf)
  expect_identical(b$source, "tailgate")
  expect_equal(b$threshold, tg)
  # no dead cells: the empty tail is accepted on both paths
  set.seed(31)
  d <- rnorm(20000, 2.7, 0.3)
  bd <- best_threshold(d, floor = 2.7, min_quality = Inf)
  expect_identical(bd$source, "tailgate")
  expect_gt(mean(d < bd$threshold), 0.999)
})

test_that("a single mode yields no valley from its sparse tails or the histogram's ends", {
  # The end bins used to be left unsmoothed, and a dip of 30% between tail bins
  # of a few events passed as a separation. Either cut the mode from below and
  # called nearly every cell positive.
  v <- vapply(1:30, function(s) {
    set.seed(s); density_valley(rnorm(10000, 2, 0.4), floor = 0)
  }, numeric(1))
  expect_true(all(is.na(v)))
  # the Vd2 case: few cells, every value above zero, no second population
  set.seed(61)
  expect_true(is.na(density_valley(rlnorm(3000, log(0.9), 0.45), floor = 0)))
  # a real minority population is still separated
  set.seed(62)
  x <- c(rnorm(2850, 0.5, 0.35), rnorm(150, 2.5, 0.35))
  cut <- density_valley(x, floor = 0)
  expect_true(cut > 1.2 && cut < 2.2)
})

test_that("a valley needs slopes beyond counting noise, not a noisy dip", {
  # single modes of 3,000 events: the uncorrected depth test found a valley in
  # about half of them
  v <- vapply(1:30, function(s) {
    set.seed(s); density_valley(rnorm(3000, 2, 0.4), floor = 0)
  }, numeric(1))
  expect_true(all(is.na(v)))
  # a flat shoulder of a few events per bin rising into a mode has no dip
  f <- vapply(1:30, function(s) {
    set.seed(s); density_valley(c(rnorm(18000, 2.8, 0.25), runif(1500, 0.5, 2.6)), floor = 0)
  }, numeric(1))
  expect_lt(sum(is.finite(f)), 4)
  # with the test off the same shoulders are split
  g <- vapply(1:30, function(s) {
    set.seed(s); density_valley(c(rnorm(18000, 2.8, 0.25), runif(1500, 0.5, 2.6)),
                                floor = 0, noise_alpha = NA)
  }, numeric(1))
  expect_gt(sum(is.finite(g)), 10)
})

test_that("the tail gate's spread is the mode's own, not a second negative hump's", {
  # a normal mode: the flank width gives its SD, kernel variance removed
  set.seed(71)
  x <- rnorm(20000, 0, 0.4)
  d <- stats::density(x, n = 1024, adjust = 2)
  expect_equal(flank_sd(d, which.max(d$y)), 0.4, tolerance = 0.05)
  expect_equal(tail_threshold(x), 3.5 * 0.4, tolerance = 0.1)
  # two negative humps (monocyte and lymphocyte background) and a 2.4%
  # positive plateau: the mirrored-half MAD took in the far hump and cut past
  # every cell
  set.seed(5)
  x <- c(rnorm(9000, -0.6, 0.35), rnorm(11000, 0.6, 0.35), runif(500, 2.5, 4))
  p <- mean(x > tail_threshold(x))
  expect_gt(p, 0.015)
  expect_lt(p, 0.035)
})

test_that("a declared lower tail cuts below a positive mode, and only when declared", {
  # CD127 within CD4 T cells: most are CD127+, a CD127-low shoulder, no trough
  set.seed(81)
  x <- c(rnorm(18000, 3, 0.3), runif(2000, 0, 2.4))
  up <- resolve_threshold("CD127", x, floor = 0)
  lo <- resolve_threshold("CD127", x, floor = 0, tail = "lower")
  expect_identical(lo$source, "tailgate")
  expect_lt(mean(x > up$threshold), 0.01)
  expect_gt(mean(x > lo$threshold), 0.85)
  expect_identical(lo$tail, "lower")
  # read from the configuration
  expect_identical(cfg_tail(list(tail = "lower")), "lower")
  expect_identical(cfg_tail(list(threshold = 1.2)), "upper")
  expect_identical(cfg_tail(1.2), "upper")
  expect_null(cfg_threshold(list(tail = "lower")))
  expect_identical(cfg_threshold(list(threshold = 1.2)), 1.2)
})

test_that("no derived cut sits below zero signal, the tail gate's and the quantile's included", {
  # an over-compensated channel: the negatives end below zero
  set.seed(91)
  x <- c(rnorm(19800, -2, 0.3), rnorm(200, 1.5, 0.4))
  r <- resolve_threshold("CD14", x, floor = 0)
  expect_identical(r$source, "tailgate")
  expect_gte(r$threshold, 0)
  # the last resort is floored too
  q <- resolve_threshold("CD14", rep(c(-3, -2.5, -2), length.out = 3000), floor = 0)
  expect_gte(q$threshold, 0)
})

test_that("the live mode is the tallest one, not a bump below it", {
  set.seed(41)
  bump <- rnorm(1500, 0.1, 0.15)
  expect_lt(abs(live_mode(c(rnorm(9000, 2.7, 0.3), bump)) - 2.7), 0.3)
  # a viability cut is never taken between the bump and the live cells
  x <- c(rnorm(9000, 2.7, 0.3), bump, rnorm(300, 5.5, 0.3))
  v <- density_valley(x, floor = max(0, live_mode(x)))
  expect_true(is.na(v) || v > 2.7)
})

test_that("linearised log data get a tail gate above their negative mode", {
  # Linearised log data start at 1 unit, so the negative mode sits well above
  # zero signal; the tail gate still goes above it rather than falling back.
  set.seed(51)
  neg <- asinh((10^(4 * pmax(rnorm(9000, 100, 40), 0) / 1024)) / 2.3)
  pos <- asinh((10^(4 * rnorm(500, 700, 60) / 1024)) / 2.3)
  x <- c(neg, pos)
  r <- resolve_threshold("CD14", x, floor = 0)
  expect_false(identical(r$source, "quantile_fallback"))
  expect_lt(mean(x > r$threshold), 0.2)
})

test_that("--adaptive-gates honours a declared lower tail and floors its quantile", {
  # A dominant positive mode with a dimmer shoulder and no trough: undeclared,
  # the last resort goes above the mode; declared lower, below it, and no
  # upper-tail candidate is offered at all.
  withr::local_seed(41)
  x <- c(rnorm(9000, 3, 0.3), rnorm(1000, 2.2, 0.5))
  expect_gt(best_threshold(x, floor = 0, min_quality = Inf)$threshold, 3)
  expect_lt(best_threshold(x, floor = 0, min_quality = Inf, tail = "lower")$threshold, 3)
  b <- best_threshold(x, floor = 0, tail = "lower")
  expect_false(identical(b$source, "tailgate_upper"))
  # resolve_threshold() carries the declaration into the adaptive path
  r <- resolve_threshold("M", x, adaptive = TRUE, floor = 0, tail = "lower")
  expect_lt(r$threshold, 3)
  # no tail gate can be placed here (lower side, under 40% above it), and the
  # 90th percentile, 1.28, sits below the floor: the quantile is floored
  y <- rnorm(5000, 0, 1)
  q <- best_threshold(y, floor = 1.5, min_quality = Inf, tail = "lower")
  expect_identical(q$source, "quantile_fallback")
  expect_gte(q$threshold, 1.5)
})
