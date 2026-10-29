# =============================================================================
# Spillover spreading
#
# The claim: a channel whose negative population widens when another channel is
# bright is identified, a channel that does not is left alone, and a marker that
# fails to resolve a cut for that reason is told apart from one that fails for
# any other.
# =============================================================================

#' Two markers plus a receiver. `spread` widens the receiver's NEGATIVE
#' population among cells positive for the source, without moving its centre,
#' which is what compensated spillover does. `shift` moves the centre instead,
#' which is what co-expression does and what must NOT be reported as spreading.
synth_spread <- function(n = 6000L, spread = 0, shift = 0, seed = 7) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  src_pos <- rep(c(FALSE, TRUE), each = n / 2)
  SRC <- ifelse(src_pos, stats::rnorm(n, 5, 0.5), stats::rnorm(n, 0, 0.5))
  # Receiver: everyone negative, but wider (and optionally displaced) where the
  # source is bright.
  sd_rec <- ifelse(src_pos, 0.5 + spread, 0.5)
  mu_rec <- ifelse(src_pos, shift, 0)
  REC <- stats::rnorm(n, mu_rec, sd_rec)
  QUIET <- stats::rnorm(n, 0, 0.5)
  cbind(SRC = SRC, REC = REC, QUIET = QUIET)
}

THR <- c(SRC = 2.5, REC = 3, QUIET = 3)

test_that("a widened negative population is detected and quantified", {
  m <- synth_spread(spread = 1.0)
  out <- spreading_pairs(m, THR, min_cells = 100L)
  expect_s3_class(out, "data.frame")
  row <- out[out$source == "SRC" & out$receiver == "REC", ]
  expect_identical(nrow(row), 1L)
  expect_gt(row$spreading, 0.5)
  expect_gt(row$spreading_ratio, 1.5)
})

test_that("a channel that spreads nothing reports no spreading", {
  m <- synth_spread(spread = 0)
  out <- spreading_pairs(m, THR, min_cells = 100L)
  row <- out[out$source == "SRC" & out$receiver == "REC", ]
  expect_lt(abs(row$spreading), 0.1)
  expect_lt(abs(row$spreading_ratio - 1), 0.15)
  # And the unrelated channel is unaffected in both directions.
  q <- out[out$receiver == "QUIET", ]
  expect_true(all(abs(q$spreading_ratio - 1) < 0.2))
})

test_that("co-expression is not reported as spreading", {
  # THE DISTINCTION THE MEASURE RESTS ON. A real biological signal in the
  # receiver moves its positive cells. Restricting to the receiver's own
  # negatives is what keeps that out of the number.
  m <- synth_spread(spread = 0, shift = 2.0)
  out <- spreading_pairs(m, THR, min_cells = 100L)
  row <- out[out$source == "SRC" & out$receiver == "REC", ]
  # The centre moved a great deal; the width of the negatives did not.
  expect_lt(abs(row$spreading_ratio - 1), 0.25)
})

test_that("groups too small to describe a distribution are skipped", {
  m <- synth_spread(spread = 1)
  expect_null(spreading_pairs(m, THR, min_cells = 10000L))
  expect_null(spreading_pairs(NULL, THR))
  # A marker with no finite threshold cannot define a negative population.
  expect_null(spreading_pairs(m, c(SRC = NA_real_, REC = NA_real_,
                                   QUIET = NA_real_), min_cells = 100L))
})

test_that("the cohort report ranks receivers and names the worst source", {
  mk <- function(seed) list(tmat = synth_spread(spread = 1.2, seed = seed),
                            thresholds = THR)
  pops <- list(S1 = mk(1), S2 = mk(2))
  gates <- list(S1 = list(masks = list(cd45_pos = rep(TRUE, 6000L))),
                S2 = list(masks = list(cd45_pos = rep(TRUE, 6000L))))
  out <- run_spreading_report(pops, gates, panel_of = c(S1 = "P1", S2 = "P1"))

  expect_named(out, c("pairs", "receivers"))
  rec <- out$receivers
  # REC takes the most spreading, and SRC is what is giving it.
  expect_identical(rec$receiver[1], "REC")
  expect_identical(rec$worst_source[1], "SRC")
  expect_gt(rec$worst_ratio[1], 1.25)
  expect_true(any(out$pairs$substantial))
})

test_that("an unresolved marker that is heavily spread is told apart from one that is not", {
  # The actionable finding. Without the join, the run can say a cut was a
  # quantile fallback and cannot say why.
  mk <- function(seed) list(tmat = synth_spread(spread = 1.2, seed = seed),
                            thresholds = THR)
  pops <- list(S1 = mk(1), S2 = mk(2))
  gates <- list(S1 = list(masks = list(cd45_pos = rep(TRUE, 6000L))),
                S2 = list(masks = list(cd45_pos = rep(TRUE, 6000L))))
  thr_all <- data.frame(
    sample_id = rep(c("S1", "S2"), each = 3), panel = "P1",
    marker = c("SRC", "REC", "QUIET"),
    # REC never resolves; QUIET always does.
    source = c("valley", "quantile_fallback", "valley",
               "valley", "quantile_fallback", "valley"),
    stringsAsFactors = FALSE)

  rec <- run_spreading_report(pops, gates, thr_all,
                              panel_of = c(S1 = "P1", S2 = "P1"))$receivers
  rownames(rec) <- rec$receiver
  expect_match(rec["REC", "verdict"], "panel design problem")
  expect_equal(rec["REC", "fallback_rate"], 1)
  expect_identical(rec["QUIET", "verdict"], "no substantial spreading detected")
})

test_that("the report consumes no randomness", {
  m <- synth_spread(spread = 1)
  set.seed(123)
  before <- .Random.seed
  expected <- stats::runif(5)
  set.seed(123)
  invisible(spreading_pairs(m, THR, min_cells = 100L))
  expect_identical(.Random.seed, before)
  expect_identical(stats::runif(5), expected)
})

test_that("the spreading matrix recovers a planted coefficient from a single stain", {
  # A PE single stain: positive and negative beads. Detector B receives
  # spreading whose variance grows with the PE signal, sigma_B^2 = 50^2 +
  # k^2 * F, so SS = k.
  withr::local_seed(12)
  n <- 20000; k <- 3
  pos <- c(rep(TRUE, n / 2), rep(FALSE, n / 2))
  pe <- ifelse(pos, rnorm(n, 40000, 2000), rnorm(n, 0, 60))
  b  <- rnorm(n, 0, sqrt(50^2 + ifelse(pos, k^2 * 40000, 0)))
  c3 <- rnorm(n, 0, 50)
  rd <- list(exprs = cbind(pe, b, c3), marker_cols = c(PE = 1L, B = 2L, C = 3L),
             channel_names = c("PE.A", "B.A", "C.A"))
  ssm <- spillover_spreading_matrix(list(pe_stain = rd))
  expect_identical(unique(ssm$source_detector), "PE.A")
  expect_equal(ssm$ss[ssm$receiver_detector == "B.A"], k, tolerance = 0.05)
  expect_lt(ssm$ss[ssm$receiver_detector == "C.A"], 0.3)
})

test_that("a residual trend between source and receiver is reported", {
  # Under-compensation: the receiver's negatives rise with the source.
  withr::local_seed(13)
  n <- 6000
  src <- c(rnorm(n / 2, 0, 0.3), rnorm(n / 2, 3, 0.5))
  rec <- rnorm(n, 0, 0.3) + ifelse(src > 1.5, 0.4 * (src - 3), 0)
  tm <- cbind(S = src, R = rec)
  sp <- spreading_pairs(tm, c(S = 1.5, R = 1.5))
  expect_gt(sp$trend_rho[sp$source == "S" & sp$receiver == "R"], 0.3)
})

test_that("a negative trend on a heavily spread pair is unresolved, not an error", {
  # Truncation at the receiver's cut pulls the trend negative wherever the
  # receiver widens with the source, so a negative trend there is ambiguous.
  # A positive trend shows against that pull, and a negative one on a receiver
  # no wider at a given source level has no such explanation; both stay
  # flagged. C's trend widens its MARGINAL spread (substantial) all the same,
  # which is why the width at a given source level decides.
  withr::local_seed(14)
  n <- 6000
  mk <- function(slope, sd_pos) {
    src <- c(rnorm(n / 2, 0, 0.3), rnorm(n / 2, 3, 0.5))
    pos <- src > 1.5
    rec <- ifelse(pos, rnorm(n, slope * (src - 3), sd_pos(src)), rnorm(n, 0, 0.3))
    list(tmat = cbind(S = src, R = rec), thresholds = c(S = 1.5, R = 1.5))
  }
  flat <- function(s) function(src) s
  pops <- list(A = mk(-0.8, flat(0.6)), B = mk(0.8, flat(0.6)),
               C = mk(-0.8, flat(0.3)),
               # correct compensation, spread growing with the source past the cut
               D = mk(0, function(src) 0.3 * (1 + 2 * pmax(src - 1.5, 0))))
  gates <- lapply(pops, function(p) list(masks = list(cd45_pos = rep(TRUE, n))))
  pr <- run_spreading_report(pops, gates,
                             panel_of = c(A = "A", B = "B", C = "C", D = "D"))$pairs
  sr <- pr[pr$source == "S" & pr$receiver == "R", ]
  rownames(sr) <- sr$panel
  expect_true(sr["C", "substantial"])
  expect_lt(sr["C", "median_ratio_given_source"], 1.25)
  expect_true(sr["A", "trend_unresolved"])
  expect_false(sr["A", "possible_compensation_error"])
  expect_true(sr["B", "possible_compensation_error"])
  expect_true(sr["C", "possible_compensation_error"])
  expect_false(any(sr[c("B", "C"), "trend_unresolved"]))
  expect_lt(sr["D", "median_trend_rho"], 0)
  expect_false(sr["D", "possible_compensation_error"])
})
