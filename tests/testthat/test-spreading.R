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
