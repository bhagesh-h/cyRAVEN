# =============================================================================
# Per-sample threshold overrides
#
# The claim: an analyst can correct one sample's cut without applying that
# correction to every other sample, and the correction is recorded with who made
# it and why. The second half is the part that makes the first half acceptable.
# =============================================================================

bimodal <- function(gap = 5, n = 4000L, seed = 3) {
  had <- exists(".Random.seed", .GlobalEnv)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  c(stats::rnorm(n, 0, 0.6), stats::rnorm(n, gap, 0.6))
}

OVR <- list(
  D07 = list(CCR7 = list(threshold = 2.15,
                         reason = "valley sat inside the negative mode",
                         set_by = "ab")),
  D08 = list(CD4 = list(threshold = 1.0))          # bare, no attribution
)

test_that("an override wins over every derived source and says so", {
  x <- bimodal()
  auto <- resolve_threshold("CCR7", x)
  man  <- resolve_threshold("CCR7", x, override = OVR$D07$CCR7)

  expect_identical(auto$source, "valley")
  expect_identical(man$source, "manual")
  expect_equal(man$threshold, 2.15)
  expect_identical(man$override_reason, "valley sat inside the negative mode")
  expect_identical(man$override_by, "ab")
  # It also beats a config value, which is the whole ordering question: config
  # applies to the cohort, an override applies to one tube.
  both <- resolve_threshold("CCR7", x, cfg_value = 9.9, override = OVR$D07$CCR7)
  expect_identical(both$source, "manual")
  expect_equal(both$threshold, 2.15)
})

test_that("`manual` is a distinct source from `config`", {
  # They are different claims. `config` says the assay declares this cut for
  # every sample; `manual` says one person moved one sample's cut. Collapsing
  # them would make the audit trail unreadable.
  x <- bimodal()
  expect_identical(resolve_threshold("m", x, cfg_value = 3)$source, "config")
  expect_identical(resolve_threshold("m", x,
                                     override = list(threshold = 3))$source, "manual")
})

test_that("the derived path is untouched when no override applies", {
  # ADDITIVITY AT THE FUNCTION LEVEL. Every field a caller read before must hold
  # the value it held before.
  x <- bimodal()
  r <- resolve_threshold("CCR7", x)
  expect_identical(r$source, "valley")
  expect_false(r$needs_review)
  expect_true(is.na(r$override_reason))
  expect_true(is.na(r$override_by))
  # And the fallback branch still flags itself for review. 20000 uniform draws,
  # not fewer: density_valley() finds a spurious minimum in uniform noise at a
  # few thousand events often enough to make a smaller fixture flaky, which is
  # the property test-uncertainty.R measures directly.
  withr::local_seed(1)
  f <- resolve_threshold("m", stats::runif(20000))
  expect_identical(f$source, "quantile_fallback")
  expect_true(f$needs_review)
})

test_that("the lookup keys on sample and marker together", {
  expect_equal(sample_override(OVR, "D07", "CCR7")$threshold, 2.15)
  # Right marker, wrong sample.
  expect_null(sample_override(OVR, "D01", "CCR7"))
  # Right sample, wrong marker.
  expect_null(sample_override(OVR, "D07", "CD4"))
  expect_null(sample_override(OVR, "nobody", "nothing"))
  expect_null(sample_override(NULL, "D07", "CCR7"))
  expect_null(sample_override(list(), "D07", "CCR7"))
})

test_that("a bare number is accepted but carries no attribution", {
  e <- sample_override(OVR, "D08", "CD4")
  expect_equal(e$threshold, 1.0)
  r <- resolve_threshold("CD4", bimodal(), override = e)
  expect_identical(r$source, "manual")
  expect_true(is.na(r$override_by))
  expect_true(is.na(r$override_reason))
})

test_that("a malformed entry is ignored rather than applied", {
  # An override that cannot be read must not become a threshold of NA, which
  # would silently drop the marker from every population that reads it.
  bad <- list(S1 = list(M = list(reason = "forgot the number")))
  expect_null(sample_override(bad, "S1", "M"))
  expect_null(sample_override(list(S1 = list(M = list(threshold = "abc"))), "S1", "M"))
  expect_null(sample_override(list(S1 = list(M = list(threshold = NA))), "S1", "M"))
  # And resolve_threshold falls through to the derived value.
  r <- resolve_threshold("M", bimodal(), override = sample_override(bad, "S1", "M"))
  expect_identical(r$source, "valley")
})

test_that("an override that matched nothing is reported, not swallowed", {
  # THE FAILURE THIS GUARDS. The analyst believes a cut was corrected, the run
  # says nothing, and the uncorrected number is published.
  msgs <- capture.output(
    unused <- report_unused_overrides(OVR, applied = "D07\rCCR7"),
    type = "message")
  expect_identical(unused, "D08\rCD4")
  expect_match(paste(msgs, collapse = " "), "matched no sample and marker")
  expect_match(paste(msgs, collapse = " "), "D08 / CD4")

  # Nothing to say when everything applied.
  quiet <- capture.output(
    none <- report_unused_overrides(OVR, applied = c("D07\rCCR7", "D08\rCD4")),
    type = "message")
  expect_length(none, 0L)
  expect_length(quiet, 0L)
})

test_that("conformance reports a hand-set marker as such, not as agreement", {
  # A marker made to match the baseline by hand did not conform. Counting it as
  # a pass is the one outcome that would make the whole check misleading.
  base_thr <- data.frame(
    sample_id = sprintf("S%02d", 1:6), panel = "P1", marker = "CD4",
    threshold = c(2.0, 2.1, 1.9, 2.05, 2.0, 1.95), source = "valley",
    stringsAsFactors = FALSE)
  spec <- list(P = list(CD4 = "above"))
  freq <- data.frame(sample_id = sprintf("S%02d", 1:6), population = "P",
                     pct_of_cd45_pos = 10, is_control = FALSE, qc_status = "pass",
                     stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".rds")
  write_spec_baseline(path, base_thr, freq, spec, list(transform = "arcsinh"),
                      list(P1 = 5))
  base <- read_spec_baseline(path)

  cur <- base_thr
  cur$source[1] <- "manual"
  out <- specification_conformance(cur, freq, spec, base, transform = "arcsinh")
  row <- out$markers[out$markers$marker == "CD4", ]
  expect_identical(row$verdict, "manually set")
  expect_identical(row$n_manual_overrides, 1L)
  expect_match(row$note, "set by hand")

  # With nothing set by hand the column is absent entirely, so a run that
  # overrides nothing writes the table it always wrote.
  clean <- specification_conformance(base_thr, freq, spec, base,
                                     transform = "arcsinh")
  expect_false("n_manual_overrides" %in% names(clean$markers))
  expect_identical(clean$markers$verdict[1], "pass")
})
