# The decision table staining_verdict() implements.
#
# WHY THIS FILE EXISTS. This function sits behind the most severe defect the
# package has recorded, documented at length in vignettes/known-limitations.Rmd:
# an absent `is_control` column was read as "this might be a control", so any
# sample whose CD45 gate found no density minimum was promoted to the unstained
# reference for its whole panel. Every threshold there then became the 99.5th
# percentile of that one sample, and where it was really a stained sample every
# population collapsed -- observed at 0.034% of CD45+ where the true answer was
# a hundred times that, with no error and no warning.
#
# It had no test at all.
#
# READ THIS BEFORE ASSUMING THESE TESTS PROTECT THE FIX. They do not, and
# cannot on their own. The promotion still happens here: with declared_control
# = NA a fallback sample is still returned with is_reference = TRUE, and that is
# deliberate, because this function only reports what the gate looks like. The
# guarantee is enforced one level up, in run_cyraven() at R/pipeline.R:520-542,
# which demotes every discovered reference unless the sheet actually declares a
# control or --discover-controls is passed. Tests below pin this function's half
# of the contract; the demotion half needs a pipeline-level test and does not
# have one.

fake_gate <- function(pct, source = "valley") {
  list(counts = data.frame(gate = c("live_cells", "cd45_pos"),
                           pct_of_parent = c(100, pct),
                           stringsAsFactors = FALSE),
       cd45_source = source)
}

test_that("a declared control is a reference and never a sample", {
  v <- staining_verdict(fake_gate(60), declared_control = TRUE)
  expect_true(v$is_control)
  expect_true(v$is_reference)
  expect_false(v$include)
  expect_equal(v$qc_status, "control")
})

test_that("is_control = FALSE is a positive assertion, so a failure is a failure", {
  # This is the distinction the defect turned on. The same gate, with the sheet
  # saying "this is a sample", must fail rather than become the reference.
  v <- staining_verdict(fake_gate(10, "quantile_fallback"),
                        declared_control = FALSE)
  expect_equal(v$qc_status, "failed")
  expect_false(v$is_control)
  expect_false(v$is_reference)
  expect_false(v$include)
  expect_match(v$verdict, "declared a sample, not a control")
})

test_that("an absent column still promotes here, and the caller must demote it", {
  # Pinning the behaviour that made the defect possible, so that if someone
  # changes it they change it knowingly. The protection is at
  # R/pipeline.R:520-542, not here.
  v <- staining_verdict(fake_gate(10, "quantile_fallback"),
                        declared_control = NA)
  expect_true(v$is_reference)
  expect_true(v$is_control)
  expect_equal(v$qc_status, "control")
  expect_match(v$verdict, "treated as unstained/control")
})

test_that("a quantile fallback is judged on the fallback, not on the percentage", {
  # The fraction above a q-quantile cut is (1 - q) by construction, so testing
  # it against a floor would pass an unstained file every time. Absence of a
  # CD45 mode is the diagnosis, whatever the percentage says.
  high <- staining_verdict(fake_gate(95, "quantile_fallback"),
                           declared_control = FALSE)
  expect_equal(high$qc_status, "failed")
  expect_false(high$include)
})

test_that("--include-qc-failed forces a declared sample through, still flagged", {
  v <- staining_verdict(fake_gate(10, "quantile_fallback"),
                        declared_control = FALSE, force_include = TRUE)
  expect_true(v$include)
  expect_equal(v$qc_status, "pass")
  expect_false(v$is_reference)
  expect_match(v$verdict, "included anyway")
  expect_match(v$verdict, "NO evidence")
})

test_that("force_include does not override a declared control", {
  v <- staining_verdict(fake_gate(10), declared_control = TRUE,
                        force_include = TRUE)
  expect_false(v$include)
  expect_true(v$is_reference)
})

test_that("a panel without CD45 is not judged on staining at all", {
  v <- staining_verdict(fake_gate(NA_real_, "skipped"), declared_control = FALSE)
  expect_true(v$include)
  expect_equal(v$qc_status, "pass")
  expect_false(v$is_reference)
  expect_match(v$verdict, "CD45 absent from panel")
})

test_that("a low CD45 fraction fails a declared sample", {
  v <- staining_verdict(fake_gate(1), declared_control = FALSE,
                        min_cd45_pct = 5)
  expect_false(v$include)
  expect_false(v$is_reference)
})

test_that("a stained or FMO control is never the unstained reference", {
  # The sepsis cohort's failure: a minus-multiple tube, declared is_control so
  # it stays out of the statistics, became the panel's unstained reference once
  # it shared the samples' panel, and set CD3, CD8 and CD56 cuts from its own
  # stained distribution.
  gates <- list(unstained = fake_gate(0.4, "quantile_fallback"),
                fmm       = fake_gate(62),
                full      = fake_gate(71),
                blank_lo  = fake_gate(1.2))
  rc <- reference_controls(names(gates), gates, fmo_tubes = "fmm")
  expect_identical(rc$keep, c("unstained", "blank_lo"))
  expect_setequal(names(rc$refused), c("fmm", "full"))
  expect_match(rc$refused[["fmm"]], "fmo_for markers only", fixed = TRUE)
  expect_match(rc$refused[["full"]], "stained", fixed = TRUE)
  # Nothing to choose from is not an error.
  expect_identical(reference_controls(character(0), list())$keep, character(0))
})

test_that("the live share is reported, and judged only against a declared minimum", {
  # den Braanker 2021: low-viability samples "should be excluded from analysis",
  # with the minimum left to the assay. No minimum, no judgement.
  g <- fake_gate(60)
  g$counts$pct_of_parent[g$counts$gate == "live_cells"] <- 62
  g$viability_source <- "valley"
  v <- staining_verdict(g, declared_control = FALSE)
  expect_equal(v$pct_live_of_singlets, 62)
  expect_equal(v$qc_status, "pass")

  v <- staining_verdict(g, declared_control = FALSE, min_live_pct = 80)
  expect_equal(v$qc_status, "failed")
  expect_false(v$include)
  expect_false(v$is_reference)
  expect_match(v$verdict, "low viability \\(62.0% of single cells live < 80%\\)")
  # not a control either when the sheet declares nothing: low viability says
  # nothing about whether the tube was stained
  v <- staining_verdict(g, declared_control = NA, min_live_pct = 80)
  expect_false(v$is_control)
  expect_equal(v$qc_status, "failed")

  v <- staining_verdict(g, declared_control = FALSE, min_live_pct = 80,
                        force_include = TRUE)
  expect_true(v$include)
  expect_match(v$verdict, "included anyway")

  expect_equal(staining_verdict(g, FALSE, min_live_pct = 50)$qc_status, "pass")
})

test_that("no viability gate means no live share, and nothing to judge", {
  g <- fake_gate(60)
  for (src in list(NULL, "skipped", "skipped_autofix")) {
    g$viability_source <- src
    v <- staining_verdict(g, declared_control = FALSE, min_live_pct = 80)
    expect_true(is.na(v$pct_live_of_singlets))
    expect_equal(v$qc_status, "pass")
  }
})

test_that("a mostly dead sample shows its second live reading", {
  # 30% live, 70% dead: the dye's tallest mode is the dead cells, so the gate
  # (a valley only above that mode, else the tail gate) keeps nearly everyone.
  # The minimum below the dominant mode gives the live share for this case.
  withr::local_seed(31)
  x <- c(rnorm(3000, 0.5, 0.3), rnorm(7000, 3.5, 0.4))
  vmode <- live_mode(x)
  expect_gt(vmode, 3)
  expect_true(is.na(density_valley(x, floor = vmode)))
  v_alt <- density_valley(x, floor = 0, ceiling = vmode)
  expect_equal(100 * mean(x < v_alt), 30, tolerance = 0.05)
  # a ceiling below every valley finds none
  expect_true(is.na(density_valley(x, floor = 0, ceiling = 1)))

  g <- fake_gate(60)
  g$counts$pct_of_parent[g$counts$gate == "live_cells"] <- 99.6
  g$viability_source <- "tailgate"
  g$live_pct_if_dominant_dead <- 30.1
  v <- staining_verdict(g, declared_control = FALSE)
  expect_equal(v$pct_live_if_dominant_dead, 30.1)
  g$viability_source <- "skipped_autofix"
  expect_true(is.na(staining_verdict(g, FALSE)$pct_live_if_dominant_dead))
})
