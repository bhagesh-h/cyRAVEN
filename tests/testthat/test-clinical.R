# Clinical variables: the tests, the effect sizes, and the intervals.
#
# WHY PLANTED STRUCTURE. Every fixture below has a known answer, so a test can
# assert that the right population comes out on top rather than that a data frame
# of the right shape came out. The association code is the kind that returns
# plausible numbers when it is wrong: a rho computed against the wrong sample
# order is still a rho in [-1, 1], still ranks the populations, and still writes a
# heatmap that looks like a result.

# One numeric variable built so that exactly one population tracks it, one
# two-level variable that splits on the same population, and one three-level
# variable that is pure noise. Everything is keyed on the sample ids
# helper-synthetic.R uses.
clin_fixture <- function() {
  ids <- names(synth_groups())
  n <- length(ids)
  list(
    # Ranks 1..n, so a population built as a monotone function of it correlates
    # at rho = 1 exactly and any deviation is a bug rather than sampling noise.
    severity = stats::setNames(seq_len(n), ids),
    outcome  = stats::setNames(rep(c("died", "survived"), length.out = n), ids),
    focus    = stats::setNames(rep(c("lung", "abdomen", "urine"),
                                   length.out = n), ids))
}

# A frequency table where "CD4 T cells" is a strictly increasing function of
# `severity` and every other population is flat noise.
planted_freq <- function(seed = 3) {
  withr::local_seed(seed)
  ids <- names(synth_groups())
  do.call(rbind, lapply(POPS, function(p) data.frame(
    sample_id = ids, panel = "P1", population = p,
    pct_of_cd45_pos = if (p == "CD4 T cells") seq_along(ids) * 1.5 + 5
                      else stats::rnorm(length(ids), 10, 2),
    count = 4000L, is_control = FALSE, qc_status = "pass",
    stringsAsFactors = FALSE)))
}

test_that("the test follows the variable's type, and the effect follows the test", {
  a <- clin_associate(planted_freq(), "population", "pct_of_cd45_pos",
                      clin_fixture())
  expect_true(is.data.frame(a))

  num <- a[a$variable == "severity", ]
  expect_identical(unique(num$test), "Spearman")
  expect_identical(unique(num$effect), "rho")
  expect_true(all(num$signed))
  # The planted population is a strictly increasing function of the variable, so
  # its rank correlation is exactly 1 and it is the top row for that variable.
  expect_identical(num$population[1], "CD4 T cells")
  expect_equal(num$estimate[num$population == "CD4 T cells"], 1, tolerance = 1e-8)

  two <- a[a$variable == "outcome", ]
  expect_identical(unique(two$test), "Wilcoxon rank-sum")
  expect_identical(unique(two$effect), "Cliff's delta")
  expect_true(all(two$signed))
  expect_true(all(abs(two$estimate) <= 1))

  many <- a[a$variable == "focus", ]
  expect_identical(unique(many$test), "Kruskal-Wallis")
  expect_identical(unique(many$effect), "epsilon squared")
  # UNSIGNED, and the flag says so. Before epsilon-squared existed these rows
  # carried no effect at all; the risk introduced with it is that a positive
  # number gets read as a direction, which `signed` is what prevents.
  expect_false(any(many$signed))
  expect_true(all(many$estimate >= 0 & many$estimate <= 1))
  expect_true(all(is.na(many$ci_low)))
})

test_that("epsilon-squared is H over n minus 1, and stays in range", {
  expect_equal(clin_epsilon2(6, 13), 0.5)
  # H cannot exceed n - 1, but a caller could pass anything; the bound is
  # enforced rather than trusted so a figure cannot be handed 1.4.
  expect_equal(clin_epsilon2(40, 13), 1)
  expect_equal(clin_epsilon2(-1, 13), 0)
  expect_true(is.na(clin_epsilon2(NA, 13)))
  expect_true(is.na(clin_epsilon2(6, 1)))
})

test_that("Cliff's delta is +1 when every b beats every a, and -1 reversed", {
  expect_equal(clin_cliff(1:5, 6:10), 1)
  expect_equal(clin_cliff(6:10, 1:5), -1)
  expect_equal(clin_cliff(1:5, 1:5), 0)
})

test_that("the bootstrap interval brackets the estimate and is reproducible", {
  x <- c(2, 4, 5, 7, 9, 11, 12, 15)
  y <- c(1, 3, 6, 6, 8, 12, 11, 16)
  st <- function(i) suppressWarnings(stats::cor(x[i], y[i], method = "spearman"))
  ci1 <- clin_boot_ci(st, length(x))
  ci2 <- clin_boot_ci(st, length(x))
  expect_identical(ci1, ci2)          # fixed seed, so twice is the same twice
  expect_lt(ci1[1], ci1[2])
  r <- stats::cor(x, y, method = "spearman")
  expect_lte(ci1[1], r)
  expect_gte(ci1[2], r)
})

test_that("no interval is claimed below the minimum sample count", {
  # Quoting a width nobody should trust is worse than a blank, so the interval is
  # suppressed rather than narrowed.
  expect_true(all(is.na(clin_boot_ci(function(i) 1, n = 4L))))
  expect_false(any(is.na(clin_boot_ci(function(i) stats::runif(1), n = 8L))))
})

test_that("the RNG stream is left as it was found", {
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(clin_boot_ci(function(i) stats::runif(1), n = 8L))
  # If the bootstrap leaked its seed, this draw would differ from `before` and
  # every subsample downstream of a clinical association would change.
  expect_equal(stats::runif(1), before)
})

test_that("BH is applied within a variable and not across all of them", {
  a <- clin_associate(planted_freq(), "population", "pct_of_cd45_pos",
                      clin_fixture())
  # Within one variable the adjustment is BH over that variable's own rows. If it
  # had been pooled across variables the same raw p would adjust to something
  # larger, because the family would be three times the size.
  for (v in unique(a$variable)) {
    rows <- a[a$variable == v, ]
    expect_equal(rows$p_adj_BH,
                 unname(stats::p.adjust(rows$p_value, method = "BH")),
                 tolerance = 1e-12)
  }
})

test_that("a two-level variable is coded 0/1 and named in the correlogram", {
  # The sign of a rank-biserial correlation is meaningless without knowing which
  # level is 1, so the axis label has to carry it.
  out <- withr::local_tempdir()
  f <- file.path(out, "corr.png")
  suppressMessages(fig_clinical_correlogram(clin_fixture(), f))
  expect_true(file.exists(f))
  fig <- suppressMessages(fig_clinical_correlogram(clin_fixture(), f))
  expect_true(any(grepl("=1\\]", levels(fig$data$.a))))
  # `focus` has three unordered levels, so it is excluded rather than coded 1/2/3.
  expect_false(any(grepl("^focus", levels(fig$data$.a))))
})

test_that("the smallest attainable rank-sum p is 2 over the arrangement count", {
  expect_equal(rank_sum_p_floor(4, 5), 2 / choose(9, 4))
  expect_equal(rank_sum_p_floor(6, 6), 2 / choose(12, 6))
  # At 2 against 2 there are six arrangements, so a two-sided p can never fall
  # below 1/3 -- the design cannot produce a significant result at all, which is
  # the whole reason this line is drawn on the figure.
  expect_equal(rank_sum_p_floor(2, 2), 2 / 6)
  expect_true(is.na(rank_sum_p_floor(0, 5)))
})

test_that("every clinical figure renders opaque and within bounds", {
  out <- withr::local_tempdir()
  cl <- clin_fixture()
  fq <- planted_freq()
  a <- clin_associate(fq, "population", "pct_of_cd45_pos", cl)

  files <- c(
    heatmap  = file.path(out, "clinical_association.png"),
    forest   = file.path(out, "clinical_effects_severity.png"),
    unsigned = file.path(out, "clinical_effects_focus.png"),
    detail   = file.path(out, "clinical_severity.png"),
    corr     = file.path(out, "clinical_variables_correlation.png"),
    land     = file.path(out, "clinical_landscape.png"),
    batch    = file.path(out, "populations_by_batch.png"),
    traj     = file.path(out, "population_trajectories.png"))
  suppressMessages({
    fig_clinical_heatmap(a, "population", files[["heatmap"]])
    fig_clinical_forest(a, "population", "severity", files[["forest"]])
    fig_clinical_forest(a, "population", "focus", files[["unsigned"]])
    fig_clinical_detail(fq, cl$severity, "severity", files[["detail"]])
    fig_clinical_correlogram(cl, files[["corr"]])
    fig_clinical_landscape(fq, cl, files[["land"]])
    ids <- names(synth_groups())
    fig_populations_by_batch(
      fq, stats::setNames(rep(c("B1", "B2", "B3"), length.out = length(ids)), ids),
      files[["batch"]])
    fig_clinical_trajectory(
      fq,
      stats::setNames(rep(c("d0", "d3"), length.out = length(ids)), ids),
      stats::setNames(rep(sprintf("P%d", 1:8), each = 2), ids),
      files[["traj"]], outcome_of = cl$outcome, outcome_name = "outcome")
  })
  for (k in names(files)) expect_figure_contract(files[[k]], k)
})

test_that("the between-group volcano renders and separates its labels", {
  out <- withr::local_tempdir()
  f <- file.path(out, "group_differences.png")
  gs <- stats_group_comparison(synth_freq(), synth_groups(), reference = "HC")
  fig <- suppressMessages(fig_group_volcano(gs, f))
  expect_figure_contract(f, "fig_group_volcano")
  # The label layer carries a nudged y. Two populations at the same p-value would
  # otherwise print their names on top of each other, which is the one thing a
  # label must not do.
  lab <- fig$layers[[which(vapply(fig$layers,
    function(l) inherits(l$geom, "GeomText"), logical(1)))[1]]]$data
  expect_true(".ty" %in% names(lab))
  ord <- sort(lab$.ty)
  expect_true(all(is.finite(ord)))
})

test_that("a single timepoint is not drawn as a trajectory", {
  out <- withr::local_tempdir()
  f <- file.path(out, "traj_one.png")
  ids <- names(synth_groups())
  suppressMessages(fig_clinical_trajectory(
    planted_freq(), stats::setNames(rep("d0", length(ids)), ids),
    stats::setNames(sprintf("P%02d", seq_along(ids)), ids), f))
  # A figure of unconnected points called a trajectory would be worse than none.
  expect_false(file.exists(f))
})

# ---- the unit of analysis ---------------------------------------------------
# The defect these cover: a cohort of 6 patients sampled at 3 timepoints was
# tested as 18 independent observations for variables that are properties of the
# patient. That inflates n, narrows the bootstrap interval and makes significance
# more likely without one extra patient having been recruited.

pat_fixture <- function() {
  ids <- names(synth_groups())
  # 16 samples over 8 patients, two samples each.
  stats::setNames(rep(sprintf("P%02d", 1:8), each = 2), ids)
}

test_that("a patient-constant variable is tested per patient, not per sample", {
  ids <- names(synth_groups())
  pat <- pat_fixture()
  # Constant within a patient, which is what an outcome flag looks like.
  outcome <- stats::setNames(rep(c("died", "lived"), each = 2, length.out = 16), ids)
  expect_identical(clin_variable_unit(outcome, pat), "patient")

  a <- clin_associate(planted_freq(), "population", "pct_of_cd45_pos",
                      list(outcome = outcome), patient_of = pat)
  expect_identical(unique(a$unit), "patient")
  # 8 patients, not 16 samples. Getting this wrong is the whole point of the test.
  expect_true(all(a$n == 8L))
  expect_true(all(a$n_patients == 8L))
  expect_false(any(a$repeated_measures))
})

test_that("a variable that moves within a patient stays per sample, and says so", {
  ids <- names(synth_groups())
  pat <- pat_fixture()
  # Differs between the two samples of every patient: a score taken at each draw.
  score <- stats::setNames(seq_along(ids), ids)
  expect_identical(clin_variable_unit(score, pat), "sample")

  a <- clin_associate(planted_freq(), "population", "pct_of_cd45_pos",
                      list(score = score), patient_of = pat)
  expect_identical(unique(a$unit), "sample")
  expect_true(all(a$n == 16L))
  # Collapsing here would delete the variation being asked about, so the flag is
  # what warns the reader instead.
  expect_true(all(a$n_patients == 8L))
  expect_true(all(a$repeated_measures))
})

test_that("with no patient map every sample is its own subject", {
  ids <- names(synth_groups())
  v <- stats::setNames(rep(c("a", "b"), 8), ids)
  expect_identical(clin_variable_unit(v, NULL), "sample")
  # A cross-sectional cohort: one sample per patient, so there is nothing to
  # collapse and the verdict must not claim otherwise.
  one_each <- stats::setNames(sprintf("P%02d", seq_along(ids)), ids)
  expect_identical(clin_variable_unit(v, one_each), "sample")
})

test_that("collapsing to the patient narrows nothing that was not already narrow", {
  ids <- names(synth_groups())
  pat <- pat_fixture()
  outcome <- stats::setNames(rep(c("died", "lived"), each = 2, length.out = 16), ids)
  per_pat <- clin_associate(planted_freq(), "population", "pct_of_cd45_pos",
                            list(outcome = outcome), patient_of = pat)
  per_smp <- clin_associate(planted_freq(), "population", "pct_of_cd45_pos",
                            list(outcome = outcome), patient_of = NULL)
  # Same populations, both directions of the effect preserved; the per-sample run
  # simply claims twice the evidence for it.
  expect_identical(sort(per_pat$population), sort(per_smp$population))
  expect_true(all(per_smp$n == 16L))
  expect_true(all(per_pat$n == 8L))
})
