test_that("write_design_feasibility reports a group too small to test", {
  out <- withr::local_tempdir()
  # Three donors in 'ctrl', one in 'case': the case arm cannot be tested.
  g <- c(s1 = "ctrl", s2 = "ctrl", s3 = "ctrl", s4 = "case")
  d <- c(s1 = "p1", s2 = "p2", s3 = "p3", s4 = "p4")

  res <- cyRAVEN:::write_design_feasibility(
    g, "cohort", min_n = 3L, reference = "ctrl", outdir = out, donor_of = d)

  expect_true(file.exists(file.path(out, "design_feasibility.csv")))
  expect_equal(nrow(res), 2L)

  case <- res[res$group == "case", ]
  ctrl <- res[res$group == "ctrl", ]
  expect_equal(case$will_be_tested, "no")
  expect_equal(ctrl$will_be_tested, "yes")
  expect_match(case$reason, "below --min-group-n 3")
  # The reference flag follows the argument, not alphabetical order.
  expect_true(ctrl$is_reference)
  expect_false(case$is_reference)
})

test_that("a donor in two groups is flagged as an invalid unpaired comparison", {
  out <- withr::local_tempdir()
  # The repeated-measures trap: three donors, each sampled at two timepoints.
  # Every group clears min_n, so the test RUNS, and running it is the problem.
  g <- c(a0 = "d0", b0 = "d0", c0 = "d0", a3 = "d3", b3 = "d3", c3 = "d3")
  d <- c(a0 = "p1", b0 = "p2", c0 = "p3", a3 = "p1", b3 = "p2", c3 = "p3")

  res <- cyRAVEN:::write_design_feasibility(
    g, "timepoint", min_n = 3L, reference = "d0", outdir = out, donor_of = d)

  expect_true(all(res$will_be_tested == "yes"))
  expect_true(all(res$valid_unpaired == "no"))
  expect_match(res$reason[1], "appear in another group")
  expect_match(res$reason[1], "--paired-column")
  expect_true(all(res$n_donors == 3L))
})

test_that("independent groups of adequate size are valid", {
  out <- withr::local_tempdir()
  g <- setNames(rep(c("ctrl", "case"), each = 4), paste0("s", 1:8))
  d <- setNames(paste0("p", 1:8), paste0("s", 1:8))

  res <- cyRAVEN:::write_design_feasibility(
    g, "cohort", min_n = 3L, reference = "ctrl", outdir = out, donor_of = d)

  expect_true(all(res$will_be_tested == "yes"))
  expect_true(all(res$valid_unpaired == "yes"))
  expect_true(all(res$reason == ""))
})

test_that("tested = FALSE records the flag as the reason and tests nothing", {
  out <- withr::local_tempdir()
  g <- setNames(rep(c("ctrl", "case"), each = 4), paste0("s", 1:8))
  d <- setNames(paste0("p", 1:8), paste0("s", 1:8))

  res <- cyRAVEN:::write_design_feasibility(
    g, "cohort", min_n = 3L, reference = "ctrl", outdir = out,
    donor_of = d, tested = FALSE)

  expect_true(all(res$will_be_tested == "no"))
  expect_match(res$reason[1], "--no-group-tests")
  # The design itself is still sound; only the run declined to test it.
  expect_true(all(res$valid_unpaired == "yes"))
})

test_that("without a donor column the repeated-measures check says so", {
  out <- withr::local_tempdir()
  g <- setNames(rep(c("ctrl", "case"), each = 4), paste0("s", 1:8))

  res <- cyRAVEN:::write_design_feasibility(
    g, "cohort", min_n = 3L, reference = "ctrl", outdir = out, donor_of = NULL)

  expect_true(all(is.na(res$n_donors)))
  expect_true(all(res$valid_unpaired == "unknown"))
  expect_match(res$reason[1], "not checked")
})

test_that("more samples than donors inside one group is reported", {
  out <- withr::local_tempdir()
  # Four samples but two donors: n is inflated by repeat sampling.
  g <- setNames(rep(c("ctrl", "case"), each = 4), paste0("s", 1:8))
  d <- setNames(c("p1", "p1", "p2", "p2", "p3", "p4", "p5", "p6"),
                paste0("s", 1:8))

  res <- cyRAVEN:::write_design_feasibility(
    g, "cohort", min_n = 3L, reference = "ctrl", outdir = out, donor_of = d)

  ctrl <- res[res$group == "ctrl", ]
  expect_equal(ctrl$n_samples, 4L)
  expect_equal(ctrl$n_donors, 2L)
  expect_equal(ctrl$valid_unpaired, "no")
  expect_match(ctrl$reason, "unit of replication is the donor")
})

test_that("min_n is honoured by stats_group_comparison itself", {
  freq <- data.frame(
    sample_id = paste0("s", 1:6),
    population = "CD4 T cells",
    pct_of_cd45_pos = c(10, 11, 12, 20, 21, 22),
    stringsAsFactors = FALSE)
  g <- setNames(rep(c("ctrl", "case"), each = 3), paste0("s", 1:6))

  ok <- cyRAVEN::stats_group_comparison(freq, g, reference = "ctrl", min_n = 3L)
  expect_equal(nrow(ok), 1L)

  # Same data, higher bar: no comparison survives.
  none <- cyRAVEN::stats_group_comparison(freq, g, reference = "ctrl", min_n = 4L)
  expect_null(none)
})
