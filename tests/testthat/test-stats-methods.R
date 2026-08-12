test_that("brown_forsythe_test detects unequal variance and is undefined at n=1", {
  set.seed(1)
  tight <- rnorm(20, 0, 1)
  wide <- rnorm(20, 0, 6)
  g <- rep(c("a", "b"), each = 20)
  bf <- cyRAVEN:::brown_forsythe_test(c(tight, wide), g)
  expect_true(is.finite(bf$p_value))
  expect_lt(bf$p_value, 0.05)

  same <- cyRAVEN:::brown_forsythe_test(c(rnorm(20), rnorm(20)), g)
  expect_gt(same$p_value, 0.05)

  # One group is not a comparison.
  expect_true(is.na(cyRAVEN:::brown_forsythe_test(rnorm(5), rep("a", 5))$p_value))
})

test_that("normality_report refuses to call a small non-significant test 'normal'", {
  freq <- data.frame(
    sample_id = paste0("s", 1:8),
    population = "CD4 T cells",
    pct_of_cd45_pos = c(10, 11, 12, 13, 20, 21, 22, 23),
    qc_status = "pass", stringsAsFactors = FALSE)
  group_of <- setNames(rep(c("hc", "case"), each = 4), paste0("s", 1:8))

  nt <- cyRAVEN:::normality_report(freq, group_of)
  expect_true(is.data.frame(nt))
  expect_true(all(c("shapiro_p", "interpretation") %in% names(nt)))

  # The trap this column exists to close: n = 4 cannot demonstrate normality,
  # so a non-significant Shapiro-Wilk must NOT read as licence for a t-test.
  small <- nt[nt$group != "(across groups)" & nt$n_donors < 15 &
                is.finite(nt$shapiro_p) & nt$shapiro_p >= 0.05, ]
  if (nrow(small))
    expect_true(all(grepl("NOT evidence of normality", small$interpretation)))

  expect_true(any(nt$group == "(across groups)"))
})

test_that("normality_report survives a population with no variance", {
  freq <- data.frame(sample_id = paste0("s", 1:6), population = "flat",
                     pct_of_cd45_pos = rep(3, 6), qc_status = "pass",
                     stringsAsFactors = FALSE)
  go <- setNames(rep(c("a", "b"), each = 3), paste0("s", 1:6))
  nt <- cyRAVEN:::normality_report(freq, go)
  expect_true(is.data.frame(nt))
  expect_true(any(grepl("no variance", nt$interpretation)))
})

test_that("the methods catalogue names the alternatives and says why not", {
  tb <- cyRAVEN:::statistical_methods_table(n_groups = 3L, n_tests = 24L)
  expect_true(all(c("method", "role", "rationale") %in% names(tb)))
  expect_true(any(tb$role == "used"))
  expect_true(any(tb$role == "not used"))

  # A reader arriving from an immunophenotyping paper looks for these by name.
  for (m in c("Student's t-test", "One-way ANOVA with Tukey", "Kruskal-Wallis",
              "Benjamini-Hochberg", "Cliff's delta", "Shapiro-Wilk"))
    expect_true(any(grepl(m, tb$method, fixed = TRUE)),
                info = paste("catalogue does not mention", m))

  # Every entry must justify itself; a bare "not used" is what this table exists
  # to avoid.
  expect_true(all(nzchar(tb$rationale)))
  expect_true(all(nchar(tb$rationale) > 30))
})

test_that("the catalogue adapts to the design", {
  two <- cyRAVEN:::statistical_methods_table(n_groups = 2L)
  kw <- two[grepl("^Kruskal-Wallis", two$method), ]
  expect_equal(kw$role, "not used")
  expect_true(grepl("three or more", kw$rationale))

  three <- cyRAVEN:::statistical_methods_table(n_groups = 3L)
  expect_equal(three$role[grepl("^Kruskal-Wallis", three$method)], "used")
})
