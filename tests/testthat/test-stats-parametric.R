# =============================================================================
# Parametric tests, post-hoc comparisons, and the assumption checks that decide
# which of them is defensible.
# =============================================================================

test_that("the arcsine transform maps percentages onto [0, pi/2] and is monotone", {
  p <- c(0, 1, 25, 50, 75, 99, 100)
  a <- cyRAVEN:::asin_sqrt_pct(p)
  expect_equal(a[1], 0)
  expect_equal(a[length(a)], pi / 2)
  expect_true(all(diff(a) > 0))
  # Out-of-range values are clamped rather than returned as NaN, which is what
  # asin(sqrt(x > 1)) would otherwise produce.
  expect_true(is.finite(cyRAVEN:::asin_sqrt_pct(101)))
  expect_true(is.finite(cyRAVEN:::asin_sqrt_pct(-1)))
})

test_that("cohens_d has the right sign, scale and undefined cases", {
  a <- c(10, 11, 12, 11, 10)
  b <- a + 5
  d <- cyRAVEN:::cohens_d(a, b)
  expect_gt(d, 0)                                  # b above a is positive
  expect_lt(cyRAVEN:::cohens_d(b, a), 0)           # and reverses
  expect_equal(cyRAVEN:::cohens_d(a, a), 0)
  # One observation cannot give a spread, so the effect size is undefined
  # rather than infinite.
  expect_true(is.na(cyRAVEN:::cohens_d(1, c(2, 3, 4))))
  expect_true(is.na(cyRAVEN:::cohens_d(rep(5, 4), rep(5, 4))))
})

test_that("games_howell finds a separated group and spares the two that overlap", {
  set.seed(42)
  v <- c(rnorm(8, 10, 1), rnorm(8, 10.2, 1), rnorm(8, 20, 1))
  g <- rep(c("a", "b", "c"), each = 8)
  gh <- cyRAVEN:::games_howell(v, g)

  expect_equal(nrow(gh), 3L)                        # three pairs from three groups
  expect_true(all(gh$test == "Games-Howell"))
  expect_true(all(gh$p_value >= 0 & gh$p_value <= 1))

  ab <- gh[gh$group_a == "a" & gh$group_b == "b", ]
  ac <- gh[gh$group_a == "a" & gh$group_b == "c", ]
  expect_gt(ab$p_value, 0.05)
  expect_lt(ac$p_value, 0.05)
})

test_that("games_howell tolerates unequal variance and unequal n, which is its point", {
  set.seed(7)
  v <- c(rnorm(20, 10, 1), rnorm(4, 10, 8))
  g <- c(rep("tight", 20), rep("wide", 4))
  gh <- cyRAVEN:::games_howell(v, g)
  expect_equal(nrow(gh), 1L)
  expect_true(is.finite(gh$p_value))
  # Welch-Satterthwaite df is far below the pooled n - 2 a Tukey test would use.
  expect_lt(gh$df, 22)
})

test_that("games_howell returns NULL when fewer than two groups are usable", {
  expect_null(cyRAVEN:::games_howell(c(1, 2, 3), rep("a", 3)))
  # A group of one has no variance, so it cannot enter the comparison.
  expect_null(cyRAVEN:::games_howell(c(1, 2, 3, 9), c("a", "a", "a", "b")))
})

test_that("dunn_test ranks, corrects for ties, and separates the right group", {
  set.seed(11)
  v <- c(rnorm(8, 10), rnorm(8, 10.1), rnorm(8, 25))
  g <- rep(c("a", "b", "c"), each = 8)
  dn <- cyRAVEN:::dunn_test(v, g)

  expect_equal(nrow(dn), 3L)
  expect_true(all(dn$test == "Dunn"))
  ac <- dn[dn$group_a == "a" & dn$group_b == "c", ]
  expect_lt(ac$p_value, 0.05)

  # Heavily tied data must still produce finite p-values.
  tied <- cyRAVEN:::dunn_test(c(rep(1, 6), rep(2, 6)), rep(c("a", "b"), each = 6))
  expect_true(all(is.finite(tied$p_value)))
})

test_that("parametric_group_tests records assumptions and does not hide a failure", {
  set.seed(5)
  freq <- data.frame(
    sample_id = paste0("s", 1:12),
    population = "CD4 T cells",
    pct_of_cd45_pos = c(rnorm(6, 20, 2), rnorm(6, 30, 2)),
    stringsAsFactors = FALSE)
  g <- setNames(rep(c("ctrl", "case"), each = 6), paste0("s", 1:12))

  res <- cyRAVEN:::parametric_group_tests(freq, g, reference = "ctrl")
  expect_equal(nrow(res), 1L)
  expect_equal(res$primary_test, "Welch's t-test")
  expect_equal(res$transform, "arcsine sqrt")
  expect_lt(res$p_value, 0.05)

  # Both assumption checks are present on the row, whatever they say.
  expect_true(all(c("shapiro_p", "brown_forsythe_p", "residuals_normal",
                    "equal_variance", "assumptions_met") %in% names(res)))
  # Means are reported on the percentage scale, not the arcsine scale, because
  # nobody can interpret the latter.
  expect_gt(res$mean_reference_pct, 10)
  expect_lt(res$mean_reference_pct, 40)
})

test_that("three groups give Welch's ANOVA with eta squared", {
  set.seed(9)
  freq <- data.frame(
    sample_id = paste0("s", 1:18),
    population = "B cells",
    pct_of_cd45_pos = c(rnorm(6, 5, 1), rnorm(6, 10, 1), rnorm(6, 20, 1)),
    stringsAsFactors = FALSE)
  g <- setNames(rep(c("g1", "g2", "g3"), each = 6), paste0("s", 1:18))

  res <- cyRAVEN:::parametric_group_tests(freq, g, reference = "g1")
  expect_equal(res$primary_test, "Welch's ANOVA")
  expect_equal(res$n_groups, 3L)
  expect_lt(res$p_value, 0.05)
  expect_true(is.finite(res$eta_squared))
  expect_gt(res$eta_squared, 0.5)          # groups are well separated
})

test_that("a group below min_n produces no parametric row rather than a bad one", {
  freq <- data.frame(
    sample_id = paste0("s", 1:5),
    population = "NK cells",
    pct_of_cd45_pos = c(5, 6, 7, 20, 21),
    stringsAsFactors = FALSE)
  g <- setNames(c("ctrl", "ctrl", "ctrl", "case", "case"), paste0("s", 1:5))
  expect_null(cyRAVEN:::parametric_group_tests(freq, g, reference = "ctrl", min_n = 3L))
})

test_that("posthoc_group_tests runs all three methods and needs three groups", {
  set.seed(13)
  freq <- data.frame(
    sample_id = paste0("s", 1:18),
    population = "B cells",
    pct_of_cd45_pos = c(rnorm(6, 5, 1), rnorm(6, 10, 1), rnorm(6, 20, 1)),
    stringsAsFactors = FALSE)
  g <- setNames(rep(c("g1", "g2", "g3"), each = 6), paste0("s", 1:18))

  ph <- cyRAVEN:::posthoc_group_tests(freq, g)
  expect_setequal(unique(ph$test), c("Games-Howell", "Tukey HSD", "Dunn"))
  expect_true(all(c("p_value", "p_adjusted_BH") %in% names(ph)))
  # Adjustment is within a method, so each family keeps its own three pairs.
  expect_equal(as.integer(table(ph$test)[["Games-Howell"]]), 3L)

  # Two groups is not a post-hoc situation.
  g2 <- setNames(rep(c("g1", "g2"), each = 9), paste0("s", 1:18))
  expect_null(cyRAVEN:::posthoc_group_tests(freq, g2))
})

test_that("ignore_channel_match matches whole names and globs, not substrings", {
  sym <- c("CD16", "CD161", "[AF color 1] - Area", "[AF color 2] - Area", "CD3")

  # Exact, so CD16 must not take CD161 with it.
  hit <- cyRAVEN:::ignore_channel_match(sym, "CD16")
  expect_equal(sym[hit], "CD16")

  # A glob with regex metacharacters in it still behaves as a glob.
  hit2 <- cyRAVEN:::ignore_channel_match(sym, "[AF color*")
  expect_equal(sum(hit2), 2L)

  # Case-insensitive, comma-separated lists, and no pattern matches nothing.
  expect_true(cyRAVEN:::ignore_channel_match("cd3", "CD3"))
  expect_equal(sum(cyRAVEN:::ignore_channel_match(sym, c("CD3", "CD16"))), 2L)
  expect_false(any(cyRAVEN:::ignore_channel_match(sym, NULL)))
  expect_false(any(cyRAVEN:::ignore_channel_match(sym, "")))
})
