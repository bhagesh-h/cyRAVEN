# Exact rank tests: identical to R's exact test without ties, exact (as coin
# computes it) with them, and never the small-sample normal approximation.

test_that("without ties the rank-sum p equals R's exact test", {
  for (nn in list(c(3, 3), c(3, 4), c(4, 4), c(2, 8), c(5, 5), c(6, 9))) {
    set.seed(sum(nn)); a <- rnorm(nn[1]); b <- rnorm(nn[2], 1)
    expect_equal(rank_sum_p(a, b), stats::wilcox.test(a, b)$p.value, tolerance = 1e-10)
  }
  # the case the approximation got wrong: 0.081 instead of 0.100
  expect_equal(rank_sum_p(1:3, 4:6), 0.1)
})

test_that("with ties the rank-sum p is the exact permutation p", {
  a <- c(0, 0, 1.2, 3); b <- c(0, 2, 4, 5, 5)
  # brute force over every split of the pooled midranks
  r <- rank(c(a, b)); obs <- sum(r[1:4])
  sums <- combn(9, 4, function(i) sum(r[i]))
  p_bf <- min(1, 2 * min(mean(sums <= obs + 1e-9), mean(sums >= obs - 1e-9)))
  expect_equal(rank_sum_p(a, b), p_bf, tolerance = 1e-12)
})

test_that("the signed-rank p is exact, matches R without ties, and drops zeros", {
  d <- c(1.5, 2.1, -0.4, 3.3, 2.8, 0.9)
  expect_equal(signed_rank_p(d), stats::wilcox.test(d)$p.value, tolerance = 1e-10)
  dz <- c(d, 0, 0)
  expect_equal(signed_rank_p(dz), signed_rank_p(d))
  # ties in |d|: brute force over every sign pattern
  dt <- c(1, 1, -2, 2, 3, 3, 3)
  r <- rank(abs(dt)); obs <- sum(r[dt > 0])
  sgn <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(dt))))
  v <- apply(sgn, 1, function(z) sum(r[z]))
  p_bf <- min(1, 2 * min(mean(v <= obs + 1e-9), mean(v >= obs - 1e-9)))
  expect_equal(signed_rank_p(dt), p_bf, tolerance = 1e-12)
})

test_that("large samples use the normal approximation", {
  set.seed(1); a <- rnorm(40); b <- rnorm(40, 0.5)
  expect_equal(rank_sum_p(a, b), stats::wilcox.test(a, b, exact = FALSE)$p.value)
})
