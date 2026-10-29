# =============================================================================
# EXACT RANK TESTS, TIES INCLUDED
# =============================================================================
#
# WHY NOT wilcox.test(exact = FALSE). The normal approximation is what R uses
# from 50 observations up, and it is poor below that, in the direction that
# matters: on two completely separated groups of three it gives p = 0.081 where
# the exact answer is 0.100, the smallest p-value three against three can
# produce. A cohort of a few patients per arm is the normal case here, so the
# approximation was reporting significance the design cannot deliver.
#
# WHY NOT wilcox.test() WITH ITS DEFAULT EITHER. R's exact distribution assumes
# no ties, and on tied data it silently falls back to the normal approximation.
# Population frequencies tie often -- several samples at 0% is a common
# reading -- so the default would be exact on some rows and approximate on the
# rest of the same table.
#
# WHAT IS COMPUTED. The exact conditional distribution of the rank statistic
# given the midranks actually observed: every way of assigning the pooled
# midranks to the two groups (or every sign pattern, for the paired test) is
# equally likely under the null, and the distribution of the statistic over
# them is counted exactly. This is the "exact" test coin's wilcox_test() and
# wilcoxsign_test() compute (Hothorn et al. 2008, the shift algorithm of
# Streitberg and Roehmel 1986), done here by dynamic programming over doubled
# midranks so every score is an integer. Without ties it is identical to R's
# exact test. Two-sided p-values follow wilcox.test(): twice the smaller tail.
# Above `max_n` observations the counts are large, the normal approximation is
# good, and it is used.

#' Exact two-sample Wilcoxon rank-sum p-value, ties included
#'
#' @param a,b numeric vectors; non-finite values are dropped.
#' @param max_n above this many observations in total, the normal
#'   approximation with continuity correction is used instead. Default `60L`.
#' @return two-sided p-value, or NA when either group is empty.
#' @keywords internal
rank_sum_p <- function(a, b, max_n = 60L) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n <- n1 + length(b)
  if (!n1 || n1 == n) return(NA_real_)
  if (n > max_n)
    return(suppressWarnings(stats::wilcox.test(a, b, exact = FALSE)$p.value))
  s <- as.integer(round(2 * rank(c(a, b))))
  tot <- sum(s)
  # dp[k + 1, w + 1]: the number of k-element subsets whose scores sum to w
  dp <- matrix(0, n1 + 1L, tot + 1L)
  dp[1L, 1L] <- 1
  for (i in seq_len(n)) {
    si <- s[i]
    for (k in min(i, n1):1L)
      dp[k + 1L, (si + 1L):(tot + 1L)] <- dp[k + 1L, (si + 1L):(tot + 1L)] +
        dp[k, seq_len(tot + 1L - si)]
  }
  two_sided_p(dp[n1 + 1L, ], sum(s[seq_len(n1)]))
}

#' Exact Wilcoxon signed-rank p-value for paired differences, ties included
#'
#' Zero differences are dropped, as [stats::wilcox.test()] drops them.
#'
#' @param d numeric vector of paired differences.
#' @param max_n above this many non-zero differences, the normal approximation
#'   with continuity correction is used instead. Default `60L`.
#' @return two-sided p-value, or NA when no difference is non-zero.
#' @keywords internal
signed_rank_p <- function(d, max_n = 60L) {
  d <- d[is.finite(d) & d != 0]
  n <- length(d)
  if (!n) return(NA_real_)
  if (n > max_n)
    return(suppressWarnings(stats::wilcox.test(d, exact = FALSE)$p.value))
  s <- as.integer(round(2 * rank(abs(d))))
  tot <- sum(s)
  # each score enters the positive sum or does not, with equal probability
  dist <- c(1, numeric(tot))
  for (si in s) dist <- dist + c(numeric(si), dist[seq_len(tot + 1L - si)])
  two_sided_p(dist, sum(s[d > 0]))
}

#' Two-sided p-value from a counted null distribution
#' @param counts counts of the statistic at 0, 1, ..., length - 1.
#' @param obs observed statistic.
#' @return twice the smaller tail probability, at most 1.
#' @keywords internal
two_sided_p <- function(counts, obs) {
  w <- seq_along(counts) - 1L
  pr <- counts / sum(counts)
  min(1, 2 * min(sum(pr[w <= obs]), sum(pr[w >= obs])))
}
