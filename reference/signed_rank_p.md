# Exact Wilcoxon signed-rank p-value for paired differences, ties included

Zero differences are dropped, as
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) drops
them.

## Usage

``` r
signed_rank_p(d, max_n = 60L)
```

## Arguments

- d:

  numeric vector of paired differences.

- max_n:

  above this many non-zero differences, the normal approximation with
  continuity correction is used instead. Default `60L`.

## Value

two-sided p-value, or NA when no difference is non-zero.
