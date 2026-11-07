# Exact two-sample Wilcoxon rank-sum p-value, ties included

Exact two-sample Wilcoxon rank-sum p-value, ties included

## Usage

``` r
rank_sum_p(a, b, max_n = 60L)
```

## Arguments

- a, b:

  numeric vectors; non-finite values are dropped.

- max_n:

  above this many observations in total, the normal approximation with
  continuity correction is used instead. Default `60L`.

## Value

two-sided p-value, or NA when either group is empty.
