# Spearman correlation of two variables with a grouping held fixed

The partial Spearman correlation controlling for group membership: both
variables are ranked, each group's mean rank is removed, and the
residuals are correlated, which is the partial correlation of the ranks
on the group indicators. Tested with t on n - 2 - (k - 1) degrees of
freedom for k groups, the inference ppcor uses for a partial
correlation.

## Usage

``` r
partial_spearman_within(x, y, g)
```

## Arguments

- x, y:

  numeric vectors.

- g:

  grouping, the same length.

## Value

list(rho, p), NA where it cannot be estimated.
