# Empirical logit for percentages

log((p + e) / (1 - p + e)), with `e` the smallest non-zero value of p or
of 1 - p among the values transformed together, so a population at 0% or
100% in some samples stays finite (Warton & Hui 2011).

## Usage

``` r
logit_pct(pct)
```

## Arguments

- pct:

  Numeric vector of percentages.

## Value

Numeric vector on the logit scale.
