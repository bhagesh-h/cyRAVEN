# Rank F test of a categorical covariate after a grouping term

`lm(rank(y) ~ g + f)` against `lm(rank(y) ~ g)`: the rank-transform test
(Conover & Iman 1981) of whether `f` explains the outcome once the
groups are accounted for.

## Usage

``` r
rank_f_within(y, f, g)
```

## Arguments

- y:

  numeric outcome.

- f:

  categorical covariate.

- g:

  grouping.

## Value

the p-value, or NA where the model cannot be fitted.
