# Games-Howell post-hoc comparisons

The post-hoc test that pairs with Welch's ANOVA: it uses Welch's
standard error and pair-specific degrees of freedom, so it does not
assume the equal variances or equal group sizes that Tukey's HSD does.

## Usage

``` r
games_howell(values, groups)
```

## Arguments

- values:

  Numeric vector.

- groups:

  Grouping vector the same length.

## Value

data.frame of pairwise comparisons, or NULL when fewer than two groups
have at least two observations.
