# Dunn's post-hoc comparisons with a tie correction

The rank-based post-hoc test that follows a Kruskal-Wallis result.
Compares mean ranks computed over the whole data, which is what makes it
the correct follow-up rather than a set of pairwise Wilcoxon tests.

## Usage

``` r
dunn_test(values, groups)
```

## Arguments

- values:

  Numeric vector.

- groups:

  Grouping vector the same length.

## Value

data.frame of pairwise comparisons, or NULL when there are fewer than
two usable groups.
