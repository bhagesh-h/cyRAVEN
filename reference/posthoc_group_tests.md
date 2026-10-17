# Post-hoc comparisons for every population with more than two groups

Runs all three: Games-Howell for the unequal-variance case, Tukey HSD
for the equal-variance case, and Dunn for the rank case. All three are
written, with the assumption columns from the parametric table saying
which one to read, so the choice is documented rather than made
silently.

## Usage

``` r
posthoc_group_tests(
  freq,
  group_of,
  value_col = NULL,
  min_n = 3L,
  transform = "asin_sqrt"
)
```

## Arguments

- freq, group_of, value_col, min_n:

  As in
  [`parametric_group_tests()`](https://bhagesh-h.github.io/cyRAVEN/reference/parametric_group_tests.md).

- transform:

  "asin_sqrt" or "none".

## Value

data.frame, or NULL.
