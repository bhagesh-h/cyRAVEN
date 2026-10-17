# Parametric group comparison with its assumptions recorded

Runs the test an immunology paper usually reports, on arcsine square
root transformed percentages, and states beside each result whether the
two assumptions behind it held.

## Usage

``` r
parametric_group_tests(
  freq,
  group_of,
  reference = NULL,
  value_col = NULL,
  min_n = 3L,
  transform = "asin_sqrt"
)
```

## Arguments

- freq:

  Long frequency table with `sample_id`, `population` and a value
  column.

- group_of:

  Named vector mapping sample_id to group.

- reference:

  Reference group, or NULL for the first alphabetically.

- value_col:

  Value column; defaults to the resolved abundance measure.

- min_n:

  Minimum samples per group.

- transform:

  "asin_sqrt" (default) or "none".

## Value

data.frame, or NULL when nothing is testable.

## Details

Two groups give Welch's t-test, which does not assume equal variances,
plus Student's t-test for reference. More than two give Welch's ANOVA
plus the classical one-way ANOVA. Which one is defensible is decided by
the Brown-Forsythe column, not by the analyst after seeing the p-values.
