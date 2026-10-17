# Normality and equal-variance evidence, per population and group

The diagnostics behind cyRAVEN's use of rank tests. One Shapiro-Wilk per
population per group, plus one Brown-Forsythe per population across
groups.

## Usage

``` r
normality_report(freq, group_of, value_col = "pct_of_cd45_pos")
```

## Arguments

- freq:

  Long data frame with `sample_id`, `population` and a value column.

- group_of:

  Named character vector, sample id -\> group.

- value_col:

  Column holding the per-sample value.

## Value

data.frame, or NULL when there is nothing testable.
