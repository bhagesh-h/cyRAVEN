# Brown-Forsythe test for equality of variances

Levene's test computed on absolute deviations from the group MEDIAN
rather than the mean, which is the variant that holds up when the data
are skewed – the usual case for population frequencies bounded at zero.

## Usage

``` r
brown_forsythe_test(x, g)
```

## Arguments

- x:

  Numeric values.

- g:

  Grouping factor.

## Value

list(statistic, p_value, df1, df2), or NULLs when undefined.
