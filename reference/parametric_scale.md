# Transform a measure for a parametric test

Transform a measure for a parametric test

## Usage

``` r
parametric_scale(x, transform = "logit", absolute = FALSE)
```

## Arguments

- x:

  Values.

- transform:

  "logit", "asin_sqrt" or "none".

- absolute:

  TRUE when `x` is a concentration rather than a percentage.

## Value

list(value, label).
