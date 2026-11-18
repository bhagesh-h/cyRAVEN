# Whether a cut was placed without a density minimum

`quantile_fallback` and the default path's `tailgate` both mean the
sample's density showed no minimum for this marker. Counted together
wherever the run reports how often that happened, so the count does not
fall just because the tail gate now stands in for the quantile.

## Usage

``` r
no_minimum(source)
```

## Arguments

- source:

  `source` strings from
  [`resolve_threshold()`](https://bhagesh-h.github.io/cyRAVEN/reference/resolve_threshold.md).

## Value

logical vector.
