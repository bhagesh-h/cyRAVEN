# Marker symbol to the files' `$PnN`, over the run's reads

Marker symbol to the files' `$PnN`, over the run's reads

## Usage

``` r
gating_ml_markers(reads)
```

## Arguments

- reads:

  the run's reads, each carrying `marker_cols` (indices into
  `channel_names`), `channel_names` and `keywords`.

## Value

named character vector, marker symbol to the first `$PnN` seen for it;
the `channel_map` of
[`write_gating_ml()`](https://bhagesh-h.github.io/cyRAVEN/reference/write_gating_ml.md).
