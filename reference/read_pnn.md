# Each channel's `$PnN` as the file writes it

`channel_names` carry the `alter.names` spelling (`FITC.A`); Gating-ML
names a dimension by the `$PnN` keyword itself (`FITC-A`).

## Usage

``` r
read_pnn(r)
```

## Arguments

- r:

  one read, carrying `channel_names` and `keywords`.

## Value

character vector parallel to `channel_names`.
