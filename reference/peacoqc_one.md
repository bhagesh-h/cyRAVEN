# Run PeacoQC on one read, prepared as its vignette prescribes

Run PeacoQC on one read, prepared as its vignette prescribes

## Usage

``` r
peacoqc_one(rd)
```

## Arguments

- rd:

  a list from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md),
  compensated.

## Value

list(good = logical per event or NULL, pct_removed, pct_it, pct_mad,
pct_consecutive, reason).
