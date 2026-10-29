# Acquisition-time QC for one sample

Runs before gating, so it uses log10 of the marker channels rather than
the fitted analysis transform. The test is a robust z on per-bin medians
and both scales are monotone, so the detection is equivalent; log10
needs no cofactor and is therefore available at the point where a
decision about the events is still possible.

## Usage

``` r
acquisition_qc_one(rd, tmat = NULL, n_bins = 40L, mad_k = 5)
```

## Arguments

- rd:

  a list from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)

- tmat:

  transformed marker matrix for the same events. When NULL, log10 of the
  marker columns is used.

- n_bins, mad_k:

  passed to
  [`detect_time_anomalies()`](https://bhagesh-h.github.io/cyRAVEN/reference/detect_time_anomalies.md)

## Value

list(summary = one-row data.frame, bins, flagged)
