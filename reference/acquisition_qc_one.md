# Acquisition-time QC for one sample

Runs before gating, so it judges the marker channels as stored,
compensated and linear, rather than on the fitted analysis transform,
which does not exist yet. The test is a robust z on per-bin medians,
each scaled by its own channel's spread, so a smooth monotone rescaling
barely changes it where bin medians vary little; what matters is that
the scale has no kink. log10 of max(x, 1), used before, pinned every
value at or below 1 to 0, and the median of a mostly negative
compensated channel sits on that kink: simulated with negatives centred
on zero it flagged 0.22 stable bins per channel and found a +80 drift in
1.9 of 3 bins, against 0 and 3.1 on the stored values.

## Usage

``` r
acquisition_qc_one(rd, tmat = NULL, n_bins = 40L, mad_k = 5)
```

## Arguments

- rd:

  a list from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)

- tmat:

  transformed marker matrix for the same events. When NULL, the marker
  columns as stored are used.

- n_bins, mad_k:

  passed to
  [`detect_time_anomalies()`](https://bhagesh-h.github.io/cyRAVEN/reference/detect_time_anomalies.md)

## Value

list(summary = one-row data.frame, bins, flagged)
