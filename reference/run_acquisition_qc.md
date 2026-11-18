# Run acquisition-time QC across a cohort

Run acquisition-time QC across a cohort

## Usage

``` r
run_acquisition_qc(reads, pops = NULL, n_bins = 40L, mad_k = 5)
```

## Arguments

- reads:

  named list of
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md)
  results

- pops:

  per-sample scoring results, used for the transformed matrix

- n_bins, mad_k:

  detector settings

## Value

list(summary, bins, flagged) where `flagged` is a named list of
per-event logical vectors
