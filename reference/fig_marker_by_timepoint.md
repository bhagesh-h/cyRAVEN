# Every marker's median intensity per population, split by timepoint

The state heatmap asks what a population expresses; this asks whether
that changed over the admission. Split rather than pooled, because
pooling the timepoints is what hides a change that happens in all
patients at once.

## Usage

``` r
fig_marker_by_timepoint(mfi, smap, outfile, max_pops = 8L)
```

## Arguments

- mfi:

  population_marker_mfi table (median_asinh per sample x population x
  marker).

- smap:

  Sample map carrying `timepoint`.

- outfile:

  Destination PNG.

- max_pops:

  Population cap, largest first.
