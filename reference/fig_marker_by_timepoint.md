# Every marker's median intensity per population, split by timepoint

The state heatmap asks what a population expresses; this asks whether
that changed over the admission. Split rather than pooled, because
pooling the timepoints is what hides a change that happens in all
patients at once.

## Usage

``` r
fig_marker_by_timepoint(
  mfi,
  smap,
  outfile,
  max_pops = Inf,
  scale_label = "asinh intensity"
)
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

  Population cap, largest first. `Inf` by default, for the reason given
  in
  [`fig_timepoint_trajectories()`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_timepoint_trajectories.md).

- scale_label:

  Name of the scale the intensities are on, the run's `transform$label`
  ("asinh intensity", "logicle intensity"). Default `"asinh intensity"`.
