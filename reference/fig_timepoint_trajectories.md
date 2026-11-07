# Per-patient trajectories across timepoints, one panel per population

The paired view: each line is one patient followed through their visits,
so a consistent within-patient direction is visible even when the cohort
spread swamps it. The cohort median is drawn over the lines, not in
place of them.

## Usage

``` r
fig_timepoint_trajectories(
  freq,
  smap,
  outfile,
  value_col = "pct_of_cd45_pos",
  log_scale = TRUE,
  max_pops = Inf
)
```

## Arguments

- freq:

  population_frequencies table.

- smap:

  Sample map carrying `timepoint` and `patient_id`.

- outfile:

  Destination PNG.

- value_col:

  Column to plot.

- log_scale:

  Use a log10 y axis. Population sizes span orders of magnitude, so a
  fold change is the comparable quantity.

- max_pops:

  Facet cap, largest populations first. `Inf` by default: every
  population a specification declares was asked for by name, and a cap
  ranked by abundance removes the rarest first, which on a sepsis panel
  is exactly the gamma delta and NK compartments. A caller drawing
  dozens of unsupervised clusters passes a finite cap.
