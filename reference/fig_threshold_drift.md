# Threshold-drift figure: per-sample cut for every marker, coloured by cohort

Threshold-drift figure: per-sample cut for every marker, coloured by
cohort

## Usage

``` r
fig_threshold_drift(
  thr,
  outfile,
  group_of,
  stats = NULL,
  reference = NULL,
  panel_label = "",
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- thr:

  The thr.

- outfile:

  Path to write the figure to.

- group_of:

  Named character vector mapping sample_id to group label.

- stats:

  Statistics table from the matching stats\_ function, used to annotate
  the figure.

- reference:

  The group every other group is compared against.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `200`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
