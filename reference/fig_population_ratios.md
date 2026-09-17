# Population-ratio figure – same layout as fig_group_comparison(), one panel per ratio defined in the config's `ratios:` block

Population-ratio figure – same layout as fig_group_comparison(), one
panel per ratio defined in the config's `ratios:` block

## Usage

``` r
fig_population_ratios(
  rt,
  outfile,
  group_of = NULL,
  stats = NULL,
  reference = NULL,
  p_source = c("raw", "BH"),
  ncol = 3L,
  panel_label = "",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- rt:

  The rt.

- outfile:

  Path to write the figure to.

- group_of:

  Named character vector mapping sample_id to group label.

- stats:

  Statistics table from the matching stats\_ function, used to annotate
  the figure.

- reference:

  The group every other group is compared against.

- p_source:

  Which p-value the figure annotates: raw or BH-adjusted. Default
  `c("raw", "BH")`.

- ncol:

  Number of panel columns; NULL computes one that keeps the canvas
  roughly square. Default `3L`.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `300`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
