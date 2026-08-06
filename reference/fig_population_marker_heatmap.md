# Population x marker phenotype heatmap

Population x marker phenotype heatmap

## Usage

``` r
fig_population_marker_heatmap(
  mfi,
  outfile,
  scale_by = c("marker", "none"),
  min_cells = 20L,
  annotate_expected = NULL,
  panel_label = "",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- mfi:

  population_marker_mfi table

- outfile:

  Path to write the figure to.

- scale_by:

  "marker" (row z-score, the annotation/validation view), "none" (raw
  arcsinh medians, comparable across populations within a marker but not
  between markers)

- min_cells:

  per-sample floor before a sample contributes to a cell

- annotate_expected:

  optional named list population -\> character vector of markers the
  gate definition REQUIRES to be positive. Drawn as an outline on those
  cells, turning the figure into an explicit gate audit.

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
