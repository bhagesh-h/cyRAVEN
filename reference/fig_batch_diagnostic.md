# Batch diagnostic figure: the embedding coloured by batch, plus the LISI distribution against its permutation null

WHY BOTH PANELS: the scatter is what a reader will look at anyway and it
shows WHERE any batch structure sits; the LISI distribution is what says
whether the pattern the eye finds is more than chance. Neither alone is
enough – a scattered-looking plot can still be significantly structured
at these n, and a significant score with no visible pattern is usually
one small batch.

## Usage

``` r
fig_batch_diagnostic(
  cells,
  report,
  outfile,
  batch_col,
  panel_label = "",
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`,
  `sample_id`, `population_label` and one column per marker.

- report:

  The report.

- outfile:

  Path to write the figure to.

- batch_col:

  Column naming the acquisition batch.

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
