# Unsupervised-clustering figure: the embedding by cluster, next to the same embedding by gate label

WHY SIDE BY SIDE AND NOT TWO FILES: the entire value is the comparison.
A cluster panel alone is a pretty picture; against the gate panel it is
an audit, and a reader can see in one glance which islands the spec
named and which it missed.

## Usage

``` r
fig_unsupervised_clusters(
  cells,
  cluster,
  outfile,
  agreement = NULL,
  panel_label = "",
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`,
  `sample_id`, `population_label` and one column per marker.

- cluster:

  Integer vector of cluster assignments, one per cell.

- outfile:

  Path to write the figure to.

- agreement:

  The agreement.

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
