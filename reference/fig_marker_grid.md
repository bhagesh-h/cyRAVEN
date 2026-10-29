# Marker-expression grid over the shared embedding

WHAT: one small UMAP per marker, coloured by that marker's asinh
intensity. WHY: this is the check that the embedding is biologically
structured – each lineage marker should light up a coherent region. If
markers are smeared uniformly, the embedding is driven by noise (which
is exactly what the ungated template produced).

## Usage

``` r
fig_marker_grid(
  cells,
  markers,
  outfile,
  panel_label = "",
  ncol = NULL,
  dpi = 300,
  panel_size = 2.5,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  data.frame with umap_1/umap_2 plus one column per marker.

- markers:

  character vector of marker column names to render.

- outfile:

  Path to write the figure to.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- ncol:

  Number of panel columns; NULL computes one that keeps the canvas
  roughly square.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `300`.

- panel_size:

  The panel size. Default `2.5`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
