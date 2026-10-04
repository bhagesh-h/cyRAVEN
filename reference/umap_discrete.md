# UMAP coloured by a discrete variable

UMAP coloured by a discrete variable

## Usage

``` r
umap_discrete(
  df,
  colour_by,
  title = NULL,
  subtitle = NULL,
  point_size = NULL,
  alpha = NULL,
  legend_rows = NULL,
  palette = NULL,
  colors = fcs_colors()
)
```

## Arguments

- df:

  The df.

- colour_by:

  The colour by.

- title:

  The title.

- subtitle:

  The subtitle.

- point_size:

  The point size.

- alpha:

  Significance threshold.

- legend_rows:

  The legend rows.

- palette:

  The palette.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
