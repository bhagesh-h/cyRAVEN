# UMAP coloured by a continuous variable (marker intensity or covariate)

UMAP coloured by a continuous variable (marker intensity or covariate)

## Usage

``` r
umap_continuous(
  df,
  colour_by,
  title = NULL,
  subtitle = NULL,
  point_size = NULL,
  alpha = NULL,
  limits = NULL,
  legend_name = "asinh",
  compact_bar = FALSE,
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

- limits:

  Length-2 numeric vector giving the scale limits.

- legend_name:

  The legend name. Default `"asinh"`.

- compact_bar:

  The compact bar. Default `FALSE`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
