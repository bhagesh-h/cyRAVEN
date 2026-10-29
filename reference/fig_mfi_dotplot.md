# Dot plot: percent positive as size, scaled intensity as colour

Dot plot: percent positive as size, scaled intensity as colour

## Usage

``` r
fig_mfi_dotplot(
  summ,
  outfile,
  label_col = "population",
  cluster_rows = FALSE,
  row_labels = NULL,
  title = "Marker positivity and intensity",
  subtitle = NULL,
  dpi = 200
)
```

## Arguments

- summ:

  Output of
  [`summarise_mfi()`](https://bhagesh-h.github.io/cyRAVEN/reference/summarise_mfi.md),
  restricted to one panel.

- outfile:

  Destination PNG.

- label_col:

  Population or cluster column in `summ`.

- cluster_rows:

  See
  [`mfi_axis_order()`](https://bhagesh-h.github.io/cyRAVEN/reference/mfi_axis_order.md).

- row_labels:

  Optional named character, label -\> text shown on the axis.

- title, subtitle:

  Figure title and subtitle.

- dpi:

  Resolution.

## Value

The plot, invisibly, or NULL when there was nothing to draw.
