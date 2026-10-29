# Heatmap of median fluorescence intensity per population and marker

Heatmap of median fluorescence intensity per population and marker

## Usage

``` r
fig_mfi_heatmap(
  summ,
  outfile,
  label_col = "population",
  by_group = FALSE,
  freq = NULL,
  freq_label = "% of cells",
  cluster_rows = FALSE,
  row_labels = NULL,
  title = "Median fluorescence intensity",
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

- by_group:

  Facet by study group instead of drawing "All samples".

- freq:

  Optional named numeric, label -\> percent, drawn as a bar beside each
  row. Ignored when `by_group`.

- freq_label:

  Axis title for that bar.

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
