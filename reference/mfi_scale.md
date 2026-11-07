# Per-marker 0-1 rescaling of analysis-scale medians, for colour

Between the 1st and 99th percentile of the values being drawn, per
marker, as CATALYST scales its expression heatmaps. Scaling over every
group at once keeps one colour meaning one intensity across the faceted
figure.

## Usage

``` r
mfi_scale(d)
```

## Arguments

- d:

  Rows of
  [`summarise_mfi()`](https://bhagesh-h.github.io/cyRAVEN/reference/summarise_mfi.md)
  output for one panel.
