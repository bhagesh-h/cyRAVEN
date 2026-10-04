# Per-population detail against one clinical variable

The heatmap says which associations are worth looking at; this shows the
data behind them, because a rank correlation of 0.8 on nine points can
be one outlier and the scatter is the only place that is visible.

## Usage

``` r
fig_clinical_detail(
  freq,
  clin_values,
  var_name,
  outfile,
  value_col = NULL,
  ncol = 4L,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- freq:

  population_frequencies table.

- clin_values:

  named vector sample_id -\> value.

- var_name:

  the variable's name.

- outfile:

  path.

- value_col:

  frequency column.

- ncol:

  panels per row.

- dpi:

  resolution.

- colors:

  palette.
