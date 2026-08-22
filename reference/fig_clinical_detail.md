# Per-population detail against one clinical variable

The data behind one column of
[`fig_clinical_heatmap`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_heatmap.md):
a scatter against a numeric variable, or a box plot against a
categorical one, one point per sample.

A rank correlation of 0.8 on nine points can be a single outlier, and
the scatter is the only place that is visible. The fitted line is a
least-squares guide for the eye; the test itself was Spearman on ranks.

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

  named vector, sample_id -\> value.

- var_name:

  the variable's name, used in the title and axis.

- outfile:

  path.

- value_col:

  frequency column; defaults to the run's abundance measure.

- ncol:

  panels per row.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.
