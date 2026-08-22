# Ordered effect sizes with bootstrap intervals for one clinical variable

One row per population (or marker), ordered by effect: the point is
Spearman's rho or Cliff's delta, the bar is its 95% percentile bootstrap
interval, and the dashed line is no effect. Rows whose interval clears
the line are the ones a follow-up cohort should be powered on.

A three-level variable is drawn on epsilon-squared instead, which is
unsigned, so its axis starts at zero and no interval is drawn: there is
no direction to put an interval around.

## Usage

``` r
fig_clinical_forest(
  assoc,
  key_col,
  var_name,
  outfile,
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- assoc:

  data frame from
  [`stats_clinical_association`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_clinical_association.md).

- key_col:

  "population" or "marker".

- var_name:

  which clinical variable to draw.

- outfile:

  path.

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.

## Details

Ordered by effect and not by p-value. At these sample sizes the two
orderings are close but not identical, and where they disagree the
effect is the one worth reading: a p-value carries how many samples
there were as much as how large the difference is.

## See also

[`fig_clinical_heatmap`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_heatmap.md)
for every variable at once,
[`fig_clinical_detail`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_clinical_detail.md)
for the points behind one variable.
