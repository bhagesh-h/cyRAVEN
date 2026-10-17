# Ordered effect sizes with bootstrap intervals for one clinical variable

WHAT IS ON IT. One row per population (or marker), ordered by effect:
the point is Spearman's rho or Cliff's delta, the bar is its 95%
percentile bootstrap interval, and the dashed line is no effect. Rows
whose interval clears the line are the ones a follow-up cohort should be
powered on.

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
  [`stats_clinical_association()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_clinical_association.md).

- key_col:

  "population" or "marker".

- var_name:

  which clinical variable to draw.

- outfile:

  path.

- dpi:

  resolution.

- colors:

  palette.

## Details

A three-level variable is drawn on epsilon-squared instead, which is
unsigned, so its axis starts at zero and no interval is drawn: there is
no direction to put an interval around.

WHY IT IS ORDERED BY EFFECT AND NOT BY p. At these sample sizes the
p-value ordering is close to the effect ordering but not identical, and
where they disagree the effect is the one worth reading: p carries how
many samples there were as much as how large the difference is.
