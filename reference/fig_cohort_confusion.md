# Cohort composition heatmap ("confusion matrix")

WHAT: for each population, the share of its cells contributed by each
cohort AFTER equalising every cohort to the same notional cell count.

## Usage

``` r
fig_cohort_confusion(
  cells,
  outfile,
  group_col = "cohort",
  pop_col = "population_label",
  norm_to = 1000,
  panel_label = "",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- cells:

  Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`,
  `sample_id`, `population_label` and one column per marker.

- outfile:

  Path to write the figure to.

- group_col:

  Name of the column holding the biological grouping (cohort). Default
  `"cohort"`.

- pop_col:

  Name of the column holding the population label. Default
  `"population_label"`.

- norm_to:

  Notional cell count each group is normalised to before shares are
  taken. Default `1000`.

- panel_label:

  Marker-panel name added to the figure title; empty for none. Default
  `""`.

- dpi:

  Resolution in dots per inch. May be reduced automatically to respect
  the raster ceiling; see
  [`safe_ggsave()`](https://bhagesh-h.github.io/cyRAVEN/reference/safe_ggsave.md).
  Default `300`.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY THE EQUALISATION IS THE WHOLE POINT: cohorts differ in sample count
and in acquisition depth, so raw contributions are dominated by
whichever cohort has the most cells and every population looks like it
belongs to that cohort. The normalisation – the same device used by the
cohort-confusion heatmaps elsewhere in the field, which equalise to a
fixed cell count per group – removes that, so a row that is not
near-uniform means a population is genuinely cohort-skewed.

WHY IT IS NOT A SUBSTITUTE FOR group_comparison.png: this pools cells
across donors, so it shows a pattern, not a tested effect. A population
that looks skewed here is a lead to check against the sample-level
abundance test, which is the figure that carries statistics.
