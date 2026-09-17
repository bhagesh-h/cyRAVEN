# Functional-marker positivity figure

Delegates to fig_group_comparison(), exactly as
fig_population_frequencies() does, so a monocyte HLA-DR shift gets the
same bar/median/whisker/points treatment and the same Wilcoxon brackets
as a population-abundance shift – before this figure existed,
functional_markers.csv was written but never plotted, so a functional
difference between cohorts was visible only as rows in a CSV.

## Usage

``` r
fig_functional_markers(
  fx,
  outfile,
  group_of = NULL,
  stats = NULL,
  reference = NULL,
  p_source = c("raw", "BH"),
  ncol = 3L,
  panel_label = "",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- fx:

  functional_markers table (sample_id, block, population, marker,
  pct_positive, is_control). Population x marker becomes the panel key,
  because a positivity threshold only means something inside the
  population it was derived for.

- outfile:

  Path to write the figure to.

- group_of:

  optional sample_id -\> group map; NULL pools every sample into "All
  samples" (mirrors fig_population_frequencies()).

- stats:

  Statistics table from the matching stats\_ function, used to annotate
  the figure.

- reference:

  The group every other group is compared against.

- p_source:

  Which p-value the figure annotates: raw or BH-adjusted. Default
  `c("raw", "BH")`.

- ncol:

  Number of panel columns; NULL computes one that keeps the canvas
  roughly square. Default `3L`.

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
