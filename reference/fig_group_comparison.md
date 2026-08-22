# Grouped abundance figure: one panel per population, bar + SD + every sample

WHAT: mean bar, SD whiskers, one point per sample, significance brackets
from the reference group to each other group. Panels are lettered.

## Usage

``` r
fig_group_comparison(
  freq,
  outfile,
  group_of,
  stats = NULL,
  reference = NULL,
  p_source = c("raw", "BH"),
  measure = "auto",
  ncol = 3L,
  panel_label = "",
  dpi = 300,
  value_col = NULL,
  value_label = NULL,
  value_caveat = NULL,
  title_noun = "Population abundance",
  colors = fcs_colors()
)
```

## Arguments

- freq:

  population_frequencies table

- outfile:

  Path to write the figure to.

- group_of:

  named vector sample_id -\> group

- stats:

  output of stats_group_comparison(); NULL to skip brackets

- reference:

  The group every other group is compared against.

- p_source:

  "raw" or "BH" – which p-value the brackets display

- measure:

  Which abundance measure to use. Default `"auto"`.

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

- value_col, value_label, value_caveat:

  override the abundance measure with an arbitrary numeric column of
  `freq` (must still be one row per sample_id x population). Used to
  draw the same bar/median/whisker/points layout, with the same
  Wilcoxon-bracket machinery, for tables that share `freq`'s shape but
  not its meaning – functional-marker positivity, or a derived
  population ratio.

- title_noun:

  the figure title's subject, e.g. "Population abundance" or "Functional
  marker positivity"; " by group" is appended automatically when 2+
  groups are plotted.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).
  Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Details

WHY POINTS AND NOT JUST BARS: at n = 6-11 a bar and a whisker hide
whether a "difference" rests on one outlier or on a consistent shift.
Every sample is drawn, so the reader can see the shape of the evidence –
which is also why the bar is left unfilled for the reference group: the
fill carries group identity, not emphasis.
