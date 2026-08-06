# Differential-state figure, drawn by the baseline's own comparison figure

WHY IT DELEGATES RATHER THAN DRAWING: fig_group_comparison() already
renders exactly this shape – one panel per unit, bar = mean, rule =
median, whiskers = SD, points = samples, Wilcoxon brackets. Reusing it
means the DS figure cannot drift from the abundance figure's layout, and
a reader who has learned to read one has learned to read the other.
Exactly the pattern fig_functional_markers() already uses.

## Usage

``` r
fig_marker_state(
  mfi,
  outfile,
  group_of = NULL,
  stats = NULL,
  reference = NULL,
  p_source = c("raw", "BH"),
  min_cells = 20L,
  ncol = NULL,
  panel_label = "",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- mfi:

  Per-sample x population x marker summary table.

- outfile:

  Path to write the figure to.

- group_of:

  Named character vector mapping sample_id to group label.

- stats:

  Statistics table from the matching stats\_ function, used to annotate
  the figure.

- reference:

  The group every other group is compared against.

- p_source:

  Which p-value the figure annotates: raw or BH-adjusted. Default
  `c("raw", "BH")`.

- min_cells:

  Minimum cells a sample must contribute before it is used. Default
  `20L`.

- ncol:

  Number of panel columns; NULL computes one that keeps the canvas
  roughly square.

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

WHY pct_positive AND NOT median_asinh IS PLOTTED: fig_group_comparison()
fixes the y axis at c(0, headroom) – deliberately, since every quantity
it was built for (percentages, concentrations, ratios) is non-negative.
Arcsinh medians run negative for dim markers, and handing them to that
scale would silently clip the bar to zero and draw a confident-looking
panel of a value that is not there. pct_positive is bounded �0 to 100,
is a genuine differential-state measure, and is what a bimodal shift
moves most. The median_asinh results are not discarded – they carry the
full test in the CSV and are the quantity the state heatmap in
fig_population_marker_heatmap() draws, where a diverging scale handles
sign correctly.
