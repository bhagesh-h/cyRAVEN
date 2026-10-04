# Population frequency figure – one lettered panel per population

Delegates to fig_group_comparison() so the frequency overview and the
between-group comparison are THE SAME LAYOUT. The previous version was
one dodged bar per sample per population on a shared axis, which showed
each sample's value but no summary at all: mean, median and spread could
not be read off it, and the rare populations were invisible next to the
abundant ones because every population shared one x-axis.

## Usage

``` r
fig_population_frequencies(
  freq,
  outfile,
  panel_label = "",
  dpi = 300,
  colors = fcs_colors()
)
```

## Arguments

- freq:

  population_frequencies table; rows flagged `is_control` are DROPPED,
  because a control's percentages are fractions of a quantile- fallback
  slice rather than of a real CD45+ population – plotting them beside
  stained samples invites a false comparison. split by group with tests;
  when absent all samples form one group and the panels show the pooled
  mean/median/SD/points.

- outfile:

  Path to write the figure to.

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

Giving each population its own panel gives each its own y-scale, so a
population at 0.3% is as legible as one at 40%.
