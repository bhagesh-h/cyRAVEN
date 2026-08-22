# Between-group differences for every population in one figure

One point per population: the effect on the x-axis, the evidence on the
y-axis. `group_comparison.png` gives each population a panel and shows
every sample, which is the right figure for reading one population
carefully and the wrong one for finding which population to read – a
dozen panels have to be compared by eye, and the effect sizes never
appear on it at all. This is the shortlist.

## Usage

``` r
fig_group_volcano(
  stats,
  outfile,
  effect = c("cliff", "log2fc"),
  p_source = c("raw", "BH"),
  label_n = 8L,
  title_noun = "Population abundance",
  dpi = 200,
  colors = fcs_colors()
)
```

## Arguments

- stats:

  output of
  [`stats_group_comparison`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md).

- outfile:

  path.

- effect:

  `"cliff"` (default) or `"log2fc"`.

- p_source:

  `"raw"` or `"BH"` – which p-value on the y-axis.

- label_n:

  how many of the strongest populations to name on the figure.

- title_noun:

  the subject of the title, e.g. "Population abundance".

- dpi:

  resolution.

- colors:

  Named list of colours; defaults to the package palette. See
  [`fcs_colors`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

## Value

The ggplot object, invisibly.

## Details

Cliff's delta is the default x-axis rather than a fold change. A fold
change of medians is unbounded and, at these group sizes, dominated by
whichever sample sits at the median; a population whose reference median
is near zero produces a fold change of 40 that means nothing. Cliff's
delta is bounded `[-1, 1]`, is the effect size the rank test corresponds
to, and is directly readable: 0.5 means the comparison group was higher
in three quarters of the cross-sample pairs.

Two horizontal lines. The dashed one is p = 0.05. The dot-dash one is
the smallest p the design can produce at all (see `rank_sum_p_floor`);
when it sits below the 0.05 line the design can reach significance, and
when it sits above, no population can, and an empty upper region says
nothing about biology.

That floor is drawn per comparison rather than once for the figure,
because each comparison group has its own size and therefore its own
attainable minimum: 3 against 6 cannot go below 0.024 while 6 against 6
reaches 0.0022. With three or more study groups the figure facets by
comparison and each panel carries its own line.

## See also

[`fig_group_comparison`](https://bhagesh-h.github.io/cyRAVEN/reference/fig_group_comparison.md)
for one panel per population with every sample drawn,
[`stats_group_comparison`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md)
for the table behind both.
