# Between-group differences for every population in one figure

WHAT IT IS. One point per population: the effect on the x-axis, the
evidence on the y-axis. `group_comparison.png` gives each population a
panel and shows every sample, which is the right figure for reading one
population carefully and the wrong one for finding which population to
read – a dozen panels have to be compared by eye, and the effect sizes
never appear on it at all. This is the shortlist: everything far from
the centre line moved, everything high up moved consistently.

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
  [`stats_group_comparison()`](https://bhagesh-h.github.io/cyRAVEN/reference/stats_group_comparison.md).

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

  palette.

## Details

WHY CLIFF'S DELTA AND NOT A FOLD CHANGE by default. A fold change of
medians is unbounded and, at these group sizes, dominated by whichever
sample sits at the median; a population whose reference median is near
zero produces a fold change of 40 that means nothing. Cliff's delta is
bounded `[-1, 1]`, is the effect size the rank test actually corresponds
to, and is directly readable: 0.5 means that in three quarters of the
cross-sample pairs the comparison group was higher. `effect = "log2fc"`
gives the conventional axis for anyone who wants it.

WHAT THE TWO HORIZONTAL LINES MEAN. The upper is p = 0.05. The lower is
the smallest p this design can produce at all (see
[`rank_sum_p_floor()`](https://bhagesh-h.github.io/cyRAVEN/reference/rank_sum_p_floor.md));
when it sits BELOW the 0.05 line the design can reach significance, and
when it sits above, no population can, and an empty upper region says
nothing about biology.
