# Well-separated colours for a small set of categories

WHY NOT `population_colours()`. That spreads the population palette by
index, and the population palette holds two greens: at three levels it
returns red, green and spring-green, whose last two pair is the closest
in the whole palette. `study_palette` is the set chosen for exactly this
job – few categories, each about a sixth of the hue wheel from the last
– so batches and clinical strata come out as distinguishable as cohorts
already do.

## Usage

``` r
clin_cat_palette(levels, colors = fcs_colors())
```

## Arguments

- levels:

  category labels.

- colors:

  palette.
