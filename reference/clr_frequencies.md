# Zero-safe centred log-ratio of each sample's population composition

WHY ZEROS NEED HANDLING: log(0) is -Inf, and a genuinely absent
population (0 cells scored) is common for rare subsets in a small
sample. Dropping those rows would drop exactly the populations most
likely to differ between cohorts.

## Usage

``` r
clr_frequencies(
  freq,
  value_col = "pct_of_cd45_pos",
  populations = NULL,
  delta_frac = 0.65
)
```

## Arguments

- freq:

  population_frequencies table

- value_col:

  the compositional column (percentages)

- populations:

  optional restriction to a closed set; defaults to every population
  present, the right closure when they partition CD45+

- delta_frac:

  Fraction of the smallest non-zero part used to replace zeros. Default
  `0.65`.

## Value

`freq` plus a `clr` column, or NULL

## Details

WHY MULTIPLICATIVE REPLACEMENT: the standard treatment for rounded zeros
in compositional data (Martin-Fernandez et al. 2003). Each zero takes a
small delta and the non-zero parts are scaled down to preserve the
total, so the result is still a composition. Delta defaults to 65% of
the smallest non-zero part observed, the usual choice when no detection
limit is declared.

WHY THE GEOMETRIC MEAN IS PER SAMPLE: the constraint is per sample – one
sample's populations sum to 100, not one population's samples. Centring
across samples instead would be a different, and wrong, operation.
