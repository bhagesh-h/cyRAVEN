# Find the deepest density valley between two modes

WHAT: smooths a histogram, locates peaks above `peak_frac` of the
maximum, and returns the valley between the peak pair maximising depth
\* min(height). Returns NA when the distribution is unimodal – the
honest answer, which then routes threshold resolution to the
unstained-control fallback. WHY: for a bimodal marker the valley is the
least-arbitrary cutoff available and adapts to each sample's own
staining intensity.

## Usage

``` r
density_valley(
  x,
  bins = 220L,
  smooth = 4,
  peak_frac = 0.02,
  min_gap_frac = 0.06,
  range_q = c(0.001, 0.999),
  min_rel_depth = 0.3,
  min_upper_frac = 0.002,
  details = FALSE
)
```

## Arguments

- x:

  A vector of values.

- bins:

  The bins. Default `220L`.

- smooth:

  The smooth. Default `4`.

- peak_frac:

  The peak frac. Default `0.02`.

- min_gap_frac:

  The min gap frac. Default `0.06`.

- range_q:

  The range q. Default `c(0.001, 0.999)`.

- min_rel_depth:

  minimum valley depth AS A FRACTION of the smaller of the two flanking
  peaks. This is the parameter that makes minority positive populations
  findable – see the note on scoring below.

- min_upper_frac:

  the upper mode must hold at least this fraction of events, so a cut
  cannot be placed out in the extreme tail where a handful of events
  form a bump.

- details:

  when TRUE return a list carrying the cut plus the relative depth of
  the valley and why no cut was found, instead of the bare cut. The
  default returns exactly what it always has: one number or `NA_real_`.
  [`threshold_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/threshold_uncertainty.md)
  uses the list form, because how deep a valley is determines how far
  the cut can move under resampling, and recomputing the smoothed
  histogram elsewhere would let the two drift apart.
