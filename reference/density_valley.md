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
  floor = -Inf,
  prefer = c("deepest", "lowest"),
  min_lower_frac = 0.01,
  details = FALSE,
  noise_alpha = 0.05,
  ceiling = Inf
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

- floor:

  no cut is accepted at or below this value. Callers pass the
  transformed value of zero signal, `transform$fn(0, marker)`: a stained
  population cannot sit below zero. Without it, the small bump that
  compensated events form below zero, next to the negative peak, scores
  as a deep separation, and the cut lands under the negative peak, so
  every cell is called positive. On one cohort that put 87% of CD45+
  cells "Vd1+" and emptied every population defined as Vd1-. Default
  `-Inf`, no floor.

- prefer:

  which qualifying valley to return. `"deepest"`, the default, takes the
  best-separated pair of modes. `"lowest"` takes the qualifying valley
  nearest the bottom of the scale, which is what a gate meant to remove
  only NEGATIVE events needs: in whole blood the CD45 distribution has a
  negative mode, a dim mode (granulocytes) and a bright one
  (lymphocytes), the dim/bright valley is often the deepest, and taking
  it removed the granulocytes from the CD45+ parent. Used for the CD45
  gate. Under `"lowest"` a valley also needs `min_lower_frac` of the
  events below it, so a ripple in the tail of the negative mode does not
  count as a separation.

- min_lower_frac:

  see `prefer`. Default `0.01`.

- details:

  when TRUE return a list carrying the cut plus the relative depth of
  the valley and why no cut was found, instead of the bare cut. The
  default returns exactly what it always has: one number or `NA_real_`.
  [`threshold_uncertainty()`](https://bhagesh-h.github.io/cyRAVEN/reference/threshold_uncertainty.md)
  uses the list form, because how deep a valley is determines how far
  the cut can move under resampling, and recomputing the smoothed
  histogram elsewhere would let the two drift apart.

- noise_alpha:

  the family-wise rate at which counting noise alone may pass for a
  slope. A valley counts only where the density falls significantly from
  the first peak to it and rises significantly from it to the second,
  each judged against Poisson noise with a two-sided Bonferroni quantile
  over the non-overlapping smoothing windows (a SiZer-style test,
  Chaudhuri and Marron 1999). Relative depth alone does not see sample
  size: on single modes of 3,000 simulated events it found a valley in
  81 of 150 samples, most of them below 95% of the events, and on a flat
  shoulder of a few events per bin it split the shoulder. Default
  `0.05`; `NA` turns the test off.

- ceiling:

  no cut is accepted at or above this value. Default `Inf`; the
  viability step passes the dye's dominant mode to look for a minimum
  below it (see
  [`apply_gate_hierarchy()`](https://bhagesh-h.github.io/cyRAVEN/reference/apply_gate_hierarchy.md)).
