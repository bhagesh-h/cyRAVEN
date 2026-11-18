# Write a conformance baseline from an accepted run

Summaries, not cells: the file holds per-marker threshold location and
spread, the fallback rate, per-population frequency location and spread,
and the specification text. It carries no event-level data and no
patient data, so it can be version-controlled next to the config it
describes.

## Usage

``` r
write_spec_baseline(
  path,
  thr_all,
  freq,
  spec,
  opt = list(),
  cofactors = list(),
  transforms = list()
)
```

## Arguments

- path:

  destination `.rds`

- thr_all:

  the thresholds table (`thresholds_used.csv` in memory)

- freq:

  the frequency table

- spec:

  population specification

- opt:

  the resolved option list

- cofactors:

  per-panel cofactors

- transforms:

  per-panel transforms from
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md);
  their method and parameters are stored so a later run's cuts can be
  put on this run's scale (see
  [`specification_conformance()`](https://bhagesh-h.github.io/cyRAVEN/reference/specification_conformance.md)).

## Value

`path`, invisibly
