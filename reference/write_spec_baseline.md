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
  cofactors = list()
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

## Value

`path`, invisibly
