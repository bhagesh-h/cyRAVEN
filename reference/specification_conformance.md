# Test this run against a baseline

Compares each marker's threshold location against where the baseline put
it, scaled by the baseline's own spread for that marker, so the
tolerance adapts to how variable the marker legitimately is. The same
comparison is made for population frequencies.

## Usage

``` r
specification_conformance(
  thr_all,
  freq,
  spec,
  baseline,
  transform = "arcsinh",
  z_qualify = 3.5,
  z_fail = 6
)
```

## Arguments

- thr_all:

  this run's thresholds table

- freq:

  this run's frequency table

- spec:

  this run's population specification

- baseline:

  the list from
  [`read_spec_baseline()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_spec_baseline.md)

- transform:

  this run's transform name

- z_qualify:

  robust z above which a marker is qualified

- z_fail:

  robust z above which a marker fails

## Value

list(markers, populations, spec_changes, summary)

## Details

THE SCALE IS FLOORED, for the reason it is floored in the within-run
check: a marker whose baseline samples happened to agree closely would
otherwise divide by near zero and fail on a difference of no
consequence.

A CHANGED TRANSFORM INVALIDATES EVERY THRESHOLD COMPARISON, because
thresholds are expressed on the analysis scale. That case is reported as
its own verdict rather than as dozens of drifted markers.
