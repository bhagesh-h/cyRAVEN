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
  z_fail = 6,
  transforms = NULL
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

- transforms:

  this run's per-panel transforms from
  [`make_transform()`](https://bhagesh-h.github.io/cyRAVEN/reference/make_transform.md),
  or NULL to compare on the analysis scale only when it is unchanged.

## Value

list(markers, populations, spec_changes, summary)

## Details

THE SCALE IS FLOORED, for the reason it is floored in the within-run
check: a marker whose baseline samples happened to agree closely would
otherwise divide by near zero and fail on a difference of no
consequence.

THE CUTS ARE COMPARED IN INSTRUMENT UNITS, on the baseline's scale.
Thresholds are expressed on the analysis scale, and that scale is
re-derived from each run's data: the arcsinh cofactor of the same 30
files moved from 261 to 403 between two versions of this package, which
moves every cut on the transformed axis without moving it on the
instrument's. Comparing the transformed values across runs therefore
measured the change of ruler. Each marker's median cut is carried back
to instrument units by this run's transform and forward by the
baseline's, then compared there against the baseline's spread. A
baseline written before transforms were stored is rebuilt from its
arcsinh cofactor; where either transform cannot be rebuilt, a changed
transform makes the comparison "not comparable", as before.
