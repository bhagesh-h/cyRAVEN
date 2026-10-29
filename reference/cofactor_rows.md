# The events the transform is estimated from: scatter-gated single cells

Debris is most of the events in many files, two thirds on one cohort,
and it scatters little and fluoresces less, so a background estimated
over every event is the debris's and is narrower than the cells'. The
scatter and singlet gates need no fluorescence transform, so they can
run first and the cofactor can be estimated on the cells it will be
applied to.

## Usage

``` r
cofactor_rows(rd)
```

## Arguments

- rd:

  A list from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md).

## Value

logical vector over the events, or NULL when the gates cannot run.
