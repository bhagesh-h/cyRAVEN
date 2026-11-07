# Which FMO file, if any, controls a given marker for a given sample

A control with no `control_group` applies everywhere. A control with one
applies only to samples sharing it, and a sample whose group has no
control for that marker falls back to a group-less control if one
exists.

## Usage

``` r
fmo_for_sample(fmo_map, sample_id, marker, group_of = NULL)
```

## Arguments

- fmo_map:

  from
  [`parse_fmo_map()`](https://bhagesh-h.github.io/cyRAVEN/reference/parse_fmo_map.md)

- sample_id:

  the sample being gated

- marker:

  the marker being thresholded

- group_of:

  named character vector mapping sample_id to control group

## Value

the controlling sample_id, or NA
