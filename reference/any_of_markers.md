# The markers an `any_of` group reads

A member of the group is either one condition (`CD56: above`) or a named
set of conditions that must all hold
(`CD16, not HLA-DR: {CD16: above, HLA-DR: below}`), so "CD56+, or CD16+
and HLA-DR-" can be written as one population.

## Usage

``` r
any_of_markers(a)
```

## Arguments

- a:

  an `any_of` list, or NULL.

## Value

character vector of marker names.
