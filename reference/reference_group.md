# The cohort every comparison is drawn against

The reference group is figure-wide state: it decides which cohort is
plotted first, which one is drawn unfilled, which one keeps the
reference colour, and which one every statistical comparison is made
against. Setting it once is what keeps those four answers consistent
across a result set.

## Usage

``` r
reference_group()

set_reference_group(group)
```

## Arguments

- group:

  The group.

## Value

`reference_group()` returns the active label or `NULL`.
`set_reference_group()` returns the previous value invisibly.

## Examples

``` r
old <- set_reference_group("Healthy controls")
reference_group()
#> [1] "Healthy controls"
set_reference_group(old)
```
