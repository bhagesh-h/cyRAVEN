# Whether a CD45 cut excludes a population that is CD45-positive

Whether a CD45 cut excludes a population that is CD45-positive

## Usage

``` r
cd45_cut_splits_positive(x, cut, neg_ref, max_frac = 0.1)
```

## Arguments

- x:

  CD45 values of the gate's parent.

- cut:

  The proposed cut.

- neg_ref:

  The negative reference from
  [`cd45_negative_reference()`](https://bhagesh-h.github.io/cyRAVEN/reference/cd45_negative_reference.md).

- max_frac:

  More than this fraction of the parent between `neg_ref` and `cut`
  means the cut is splitting positive cells. Default `0.10`.

## Value

`TRUE`, `FALSE`, or `NA` when either value is missing.
