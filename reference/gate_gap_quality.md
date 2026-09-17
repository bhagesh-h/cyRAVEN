# How deep a density gap does a candidate cut sit in

How deep a density gap does a candidate cut sit in

## Usage

``` r
gate_gap_quality(x, cut, adjust = 2)
```

## Arguments

- x:

  transformed values for one marker within the parent gate.

- cut:

  candidate threshold.

- adjust:

  kernel bandwidth multiplier. Default `2`.

## Value

a number in `(0, 1)`: 1 is a cut at the floor of a clean trough, 0 a cut
through the middle of a single mode. `-Inf` for a degenerate split.
