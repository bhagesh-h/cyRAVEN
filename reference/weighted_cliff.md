# Cliff's delta of `b` against `a` with per-cell weights

The sum over pairs of w_a w_b sign(b - a), over the sum of the weights'
products: P(b \> a) - P(b \< a) for weighted samples, as the unweighted
delta is for equal weights. Sorted once, so no pair is enumerated.

## Usage

``` r
weighted_cliff(a, b, wa, wb)
```

## Arguments

- a, b:

  numeric vectors.

- wa, wb:

  their non-negative weights.

## Value

numeric in `[-1, 1]`.
