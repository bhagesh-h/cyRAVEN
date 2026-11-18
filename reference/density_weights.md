# Inverse-density sampling weights in marker space

The weight is the distance to the k-th nearest neighbour, which is a
standard non-parametric density estimate: large in sparse regions, small
in dense ones. Distances are computed on a bounded random reference
subset rather than between all pairs, which makes the cost linear in
cells rather than quadratic.

## Usage

``` r
density_weights(X, k = 20L, n_ref = 2000L, seed = 1L)
```

## Arguments

- X:

  numeric matrix of features, one row per cell

- k:

  neighbour rank used for the density estimate

- n_ref:

  reference cells the distances are measured against

- seed:

  seed for the local RNG stream, which is restored on exit

## Value

numeric vector of non-negative weights, one per row of `X`
