# Draw a subsample that preserves sparse regions of marker space

Draw a subsample that preserves sparse regions of marker space

## Usage

``` r
draw_subsample_rare(gated_idx, X, n_take, seed = 1L, k = 20L)
```

## Arguments

- gated_idx:

  Integer indices of gated cells.

- X:

  feature matrix for those cells, in the same row order

- n_take:

  Number of cells to draw.

- seed:

  Random seed. The stream is restored afterwards.

- k:

  neighbour rank for
  [`density_weights()`](https://bhagesh-h.github.io/cyRAVEN/reference/density_weights.md)

## Value

list(idx = drawn indices, weight = sampling weight per drawn cell)
