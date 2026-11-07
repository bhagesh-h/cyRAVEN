# Sample row indices reproducibly

Sample row indices reproducibly

## Usage

``` r
draw_subsample(gated_idx, n_take, seed = 1L)
```

## Arguments

- gated_idx:

  Integer indices of gated cells.

- n_take:

  Number of cells to draw.

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `1L`.

## Value

integer vector of selected indices into the per-sample gated index set
