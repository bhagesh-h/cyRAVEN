# Train a self-organising map on the marker matrix

WHAT: batch-SOM training. Each epoch assigns every cell to its
best-matching node, then moves every node toward the weighted mean of
the cells assigned near it, with the neighbourhood radius shrinking over
epochs.

## Usage

``` r
som_train(X, xdim = 10L, ydim = 10L, epochs = 10L, seed = 42L)
```

## Arguments

- X:

  numeric matrix, cells x markers, already scaled

- xdim, ydim:

  SOM grid dimensions

- epochs:

  training epochs

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `42L`.

## Value

list(codes = nodes x markers prototypes, mapping = node per cell)

## Details

WHY BATCH AND NOT ONLINE: batch SOM is deterministic given the
initialisation, order-independent (an online SOM's result depends on the
order cells arrive, which here is the order files were read), and
vectorises. Reproducibility of cluster identity between runs is a hard
requirement in this pipeline – the whole subcluster-lettering machinery
exists for that reason.
