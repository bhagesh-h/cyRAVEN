# Split per-cluster medians into a low and a high mode

2-means on one column of per-cluster medians. Used only when no
per-sample threshold is available, i.e. standalone explore. Returns the
indices of the HIGH group and the boundary between them, or NULL when
the split is not meaningful (fewer than 3 clusters, or no separation).

## Usage

``` r
two_mode_split(v, seed = 42L)
```

## Arguments

- v:

  Numeric vector, one median per cluster.

- seed:

  Integer seed.

## Value

list(high = integer, cut = numeric) or NULL.
