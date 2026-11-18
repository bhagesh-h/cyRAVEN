# Mean silhouette width for a k-means partition

Computed on a subsample: silhouette is O(n^2) in distances, and the mean
over a few thousand cells is stable to well within the differences
between adjacent k that this is used to resolve.

## Usage

``` r
mean_silhouette(X, k, sample_n = 1500L, seed = 42L, nstart = 5L)
```

## Arguments

- X:

  Numeric matrix, cells x features, already scaled.

- k:

  Neighbourhood size, or number of clusters, depending on the function.

- sample_n:

  Number of cells to subsample. Default `1500L`.

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `42L`.

- nstart:

  Number of k-means restarts. Default `5L`.
