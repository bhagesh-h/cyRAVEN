# Scale each feature before embedding or clustering

One implementation for
[`run_umap()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_umap.md)
and
[`run_unsupervised_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_unsupervised_clusters.md),
so clusters and coordinates are computed in the same space.

## Usage

``` r
scale_features(mat, method = "range")
```

## Arguments

- mat:

  numeric matrix, cells x features.

- method:

  "none", "range" (0.5th-99.5th percentile range), "robust" (median/MAD)
  or "zscore".

## Value

list(mat, params); `params` carries `center` and `scale` so a saved
model can be applied to later data in the same space, or NULL for
"none".
