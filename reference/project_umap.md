# Project new cells into a saved embedding

Refuses rather than improvises when the feature sets differ: a model
trained on 11 markers cannot place cells described by 9, and quietly
filling the gap with zeros would produce coordinates that look plausible
and mean nothing.

## Usage

``` r
project_umap(saved, mat, n_threads = 1L)
```

## Arguments

- saved:

  A model loaded by
  [`load_umap_model()`](https://bhagesh-h.github.io/cyRAVEN/reference/load_umap_model.md).

- mat:

  Numeric matrix, cells x features.

- n_threads:

  Number of threads. Default `1L`.
