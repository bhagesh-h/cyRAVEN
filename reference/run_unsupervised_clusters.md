# Unsupervised clustering of the embedded cells

WHY IT CLUSTERS THE MARKER MATRIX AND NOT THE UMAP COORDINATES: UMAP is
a visualisation. Its distances are not metric, it does not preserve
density, and the gaps between its islands are partly an artefact of
min_dist. Clustering it would give clusters of the PICTURE. Every
established tool – FlowSOM and Phenograph among them – clusters the
high-dimensional space and uses the embedding only to display the
result, which is what this does.

## Usage

``` r
run_unsupervised_clusters(
  cells,
  markers,
  n_clusters = 12L,
  grid = 10L,
  epochs = 10L,
  scale_method = "robust",
  seed = 42L,
  max_cells = 200000L
)
```

## Arguments

- cells:

  embedding cell table

- markers:

  marker columns to cluster on (normally the UMAP feature set)

- n_clusters:

  number of metaclusters

- grid:

  SOM grid side length; grid^2 nodes

- epochs:

  Number of SOM training epochs. Default `10L`.

- scale_method:

  matched to run_umap()'s default so clusters and the embedding see the
  same geometry

- seed:

  Random seed. The RNG stream is restored afterwards, so this cannot
  perturb later sampling. Default `42L`.

- max_cells:

  Ceiling on the number of cells used. Default `200000L`.

## Value

list(cluster = integer per cell, method, codes, node)
