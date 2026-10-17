# Align quantiles within each cell type rather than over the whole file

WHY THIS EXISTS. Whole-file alignment moves every cell of a sample by
the same map. That is only right when the batch effect is also the same
for every cell, and it usually is not: a shift in a detector moves a
bright population and a dim one by different amounts, so one map fitted
to the pooled distribution over-corrects one and under-corrects the
other. Worse, a file's pooled density peaks move with the biology as
well as with the batch, so a per-file map can remove the difference it
was meant to preserve.

## Usage

``` r
align_quantiles_by_cluster(tmat, batch, markers, k = 10L, seed = 42L)
```

## Arguments

- tmat:

  numeric matrix, cells x markers, already transformed

- batch:

  batch label per cell

- markers:

  markers to align

- k:

  number of cell-type clusters to fit

- seed:

  RNG seed, passed through so a run reproduces

## Value

list(tmat, clusters, method)

## Details

Fitting one map per cell type is the published answer to that, and is
what CytoNorm does (Van Gassen et al. 2020; CytoNorm 2.0, Quintelier et
al., Cytometry A, 2025). It is implemented here rather than by calling
that package for two reasons, both recorded so nobody has to rediscover
them: CytoNorm installs only from GitHub, which breaks the
dated-snapshot guarantee the Docker image rests on, and its API is
file-based (`QuantileNorm.train()` takes FCS paths and
`QuantileNorm.normalize()` writes new FCS to disk), while correction
here happens on an in-memory matrix that has already been read,
transformed and gated.

The clustering comes from
[`run_unsupervised_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_unsupervised_clusters.md),
which is already seeded, already stream-safe, and already falls back to
a built-in SOM when FlowSOM is absent. So this path adds no dependency
the package did not already have.

CELLS WITHOUT A CLUSTER. Incomplete cases get no label. They are aligned
whole-file rather than left alone, because a matrix where some rows are
corrected and others are not is worse than either choice applied
consistently. The count is logged.
