# Assign each unsupervised cluster the declared population it corresponds to

Assign each unsupervised cluster the declared population it corresponds
to

## Usage

``` r
cluster_identity_table(ct)
```

## Arguments

- ct:

  Long cross-tabulation with columns `cluster`, `population`, `cells`,
  and optionally `phenotype`: one row per cluster x population pair.

## Value

data.frame, one row per cluster, carrying the called identity, the
precision/recall/F1 behind it, the runner-up, and the catch-all share.
