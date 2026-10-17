# Figures for explore mode

Figures for explore mode

## Usage

``` r
explore_figures(cells, feats, prof, ex_dir, tag = "", group_col = NULL)
```

## Arguments

- cells:

  data.frame with umap_1, umap_2, cluster, sample_id and the feature
  columns.

- feats:

  Feature names.

- prof:

  Cluster profile from the orchestrator, carrying `phenotype`.

- ex_dir:

  Output directory.

- tag:

  Panel suffix.

- group_col:

  Name of the group column, or NULL.

## Value

Character vector of file names written.
