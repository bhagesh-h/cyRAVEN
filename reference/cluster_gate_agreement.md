# Cross-tabulate unsupervised clusters against gate labels

Cross-tabulate unsupervised clusters against gate labels

## Usage

``` r
cluster_gate_agreement(cells, cluster, other_pattern = "^Other|unclassified")
```

## Arguments

- cells:

  embedding cell table with population_label

- cluster:

  integer vector from run_unsupervised_clusters()

- other_pattern:

  regex identifying the catch-all label

## Value

list(per_cluster, per_population)
