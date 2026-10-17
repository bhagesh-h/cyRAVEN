# Annotate clusters with the immune subset their marker profile matches

Annotate clusters with the immune subset their marker profile matches

## Usage

``` r
annotate_clusters_with_subsets(prof, subsets = NULL, pos_cut = 0.5)
```

## Arguments

- prof:

  Cluster profile: one row per cluster, with a `cluster` column and one
  `frac_pos.<marker>` column per marker, as written to
  explore_cluster_profile.csv.

- subsets:

  Named list of subset definitions; defaults to the ones the profile's
  own markers can resolve.

- pos_cut:

  Fraction of a cluster's cells that must be positive for a marker
  before the cluster counts as positive for it.

## Value

data.frame with one row per cluster: `cluster`, `subset_label`,
`subset_markers` (the requirements that were checked), `n_requirements`,
`margin` (the smallest distance from `pos_cut` across those
requirements) and `subset_alternatives`.
