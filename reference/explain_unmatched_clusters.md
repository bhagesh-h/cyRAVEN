# Propose gating strategies for clusters the spec does not describe

WHAT IT IS FOR.
[`cluster_gate_agreement()`](https://bhagesh-h.github.io/cyRAVEN/reference/cluster_gate_agreement.md)
identifies clusters whose cells carry no population label, or whose
dominant label covers only part of them. Those are the clusters worth
explaining: the spec has nothing to say about them, so there is no gate
to inspect and no threshold to blame. This runs
[`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)
on each and writes out a strategy a cytometrist can draw.

## Usage

``` r
explain_unmatched_clusters(
  cells,
  cluster,
  features,
  agreement,
  max_clusters = 4L,
  min_cells = 200L,
  purity_max = 80,
  ...
)
```

## Arguments

- cells:

  embedding cell table carrying the marker columns

- cluster:

  integer vector from
  [`run_unsupervised_clusters()`](https://bhagesh-h.github.io/cyRAVEN/reference/run_unsupervised_clusters.md)

- features:

  marker columns to gate on

- agreement:

  the list returned by
  [`cluster_gate_agreement()`](https://bhagesh-h.github.io/cyRAVEN/reference/cluster_gate_agreement.md)

- max_clusters:

  ceiling on how many strategies to derive

- min_cells:

  smallest cluster worth explaining

- purity_max:

  only explain clusters no cleaner than this (percent)

- ...:

  passed to
  [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)

## Value

list(summary, polygons, strategies), or NULL

## Details

WHY IT ONLY TAKES THE FIRST FEW. Each strategy costs a handful of small
optimisations, which is cheap, but a report proposing fifteen new gates
is not a finding, it is a second problem. `max_clusters` keeps the
output to the clusters with the strongest claim to being real
populations – the largest undescribed ones – and the count of what was
left out is logged rather than passed over in silence.
