# Cluster-level QC gate

cyCONDOR's idea: cluster coarsely first, then judge WHOLE CLUSTERS on
their marker profile, so no per-event cutoff is invented. Three calls,
all derived rather than hard-coded:

## Usage

``` r
explore_qc_gate(
  pos,
  med,
  sizes,
  leukocyte = NULL,
  viability = NULL,
  have_thresholds = FALSE,
  drop_dead = TRUE,
  sat_quantile = 0.99,
  sat_min_fraction = 0.7,
  seed = 42L
)
```

## Arguments

- pos:

  Positivity matrix from
  [`explore_positivity()`](https://bhagesh-h.github.io/cyRAVEN/reference/explore_positivity.md).

- med:

  Per-cluster medians, clusters x features.

- sizes:

  Cells per cluster.

- leukocyte, viability:

  Channel names, or NULL.

- have_thresholds:

  Whether `pos` was built from per-sample thresholds.

- drop_dead:

  Whether to drop the dead call.

- sat_quantile, sat_min_fraction:

  Saturation rule.

- seed:

  Integer seed.

## Value

data.frame with one row per cluster and a `call` column.

## Details

debris leukocyte marker low dead viability marker high saturated
top-percentile in most channels at once – no real cell type is, so these
are aggregates and doublets

The improvement over doing it blind: when per-sample thresholds exist,
the leukocyte and viability calls use them, so "CD45 low" means "below
this sample's own CD45 cut" rather than "in the lower of two pooled
modes".
