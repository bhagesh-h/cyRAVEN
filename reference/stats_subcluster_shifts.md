# Per (cluster, subcluster, marker, study): how far is this study from reference

The figure shows the shift; this makes it sortable. Reports BOTH
deviation modes named above – the occupancy of the subcluster and the
intensity of the marker within it – so "which markers and which cells
differ" is a sort on a column rather than an eyeball over 288 panels.

## Usage

``` r
stats_subcluster_shifts(
  cells,
  markers,
  subcluster,
  group_col = "cohort",
  reference = NULL,
  min_cells = 20L
)
```

## Arguments

- cells:

  Embedded cell table: one row per cell, carrying `umap_1`/`umap_2`,
  `sample_id`, `population_label` and one column per marker.

- markers:

  Character vector of marker names to use.

- subcluster:

  Integer vector of subcluster assignments, one per cell.

- group_col:

  Name of the column holding the biological grouping (cohort). Default
  `"cohort"`.

- reference:

  The group every other group is compared against.

- min_cells:

  Minimum cells a sample must contribute before it is used. Default
  `20L`.

## Details

Effect size is the median difference in asinh units plus Cliff's delta,
which is bounded �-1 to 1 and rank-based: at these per-subcluster n a
difference in means is dominated by outliers, and these distributions
are not normal.
