# Cluster abundance as a cell number, beside the same cluster as a share

WHY BOTH PANELS AND NOT JUST THE ABSOLUTE ONE. The pair is the argument.
A cluster that moves on the right panel and not the left is a
compartment that changed size, which no share can show; one that moves
on the left and not the right is a redistribution at constant size.
Shown apart, either panel invites the wrong reading of the other.

## Usage

``` r
fig_absolute_vs_share(
  ab,
  outfile,
  group_of = NULL,
  max_clusters = 12L,
  share_col = NULL,
  panel_label = ""
)
```

## Arguments

- ab:

  Abundance table carrying a share column and cells_absolute.

- outfile:

  Destination PNG.

- group_of:

  Named vector, sample_id -\> group.

- max_clusters:

  Facet cap, largest populations first.

- share_col:

  The share column. Explore writes `pct_of_gated`, the declared run
  writes `pct_of_cd45_pos`; resolved rather than hard-coded so one
  figure serves both, and so a table carrying neither is skipped rather
  than drawn against a column of NULL.
