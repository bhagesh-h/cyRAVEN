# Figure: which declared population each unsupervised cluster corresponds to

The correspondence matrix of (2) above: clusters as rows, declared
populations as columns, fill is the share of the cluster. Rows are
grouped by the called identity so the answer is readable off the row
labels alone, and the catch-all keeps its own column at the right so an
undescribed cluster is visibly undescribed.

## Usage

``` r
fig_cluster_identity(
  ct,
  ident,
  outfile,
  panel_label = "",
  colors = fcs_colors()
)
```

## Arguments

- ct:

  Long cross-tabulation; see
  [`cluster_identity_table()`](https://bhagesh-h.github.io/cyRAVEN/reference/cluster_identity_table.md).

- ident:

  Output of
  [`cluster_identity_table()`](https://bhagesh-h.github.io/cyRAVEN/reference/cluster_identity_table.md);
  computed if NULL.

- outfile:

  Destination PNG.

- panel_label:

  Panel name for the title; empty for none.

- colors:

  Palette.
