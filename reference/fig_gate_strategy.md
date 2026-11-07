# Draw a proposed gating strategy, one panel per gate

WHY TARGETS ARE DRAWN LAST AND NON-TARGETS IN GREY: the question the
reader has is "does this polygon contain the population and exclude the
rest", and that is answered by seeing the population against the
background it was separated from. Colouring both categories equally
makes a dense non-target cloud hide the very cells the gate is about.

## Usage

``` r
fig_gate_strategy(
  X,
  y,
  strategy,
  outfile,
  label = "cluster",
  colors = fcs_colors(),
  dpi = 300
)
```

## Arguments

- X:

  marker matrix used to derive the strategy

- y:

  0/1 target indicator

- strategy:

  the list returned by
  [`explain_cluster()`](https://bhagesh-h.github.io/cyRAVEN/reference/explain_cluster.md)

- outfile:

  path to write to

- label:

  name of the population being explained

- colors:

  Named list of colours; defaults to the package palette. Default
  [`fcs_colors()`](https://bhagesh-h.github.io/cyRAVEN/reference/fcs_colors.md).

- dpi:

  Raster resolution. Default `300`.

## Value

the assembled plot, invisibly

## Details

The subtitle carries the held-out metrics, not the in-sample ones,
because a figure is where an optimistic number does the most damage.
