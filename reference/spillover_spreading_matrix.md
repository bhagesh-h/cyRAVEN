# Spillover spreading matrix from single-stain controls

The published measure (Nguyen et al. 2013): for a single-stain control
of fluorochrome i, split into its positive population P and negative
population N on its primary detector, the spreading coefficient into
detector j is

## Usage

``` r
spillover_spreading_matrix(stains, min_events = 200L)
```

## Arguments

- stains:

  named list of reads from
  [`read_fcs_resolved()`](https://bhagesh-h.github.io/cyRAVEN/reference/read_fcs_resolved.md),
  compensated.

- min_events:

  fewest events each of P and N must hold. Default 200.

## Value

data.frame, one row per control and receiving detector: `control`,
`source_detector`, `receiver_detector`, `n_positive`, `n_negative`,
`delta_median`, `spread_positive`, `spread_negative`, `ss`. NULL when no
control can be used.

## Details

    SS_ij = sqrt(max(0, sigma(P, j)^2 - sigma(N, j)^2)) / sqrt(F(P, i) - F(N, i))

where sigma is the robust standard deviation (84.13th minus 50th
percentile) of the compensated values and F the median. The square root
of the intensity difference normalises out how bright the control was,
so the coefficient is a property of the fluorochrome-detector pair on
this instrument, comparable with published values and between panels.

The primary detector of each control is the one where it is brightest
(the highest 99th percentile); the positive and negative populations are
split at the density minimum on that detector, on an arcsinh scale with
cofactor 150, or at the median where none is found.
